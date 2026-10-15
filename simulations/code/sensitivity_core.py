#!/usr/bin/env python3
"""Self-contained core for the per-hemisphere sensitivity simulation with the EXACT SMOOTHstat
null (validated bit-exact vs MATLAB). One hemisphere is modelled independently: build geometry
operators once (column-normalized sphere smoother P, graph-Laplacian eigenbasis Q, noise
Cholesky), then for each seed vertex run NITER iterations, each a full NPERM-permutation
eigenmode sign-flip SMOOTH test, scoring TP / FP / Dice.

Numerical primitives (graph_eigenbasis, column-normalized smoother, group-t, cluster stats) are
inlined here; fsaverage5 meshes are read from config.FS_SURF with a minimal reader (read_fs_surface).
Requires ONLY numpy + scipy."""
import os, re, warnings, numpy as np
# empty/degenerate slices are expected (vertices/parcels with no covering subject) and yield NaN, which
# the callers mask -- silence the benign numpy warnings so cluster logs stay readable.
warnings.filterwarnings("ignore", message="Mean of empty slice")
warnings.filterwarnings("ignore", message="Degrees of freedom <= 0 for slice.")
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from scipy.sparse import csr_matrix, block_diag
from scipy.sparse.csgraph import connected_components
from scipy.stats import norm
from scipy.linalg import cholesky
import config as C
MM2SPHERE = 1.25


def geodesic_mm(A_sph, B_sph):
    """Great-circle (geodesic) distance in mm between points on the registration sphere (sphere.reg).
    The eigenbasis/smoother work on a radius-100 sphere in 'sphere units' = MM2SPHERE*mm, so
    arc_mm = (100/MM2SPHERE)*angle. This is the same surface metric SMOOTH smooths in (no sulcus-jumping)."""
    Ua = A_sph/np.linalg.norm(A_sph, axis=1, keepdims=True)
    Ub = B_sph/np.linalg.norm(B_sph, axis=1, keepdims=True)
    return (100.0/MM2SPHERE)*np.arccos(np.clip(Ua @ Ub.T, -1.0, 1.0))


def chordal_mm(A_sph, B_sph=None):
    """Chordal (Euclidean-in-3D-embedding) distance in mm on the registration sphere:
    (100/MM2SPHERE)*||Ua-Ub||. A guaranteed-PSD RBF metric (great-circle substitution is not always PD);
    chordal ~= geodesic for separations << sphere radius (<0.2% at ~13 mm)."""
    Ua = A_sph/np.linalg.norm(A_sph, axis=1, keepdims=True)
    Ub = Ua if B_sph is None else B_sph/np.linalg.norm(B_sph, axis=1, keepdims=True)
    return (100.0/MM2SPHERE)*cdist(Ua, Ub)


def aniso_dist_mm(esph, rho):
    """Transformed-chordal electrode distances for the anisotropic FWER noise: unit sphere vectors scaled
    by the area-preserving B=diag(rho^-1/2, rho^1/2, 1) (R=I, matched axis), then chordal distance.
    L_par/L_perp = rho with geometric-mean correlation length unchanged; rho=1 -> isotropic chordal."""
    Vt = (esph/np.linalg.norm(esph, axis=1, keepdims=True)) * np.array([rho**-0.5, rho**0.5, 1.0])[None, :]
    return (100.0/MM2SPHERE)*cdist(Vt, Vt)


def matern32(d, L):
    """Matern nu=3/2 correlation: (1+r)exp(-r), r=sqrt(3) d/L. Intermediate roughness between the
    exponential (nu=1/2) and the Gaussian/SE (nu=inf) -- a once-differentiable field."""
    r = np.sqrt(3.0)*np.asarray(d)/L
    return (1.0 + r)*np.exp(-r)


_MESH_EIG_CACHE = {}
def load_mesh_eig(hemi):
    """fs5 registration-sphere Laplacian eigenmodes (lam, Q) cached by autocorr/build_eigmodes.py,
    for the 1/f noise model Kee_1f = Q_e diag((lam+c)^-beta) Q_e^T (electrode rows Q_e). The env var
    SENS_MESH_EIG_TAG selects the cache mesh_eig<tag>_<hemi>.npz (config sets tag=5000, the 5000-mode
    basis used throughout)."""
    tag = os.environ.get("SENS_MESH_EIG_TAG", ""); key = (hemi, tag)
    if key not in _MESH_EIG_CACHE:
        import config as _C
        d = np.load(os.path.join(_C.EIG_DIR, f"mesh_eig{tag}_{hemi}.npz"))
        _MESH_EIG_CACHE[key] = (d["lam"].astype(float), d["Q"].astype(float))
    return _MESH_EIG_CACHE[key]


def load_mesh_eig_aniso(hemi, rho):
    """Anisotropic fs5 Laplacian eigenmodes (Gaussian arc weights on B=diag(rho^-1/2,rho^1/2,1)-warped
    sphere coords) cached by autocorr/build_eigmodes.py (aniso), for the anisotropic 1/f noise model.
    rho=1 recovers the isotropic mesh eigenmodes."""
    tag = os.environ.get("SENS_MESH_EIG_TAG", ""); key = (hemi, round(float(rho), 4), tag)
    if key not in _MESH_EIG_CACHE:
        import config as _C
        d = np.load(os.path.join(_C.EIG_DIR, f"mesh_eig{tag}_aniso_{hemi}_rho{key[1]:g}.npz"))
        _MESH_EIG_CACHE[key] = (d["lam"].astype(float), d["Q"].astype(float))
    return _MESH_EIG_CACHE[key]


def load_elec():
    """Per-subject electrode structs from config.ELEC_FILE, each with .chanpos (Nx3, fsaverage surface RAS) and .stat
    (N). Accepts the repository's demo/source.mat layout (`source` struct array: source(i).elec.chanpos,
    source(i).stat) or a flat `E` struct array (E(i).chanpos, E(i).stat)."""
    from types import SimpleNamespace
    from scipy.io import loadmat
    M = loadmat(C.ELEC_FILE, squeeze_me=True, struct_as_record=False)
    if "E" in M: return list(np.atleast_1d(M["E"]))
    return [SimpleNamespace(chanpos=s.elec.chanpos, stat=s.stat) for s in np.atleast_1d(M["source"])]


def read_fs_surface(fname):
    """Minimal FreeSurfer triangular-surface reader (matches mne.read_surface exactly) so the
    cluster code depends only on numpy + scipy (no mne)."""
    with open(fname, "rb") as f:
        f.read(3)                              # magic
        f.readline(); f.readline()             # comment line + blank line
        vnum = int(np.fromfile(f, ">i4", 1)[0]); fnum = int(np.fromfile(f, ">i4", 1)[0])
        rr = np.fromfile(f, ">f4", vnum*3).reshape(vnum, 3).astype(np.float64)
        tri = np.fromfile(f, ">i4", fnum*3).reshape(fnum, 3).astype(np.int64)
    return rr, tri


# ---------- numerical primitives (match SMOOTHstat.m / SMOOTHsub.m) ----------
def graph_eigenbasis(elec_sph, graphsigma_mm):
    n = len(elec_sph)
    if n == 0: return np.zeros((0, 0)), np.zeros(0)
    if n == 1: return np.ones((1, 1)), np.zeros(1)
    U = elec_sph / np.linalg.norm(elec_sph, axis=1, keepdims=True)
    D = 100.0 * np.arccos(np.clip(U @ U.T, -1, 1))
    sig = MM2SPHERE * graphsigma_mm
    W = np.exp(-D**2/(2*sig**2)); np.fill_diagonal(W, 0.0)
    d = W.sum(1); d[d == 0] = 1.0; Dm = 1.0/np.sqrt(d)
    L = np.eye(n) - (Dm[:, None]*W*Dm[None, :]); L = (L+L.T)/2
    lam, Q = np.linalg.eigh(L)
    return Q, lam


def hemi_smoother(elec, pial, sph, kernelwidth):
    """Column-normalized Gaussian projector P (nVert x nElec) in sphere coords; + support mask."""
    idx = cKDTree(pial).query(elec, k=1)[1]; esph = sph[idx]
    sig_u = MM2SPHERE*(kernelwidth/2.3548); R = 3*sig_u
    stree = cKDTree(sph); nV = len(sph); rows, cols, vals = [], [], []
    for e, es in enumerate(esph):
        nb = stree.query_ball_point(es, R); dd = np.linalg.norm(sph[nb]-es, axis=1)
        rows.append(nb); cols.append(np.full(len(nb), e)); vals.append(np.exp(-dd**2/(2*sig_u**2)))
    W = csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(nV, len(elec)))
    support = np.asarray((W > 0).sum(1)).ravel() > 0
    col = np.asarray(W.sum(0)).ravel(); col[col == 0] = 1e-12
    return W.multiply(1.0/col[None, :]).tocsr(), esph, support


def group_t(maps, mask, minnb):
    n = np.sum(~np.isnan(maps), axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        mu = np.nanmean(maps, axis=1); sd = np.sqrt(np.nanvar(maps, axis=1, ddof=1))
        t = mu/(sd/np.sqrt(n))
    t[~np.isfinite(t)] = np.nan; t[(n < minnb) | mask] = np.nan
    return t


def pos_clusters(t, adj, zpos):
    sel = t >= zpos; idx = np.where(sel)[0]
    if idx.size == 0: return []
    nc, lab = connected_components(adj[idx][:, idx], directed=False)
    mass = np.zeros(nc); np.add.at(mass, lab, t[idx])
    return [(mass[c], idx[lab == c]) for c in range(nc)]


def _chol_psd(K):
    """Lower Cholesky, adding escalating jitter if needed. The a=1 (no-nugget) SE covariance used by
    the FWER sweep can be numerically singular when electrodes are closer than L; the nugget (1-a)
    regularizes for a<1 so this is a no-op there."""
    for j in (0.0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2):
        try:
            return cholesky(K if j == 0.0 else K + j*np.eye(len(K)), lower=True)
        except np.linalg.LinAlgError:
            continue
    # rank-deficient (pure-SE a=1 at large L): fall back to an eigh PSD square root (valid sample generator)
    w, V = np.linalg.eigh(K); return V * np.sqrt(np.clip(w, 0.0, None))


def _moran(y, W, Wsum):
    """Moran's I spatial autocorrelation of y under (zero-diagonal) weight matrix W. Used to summarize
    the observed electrode-level SA of a null dataset (the noise-parameter -> SA mapping)."""
    yc = y - y.mean(); denom = float(yc @ yc)
    if denom <= 0 or Wsum <= 0: return np.nan
    return (len(y)/Wsum) * float(yc @ (W @ yc)) / denom


def _reduce(acc):
    """Average a per-iteration accumulator (keys tp/fp/dice/recall/precision/peak_err/n_sig/det_mass)
    into the reported metric dict. precision/peak_err use nanmean (only iters with a detection/TP)."""
    return dict(sens=float(np.mean(acc["tp"])), fpr=float(np.mean(acc["fp"])),
                dice=float(np.mean(acc["dice"])), recall=float(np.mean(acc["recall"])),
                precision=float(np.nanmean(acc["precision"])), peak_err=float(np.nanmean(acc["peak_err"])),
                n_sig=float(np.mean(acc["n_sig"])), det_mass=float(np.mean(acc["det_mass"])))


# ---------- one hemisphere ----------
class HemiSim:
    def __init__(self, hemi, subs_elec=None, sa_a=None, sa_l=None, noise="se", beta=8.0):
        """subs_elec: optional list of per-subject Nx3 electrode arrays (fsaverage surface RAS) to use instead of
        config.ELEC_FILE. sa_a/sa_l: optional noise-SA overrides (default config.SA_A/SA_L). noise: 'se'
        (geodesic squared-exponential) or '1f' (graph-Laplacian power-law, exponent beta).
        The manuscript's sensitivity and specificity runs use noise='1f' with beta=config.BETA (=29), which
        every driver passes explicitly; the 'se' default is used only by geometry-only callers (e.g.
        .coverage)."""
        self.hemi = hemi
        a = C.SA_A if sa_a is None else float(sa_a)
        L = C.SA_L if sa_l is None else float(sa_l)
        self.sa_a, self.sa_l = a, L; self.noise = noise; self.beta = float(beta)
        base = C.FS_SURF
        off = np.array(C.FT_OFFSET)
        self.vpos, faces = read_fs_surface(f"{base}/{hemi}.pial"); self.vpos = self.vpos + off
        sph, _ = read_fs_surface(f"{base}/{hemi}.sphere.reg"); sph = sph + off; self.sph = sph
        nV = len(self.vpos); self.nV = nV
        Ei = faces[:, [0, 1, 2]].ravel(); Ej = faces[:, [1, 2, 0]].ravel()
        self.adj = (csr_matrix((np.ones(len(Ei)*2), (np.r_[Ei, Ej], np.r_[Ej, Ei])), shape=(nV, nV)) > 0).astype(np.int8)
        self.zpos = norm.ppf(1 - C.CLUSTERALPHA)
        self.sig_b = C.BUMP_FWHM/2.3548
        self.sig_b_sph = MM2SPHERE*(C.BUMP_FWHM/2.3548)   # bump width in sphere units (geodesic; matches smoother sig_u)
        self.sig_thresh = C.ALPHA/2 if C.TAIL == "both" else C.ALPHA
        # per-subject electrode sets: config.ELEC_FILE (default) or a caller-supplied list
        if subs_elec is None:
            cps = [np.asarray(x.chanpos, float) for x in load_elec()]
        else:
            cps = [np.asarray(cp, float) for cp in subs_elec]
        self.ops = []
        for cp in cps:
            if cp.ndim != 2 or len(cp) == 0: continue
            m = (cp[:, 0] < 0) if hemi == "lh" else (cp[:, 0] >= 0)
            elec = cp[m]
            if len(elec) == 0: continue
            P, esph, sup = hemi_smoother(elec, self.vpos, sph, C.KERNELWIDTH)
            Q, _ = graph_eigenbasis(esph, C.GRAPHSIGMA)
            if self.noise == "1f":                                                          # graph-Laplacian power-law noise SA
                vidx = cKDTree(self.vpos).query(elec, k=1)[1]; lam_m, Qm = load_mesh_eig(hemi); Qe = Qm[vidx]
                Pw = 1.0/np.power(lam_m + np.median(lam_m), self.beta); Pw[0] = 0.0
                C1f = (Qe*Pw[None, :]) @ Qe.T; dn = np.sqrt(np.clip(np.diag(C1f), 1e-30, None))
                Kee = a*(C1f/(dn[:, None]*dn[None, :])) + (1-a)*np.eye(len(elec))
            else:
                Kee = a*np.exp(-geodesic_mm(esph, esph)**2/(2*L**2)) + (1-a)*np.eye(len(elec))   # geodesic SE noise SA
            self.ops.append(dict(elec=elec, esph=esph, P=P, Q=Q, sup=sup, chol=_chol_psd(Kee), n=len(elec)))
        # coverage + seeds (geometry-only, fixed across iterations)
        supp = np.column_stack([o["sup"] for o in self.ops])
        self.coverage = supp.sum(1); self.mask = self.coverage < C.MINNBSUB
        self.seeds = np.where(self.coverage >= C.MINNBSUB)[0]
        self._init_parc()

    _PARC_BAD = re.compile(r"unknown|corpuscallosum|background|medial.?wall|\?\?\?|^$", re.I)  # NB: not bare 'medial' (medialorbitofrontal is a real DK parcel)

    def _init_parc(self):
        """Load the ROI/atlas parcellation (config.ATLAS: aparc | schaefer300) for the ROI comparison
        arm: per-vertex parcel id (medial-wall/background/unknown -> -1), per-electrode parcel id
        (nearest pial vertex), and one representative on-surface vertex per parcel (parcel-centroid
        vertex, for localization error)."""
        z = np.load(os.path.join(C.ATLAS_DIR, f"{C.ATLAS}_fs5.npz"), allow_pickle=True)
        lab = z[f"{self.hemi}_labels"].astype(int); names = [str(x) for x in z[f"{self.hemi}_names"]]
        bad = [i for i, nm in enumerate(names) if self._PARC_BAD.search(nm)]
        pv = lab.copy(); pv[np.isin(lab, bad) | (lab < 0)] = -1
        self.parc_vert = pv
        self.parc_ids = np.array(sorted(set(pv[pv >= 0].tolist())))
        self.parc_col = {int(p): i for i, p in enumerate(self.parc_ids)}
        self.parc_rep = {}
        for pid in self.parc_ids:
            vi = np.where(pv == pid)[0]; c = self.sph[vi].mean(0)
            self.parc_rep[int(pid)] = int(vi[np.argmin(np.linalg.norm(self.sph[vi] - c, axis=1))])
        vt = cKDTree(self.vpos)
        for o in self.ops:
            o["eparc"] = pv[vt.query(o["elec"])[1]]                           # parcel id per electrode

    def _detect(self, seed, rng, nperm, rroi=None):
        sp = self.vpos[seed]; sp_sph = self.sph[seed]; emp = []; mc = []; tgt = []; ys = []
        for o in self.ops:
            bump = C.A0*np.exp(-np.sum((o["esph"]-sp_sph)**2, 1)/(2*self.sig_b_sph**2))  # geodesic (sphere) bump, not Euclidean-3D
            y = bump + o["chol"] @ rng.standard_normal(o["n"]); ys.append(y)             # electrode DV (shared by SMOOTH + ROI arms)
            m = np.asarray(o["P"] @ y).ravel(); m[~o["sup"]] = np.nan
            emp.append(m); mc.append(o["Q"].T @ y); tgt.append(np.sort(m[o["sup"]]))
        t = group_t(np.column_stack(emp), self.mask, C.MINNBSUB)
        obs = pos_clusters(t, self.adj, self.zpos)
        posdist = np.empty(nperm)
        for p in range(nperm):
            cols = []
            for i, o in enumerate(self.ops):
                sf = np.sign(rng.standard_normal(len(mc[i]))); sf[sf == 0] = 1
                m = np.asarray(o["P"] @ (o["Q"] @ (mc[i]*sf))).ravel(); m[~o["sup"]] = np.nan
                fin = o["sup"]; order = np.argsort(m[fin]); tmp = m[fin].copy(); tmp[order] = tgt[i]; m[fin] = tmp
                cols.append(m)
            pt = group_t(np.column_stack(cols), self.mask, C.MINNBSUB)
            pc = pos_clusters(pt, self.adj, self.zpos)
            posdist[p] = max((mm for mm, _ in pc), default=0.0)
        # score
        dvert = np.linalg.norm(self.sph - sp_sph, axis=1)                     # sphere (geodesic) distance, matches the bump
        hit = dvert < MM2SPHERE*C.R_HIT; bsup = dvert < MM2SPHERE*(C.BUMP_FWHM/2)
        sig = np.zeros(self.nV, bool); is_tp = is_fp = False; nsig = 0
        seed_mass = 0.0; peak_err = np.nan
        for mass, mem in obs:
            if (np.sum(posdist >= mass)+1)/(nperm+1) <= self.sig_thresh:
                sig[mem] = True; nsig += 1
                if hit[mem].any():
                    is_tp = True
                    if mass > seed_mass:                        # the seed-overlapping cluster
                        seed_mass = mass; peak_err = dvert[mem[np.argmax(t[mem])]]/MM2SPHERE  # back to mm-equiv
                else:
                    is_fp = True
        inter = float(np.sum(sig & bsup)); ns = float(sig.sum()); nb = float(bsup.sum())
        out = dict(tp=float(is_tp), fp=float(is_fp), dice=2*inter/max(ns+nb, 1),
                   recall=inter/max(nb, 1), precision=(inter/ns if ns > 0 else np.nan),
                   peak_err=peak_err, n_sig=float(nsig), det_mass=seed_mass)
        roi = self._roi_score(ys, sp_sph, hit, bsup, rroi if rroi is not None else rng, nperm)   # ROI arm, independent RNG (SMOOTH stream untouched)
        out.update({"roi_"+k: v for k, v in roi.items()})
        return out

    def _roi_score(self, ys, sp_sph, hit, bsup, rng, nperm):
        """ROI/atlas (Desikan-Killiany) arm on the SAME planted electrode DVs `ys`: average each
        subject's electrodes within each parcel -> per-parcel one-sample t across subjects -> sign-flip
        max-statistic permutation FWER across parcels (positive tail, same alpha as SMOOTH). Significant
        parcels -> vertex mask, scored identically to the SMOOTH arm (TP/FP/dice/recall/precision/
        localization). Standard parcel-then-test alternative; the head-to-head is fair (identical
        bump+noise, matched alpha)."""
        PM = np.full((len(self.ops), len(self.parc_ids)), np.nan)           # per-subject x per-parcel mean DV
        for si, o in enumerate(self.ops):
            ep = o["eparc"]; yv = ys[si]
            for pid in np.unique(ep[ep >= 0]):
                PM[si, self.parc_col[int(pid)]] = yv[ep == pid].mean()
        def ptvec(M):
            n = np.sum(~np.isnan(M), 0)
            with np.errstate(invalid="ignore", divide="ignore"):
                mu = np.nanmean(M, 0); sd = np.sqrt(np.nanvar(M, 0, ddof=1)); t = mu/(sd/np.sqrt(n))
            t[(n < C.MINNBSUB) | ~np.isfinite(t)] = np.nan
            return t
        tobs = ptvec(PM)
        nullmax = np.empty(nperm)
        for p in range(nperm):
            sf = np.sign(rng.standard_normal(len(self.ops))); sf[sf == 0] = 1
            tp = ptvec(PM * sf[:, None]); nullmax[p] = np.nanmax(tp) if np.any(np.isfinite(tp)) else 0.0
        sigv = np.zeros(self.nV, bool); is_tp = is_fp = False; nsig = 0; best_t = 0.0; peak_err = np.nan
        for i, pid in enumerate(self.parc_ids):
            if not np.isfinite(tobs[i]): continue
            if (np.sum(nullmax >= tobs[i]) + 1)/(nperm + 1) > self.sig_thresh: continue
            nsig += 1; pvv = self.parc_vert == pid; sigv |= pvv
            if (pvv & hit).any():
                is_tp = True
                if tobs[i] > best_t:                                         # localize to the seed-overlapping parcel's centroid
                    best_t = tobs[i]
                    peak_err = float(geodesic_mm(sp_sph[None, :], self.sph[[self.parc_rep[int(pid)]]])[0, 0])
            else:
                is_fp = True
        inter = float(np.sum(sigv & bsup)); ns = float(sigv.sum()); nb = float(bsup.sum())
        return dict(tp=float(is_tp), fp=float(is_fp), dice=2*inter/max(ns+nb, 1),
                    recall=inter/max(nb, 1), precision=(inter/ns if ns > 0 else np.nan),
                    peak_err=peak_err, n_sig=float(nsig), det_mass=best_t)

    def run_null(self, rng, nperm):
        """One FWER (specificity) trial: draw INDEPENDENT per-subject SA-matched noise (NO bump; each
        subject an independent GRF ⇒ no consistent effect), run the group SMOOTH test, and report
        whether ANY positive cluster is significant (a false positive) + the smallest cluster p."""
        emp = []; mc = []; tgt = []
        for o in self.ops:
            y = o["chol"] @ rng.standard_normal(o["n"])                       # noise only
            m = np.asarray(o["P"] @ y).ravel(); m[~o["sup"]] = np.nan
            emp.append(m); mc.append(o["Q"].T @ y); tgt.append(np.sort(m[o["sup"]]))
        t = group_t(np.column_stack(emp), self.mask, C.MINNBSUB)
        obs = pos_clusters(t, self.adj, self.zpos)
        posdist = np.empty(nperm)
        for p in range(nperm):
            cols = []
            for i, o in enumerate(self.ops):
                sf = np.sign(rng.standard_normal(len(mc[i]))); sf[sf == 0] = 1
                m = np.asarray(o["P"] @ (o["Q"] @ (mc[i]*sf))).ravel(); m[~o["sup"]] = np.nan
                fin = o["sup"]; order = np.argsort(m[fin]); tmp = m[fin].copy(); tmp[order] = tgt[i]; m[fin] = tmp
                cols.append(m)
            pt = group_t(np.column_stack(cols), self.mask, C.MINNBSUB)
            posdist[p] = max((mm for mm, _ in pos_clusters(pt, self.adj, self.zpos)), default=0.0)
        any_sig = False; minp = 1.0
        for mass, mem in obs:
            pval = (np.sum(posdist >= mass) + 1) / (nperm + 1)
            minp = min(minp, pval)
            if pval <= self.sig_thresh: any_sig = True
        return dict(any_sig=any_sig, min_p=minp, n_clusters=len(obs))

    METRICS = ["sens", "fpr", "dice", "recall", "precision", "peak_err", "n_sig", "det_mass",
               "roi_sens", "roi_fpr", "roi_dice", "roi_recall", "roi_precision", "roi_peak_err",
               "roi_n_sig", "roi_det_mass"]                                   # SMOOTH arm + ROI/atlas arm

    def run_seed(self, seed, niter, nperm):
        hid = 0 if self.hemi == "lh" else 1
        rng = np.random.default_rng([C.BASE_SEED, hid, int(seed)])            # SMOOTH stream
        rroi = np.random.default_rng([C.BASE_SEED, 10 + hid, int(seed)])      # independent ROI-null stream
        K = ["tp", "fp", "dice", "recall", "precision", "peak_err", "n_sig", "det_mass"]
        acc = {k: [] for k in K}; racc = {k: [] for k in K}                   # SMOOTH + ROI accumulators
        for _ in range(niter):
            d = self._detect(seed, rng, nperm, rroi)
            for k in K: acc[k].append(d[k]); racc[k].append(d["roi_"+k])
        with np.errstate(invalid="ignore"):
            return {**_reduce(acc), **{"roi_"+m: v for m, v in _reduce(racc).items()}}

    def run_seed_roi_only(self, seed, niter, nperm):
        """Fast preview: only the ROI/atlas arm (draw the same planted bump+noise, skip the expensive
        SMOOTH permutation null). Valid MC estimate of the ROI sensitivity map; NOT per-iteration paired
        with the full run's SMOOTH stream (use the full run for the official paired maps)."""
        hid = 0 if self.hemi == "lh" else 1
        rng = np.random.default_rng([C.BASE_SEED, hid, int(seed)])            # same noise seed as the full run
        rroi = np.random.default_rng([C.BASE_SEED, 10 + hid, int(seed)])
        sp_sph = self.sph[seed]; dvert = np.linalg.norm(self.sph - sp_sph, axis=1)
        hit = dvert < MM2SPHERE*C.R_HIT; bsup = dvert < MM2SPHERE*(C.BUMP_FWHM/2)
        K = ["tp", "fp", "dice", "recall", "precision", "peak_err", "n_sig", "det_mass"]
        racc = {k: [] for k in K}
        for _ in range(niter):
            ys = [C.A0*np.exp(-np.sum((o["esph"]-sp_sph)**2, 1)/(2*self.sig_b_sph**2))
                  + o["chol"] @ rng.standard_normal(o["n"]) for o in self.ops]
            d = self._roi_score(ys, sp_sph, hit, bsup, rroi, nperm)
            for k in K: racc[k].append(d[k])
        with np.errstate(invalid="ignore"):
            return {"roi_"+m: v for m, v in _reduce(racc).items()}


# ---------- whole brain (both hemispheres in ONE test) ----------
class BrainSim:
    """Whole-brain twin of HemiSim: both hemispheres go through ONE SMOOTH cluster test with a SINGLE
    permutation null (max cluster mass over the WHOLE brain), matching how SMOOTHstat is actually run.
    Clusters cannot cross hemispheres (block-diagonal adjacency), but the null is pooled -- so the
    whole-brain FWER is correct (per-hemi-then-OR would inflate it ~2x). Same numerics/geodesic model
    as HemiSim; ~2x the per-call cost of one hemisphere. (No ROI arm -- ROI comparison is per-hemi.)"""
    METRICS = ["sens", "fpr", "dice", "recall", "precision", "peak_err", "n_sig", "det_mass"]

    def __init__(self, subs_elec=None, sa_a=None, sa_l=None, noise="se", rho=1.0, q=1.0, beta=8.0, roi_atlases=None):
        a = C.SA_A if sa_a is None else float(sa_a); L = C.SA_L if sa_l is None else float(sa_l)
        self.noise = noise; self.rho = float(rho); self.q = float(q); self.beta = float(beta)   # noise model + aniso/nonstat/1f-exponent params
        self.roi_atlases = list(roi_atlases) if roi_atlases else []                             # ROI-comparator atlases (aparc|schaefer300) for run_null
        self._sqrt_a = np.sqrt(a); self._sqrt_nug = np.sqrt(max(1.0-a, 0.0))          # empirical-nugget split (nonstat mixture draw)
        base = C.FS_SURF; off = np.array(C.FT_OFFSET)
        self.geom = {}; adjs = []
        for h in ("lh", "rh"):
            vpos, faces = read_fs_surface(f"{base}/{h}.pial"); vpos = vpos + off
            sph, _ = read_fs_surface(f"{base}/{h}.sphere.reg"); sph = sph + off
            nV = len(vpos); Ei = faces[:, [0, 1, 2]].ravel(); Ej = faces[:, [1, 2, 0]].ravel()
            adjs.append((csr_matrix((np.ones(len(Ei)*2), (np.r_[Ei, Ej], np.r_[Ej, Ei])), shape=(nV, nV)) > 0).astype(np.int8))
            self.geom[h] = dict(vpos=vpos, sph=sph, nV=nV)
        self.nL = self.geom["lh"]["nV"]; self.nR = self.geom["rh"]["nV"]; self.nV = self.nL + self.nR
        self.adj = block_diag(adjs).tocsr()                                  # no cross-hemisphere edges
        self.sph = np.vstack([self.geom["lh"]["sph"], self.geom["rh"]["sph"]])
        ap = np.concatenate([self.geom[h]["vpos"][:, 1] for h in ("lh", "rh")])       # fixed anatomical A-P bounds (nonstat weight)
        self._ap_lo, self._ap_hi = np.percentile(ap, [5, 95])
        self.zpos = norm.ppf(1 - C.CLUSTERALPHA); self.sig_thresh = C.ALPHA/2 if C.TAIL == "both" else C.ALPHA
        self.sig_b_sph = MM2SPHERE*(C.BUMP_FWHM/2.3548)
        if subs_elec is None:
            cps = [np.asarray(x.chanpos, float) for x in load_elec()]
        else:
            cps = [np.asarray(cp, float) for cp in subs_elec]
        self.subs = []
        for cp in cps:
            if cp.ndim != 2 or len(cp) == 0: continue
            sub = {}
            for h in ("lh", "rh"):
                m = (cp[:, 0] < 0) if h == "lh" else (cp[:, 0] >= 0); elec = cp[m]
                if len(elec) == 0: sub[h] = None; continue
                g = self.geom[h]
                P, esph, sup = hemi_smoother(elec, g["vpos"], g["sph"], C.KERNELWIDTH)
                Q, _ = graph_eigenbasis(esph, C.GRAPHSIGMA)
                dee = geodesic_mm(esph, esph)
                We = np.exp(-dee**2/(2*15.0**2)); np.fill_diagonal(We, 0.0)             # electrode SA weights (Moran's I)
                o = dict(esph=esph, P=P, Q=Q, sup=sup, n=len(elec), We=We, Wesum=float(We.sum()))
                if self.noise == "aniso":                                              # SE on transformed (anisotropic) chordal coords
                    o["chol"] = _chol_psd(a*np.exp(-aniso_dist_mm(esph, self.rho)**2/(2*L**2)) + (1-a)*np.eye(len(elec)))
                elif self.noise == "exp":                                              # exponential (Matern nu=1/2): heavier tail than SE
                    o["chol"] = _chol_psd(a*np.exp(-dee/L) + (1-a)*np.eye(len(elec)))
                elif self.noise == "matern":                                          # Matern nu=3/2: intermediate roughness
                    o["chol"] = _chol_psd(a*matern32(dee, L) + (1-a)*np.eye(len(elec)))
                elif self.noise in ("1f", "aniso1f"):                                 # power-law (1/f) field on the (an)isotropic fs5 mesh Laplacian
                    lam_m, Qm = (load_mesh_eig_aniso(h, self.rho) if self.noise == "aniso1f" else load_mesh_eig(h))
                    vidx = cKDTree(g["vpos"]).query(elec, k=1)[1]
                    Qe = Qm[vidx]; Pw = 1.0/np.power(lam_m + np.median(lam_m), self.beta); Pw[0] = 0.0   # drop DC mode
                    C1f = (Qe*Pw[None, :]) @ Qe.T; dn = np.sqrt(np.clip(np.diag(C1f), 1e-30, None))
                    o["chol"] = _chol_psd(a*(C1f/(dn[:, None]*dn[None, :])) + (1-a)*np.eye(len(elec)))
                elif self.noise == "nonstat1f":                                        # nonstationary 1/f, SA-PRESERVING symmetric mix:
                    # posterior slope beta/sqrt(q) (rough), anterior slope beta*sqrt(q) (smooth), geometric-symmetric
                    # about beta so the mixture's overall SA stays ~fixed while the A-P contrast (ratio q) grows. Cap q so
                    # the smooth end (beta*sqrt(q)) stays out of the single-mode degenerate regime (grid caps at q~4).
                    lam_m, Qm = load_mesh_eig(h); vidx = cKDTree(g["vpos"]).query(elec, k=1)[1]; Qe = Qm[vidx]; cc = np.median(lam_m)
                    def _c1f(b):
                        Pw = 1.0/np.power(lam_m + cc, b); Pw[0] = 0.0; Cf = (Qe*Pw[None, :]) @ Qe.T
                        dn = np.sqrt(np.clip(np.diag(Cf), 1e-30, None)); return Cf/(dn[:, None]*dn[None, :])
                    sq = np.sqrt(self.q)
                    wv = np.clip((self._ap_hi - elec[:, 1])/(self._ap_hi - self._ap_lo), 0.0, 1.0)   # 1 posterior -> 0 anterior
                    o.update(chol_s=_chol_psd(_c1f(self.beta/sq)), chol_l=_chol_psd(_c1f(self.beta*sq)),  # posterior rough (beta/sqrt q), anterior smooth (beta*sqrt q)
                             sw=np.sqrt(wv), sw1=np.sqrt(1.0-wv))                                          # unit-variance, ~SA-preserving mixture (report obs-SA)
                elif self.noise == "nonstat":                                          # two-scale A-P mixture (empirical nugget added at draw)
                    Ls, Ll = L/np.sqrt(self.q), L*np.sqrt(self.q); dc = chordal_mm(esph)
                    wv = np.clip((self._ap_hi - elec[:, 1])/(self._ap_hi - self._ap_lo), 0.0, 1.0)  # 1 posterior(short) -> 0 anterior(long)
                    o.update(chol_s=_chol_psd(np.exp(-dc**2/(2*Ls**2))), chol_l=_chol_psd(np.exp(-dc**2/(2*Ll**2))),
                             sw=np.sqrt(wv), sw1=np.sqrt(1.0-wv))
                else:                                                                   # isotropic geodesic SE (default; empirical-matched)
                    o["chol"] = _chol_psd(a*np.exp(-dee**2/(2*L**2)) + (1-a)*np.eye(len(elec)))
                sub[h] = o
            if sub["lh"] is not None or sub["rh"] is not None: self.subs.append(sub)
        supp = np.zeros(self.nV)
        for sub in self.subs:
            wb = np.zeros(self.nV, bool)
            for h, sl in self._slices():
                if sub[h] is not None: wb[sl] |= sub[h]["sup"]
            supp += wb
        self.coverage = supp; self.mask = supp < C.MINNBSUB; self.seeds = np.where(~self.mask)[0]
        self._parc = {atl: self._init_parc_atlas(atl) for atl in self.roi_atlases}

    _PARC_BAD = re.compile(r"unknown|corpuscallosum|background|medial.?wall|\?\?\?|^$", re.I)

    def _init_parc_atlas(self, atlas):
        """Whole-brain parcellation for the ROI comparator: build a global parcel list over both hemis
        (medial-wall/unknown dropped), and map each subject's electrodes -> global parcel column (nearest
        fs5 vertex's label). Returns dict(ncol, ecols[subject][hemi] = per-electrode global col, -1=drop)."""
        z = np.load(os.path.join(C.ATLAS_DIR, f"{atlas}_fs5.npz"), allow_pickle=True)
        pv = {}; col = 0                                                        # per-hemi per-vertex global col (-1 = bad)
        for h in ("lh", "rh"):
            lab = z[f"{h}_labels"].astype(int); names = [str(x) for x in z[f"{h}_names"]]
            bad = {i for i, nm in enumerate(names) if self._PARC_BAD.search(nm)}
            good = sorted({int(p) for p in np.unique(lab) if p >= 0 and p not in bad})
            m = {p: col + i for i, p in enumerate(good)}; col += len(good)
            vcol = np.full(len(lab), -1)
            for p, c in m.items(): vcol[lab == p] = c
            pv[h] = vcol
        vcol_wb = np.concatenate([pv["lh"], pv["rh"]])                          # whole-brain per-vertex global parcel col (ROI detection)
        ecols = []
        for sub in self.subs:
            ec = {}
            for h in ("lh", "rh"):
                o = sub[h]
                if o is None: ec[h] = None; continue
                vidx = cKDTree(self.geom[h]["sph"]).query(o["esph"])[1]         # electrode -> its fs5 vertex
                ec[h] = pv[h][vidx]
            ecols.append(ec)
        return dict(ncol=col, ecols=ecols, vcol=vcol_wb)

    def _roi_fwer(self, ys, atlas, rng, nperm):
        """ROI FWER on the null electrode DVs: average each subject's electrodes within each parcel ->
        per-parcel one-sample t across subjects -> sign-flip max-statistic permutation (positive tail,
        same sig_thresh as SMOOTH). Returns True if ANY parcel survives (a false positive)."""
        P = self._parc[atlas]; ncol = P["ncol"]; nsub = len(self.subs)
        PM = np.full((nsub, ncol), np.nan)
        for si in range(nsub):
            ssum = np.zeros(ncol); scnt = np.zeros(ncol)
            for h in ("lh", "rh"):
                ec = P["ecols"][si][h]
                if ec is None: continue
                yv = ys[si][h]
                for c in np.unique(ec[ec >= 0]):
                    mm = ec == c; ssum[c] += yv[mm].sum(); scnt[c] += mm.sum()
            g = scnt > 0; PM[si, g] = ssum[g]/scnt[g]
        def ptvec(M):
            n = np.sum(~np.isnan(M), 0)
            with np.errstate(invalid="ignore", divide="ignore"):
                t = np.nanmean(M, 0)/(np.sqrt(np.nanvar(M, 0, ddof=1))/np.sqrt(n))
            t[(n < C.MINNBSUB) | ~np.isfinite(t)] = np.nan; return t
        tobs = ptvec(PM); nullmax = np.empty(nperm)
        for p in range(nperm):
            sf = np.sign(rng.standard_normal(nsub)); sf[sf == 0] = 1
            tp = ptvec(PM*sf[:, None]); nullmax[p] = np.nanmax(tp) if np.any(np.isfinite(tp)) else 0.0
        fin = np.isfinite(tobs)
        if not fin.any(): return False
        pv = (np.array([(np.sum(nullmax >= tobs[i]) + 1)/(nperm + 1) for i in np.where(fin)[0]]))
        return bool(np.min(pv) <= self.sig_thresh)

    def _roi_detect(self, atlas, ys, hit_mask, rng, nperm):
        """Whole-brain ROI detection at a seed: parcel-mean each subject's electrodes -> per-parcel
        one-sample t across subjects -> sign-flip max-statistic permutation (same threshold as SMOOTH).
        A surviving parcel that overlaps the hit region (within R_HIT of the seed) is a true positive;
        a surviving parcel anywhere off-target is a false positive. Returns (tp, fp) as 0/1 floats."""
        P = self._parc[atlas]; ncol = P["ncol"]; nsub = len(self.subs)
        PM = np.full((nsub, ncol), np.nan)
        for si in range(nsub):
            ssum = np.zeros(ncol); scnt = np.zeros(ncol)
            for h in ("lh", "rh"):
                ec = P["ecols"][si][h]
                if ec is None: continue
                yv = ys[si][h]
                for c in np.unique(ec[ec >= 0]):
                    mm = ec == c; ssum[c] += yv[mm].sum(); scnt[c] += mm.sum()
            g = scnt > 0; PM[si, g] = ssum[g]/scnt[g]
        def ptvec(M):
            n = np.sum(~np.isnan(M), 0)
            with np.errstate(invalid="ignore", divide="ignore"):
                t = np.nanmean(M, 0)/(np.sqrt(np.nanvar(M, 0, ddof=1))/np.sqrt(n))
            t[(n < C.MINNBSUB) | ~np.isfinite(t)] = np.nan; return t
        tobs = ptvec(PM); nullmax = np.empty(nperm)
        for p in range(nperm):
            sf = np.sign(rng.standard_normal(nsub)); sf[sf == 0] = 1
            tp_ = ptvec(PM*sf[:, None]); nullmax[p] = np.nanmax(tp_) if np.any(np.isfinite(tp_)) else 0.0
        on = set(int(c) for c in np.unique(P["vcol"][hit_mask & (P["vcol"] >= 0)]))   # parcels overlapping the hit region
        is_tp = is_fp = False
        fin = np.where(np.isfinite(tobs))[0]
        for i in fin:
            pv = (np.sum(nullmax >= tobs[i]) + 1)/(nperm + 1)
            if pv <= self.sig_thresh:
                if int(i) in on: is_tp = True
                else: is_fp = True
        return float(is_tp), float(is_fp)

    def _slices(self):
        return (("lh", slice(0, self.nL)), ("rh", slice(self.nL, self.nV)))

    def _draw_noise(self, o, rng):
        """Per-subject-hemi electrode noise for the active model. 'se'/'aniso': single Cholesky factor.
        'nonstat': mix two SE fields by the A-P weight, then add the empirical nugget (unit total variance)."""
        n = o["n"]
        if self.noise in ("nonstat", "nonstat1f"):
            corr = o["sw"]*(o["chol_s"] @ rng.standard_normal(n)) + o["sw1"]*(o["chol_l"] @ rng.standard_normal(n))
            return self._sqrt_a*corr + self._sqrt_nug*rng.standard_normal(n)
        return o["chol"] @ rng.standard_normal(n)

    def _emp_maps(self, rng, seed=None):
        """Per-subject whole-brain empirical maps: noise everywhere (+ geodesic bump at `seed`, in the
        seed's hemisphere only). Returns (list of V-vectors, list of {hemi:(mc, tgt)})."""
        seed_h = None if seed is None else ("lh" if seed < self.nL else "rh")
        sp_sph = None if seed is None else self.sph[seed]
        emp = []; per = []
        for sub in self.subs:
            m = np.full(self.nV, np.nan); store = {}
            for h, sl in self._slices():
                o = sub[h]
                if o is None: continue
                y = self._draw_noise(o, rng)
                if h == seed_h:
                    y = y + C.A0*np.exp(-np.sum((o["esph"]-sp_sph)**2, 1)/(2*self.sig_b_sph**2))
                mm = np.asarray(o["P"] @ y).ravel(); mm[~o["sup"]] = np.nan
                m[sl] = mm; store[h] = (o["Q"].T @ y, np.sort(mm[o["sup"]]))
            emp.append(m); per.append(store)
        return emp, per

    def _null_dist(self, per, nperm, rng):
        posdist = np.empty(nperm)
        for p in range(nperm):
            cols = []
            for si, sub in enumerate(self.subs):
                m = np.full(self.nV, np.nan)
                for h, sl in self._slices():
                    o = sub[h]
                    if o is None: continue
                    mc, tgt = per[si][h]
                    sf = np.sign(rng.standard_normal(len(mc))); sf[sf == 0] = 1
                    mm = np.asarray(o["P"] @ (o["Q"] @ (mc*sf))).ravel(); mm[~o["sup"]] = np.nan
                    fin = o["sup"]; order = np.argsort(mm[fin]); tmp = mm[fin].copy(); tmp[order] = tgt; mm[fin] = tmp
                    m[sl] = mm
                cols.append(m)
            pt = group_t(np.column_stack(cols), self.mask, C.MINNBSUB)
            posdist[p] = max((mm for mm, _ in pos_clusters(pt, self.adj, self.zpos)), default=0.0)
        return posdist

    def _naive_dist(self, ys, nperm, rng):
        """Naive-shuffle null: destroy spatial structure by permuting each subject's DV<->electrode
        assignment, re-smooth, and take the whole-brain max positive-cluster mass. With autocorrelated
        data this null is too narrow -> anti-conservative (inflated FWER), the contrast to SMOOTH in the FWER figure."""
        nd = np.empty(nperm)
        for p in range(nperm):
            cols = []
            for si, sub in enumerate(self.subs):
                m = np.full(self.nV, np.nan)
                for h, sl in self._slices():
                    o = sub[h]
                    if o is None: continue
                    yv = ys[si][h]; yshuf = yv[rng.permutation(len(yv))]        # naive DV<->electrode shuffle
                    mm = np.asarray(o["P"] @ yshuf).ravel(); mm[~o["sup"]] = np.nan
                    m[sl] = mm
                cols.append(m)
            pt = group_t(np.column_stack(cols), self.mask, C.MINNBSUB)
            nd[p] = max((mm for mm, _ in pos_clusters(pt, self.adj, self.zpos)), default=0.0)
        return nd

    def run_null(self, rng, nperm):
        """One whole-brain FWER trial: independent per-subject SA-matched noise (no bump). Reports the
        observed max positive-cluster mass and its p-value under BOTH the SMOOTH (eigenmode sign-flip)
        null and the NAIVE (DV<->electrode shuffle) null, the 95th-pctile mass thresholds, and the mean
        observed electrode SA (Moran's I) -- everything the FWER figure (QQ / distributions / beta->SA) needs."""
        emp, per = self._emp_maps(rng, seed=None)
        obs = pos_clusters(group_t(np.column_stack(emp), self.mask, C.MINNBSUB), self.adj, self.zpos)
        obs_max = max((mass for mass, _ in obs), default=0.0)
        ys = []                                                                # reconstruct per-subject electrode DV once
        for si, sub in enumerate(self.subs):
            ys.append({h: sub[h]["Q"] @ per[si][h][0] for h, _ in self._slices() if sub[h] is not None})
        sd = self._null_dist(per, nperm, rng)                                  # SMOOTH null max-mass distribution
        nd = self._naive_dist(ys, nperm, rng)                                  # naive-shuffle null max-mass distribution
        sp  = (np.sum(sd >= obs_max) + 1)/(nperm + 1)
        npv = (np.sum(nd >= obs_max) + 1)/(nperm + 1)
        sas = [_moran(ys[si][h], self.subs[si][h]["We"], self.subs[si][h]["Wesum"])
               for si in range(len(self.subs)) for h in ys[si]]
        out = dict(obs_max=float(obs_max), smooth_p=float(sp), naive_p=float(npv),
                   smooth_p95=float(np.percentile(sd, 95)), naive_p95=float(np.percentile(nd, 95)),
                   observed_SA=float(np.nanmean(sas)),
                   smooth_sig=bool(sp <= self.sig_thresh), naive_sig=bool(npv <= self.sig_thresh),
                   any_sig=bool(sp <= self.sig_thresh), min_p=float(sp))        # aliases of smooth_sig / smooth_p
        for atl in self.roi_atlases:                                            # ROI comparators (parcel-then-test max-stat FWER)
            out[f"roi_{atl}_sig"] = self._roi_fwer(ys, atl, rng, nperm)
        return out

    def _detect(self, seed, rng, nperm):
        sp_sph = self.sph[seed]
        emp, per = self._emp_maps(rng, seed=seed)
        t = group_t(np.column_stack(emp), self.mask, C.MINNBSUB)
        obs = pos_clusters(t, self.adj, self.zpos)
        posdist = self._null_dist(per, nperm, rng)
        dvert = np.full(self.nV, np.inf)                                     # geodesic dist within seed's hemi
        sl = slice(0, self.nL) if seed < self.nL else slice(self.nL, self.nV)
        dvert[sl] = np.linalg.norm(self.sph[sl] - sp_sph, axis=1)
        hit = dvert < MM2SPHERE*C.R_HIT; bsup = dvert < MM2SPHERE*(C.BUMP_FWHM/2)
        sig = np.zeros(self.nV, bool); is_tp = is_fp = False; nsig = 0; seed_mass = 0.0; peak_err = np.nan
        for mass, mem in obs:
            if (np.sum(posdist >= mass) + 1)/(nperm + 1) <= self.sig_thresh:
                sig[mem] = True; nsig += 1
                if hit[mem].any():
                    is_tp = True
                    if mass > seed_mass: seed_mass = mass; peak_err = dvert[mem[np.argmax(t[mem])]]/MM2SPHERE
                else:
                    is_fp = True                                             # cluster anywhere off-target (either hemi)
        inter = float(np.sum(sig & bsup)); ns = float(sig.sum()); nb = float(bsup.sum())
        out = dict(tp=float(is_tp), fp=float(is_fp), dice=2*inter/max(ns+nb, 1),
                   recall=inter/max(nb, 1), precision=(inter/ns if ns > 0 else np.nan),
                   peak_err=peak_err, n_sig=float(nsig), det_mass=seed_mass)
        if self.roi_atlases:                                                   # ROI-comparator detection (DK / Schaefer300) at this seed
            ys = {si: {h: self.subs[si][h]["Q"] @ per[si][h][0]
                       for h, _ in self._slices() if self.subs[si][h] is not None} for si in range(len(self.subs))}
            for atl in self.roi_atlases:
                rtp, rfp = self._roi_detect(atl, ys, hit, rng, nperm)
                out[f"roi_{atl}_tp"] = rtp; out[f"roi_{atl}_fp"] = rfp
        return out

    def run_seed(self, seed, niter, nperm):
        rng = np.random.default_rng([C.BASE_SEED, 2, int(seed)])
        keys = ["tp", "fp", "dice", "recall", "precision", "peak_err", "n_sig", "det_mass"]
        for atl in self.roi_atlases: keys += [f"roi_{atl}_tp", f"roi_{atl}_fp"]
        acc = {k: [] for k in keys}
        for _ in range(niter):
            d = self._detect(seed, rng, nperm)
            for k in acc: acc[k].append(d[k])
        with np.errstate(invalid="ignore"):
            out = dict(sens=float(np.mean(acc["tp"])), fpr=float(np.mean(acc["fp"])),
                       dice=float(np.mean(acc["dice"])), recall=float(np.mean(acc["recall"])),
                       precision=float(np.nanmean(acc["precision"])), peak_err=float(np.nanmean(acc["peak_err"])),
                       n_sig=float(np.mean(acc["n_sig"])), det_mass=float(np.mean(acc["det_mass"])))
            for atl in self.roi_atlases:                                        # ROI sensitivity / FPR at this seed
                out[f"roi_{atl}_sens"] = float(np.mean(acc[f"roi_{atl}_tp"]))
                out[f"roi_{atl}_fpr"] = float(np.mean(acc[f"roi_{atl}_fp"]))
        return out
