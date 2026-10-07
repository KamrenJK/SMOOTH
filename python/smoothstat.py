#!/usr/bin/env python3
"""
SMOOTHstat — Python port
========================
Group-level cluster-based permutation inference on human intracranial EEG (ECoG/sEEG) with the
exact eigenmode sign-flip surrogate-map null, ported from the MATLAB reference `smooth/SMOOTHstat.m`.

This is the Python version of SMOOTH. It is **self-contained** -- the only dependencies are numpy
and scipy (a minimal FreeSurfer surface reader is included, so neither mne nor FieldTrip is required).

Per subject, per hemisphere:
  * column-normalized sphere-coordinate Gaussian smoother P (each electrode spreads unit mass),
  * graph-Laplacian eigenbasis Q on the electrode sphere directions,
  * empirical map = P y (NaN where no electrode support).
Group:
  * coverage = #subjects with a finite map at each vertex; mask = coverage < minnbsub,
  * group statistic = vertex-wise one-sample t across subjects,
  * null = per permutation, an INDEPENDENT per-subject eigenmode sign-flip surrogate, rank-rescaled to
    the subject's empirical marginal; cluster mass = sum of supra-threshold t; the null is the max
    positive (and min negative) cluster mass over the whole brain,
  * p = (#{null >= mass} + 1) / (nperm + 1).

Electrode coordinates must be in **fsaverage** space (same frame the MATLAB pipeline uses); they are
snapped to the nearest surface vertex internally. Coordinates are shifted by +[1,1,1] to match
FieldTrip's `ft_read_headshape` frame (the frame the electrodes were localized in).

Usage
-----
    import smoothstat as ss
    surf, nL, nR, adj = ss.load_mesh("fsaverage")            # needs $SUBJECTS_DIR/fsaverage
    subs = [(chanpos_i, stat_i), ...]                        # per subject: (N_i x 3, N_i)
    res  = ss.run_smoothstat(subs, surf, nL, nR, adj, numrand=1000, tail="positive", seed=42)
    print(res["posclusters"])                                # [(mass, p, member_vertices), ...]

FreeSurfer surfaces are found via $SUBJECTS_DIR (or $FREESURFER_HOME/subjects). fsaverage matches the
manuscript; fsaverage5 is a fast, lower-resolution option (used by the simulations in Figures 3-4).
"""
import os, warnings
import numpy as np
# uncovered vertices are NaN by design (no electrode support); the group stat masks them afterwards
warnings.filterwarnings("ignore", message="Mean of empty slice")
warnings.filterwarnings("ignore", message="Degrees of freedom <= 0 for slice.")
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.stats import norm

MM2SPHERE = 1.25                      # sphere-units per mm (measured median pial/spherical edge ratio)
FT_OFFSET = np.array([1.0, 1.0, 1.0])  # ft_read_headshape returns FS coords +[1,1,1]; electrodes live in that frame


# ----------------------------------------------------------------------------- geometry / IO
def default_subjects_dir():
    if os.environ.get("SUBJECTS_DIR"):
        return os.environ["SUBJECTS_DIR"]
    if os.environ.get("FREESURFER_HOME"):
        return os.path.join(os.environ["FREESURFER_HOME"], "subjects")
    return "/Applications/freesurfer/7.3.2/subjects"


def read_surface(fname):
    """Minimal FreeSurfer triangular-surface reader (matches mne.read_surface / FieldTrip's
    ft_read_headshape up to the +[1,1,1] offset the callers add). Returns (vertices, faces)."""
    with open(fname, "rb") as f:
        f.read(3)                                  # magic
        f.readline(); f.readline()                 # comment + blank line
        vnum = int(np.fromfile(f, ">i4", 1)[0]); fnum = int(np.fromfile(f, ">i4", 1)[0])
        rr = np.fromfile(f, ">f4", vnum * 3).reshape(vnum, 3).astype(np.float64)
        tri = np.fromfile(f, ">i4", fnum * 3).reshape(fnum, 3).astype(np.int64)
    return rr, tri


def load_mesh(mesh="fsaverage", subjects_dir=None):
    """Load both hemispheres of a FreeSurfer subject (pial + sphere.reg) and the whole-brain
    (block-diagonal, no cross-hemisphere edges) vertex adjacency. Returns (surf, nL, nR, adj)."""
    sd = subjects_dir or default_subjects_dir()
    surf = {}
    for h in ("lh", "rh"):
        pial, _ = read_surface(f"{sd}/{mesh}/surf/{h}.pial")
        sph, tri = read_surface(f"{sd}/{mesh}/surf/{h}.sphere.reg")
        surf[h] = dict(pial=pial + FT_OFFSET, sph=sph + FT_OFFSET, tri=tri)
    nL, nR = len(surf["lh"]["pial"]), len(surf["rh"]["pial"])
    tri = np.vstack([surf["lh"]["tri"], surf["rh"]["tri"] + nL])
    ii = np.concatenate([tri[:, 0], tri[:, 0], tri[:, 1], tri[:, 1], tri[:, 2], tri[:, 2]])
    jj = np.concatenate([tri[:, 1], tri[:, 2], tri[:, 0], tri[:, 2], tri[:, 0], tri[:, 1]])
    nV = nL + nR
    adj = (csr_matrix((np.ones(len(ii)), (ii, jj)), shape=(nV, nV)) > 0).astype(np.int8)
    return surf, nL, nR, adj


# ----------------------------------------------------------------------------- operators
def graph_eigenbasis(elec_sph, graphsigma_mm=8.0):
    """Normalized graph-Laplacian eigenbasis (Q, lam) on the electrode sphere directions, with
    Gaussian arc-distance weights (sigma = MM2SPHERE * graphsigma_mm). Eigenvalues ascending
    (smooth -> rough), matching SMOOTHsub.m."""
    n = len(elec_sph)
    if n == 0: return np.zeros((0, 0)), np.zeros(0)
    if n == 1: return np.ones((1, 1)), np.zeros(1)
    U = elec_sph / np.linalg.norm(elec_sph, axis=1, keepdims=True)
    D = 100.0 * np.arccos(np.clip(U @ U.T, -1.0, 1.0))
    sig = MM2SPHERE * graphsigma_mm
    W = np.exp(-D**2 / (2 * sig**2)); np.fill_diagonal(W, 0.0)
    d = W.sum(1); d[d == 0] = 1.0; Dm = 1.0 / np.sqrt(d)
    L = np.eye(n) - (Dm[:, None] * W * Dm[None, :]); L = (L + L.T) / 2
    lam, Q = np.linalg.eigh(L)
    return Q, lam


def hemi_smoother(elec, pial, sph, kernelwidth=12.0):
    """Column-normalized Gaussian projector P (nVert x nElec) in sphere coords (each electrode spreads
    unit mass), + the electrodes' sphere positions and a support mask. Matches SMOOTHsub.m."""
    idx = cKDTree(pial).query(elec, k=1)[1]; esph = sph[idx]
    sig_u = MM2SPHERE * (kernelwidth / 2.3548); R = 3 * sig_u
    stree = cKDTree(sph); nV = len(sph); rows, cols, vals = [], [], []
    for e, es in enumerate(esph):
        nb = stree.query_ball_point(es, R); dd = np.linalg.norm(sph[nb] - es, axis=1)
        rows.append(nb); cols.append(np.full(len(nb), e)); vals.append(np.exp(-dd**2 / (2 * sig_u**2)))
    W = csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(nV, len(elec)))
    support = np.asarray((W > 0).sum(1)).ravel() > 0
    col = np.asarray(W.sum(0)).ravel(); col[col == 0] = 1e-12
    P = W.multiply(1.0 / col[None, :]).tocsr()
    return P, esph, support


class Subject:
    """Per-subject geometry-only operators (P, Q) + empirical map, reused across permutations."""
    def __init__(self, chanpos, y, surf, nL, nR, kernelwidth=12.0, graphsigma=8.0):
        self.nL, self.nR = nL, nR; self.H = {}
        for h, m in (("lh", chanpos[:, 0] < 0), ("rh", chanpos[:, 0] >= 0)):
            E, yy = chanpos[m], y[m]
            if len(E) == 0: self.H[h] = None; continue
            P, esph, sup = hemi_smoother(E, surf[h]["pial"], surf[h]["sph"], kernelwidth)
            Q, _ = graph_eigenbasis(esph, graphsigma)
            emp = np.asarray(P @ yy).ravel(); emp[~sup] = np.nan
            self.H[h] = dict(P=P, Q=Q, sup=sup, emp=emp, mc=Q.T @ yy, tgt=np.sort(emp[sup]))

    def empirical(self):
        out = np.full(self.nL + self.nR, np.nan)
        if self.H["lh"]: out[:self.nL] = self.H["lh"]["emp"]
        if self.H["rh"]: out[self.nL:] = self.H["rh"]["emp"]
        return out

    def surrogate(self, rng):
        """One eigenmode sign-flip surrogate map, rank-rescaled to the empirical marginal."""
        out = np.full(self.nL + self.nR, np.nan)
        for h, off in (("lh", 0), ("rh", self.nL)):
            d = self.H[h]
            if d is None: continue
            sf = np.sign(rng.standard_normal(len(d["mc"]))); sf[sf == 0] = 1
            m = np.asarray(d["P"] @ (d["Q"] @ (d["mc"] * sf))).ravel(); m[~d["sup"]] = np.nan
            fin = d["sup"]; order = np.argsort(m[fin]); tmp = m[fin].copy(); tmp[order] = d["tgt"]; m[fin] = tmp
            out[off:off + len(m)] = m
        return out


# ----------------------------------------------------------------------------- group stat + clusters
def group_t(maps, mask, minnb):
    """maps: V x S (NaN where uncovered) -> vertex-wise one-sample t across subjects (NaN where the
    coverage is < minnb or the vertex is masked). Variance uses N-1 (matches MATLAB nanvar(x,0,2))."""
    n = np.sum(~np.isnan(maps), axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        mu = np.nanmean(maps, axis=1); sd = np.sqrt(np.nanvar(maps, axis=1, ddof=1))
        t = mu / (sd / np.sqrt(n))
    t[~np.isfinite(t)] = np.nan
    t[(n < minnb) | mask] = np.nan
    return t


def cluster_stats(t, adj, zpos, zneg):
    """Connected supra-threshold clusters.
    Returns (max pos mass, min neg mass, [(mass, members)] pos, [(mass, members)] neg)."""
    def clusters(sel):
        idx = np.where(sel)[0]
        if idx.size == 0: return []
        nc, lab = connected_components(adj[idx][:, idx], directed=False)
        mass = np.zeros(nc); np.add.at(mass, lab, t[idx])
        return [(mass[c], idx[lab == c]) for c in range(nc)]
    pos = clusters(t >= zpos); neg = clusters(t <= zneg)
    return max((m for m, _ in pos), default=0.0), min((m for m, _ in neg), default=0.0), pos, neg


def run_smoothstat(subs, surf, nL, nR, adj, numrand=1000, clusteralpha=0.05, alpha=0.05,
                   minnbsub=3, tail="positive", seed=42, kernelwidth=12.0, graphsigma=8.0, verbose=True):
    """Run the group-level SMOOTH test.

    subs      : list of (chanpos [N_i x 3, fsaverage coords], y [N_i]) per subject
    surf,nL,nR,adj : from load_mesh(...)
    tail      : 'positive' | 'negative' | 'both'. Cluster p-values are one-tailed within each tail;
                `sig` applies alpha for a one-tailed test, or alpha/2 per tail for 'both'.

    Returns a dict with: stat (group t-map), coverage, mask, posdist/negdist (null cluster-mass
    distributions), posclusters / negclusters [(mass, p, members), ...] (strongest first),
    posclusterslabelmat / negclusterslabelmat, sig (mask at alpha), and the per-subject operators
    (subjects) for downstream SA diagnostics.
    """
    if tail not in ("positive", "negative", "both"):
        raise ValueError(f"tail must be 'positive', 'negative' or 'both' (got {tail!r})")
    S = [Subject(cp, np.asarray(y, float).ravel(), surf, nL, nR, kernelwidth, graphsigma) for cp, y in subs]
    allmaps = np.column_stack([s.empirical() for s in S])
    coverage = np.sum(~np.isnan(allmaps), axis=1); mask = coverage < minnbsub
    zpos, zneg = norm.ppf(1 - clusteralpha), norm.ppf(clusteralpha)
    rng = np.random.default_rng(seed)
    posdist = np.zeros(numrand); negdist = np.zeros(numrand)
    for it in range(numrand):
        if verbose and it % 100 == 0: print(f"  perm {it}/{numrand}", flush=True)
        sm = np.column_stack([s.surrogate(rng) for s in S])
        posdist[it], negdist[it], _, _ = cluster_stats(group_t(sm, mask, minnbsub), adj, zpos, zneg)
    temp = group_t(allmaps, mask, minnbsub)
    _, _, posclus, negclus = cluster_stats(temp, adj, zpos, zneg)
    posclus.sort(key=lambda x: -x[0]); negclus.sort(key=lambda x: x[0])
    sig_thresh = alpha / 2 if tail == "both" else alpha
    sig = np.zeros(len(temp), bool)

    def score(clus, null_count, mark):
        out, labelmat = [], np.zeros(len(temp), int)
        for r, (m, mem) in enumerate(clus, 1):
            p = (null_count(m) + 1) / (numrand + 1)
            out.append((float(m), float(p), mem)); labelmat[mem] = r
            if mark and p <= sig_thresh: sig[mem] = True
        return out, labelmat

    posprobs, poslabel = score(posclus, lambda m: np.sum(posdist >= m), tail in ("positive", "both"))
    negprobs, neglabel = score(negclus, lambda m: np.sum(negdist <= m), tail in ("negative", "both"))
    return dict(stat=temp, coverage=coverage, mask=mask, posdist=posdist, negdist=negdist,
                posclusters=posprobs, posclusterslabelmat=poslabel,
                negclusters=negprobs, negclusterslabelmat=neglabel, sig=sig, subjects=S,
                alpha=alpha, clusteralpha=clusteralpha, tail=tail)


def moran_sa(res, adj, nsurr=10, seed=7):
    """Per-subject Moran's I of the empirical vs surrogate maps (restricted to each subject's covered
    vertices), for the surrogate-SA validation plot (cf. SMOOTHdiag.m). Returns (emp, surr_mean, surr_sd)."""
    rng = np.random.default_rng(seed); S = res["subjects"]
    def moran(x):
        vidx = np.where(np.isfinite(x))[0]
        if len(vidx) < 2: return np.nan
        xr = x[vidx] - x[vidx].mean(); W = adj[vidx][:, vidx]
        S0 = float(W.sum()); den = float((xr**2).sum())
        if S0 == 0 or den <= 0: return np.nan
        return (len(vidx) / S0) * float(xr @ (W @ xr)) / den
    emp = np.array([moran(s.empirical()) for s in S])
    null = np.array([[moran(s.surrogate(rng)) for _ in range(nsurr)] for s in S])
    return emp, np.nanmean(null, 1), np.nanstd(null, 1)


if __name__ == "__main__":
    # smoke test on the bundled demo dataset (fsaverage5 for speed)
    import sys
    from scipy.io import loadmat
    demo = os.path.join(os.path.dirname(__file__), "..", "demo", "source.mat")
    src = loadmat(demo, squeeze_me=True, struct_as_record=False)["source"]
    subs = [(np.asarray(s.elec.chanpos, float), np.asarray(s.stat, float).ravel()) for s in np.atleast_1d(src)]
    mesh = sys.argv[1] if len(sys.argv) > 1 else "fsaverage5"
    print(f"loaded {len(subs)} subjects; running SMOOTHstat on {mesh} ...")
    surf, nL, nR, adj = load_mesh(mesh)
    res = run_smoothstat(subs, surf, nL, nR, adj, numrand=200, tail="positive", seed=42)
    print(f"coverage: max {int(res['coverage'].max())} subjects | testable vertices {int((~res['mask']).sum())}")
    for m, p, mem in res["posclusters"][:5]:
        print(f"  cluster: mass={m:8.1f}  p={p:.4f}  vertices={len(mem)}")
