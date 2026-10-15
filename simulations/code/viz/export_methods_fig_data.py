#!/usr/bin/env python3
"""Export ALL data for the sensitivity-simulation methods figure (Supplementary Figure 4,
paper/code/genFigS4.m), generated with the EXACT 1/f-GRF sim: geodesic Gaussian bump (FWHM 20 mm, peak
SNR A0) + power-law graph-Laplacian SURFACE noise ((lam+c)^-beta on the fs5 Laplacian eigenbasis) +
empirical nugget. The noise panels and the rho(d) fit panel use the 1/f model (beta, a). Renders in
MATLAB from methods_fig_data.mat.

beta and the nugget a are REFIT inline to the empirical geodesic rho(d) on the ACTIVE eigenbasis (the same
electrode-correlation model as autocorr/fit_beta.py), so the illustration is internally consistent at any
mode count; at SENS_MESH_EIG_TAG=5000 this reproduces the config-locked BETA=29, SA_A=0.53. The surface
noise, electrode DVs, group t-map/cluster, and fit panel all use that refit beta.

Reads OUT_DIR/per_vertex_maps.npz (merge_sensitivity.py) for the sensitivity map + seed metrics; inflated
surface + curv are display-only (FreeSurfer install, config.FS_SURF; needs nibabel). Run where
per_vertex_maps.npz lives:
  SENS_MESH_EIG_TAG=5000 SENS_OUT=/path/to/out/socialcog python code/viz/export_methods_fig_data.py
  -> <repo>/paper/data/methods_fig_data.mat   (override with SENS_METHODS_OUT)
"""
import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _sys.path.insert(0, _os.path.join(_HERE, ".."))
import os, numpy as np, nibabel as nib
from scipy.spatial import cKDTree
from scipy.io import savemat
import config as C, cs_common as CC, sensitivity_core as SC

HEMI = "lh"; NSHOW = 5; NPERM = int(os.environ.get("SENS_METHODS_NPERM", "300"))
FS = C.FS_SURF                                                                     # inflated + curv (display only)
rng = np.random.default_rng(7); A0 = C.A0
base = C.FS_SURF; off = np.array(C.FT_OFFSET)
TAG = os.environ.get("SENS_MESH_EIG_TAG", "")

# ---------------- geometry (LH = method hemi; RH carries the second sensitivity map) ----------------
vpial, faces = SC.read_fs_surface(f"{base}/{HEMI}.pial"); vpial = vpial + off
vinfl, _     = SC.read_fs_surface(f"{FS}/{HEMI}.inflated")                          # seed pick + display
curv = nib.freesurfer.read_morph_data(f"{FS}/{HEMI}.curv").astype(float)

# ---------------- sim operators (smoother P, eigenbasis Q, support, coverage, adjacency) ----------------
sim = SC.HemiSim(HEMI, noise="1f", beta=C.BETA)                                     # ops + coverage + adj/mask/zpos
sph = sim.sph; nV = sim.nV; nsub = len(sim.ops)

# ---------------- refit the 1/f slope beta + nugget a to empirical rho(d) on the ACTIVE eigenbasis ----------------
subs = CC.load_socialcog(None)
cen, rho_emp, ccnt, _ = CC.empirical_rho(subs); ok = ccnt > 200
dv, rv = cen[ok], rho_emp[ok]; SST = float(np.sum((rv - rv.mean())**2))
edges = np.arange(0, 63, 3); ecen = edges[:-1] + 1.5
egeo = {}
for h in ("lh", "rh"):
    vh, _ = SC.read_fs_surface(f"{base}/{h}.pial"); vh = vh + off
    sh, _ = SC.read_fs_surface(f"{base}/{h}.sphere.reg"); sh = sh + off
    lh_, Qh = SC.load_mesh_eig(h); egeo[h] = (cKDTree(vh), sh, lh_, Qh)
units = []
for cp, _z in subs:
    for h in ("lh", "rh"):
        m = (cp[:, 0] < 0) if h == "lh" else (cp[:, 0] >= 0); elec = cp[m]
        if len(elec) < 3: continue
        tree, sh, lh_, Qh = egeo[h]; idx = tree.query(elec)[1]; esph = sh[idx]
        d = SC.geodesic_mm(esph, esph); iu, ju = np.triu_indices(len(elec), 1)
        units.append((Qh[idx], lh_, d[iu, ju]))
def rho_model_curve(beta):
    cs = np.zeros(len(ecen)); ct = np.zeros(len(ecen))
    for Qe, lam, dd in units:
        Pw = 1.0/np.power(lam + np.median(lam), beta); Pw[0] = 0.0
        Cf = (Qe*Pw[None, :]) @ Qe.T; dn = np.sqrt(np.clip(np.diag(Cf), 1e-30, None)); Cf = Cf/(dn[:, None]*dn[None, :])
        iu, ju = np.triu_indices(Cf.shape[0], 1); bi = np.digitize(dd, edges) - 1; g = (bi >= 0) & (bi < len(ecen))
        np.add.at(cs, bi[g], Cf[iu, ju][g]); np.add.at(ct, bi[g], 1.0)
    return cs/np.maximum(ct, 1)
best = None
for b in np.arange(1, 40.1, 0.5):
    mm = np.interp(dv, ecen, rho_model_curve(b)); af = float(np.clip(np.sum(rv*mm)/np.sum(mm*mm), 0.1, 1.0))
    sse = float(np.sum((rv - af*mm)**2))
    if best is None or sse < best[1]: best = (b, sse, af)
beta, _, a = best; R2 = 1 - best[1]/SST; rho_curve = a*rho_model_curve(beta)
print(f"[1f fit] modes(tag={TAG or 'default'})={egeo['lh'][3].shape[1]}  beta={beta:.1f} a={a:.3f} R2={R2:.3f}"
      f"   (config BETA={C.BETA}, SA_A={C.SA_A})", flush=True)

# ---------------- surface 1/f noise: per-subject INDEPENDENT power-law field + empirical nugget ----------------
lam_m, Qm = SC.load_mesh_eig(HEMI); Pw = 1.0/np.power(lam_m + np.median(lam_m), beta); Pw[0] = 0.0; sqrtP = np.sqrt(Pw)
def draw_noise():
    s = Qm @ (sqrtP*rng.standard_normal(len(lam_m))); s = (s - s.mean())/(s.std() + 1e-12)   # unit-variance 1/f structure
    return np.sqrt(a)*s + np.sqrt(max(1.0-a, 0.0))*rng.standard_normal(nV)                    # + empirical nugget
noise_surf = np.array([draw_noise() for _ in range(nsub)])                                    # nsub x nV

# ---------------- sensitivity map + seed pick (lateral, high coverage / sensitivity / dice) ----------------
pv = np.load(os.path.join(C.OUT_DIR, "per_vertex_maps.npz")); nL = int(pv["nL"])
sens_wb = pv["sens"].astype(float); dice_wb = pv["dice"].astype(float)
sens = sens_wb[:nL]; dice = dice_wb[:nL]; cov = sim.coverage.astype(float); lat = -vinfl[:, 0]
cand = np.where((cov >= 10) & (cov <= 15) & (sens > 0.85) & (dice > 0.30) & (~sim.mask) &
                (lat > np.percentile(lat, 75)) & (np.abs(vinfl[:, 2]) < 40))[0]
if len(cand) == 0: cand = np.where((cov >= 10) & (cov <= 16) & (sens > 0.7) & (dice > 0.2) & (~sim.mask))[0]
seed = int(cand[np.argmax(dice[cand] + 0.5*sens[cand] + 0.008*lat[cand])])
print(f"seed={seed} coverage={cov[seed]:.0f} sens={sens[seed]:.2f} dice={dice[seed]:.2f} curv={curv[seed]:.2f}", flush=True)

# ---------------- bump + per-subject field (bump+noise) + electrode DVs ----------------
vt = cKDTree(vpial); vi = [vt.query(o["elec"])[1] for o in sim.ops]               # electrode -> nearest pial vertex
sp_sph = sph[seed]
bump_surf = A0*np.exp(-np.sum((sph - sp_sph)**2, 1)/(2*sim.sig_b_sph**2))          # geodesic (sphere) bump
geo_from_seed = SC.geodesic_mm(sph, sph[[seed]]).ravel()
field_surf = bump_surf[None, :] + noise_surf
elec_dv = [bump_surf[vi[i]] + noise_surf[i][vi[i]] for i in range(nsub)]           # surface field sampled at electrodes
nearest = np.array([SC.geodesic_mm(sph[[seed]], o["esph"]).min() for o in sim.ops])
show = np.argsort(nearest)[:NSHOW]
print(f"shown subjects (nearest-electrode geodesic mm): {[(int(s), round(float(nearest[s]),1)) for s in show]}", flush=True)

# ---------------- group SMOOTH test over ALL covering subjects (t-map + significant cluster) ----------------
def group_test(dvs, nperm):
    emp = []; mc = []; tgt = []
    for o, dv_ in zip(sim.ops, dvs):
        mmap = np.asarray(o["P"] @ dv_).ravel(); mmap[~o["sup"]] = np.nan
        emp.append(mmap); mc.append(o["Q"].T @ dv_); tgt.append(np.sort(mmap[o["sup"]]))
    t = SC.group_t(np.column_stack(emp), sim.mask, C.MINNBSUB)
    obs = SC.pos_clusters(t, sim.adj, sim.zpos)
    posdist = np.empty(nperm)
    for p in range(nperm):
        cols = []
        for i, o in enumerate(sim.ops):
            sf = np.sign(rng.standard_normal(len(mc[i]))); sf[sf == 0] = 1
            mmap = np.asarray(o["P"] @ (o["Q"] @ (mc[i]*sf))).ravel(); mmap[~o["sup"]] = np.nan
            fin = o["sup"]; order = np.argsort(mmap[fin]); tmp = mmap[fin].copy(); tmp[order] = tgt[i]; mmap[fin] = tmp
            cols.append(mmap)
        pt = SC.group_t(np.column_stack(cols), sim.mask, C.MINNBSUB)
        posdist[p] = max((mm for mm, _ in SC.pos_clusters(pt, sim.adj, sim.zpos)), default=0.0)
    clmask = np.zeros(nV, bool); sig_p = []
    for mass, mem in obs:
        pval = (np.sum(posdist >= mass) + 1)/(nperm + 1)
        if pval <= C.ALPHA: clmask[mem] = True; sig_p.append(pval)
    return t, clmask, sig_p
print("running group SMOOTH test (all subjects) ...", flush=True)
tmap, cluster_mask, sig_p = group_test(elec_dv, NPERM)
print(f"  t-map max={np.nanmax(tmap):.1f} | sig cluster vtx={int(cluster_mask.sum())} | sig cluster p={sig_p}", flush=True)

# ---------------- coverage/#-electrodes -> sensitivity curves (pooled both hemispheres) ----------------
simR = SC.HemiSim("rh", noise="1f", beta=C.BETA)
allcov = np.concatenate([sim.coverage, simR.coverage]); alls = sens_wb            # sens_wb is [lh|rh], aligned to allcov
covbins = np.arange(3, 22); cbc = 0.5*(covbins[:-1] + covbins[1:])
csm = np.full(len(covbins)-1, np.nan); cse = np.full(len(covbins)-1, np.nan)
for i in range(len(covbins)-1):
    k = (allcov >= covbins[i]) & (allcov < covbins[i+1]) & np.isfinite(alls)
    if k.sum() > 10: csm[i] = alls[k].mean(); cse[i] = 1.96*alls[k].std(ddof=1)/np.sqrt(k.sum())
Rsup = 3*(C.KERNELWIDTH/2.3548)
def nelec_map(s):
    aE = np.vstack([o["esph"] for o in s.ops]); return (SC.geodesic_mm(s.sph, aE) < Rsup).sum(1).astype(float)
coverage_elec = nelec_map(sim); ne_all = np.concatenate([coverage_elec, nelec_map(simR)])
nemax = np.percentile(ne_all[ne_all > 0], 98)
nebins = np.linspace(0, nemax, 13); nec = 0.5*(nebins[:-1] + nebins[1:])
nem = np.full(len(nebins)-1, np.nan); nee = np.full(len(nebins)-1, np.nan)
for i in range(len(nebins)-1):
    k = (ne_all >= nebins[i]) & (ne_all < nebins[i+1]) & np.isfinite(alls)
    if k.sum() > 15: nem[i] = alls[k].mean(); nee[i] = 1.96*alls[k].std(ddof=1)/np.sqrt(k.sum())

# ---------------- coverage berry (all electrodes + subject id) + RH surface/sens ----------------
allE = np.vstack([o["elec"] for o in sim.ops]); allEsub = np.concatenate([[i]*o["n"] for i, o in enumerate(sim.ops)])
vpial_rh, faces_rh = SC.read_fs_surface(f"{base}/rh.pial"); vpial_rh = vpial_rh + off
curv_rh = nib.freesurfer.read_morph_data(f"{FS}/rh.curv").astype(float)

dst = os.environ.get("SENS_METHODS_OUT") or os.path.join(C.PAPER_DATA, "methods_fig_data.mat")
os.makedirs(os.path.dirname(dst), exist_ok=True)
savemat(dst, dict(
    rh_pial_pos=vpial_rh, rh_pial_tri=(faces_rh+1).astype(float), rh_curv=curv_rh, rh_sens_map=sens_wb[nL:],
    hemi=HEMI, pial_pos=vpial, pial_tri=(faces+1).astype(float), infl_pos=vinfl,      # +1 -> MATLAB 1-based tri
    seed=float(seed+1), seed_pos=vpial[seed], seed_pos_infl=vinfl[seed],
    params=dict(A0=A0, FWHM=C.BUMP_FWHM, a=a, beta=beta, sig_b_sph=sim.sig_b_sph,
                kernelwidth=C.KERNELWIDTH, graphsigma=C.GRAPHSIGMA, clusteralpha=C.CLUSTERALPHA,
                alpha=C.ALPHA, minnbsub=C.MINNBSUB, mm2sphere=SC.MM2SPHERE, modes=int(len(lam_m))),
    bump_surf=bump_surf, noise_surf_show=noise_surf[show], field_surf_show=field_surf[show], show_subj=(show+1).astype(float),
    elec_show=np.array([vpial[vi[s]] for s in show], dtype=object),                    # SNAPPED to nearest pial vertex
    elec_show_dv=np.array([elec_dv[s] for s in show], dtype=object),
    elec_all=vpial[np.concatenate(vi)], elec_all_subj=(allEsub+1).astype(float), elec_all_dv=np.concatenate(elec_dv),
    curv=curv, geo_from_seed=geo_from_seed, crop_radius=float(1.5*C.BUMP_FWHM),
    coverage=cov, coverage_elec=coverage_elec, sens_map=sens, tmap=tmap, cluster_mask=cluster_mask.astype(float),
    cov_curve_x=cbc, cov_curve_y=csm, cov_curve_ci=cse, ne_curve_x=nec, ne_curve_y=nem, ne_curve_ci=nee,
    cov_pts=allcov[np.isfinite(alls)], sens_pts=alls[np.isfinite(alls)], cov_max=float(np.nanmax(allcov)),
    ne_pts=ne_all[np.isfinite(alls)], ne_max=float(nemax),
    rho_cen=cen, rho_emp=rho_emp, rho_ok=ok.astype(float), rho_dd=ecen, rho_model=rho_curve,
    spec_lam=lam_m, spec_P=Pw, sa_fit=dict(a=a, beta=beta, R2=R2), noise_model="1f",
))
print("wrote", dst, flush=True)
