#!/usr/bin/env python3
"""SA-match the anisotropic 1/f family the SAME way beta=29 was derived for the isotropic field: fit the
1/f slope beta to the empirical geodesic correlation curve rho(d), per rho, on the aniso (radius-graph)
basis. This is needed because beta=29 was fit on the *mesh-face* basis; the aniso field lives on a
different basis whose spectrum makes beta=29 far too smooth (rho=1 -> strongly over-correlated). Fitting
rho(d) per rho holds the correlation *structure* matched to the data while only the anisotropy shape varies.

Mirrors fit_beta.py: for each candidate beta, build the model electrode correlation vs geodesic distance,
jointly fit the amplitude a (correlated fraction / nugget), and pick the beta minimizing SSE to empirical
rho(d). Writes OUT_DIR/aniso_beta_match.npz {rho, beta, a, R2, I_emp}, read by fwer.py's aniso grid.
Cheap: electrode covariances only, NO permutation nulls. Needs the rho=1..10 aniso caches + the dataset DVs.
  SENS_MESH_EIG_TAG=5000 python calibrate_aniso_beta.py
"""
import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _sys.path.insert(0, _os.path.join(_HERE, ".."))
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
import cs_common as CC, sensitivity_core as SC, config as C

RHOS = [int(x) for x in _os.environ.get("SENS_ANISO_RHOS", "1,2,3,4,5,6,7,8,9,10").split(",") if x.strip()]  # env-configurable (only the cached rho are needed)
BETAS = np.arange(1.0, 55.01, 0.5)     # fit_beta.py beta grid, extended to 55 (fixed-a fits need smoother fields)
FIXA = float(_os.environ["SENS_ANISO_FIXA"]) if "SENS_ANISO_FIXA" in _os.environ else None  # hold correlated fraction (e.g. 0.53) and fit beta only
off = np.array(C.FT_OFFSET); base = C.FS_SURF

# ---- DVs from config.ELEC_FILE; empirical rho(d) via the shared (tested) estimator ----
E = SC.load_elec(); subs = []
for x in E:
    cp = np.asarray(x.chanpos, float); st = np.asarray(getattr(x, "stat", np.full(len(cp), np.nan)), float).ravel()
    ok = np.isfinite(st) & np.all(np.isfinite(cp), 1)
    if ok.sum() >= 10: subs.append((cp[ok], st[ok]))
cen, rho, ccnt, _ = CC.empirical_rho(subs); valid = ccnt > 200
dv, rv = cen[valid], rho[valid]; SST = float(np.sum((rv - rv.mean())**2))
edges = np.arange(0, 63, 3); ecen = edges[:-1] + 1.5
print(f"[{C.DATASET}] empirical rho(d): {len(dv)} bins (ccnt>200), rho(3mm)~{rv[0]:.3f}  modes={_os.environ.get('SENS_MESH_EIG_TAG','')}"
      f"  | fit={'FIXED a='+format(FIXA,'g')+' (fit beta only)' if FIXA is not None else 'joint (a floats)'}", flush=True)

# ---- per subject-hemi geometry: electrode -> fs5 vertex idx, and pairwise geodesic distances ----
geo = {}
for h in ("lh", "rh"):
    vpos, _ = SC.read_fs_surface(f"{base}/{h}.pial"); vpos = vpos + off
    sph, _ = SC.read_fs_surface(f"{base}/{h}.sphere.reg"); sph = sph + off
    geo[h] = dict(sph=sph, tree=cKDTree(vpos))
base_units = []                                                  # (h, idx, dd_upper) per subject-hemi
for cp, _st in subs:
    for h in ("lh", "rh"):
        m = (cp[:, 0] < 0) if h == "lh" else (cp[:, 0] >= 0); elec = cp[m]
        if len(elec) < 3: continue
        g = geo[h]; idx = g["tree"].query(elec)[1]; esph = g["sph"][idx]
        d = SC.geodesic_mm(esph, esph); iu, ju = np.triu_indices(len(elec), 1)
        base_units.append((h, idx, d[iu, ju]))

def model_rho(units, beta):
    """Pooled model electrode correlation vs geodesic distance on the aniso basis (units = (lam, Qe, dd))."""
    cs = np.zeros(len(ecen)); ct = np.zeros(len(ecen))
    for lam, Qe, dd in units:
        Pw = 1.0/np.power(lam + np.median(lam), beta); Pw[0] = 0.0
        Cf = (Qe*Pw[None, :]) @ Qe.T; dn = np.sqrt(np.clip(np.diag(Cf), 1e-30, None)); Cf = Cf/(dn[:, None]*dn[None, :])
        iu, ju = np.triu_indices(Cf.shape[0], 1); bi = np.digitize(dd, edges) - 1; ok = (bi >= 0) & (bi < len(ecen))
        np.add.at(cs, bi[ok], Cf[iu, ju][ok]); np.add.at(ct, bi[ok], 1.0)
    return cs/np.maximum(ct, 1)

rho_out, beta_out, a_out, r2_out = [], [], [], []
for rval in RHOS:
    caches = {h: SC.load_mesh_eig_aniso(h, rval) for h in ("lh", "rh")}   # (lam, Q) per hemi, loaded once
    units = [(caches[h][0], caches[h][1][idx], dd) for (h, idx, dd) in base_units]
    best = None
    for b in BETAS:
        m = np.interp(dv, ecen, model_rho(units, b))
        if FIXA is not None:                                                                        # fix nugget at data value, fit beta only
            a = FIXA; sse = float(np.sum((rv - a*m)**2))
        else:                                                                                       # joint amplitude fit (a floats -> basis-dependent)
            a = float(np.clip(np.sum(rv*m)/np.sum(m*m), 0.1, 1.0)); sse = float(np.sum((rv - a*m)**2))
        if best is None or sse < best[1]: best = (b, sse, a)
    b, sse, a = best; r2 = 1 - sse/SST
    rho_out.append(rval); beta_out.append(float(b)); a_out.append(a); r2_out.append(r2)
    print(f"  rho={rval:2d}: beta={b:5.1f}  a={a:.3f}  R2={r2:.3f}", flush=True)

_os.makedirs(C.OUT_DIR, exist_ok=True); dst = _os.path.join(C.OUT_DIR, "aniso_beta_match.npz")
np.savez(dst, rho=np.array(rho_out, float), beta=np.array(beta_out, float), a=np.array(a_out, float), R2=np.array(r2_out, float))
print(f"wrote {dst}\n  rho  : {rho_out}\n  beta : {[round(b,1) for b in beta_out]}\n  a    : {[round(x,3) for x in a_out]}", flush=True)
