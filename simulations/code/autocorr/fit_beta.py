#!/usr/bin/env python3
"""Reproduce the empirical 1/f-GRF slope: fit the model electrode correlation rho(d) to the pooled
within-subject empirical rho(d) at SENS_MESH_EIG_TAG modes, WITH a nugget (fit a, beta -> beta=29,
a=0.53, the config-locked values) and WITHOUT (a=1 -> the adversarial pure-1/f slope). Reports R^2 and
plots both over the empirical curve.

Needs the empirical DVs (config.ELEC_FILE, default <repo>/demo/source.mat). It is a provenance/repro
script -- BETA is already locked in config; the sweeps do not depend on it.
  SENS_MESH_EIG_TAG=5000 python fit_beta.py"""
import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _sys.path.insert(0, _os.path.join(_HERE, ".."))
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import cs_common as CC, sensitivity_core as SC, config as C

try:
    subs = CC.load_socialcog(None)
except Exception as e:
    raise SystemExit(f"empirical DVs not reachable ({e}); set SENS_ELEC to the electrode file (default <repo>/demo/source.mat). "
                     f"(BETA={C.BETA} is already locked in config; the sweeps do not need this fit.)")

cen, rho, ccnt, _ = CC.empirical_rho(subs); valid = ccnt > 200
dv, rv = cen[valid], rho[valid]; SST = float(np.sum((rv-rv.mean())**2))
off = np.array(C.FT_OFFSET); base = C.FS_SURF; geo = {}
for h in ("lh", "rh"):
    v, _ = SC.read_fs_surface(f"{base}/{h}.pial"); v = v+off
    s, _ = SC.read_fs_surface(f"{base}/{h}.sphere.reg"); s = s+off
    lam, Q = SC.load_mesh_eig(h); geo[h] = dict(tree=cKDTree(v), s=s, lam=lam, Q=Q)
KM = geo["lh"]["Q"].shape[1]; print(f"modes = {KM} (SENS_MESH_EIG_TAG={_os.environ.get('SENS_MESH_EIG_TAG','')})")
edges = np.arange(0, 63, 3); ecen = edges[:-1]+1.5; elec_geo = []
for cp, _z in subs:
    for h in ("lh", "rh"):
        m = (cp[:, 0] < 0) if h == "lh" else (cp[:, 0] >= 0); elec = cp[m]
        if len(elec) < 3: continue
        g = geo[h]; idx = g["tree"].query(elec, k=1)[1]; esph = g["s"][idx]
        d = SC.geodesic_mm(esph, esph); iu, ju = np.triu_indices(len(elec), 1); elec_geo.append((g["Q"][idx], g["lam"], d[iu, ju]))

def model(beta):
    cs = np.zeros(len(ecen)); ct = np.zeros(len(ecen))
    for Qe, lam, dd in elec_geo:
        Pw = 1.0/np.power(lam+np.median(lam), beta); Pw[0] = 0.0; Cf = (Qe*Pw[None, :])@Qe.T
        dn = np.sqrt(np.clip(np.diag(Cf), 1e-30, None)); Cf = Cf/(dn[:, None]*dn[None, :])
        iu, ju = np.triu_indices(Cf.shape[0], 1); bi = np.digitize(dd, edges)-1; ok = (bi >= 0) & (bi < len(ecen))
        np.add.at(cs, bi[ok], Cf[iu, ju][ok]); np.add.at(ct, bi[ok], 1.0)
    return cs/np.maximum(ct, 1)

betas = np.arange(1, 40.1, 0.5); wn = nn = None
for b in betas:
    m = np.interp(dv, ecen, model(b))
    a = float(np.clip(np.sum(rv*m)/np.sum(m*m), 0.1, 1.0)); sse_w = float(np.sum((rv-a*m)**2)); sse_n = float(np.sum((rv-m)**2))
    if wn is None or sse_w < wn[1]: wn = (b, sse_w, a)
    if nn is None or sse_n < nn[1]: nn = (b, sse_n)
bW, _, aW = wn; bN = nn[0]; R2W = 1-wn[1]/SST; R2N = 1-nn[1]/SST
print(f"WITH nugget : beta={bW:.1f}  a={aW:.3f} (nugget={1-aW:.3f})  R2={R2W:.3f}   [config locks BETA={C.BETA}, SA_A={C.SA_A}]")
print(f"NO   nugget : beta={bN:.1f}  a=1.000                        R2={R2N:.3f}")

fig, ax = plt.subplots(figsize=(7.4, 5.2))
ax.plot(dv, rv, "o", color="k", ms=6, label=f"empirical ({C.DATASET}, geodesic)", zorder=5)
ax.plot(ecen, aW*model(bW), "-", color="#c0392b", lw=2.2, label=f"1/f GRF + nugget:  β={bW:.1f}, a={aW:.2f}  (R²={R2W:.3f})")
ax.plot(ecen, model(bN), "--", color="#2f6db0", lw=2.2, label=f"1/f GRF, no nugget (a=1):  β={bN:.1f}  (R²={R2N:.3f})")
ax.axhline(0, color="0.6", lw=0.8); ax.set_xlabel("geodesic distance (mm)"); ax.set_ylabel("spatial correlation  ρ(d)")
ax.set_title(f"1/f GRF fit to empirical SA — {KM} modes, with vs without nugget")
ax.legend(frameon=False, fontsize=9); ax.spines[["top", "right"]].set_visible(False); fig.tight_layout()
figdir = _os.path.join(C.OUT_DIR, "figs"); _os.makedirs(figdir, exist_ok=True)
for ext in ("png", "pdf"): fig.savefig(_os.path.join(figdir, f"beta_fit_{KM}.{ext}"), dpi=150)
print("wrote", _os.path.join(figdir, f"beta_fit_{KM}.png"))
