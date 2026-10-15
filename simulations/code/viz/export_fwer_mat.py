#!/usr/bin/env python3
"""Convert the fwer.py beta-sweep cells (a=1, no nugget) into per-beta `results` .mat files that
paper/code/genFig3.m consumes (one file per beta in paper/data/fwer_sims_1f/, fwer_beta_<b>_a1.mat).

Field mapping (npz cell -> results struct genFig3 reads):
  smooth_p     -> results.smoothPosPvals        naive_p    -> results.naivePosPvals
  obs_max      -> results.observedMaxMass        smooth_p95 -> results.smoothNullMassP95
  naive_p95    -> results.naiveNullMassP95        observed_SA-> results.observedSA (nSim x 1)
  k_smooth/n   -> results.fwerPosSmooth           k_naive/n  -> results.fwerPosNaive
  pval         -> results.simcfg.beta             (+ results.n, kSmooth, kNaive for binomial CIs)
Per-sim arrays are written for cells whose partial file carries them (in the cached outputs: the
exemplar betas 1, 29, 50); counts-only cells carry the scalar FWER + median SA (enough for genFig3's
FWER-vs-beta + beta->SA panels). empiricalSA is the per-subject-hemi Moran's I of the real DVs, computed
with the SAME 15 mm geodesic-Gaussian weights as the sim, so the beta-hat crossing in genFig3 is
consistent with observedSA.

  python export_fwer_mat.py            # -> <repo>/paper/data/fwer_sims_1f/  (override with SENS_MAT_OUT)
"""
import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _sys.path.insert(0, _os.path.join(_HERE, ".."))
import glob, numpy as np
from scipy.io import loadmat, savemat
from scipy.spatial import cKDTree
import config as C, sensitivity_core as SC

OUT = _os.environ.get("SENS_MAT_OUT") or _os.path.join(C.PAPER_DATA, "fwer_sims_1f")
_os.makedirs(OUT, exist_ok=True); SA_W = 15.0

# ---- empiricalSA: per subject-hemi Moran's I of the real DVs (geodesic 15 mm, matches the sim) ----
off = np.array(C.FT_OFFSET); base = C.FS_SURF; geo = {}
for h in ("lh", "rh"):
    vp, _ = SC.read_fs_surface(f"{base}/{h}.pial"); vp = vp + off
    sp, _ = SC.read_fs_surface(f"{base}/{h}.sphere.reg"); sp = sp + off
    geo[h] = dict(sph=sp, tree=cKDTree(vp))
E = SC.load_elec(); empSA = []
for x in E:
    cp = np.asarray(x.chanpos, float); st = np.asarray(getattr(x, "stat", np.full(len(cp), np.nan)), float).ravel()
    ok = np.isfinite(st) & np.all(np.isfinite(cp), 1); cp, st = cp[ok], st[ok]
    for h in ("lh", "rh"):
        m = (cp[:, 0] < 0) if h == "lh" else (cp[:, 0] >= 0); elec = cp[m]; y = st[m]
        if len(elec) < 3: continue
        g = geo[h]; esph = g["sph"][g["tree"].query(elec)[1]]
        dee = SC.geodesic_mm(esph, esph); We = np.exp(-dee**2/(2*SA_W**2)); np.fill_diagonal(We, 0.0)
        empSA.append(SC._moran(y, We, float(We.sum())))
empSA = np.array(empSA, float)
print(f"empiricalSA: {np.sum(np.isfinite(empSA))} subj-hemis, median={np.nanmedian(empSA):.4f}", flush=True)

# ---- per-cell -> per-beta results .mat ----
SRC = _os.environ.get("SENS_FWER_BETA_DIR", _os.path.join(C.OUT_DIR, "fwer_beta"))   # cells dir (override for a downloaded folder)
cells = sorted(glob.glob(_os.path.join(SRC, "a1_beta*_chunk00.npz")))
if not cells: raise SystemExit(f"no a=1 beta cells in {SRC}")
PS = {"smooth_p": "smoothPosPvals", "naive_p": "naivePosPvals", "obs_max": "observedMaxMass",
      "smooth_p95": "smoothNullMassP95", "naive_p95": "naiveNullMassP95"}
n_persim = 0
for f in cells:
    d = np.load(f); beta = int(d["pval"]); n = int(d["n"])
    res = {"simcfg": {"beta": float(beta), "kModes": float(_os.environ.get("SENS_MESH_EIG_TAG", "5000"))},
           "fwerPosSmooth": float(d["k_smooth"])/n, "fwerPosNaive": float(d["k_naive"])/n,
           "n": float(n), "kSmooth": float(d["k_smooth"]), "kNaive": float(d["k_naive"]),
           "empiricalSA": empSA}
    if "smooth_p" in d.files:                                     # cell with full per-sim arrays
        for k, v in PS.items(): res[v] = np.asarray(d[k], float).reshape(-1, 1)
        res["observedSA"] = np.asarray(d["observed_SA"], float).reshape(-1, 1)
        n_persim += 1
    else:                                                         # counts-only cell: scalar SA is enough for FWER/betaSA panels
        res["observedSA"] = np.array([[float(d["SA"])]])
    dst = _os.path.join(OUT, f"fwer_beta_{beta}_a1.mat")
    savemat(dst, {"results": res})
print(f"wrote {len(cells)} results .mat ({n_persim} with per-sim arrays) -> {OUT}", flush=True)
