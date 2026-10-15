#!/usr/bin/env python3
"""Merge the roi_fdr_full_sweep.py chunk files and export benchmark_fig_data_fdr.mat for
paper/code/genFigS5.m (Supplementary Figure 5: atlas benchmark with BH-FDR-corrected parcels).

Same field layout as benchmark_fig_data.mat (export_benchmark_mat.py), whose mesh, curvature, SMOOTH and
coverage fields are reused unchanged; the DK / Schaefer-300 sensitivity fields are replaced by the FDR
sweep's roi_{atlas}_fdr_sens, and per-vertex specificity (1 - FPR) is added for every method:
  <h>_dk_roi = <h>_dk_atlas = roi_aparc_fdr_sens        <h>_sch_roi = <h>_sch_atlas = roi_schaefer300_fdr_sens
  <h>_smooth_spec = 1 - fpr (per_vertex_maps.npz)       <h>_dk_spec / <h>_sch_spec = 1 - roi_<atlas>_fdr_fpr
  seed_* = the same quantities over testable seeds (finite SMOOTH, DK and Schaefer sensitivity)

  python code/viz/export_benchmark_fdr_mat.py     # -> <repo>/paper/data/benchmark_fig_data_fdr.mat
Overrides: SENS_PER_VERTEX_MAPS, SENS_BENCH_IN (benchmark_fig_data.mat), SENS_BENCH_FDR_OUT.
"""
import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _sys.path.insert(0, _os.path.join(_HERE, ".."))
import glob, numpy as np
from scipy.io import loadmat, savemat
import config as C

ARR = ["roi_aparc_fdr_sens", "roi_aparc_fdr_fpr", "roi_schaefer300_fdr_sens", "roi_schaefer300_fdr_fpr"]
pv = np.load(_os.environ.get("SENS_PER_VERTEX_MAPS") or _os.path.join(C.PAPER_DATA, "per_vertex_maps.npz"))
nL, nV = int(pv["nL"]), int(pv["nV"])

# ---- merge per-worker chunks -> whole-brain per-vertex arrays ----
src = _os.path.join(C.PROJ, "out", "roi_fdr_sweep")
files = sorted(glob.glob(_os.path.join(src, "chunk_*.npz")))
if not files: raise SystemExit(f"no chunk files in {src}")
M = {k: np.full(nV, np.nan) for k in ARR}
for f in files:
    d = np.load(f)
    for k in ARR: M[k][d["seeds"]] = d[k]
np.savez(_os.path.join(C.PROJ, "out", "roi_fdr_merged.npz"), seeds=pv["seeds"], nL=nL, nV=nV, **M)

# ---- benchmark_fig_data.mat layout with FDR atlas arms + specificity ----
base = loadmat(_os.environ.get("SENS_BENCH_IN") or _os.path.join(C.PAPER_DATA, "benchmark_fig_data.mat"))
out = {k: v for k, v in base.items() if not k.startswith("__")}
S = pv["sens"].astype(float); DK = M["roi_aparc_fdr_sens"]; SCH = M["roi_schaefer300_fdr_sens"]
spec = {"smooth": 1 - pv["fpr"].astype(float), "dk": 1 - M["roi_aparc_fdr_fpr"], "sch": 1 - M["roi_schaefer300_fdr_fpr"]}
for h in ("lh", "rh"):
    sl = slice(0, nL) if h == "lh" else slice(nL, None)
    out[f"{h}_dk_roi"] = DK[sl];  out[f"{h}_dk_atlas"] = DK[sl]
    out[f"{h}_sch_roi"] = SCH[sl]; out[f"{h}_sch_atlas"] = SCH[sl]
    for m, v in spec.items(): out[f"{h}_{m}_spec"] = v[sl]
ok = np.isfinite(S) & np.isfinite(DK) & np.isfinite(SCH)
if int(ok.sum()) != np.size(out["seed_cov"]): raise SystemExit("testable-seed set differs from benchmark_fig_data.mat")
out.update(seed_dk=DK[ok], seed_sch=SCH[ok], **{f"seed_{m}_spec": v[ok] for m, v in spec.items()})

dst = _os.environ.get("SENS_BENCH_FDR_OUT") or _os.path.join(C.PAPER_DATA, "benchmark_fig_data_fdr.mat")
savemat(dst, out); print(f"wrote {dst} | {len(files)} chunks | testable seeds = {int(ok.sum())}")
