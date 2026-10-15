#!/usr/bin/env python3
"""Export the 1/f SMOOTH/DK/Schaefer sensitivity benchmark (per_vertex_maps.npz from merge_sensitivity.py)
to benchmark_fig_data.mat for paper/code/genFig4.m (Figure 4: maps + coverage curve + distributions +
paired 2D-histograms). Needs nibabel (curvature) and the FreeSurfer fsaverage5 surfaces.

  SENS_OUT=/path/to/out/socialcog python code/viz/export_benchmark_mat.py   # -> <repo>/paper/data/benchmark_fig_data.mat
"""
import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _sys.path.insert(0, _os.path.join(_HERE, ".."))
import numpy as np, nibabel as nib
from scipy.io import savemat
import sensitivity_core as SC, config as C

FS = C.FS_SURF
pv = np.load(_os.path.join(C.OUT_DIR, "per_vertex_maps.npz"))
nL = int(pv["nL"])
S = pv["sens"].astype(float); DK = pv["roi_aparc_sens"].astype(float); SCH = pv["roi_schaefer300_sens"].astype(float)
base = C.FS_SURF; off = np.array(C.FT_OFFSET)

out = {}; covs = {}
for h in ("lh", "rh"):
    sl = slice(0, nL) if h == "lh" else slice(nL, None)
    vpial, faces = SC.read_fs_surface(f"{base}/{h}.pial"); vpial = vpial + off
    out[f"{h}_pial_pos"] = vpial; out[f"{h}_pial_tri"] = (faces+1).astype(float)
    out[f"{h}_curv"] = nib.freesurfer.read_morph_data(f"{FS}/{h}.curv").astype(float)
    out[f"{h}_smooth"]  = S[sl]
    out[f"{h}_dk_roi"]  = DK[sl];  out[f"{h}_dk_atlas"]  = DK[sl]     # data is already ROI-level per-vertex; atlas==roi here
    out[f"{h}_sch_roi"] = SCH[sl]; out[f"{h}_sch_atlas"] = SCH[sl]
    covs[h] = SC.HemiSim(h).coverage.astype(float)                    # per-vertex coverage (# subjects), geometry-only

cov = np.concatenate([covs["lh"], covs["rh"]])
ok = np.isfinite(S) & np.isfinite(DK) & np.isfinite(SCH)
out.update(seed_cov=cov[ok], seed_smooth=S[ok], seed_dk=DK[ok], seed_sch=SCH[ok], dataset=C.DATASET)

dst = _os.environ.get("SENS_BENCH_OUT") or _os.path.join(C.PAPER_DATA, "benchmark_fig_data.mat")
_os.makedirs(_os.path.dirname(dst), exist_ok=True)
savemat(dst, out); print(f"wrote {dst} | testable seeds = {int(ok.sum())} | cov range {cov[ok].min():.0f}-{cov[ok].max():.0f}")
