#!/usr/bin/env python3
"""Stitch the whole-brain run_sensitivity.py chunk files (sens_wb_c*of*.npz) for the active dataset
into whole-brain per-vertex maps, write per_vertex_maps.npz, print summary stats, and render a preview
comparison (SMOOTH vs Desikan-Killiany vs Schaefer300 sensitivity on the inflated surface).
  SENS_DATASET=socialcog python merge_sensitivity.py"""
import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _sys.path.insert(0, _os.path.join(_HERE, ".."))
import glob, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import config as C, cs_common as CC

MAPS = [("sens", "SMOOTH sensitivity"), ("fpr", "SMOOTH FPR"), ("peak_err", "SMOOTH loc. error (mm)"),
        ("roi_aparc_sens", "Desikan-Killiany sensitivity"), ("roi_aparc_fpr", "Desikan-Killiany FPR"),
        ("roi_schaefer300_sens", "Schaefer300 sensitivity"), ("roi_schaefer300_fpr", "Schaefer300 FPR")]

def main():
    files = sorted(glob.glob(_os.path.join(C.OUT_DIR, "sens_wb_c*of*.npz")))
    if not files: raise SystemExit(f"no sens_wb chunk files in {C.OUT_DIR}")
    res = {}
    for f in files:
        d = np.load(f); ms = [k for k in d.files if k != "seed"]
        for i, s in enumerate(d["seed"]): res[int(s)] = {m: float(d[m][i]) for m in ms}
    lhv, lhf, lhi = CC.load_fs5("lh"); rhv, rhf, rhi = CC.load_fs5("rh")
    nL = len(lhv); nV = nL + len(rhv)
    keys = sorted(next(iter(res.values())).keys())
    M = {k: np.full(nV, np.nan) for k in keys}
    for s, m in res.items():
        for k in keys: M[k][s] = m[k]
    np.savez(_os.path.join(C.OUT_DIR, "per_vertex_maps.npz"), nL=nL, nV=nV, seeds=np.array(sorted(res), int),
             **{k: M[k] for k in keys})
    seeds = np.array(sorted(res), int)
    print(f"[{C.DATASET}] {len(files)} chunks -> {len(res)} seeds ({len(res)/nV*100:.0f}% of {nV} vertices)")
    for k, lbl in MAPS:
        if k in M:
            v = M[k][seeds]; v = v[np.isfinite(v)]
            print(f"  {lbl:34s} median={np.nanmedian(v):.3f}  mean={np.nanmean(v):.3f}")

    # ---- preview: three sensitivity maps on the inflated surface (lateral, both hemis) ----
    def draw(ax, tris, sh, mn, mx, fv, faces, azim, cmap, vmn, vmx):
        vf = fv[faces].mean(1); fc = plt.get_cmap(cmap)(plt.Normalize(vmn, vmx)(vf)).copy()
        fc[~np.isfinite(vf), :3] = 0.82; fc[:, :3] *= sh
        pc = Poly3DCollection(tris, facecolors=fc, linewidths=0, shade=False); pc.set_rasterized(True); ax.add_collection3d(pc)
        ax.set_xlim(mn[0], mx[0]); ax.set_ylim(mn[1], mx[1]); ax.set_zlim(mn[2], mx[2]); ax.set_box_aspect(mx-mn)
        ax.view_init(elev=0, azim=azim); ax.axis("off")
    prep = {}
    for h, iv, fc in (("lh", lhi, lhf), ("rh", rhi, rhf)):
        tris = iv[fc]; tn = np.cross(tris[:, 1]-tris[:, 0], tris[:, 2]-tris[:, 0]); tn /= (np.linalg.norm(tn, axis=1, keepdims=True)+1e-9)
        sh = np.clip(0.55+0.45*(-tn[:, 0] if h == "lh" else tn[:, 0]), 0.4, 1.0)[:, None]
        prep[h] = (tris, sh, iv.min(0), iv.max(0), fc)
    arms = [("sens", "SMOOTH", "viridis"), ("roi_aparc_sens", "Desikan-Killiany", "viridis"),
            ("roi_schaefer300_sens", "Schaefer300", "viridis")]
    fig = plt.figure(figsize=(11, 6.4))
    for r, (k, lbl, cmap) in enumerate(arms):
        for c, (h, azim) in enumerate((("lh", 180), ("lh", 0), ("rh", 180), ("rh", 0))):
            tris, sh, mn, mx, fc = prep[h]
            fv = (M[k][:nL] if h == "lh" else M[k][nL:])
            ax = fig.add_subplot(3, 4, r*4+c+1, projection="3d")
            draw(ax, tris, sh, mn, mx, fv, fc, azim, cmap, 0.0, 1.0)
            if r == 0: ax.set_title(["LH lat", "LH med", "RH lat", "RH med"][c], fontsize=9)
            if c == 0: ax.text2D(-0.05, 0.5, lbl, transform=ax.transAxes, rotation=90, va="center", fontsize=11)
    sm = plt.cm.ScalarMappable(norm=plt.Normalize(0, 1), cmap="viridis")
    fig.colorbar(sm, ax=fig.axes, fraction=0.02, pad=0.01, label="sensitivity")
    fig.suptitle(f"Detection sensitivity — {C.DATASET} (1/f GRF β={C.BETA:g}, a={C.SA_A:g}, whole-brain)", fontsize=12)
    figdir = _os.path.join(C.OUT_DIR, "figs"); _os.makedirs(figdir, exist_ok=True)
    for ext in ("png", "pdf"): fig.savefig(_os.path.join(figdir, "sensitivity_maps.%s" % ext), dpi=150, bbox_inches="tight")
    print("wrote", _os.path.join(figdir, "sensitivity_maps.png"))

if __name__ == "__main__":
    main()
