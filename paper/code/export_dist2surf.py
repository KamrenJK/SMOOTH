#!/usr/bin/env python3
"""Recompute corrsource{i}.elec.dist2surf (electrode-to-native-pial distance, mm) in paper/data/corrsource.mat.

genFig2.m (Figure 2b) keeps contacts within 7.5 mm of the cortical surface (all ECoG contacts are kept)
before fitting the distance-correlation curves. The dist2surf values originally stored were computed on a
machine whose read_surf.m added 1 mm to every vertex coordinate, i.e. against the native pial shifted by
+[1,1,1]. This script replaces them with the distance to the surface as stored (as nibabel and
smooth_read_surf.m read it). The filter only acts on depth (sEEG) contacts, which are localized from
the CT and so are in that frame. (For subdural grids, a signed-height test cannot tell the two frames
apart at the 1 mm level; it does not matter here, because ECoG contacts are always kept.)

corrsource.mat does not record which FreeSurfer recon belongs to which participant, so each participant's
recon is identified as the one that reproduces the stored dist2surf to < 1e-3 mm, in either the original
+[1,1,1] frame or the corrected one (so the script can be re-run on its own output). The recon trees
($SMOOTH_RAWDATA/data/<id>/recon/freesurfer/surf) are not part of the published repository.

Usage: SMOOTH_RAWDATA=<root> python3 export_dist2surf.py
"""
import os, sys, glob
import numpy as np
import nibabel as nib
import h5py
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
SRC  = os.path.abspath(os.path.join(HERE, "..", "data", "corrsource.mat"))
if not os.environ.get("SMOOTH_RAWDATA"): sys.exit("Set SMOOTH_RAWDATA to the raw-data root (not distributed).")
LEGACY_OFF = np.array([1.0, 1.0, 1.0])                       # frame of the originally stored dist2surf
LIM = 7.5                                                    # genFig2.m inclusion limit (reported only)

trees = {}
for f in sorted(glob.glob(os.path.join(os.environ["SMOOTH_RAWDATA"], "data", "*", "recon", "freesurfer", "surf", "lh.pial"))):
    surf = os.path.dirname(f)
    trees[f] = cKDTree(np.vstack([nib.freesurfer.read_geometry(os.path.join(surf, f"{h}.pial"))[0] for h in ("lh", "rh")]))
print(f"{len(trees)} candidate recons")

with h5py.File(SRC, "r+") as h:                              # in place: same shape and dtype
    refs = h["corrsource"][:, 0]; nchange = 0
    for i, r in enumerate(refs, 1):
        e = h[r]["elec"]; cp = np.asarray(e["nativechanpos"]).T
        old = np.asarray(e["dist2surf"]).ravel(); ecog = np.asarray(e["ecog"]).ravel().astype(bool)
        hits = [(t, off) for t in trees.values() for off in (LEGACY_OFF, 0.0)
                if np.abs(t.query(cp - off)[0] - old).max() < 1e-3]
        if len({id(t) for t, _ in hits}) != 1:
            raise SystemExit(f"corrsource {i}: {len(hits)} matching recons; cannot identify the surface")
        new = hits[0][0].query(cp)[0]
        changed = int((((old <= LIM) | ecog) != ((new <= LIM) | ecog)).sum()); nchange += changed
        e["dist2surf"][...] = new.reshape(e["dist2surf"].shape)
        print(f"corrsource {i:2d}: {len(cp):3d} contacts | median dist2surf {np.median(old):.2f} -> {np.median(new):.2f} mm"
              f" | inclusion at {LIM} mm changes for {changed}")
print(f"wrote {SRC} | contacts whose inclusion changes: {nchange}")
