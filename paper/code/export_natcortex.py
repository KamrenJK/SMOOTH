#!/usr/bin/env python3
"""Cache participant P19's native pial surface + curvature for genFig2_seedcorr (Figure 2, panel a).

The three files that describe this participant use three different orderings:
`corrsource.mat` (correlation matrices + native electrode positions), `demo/source.mat` (fsaverage
electrode positions, used to look up DK parcels), and the FreeSurfer recon trees under $SMOOTH_RAWDATA/data/.
corrsource <-> source was recovered by matching channel labels; source <-> recon by finding the
recon that reproduces the stored `elec.dist2surf` -- which it does to 0.0000 mm, so the
identification is verified rather than assumed.

nibabel returns raw FreeSurfer coordinates; ft_read_headshape adds [1,1,1] and
`corrsource{i}.elec.nativechanpos` lives in that FieldTrip frame, hence FT_OFF. The exact
dist2surf agreement is what confirms the offset.

The recon trees are not part of the published repository, so the surface is cached to a .mat the
figure script loads directly.

P19 (its demo/source.mat index) was chosen from the cohort's ECoG grid participants by rendering
each candidate's native pial and comparing grid regularity and mesh quality; only P19 is exported.
    P19  18 19  lh    43 TPJ  TG:112 AT:16 PST:16 AST:10
(columns: corrsource index, demo/source.mat index, hemisphere, n DK supramarginal/IP contacts, arrays)

The raw-data subject ID that locates P19's recon under $SMOOTH_RAWDATA/data/ is not distributed; it
is read from $SMOOTH_RAWSUBJECT. Everything written uses the pseudonym P19.

Usage: SMOOTH_RAWDATA=<root> SMOOTH_RAWSUBJECT=<id> python3 export_natcortex.py
"""
import os, sys
import numpy as np
import nibabel as nib
import h5py
from scipy.io import savemat
from scipy.spatial import cKDTree

# pseudonym -> (corrsource index, demo/source.mat index, hemisphere, n DK TPJ contacts, arrays)
SUBJECTS = {
    "P19": (18, 19, "lh", 43, "TG:112 AT:16 PST:16 AST:10"),
}

HERE   = os.path.dirname(os.path.abspath(__file__))
DATA   = os.path.abspath(os.path.join(HERE, "..", "data"))
if not os.environ.get("SMOOTH_RAWDATA"): sys.exit("Set SMOOTH_RAWDATA to the raw-data root (not distributed).")
ROOT   = os.path.join(os.environ["SMOOTH_RAWDATA"], "data")
if not os.environ.get("SMOOTH_RAWSUBJECT"): sys.exit("Set SMOOTH_RAWSUBJECT to P19's raw-data subject ID (not distributed).")
RAWSUB = os.environ["SMOOTH_RAWSUBJECT"]                      # raw-tree folder name; outputs use "P19"
FT_OFF = np.array([1.0, 1.0, 1.0])                            # nibabel -> FieldTrip frame

h = h5py.File(os.path.join(DATA, "corrsource.mat"))
refs = h["corrsource"][:, 0]

for sub in SUBJECTS:
    csub, ssub, hemi, ntpj, arrays = SUBJECTS[sub]
    recon = os.path.join(ROOT, RAWSUB, "recon", "freesurfer", "surf")

    pos, tri, curv, nv, nL = [], [], [], 0, 0
    for hh in ("lh", "rh"):
        v, f = nib.freesurfer.read_geometry(os.path.join(recon, f"{hh}.pial"))
        c = nib.freesurfer.read_morph_data(os.path.join(recon, f"{hh}.curv")).astype(float)
        pos.append(v + FT_OFF); tri.append(f + nv); curv.append(c); nv += len(v)
        if hh == "lh": nL = nv                                # lh occupies vertices 1..nL
    pos = np.vstack(pos); tri = np.vstack(tri); curv = np.concatenate(curv)

    g    = h[refs[csub - 1]]
    cp   = np.asarray(g["elec"]["nativechanpos"]).T
    d2s  = np.asarray(g["elec"]["dist2surf"]).ravel()
    d, _ = cKDTree(pos).query(cp)
    err  = np.abs(d - d2s).max()
    print(f"{sub}: {len(pos)} vertices, {len(cp)} electrodes, {hemi}, {ntpj} TPJ | {arrays}")
    print(f"  corrsource {csub} / source {ssub} | dist2surf err {err:.4f} mm | n_lh {nL}")
    if err > 0.05:
        raise SystemExit(f"{sub} does not match corrsource{{{csub}}} (err {err:.3f} mm)")

    dst = os.path.join(DATA, f"natcortex_{sub}.mat")
    savemat(dst, dict(natcortex=dict(pos=pos, tri=(tri + 1).astype(float), curv=curv, unit="mm"),
                      subject=sub, corrsource_idx=float(csub), source_idx=float(ssub),
                      n_lh=float(nL), hemi=hemi, arrays=arrays))
    print(f"  wrote {dst}")
