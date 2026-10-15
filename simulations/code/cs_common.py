#!/usr/bin/env python3
"""Shared helpers for the simulation pipeline: load the primary-dataset (socialcog) electrodes and DVs,
estimate empirical spatial autocorrelation, build a Gaussian surface smoother, and render scalar fields
on the fsaverage5 inflated surface (headless matplotlib)."""
import os, numpy as np
from scipy.io import loadmat
from scipy.spatial.distance import pdist
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from sensitivity_core import read_fs_surface, geodesic_mm, load_elec   # numpy-only reader + registration-sphere geodesic

import config as C
FS5 = C.FS_SURF                                                    # fsaverage5 surf/ directory (FreeSurfer install)


def load_socialcog(hemi=None):
    """Return list of (coords Nx3, stat N) per subject; hemi in {None,'lh','rh'} filters by x."""
    out = []
    for s in load_elec():                                          # config.ELEC_FILE (<repo>/demo/source.mat)
        cp = np.asarray(s.chanpos, float); z = np.asarray(s.stat, float).ravel()
        ok = np.isfinite(z) & np.all(np.isfinite(cp), 1); cp, z = cp[ok], z[ok]
        if hemi == "lh": m = cp[:, 0] < 0;  cp, z = cp[m], z[m]
        elif hemi == "rh": m = cp[:, 0] >= 0; cp, z = cp[m], z[m]
        out.append((cp, z))
    return out


_GEO = None
def _geo_meshes():
    """fs5 pial+sphere per hemi (cached) for snapping electrodes -> registration sphere."""
    global _GEO
    if _GEO is None:
        _GEO = {}
        for h in ("lh", "rh"):
            v, _ = read_fs_surface(f"{FS5}/{h}.pial"); s, _ = read_fs_surface(f"{FS5}/{h}.sphere.reg")
            _GEO[h] = (v, s, cKDTree(v))
    return _GEO


def empirical_rho(subs, binw=3.0, maxd=60.0, geodesic=True):
    """Pooled within-subject correlation rho(d) (values demeaned per subject). geodesic=True measures
    electrode separation as registration-sphere geodesic distance (electrodes snapped to fs5, closer
    hemi), matching how SMOOTH smooths; cross-hemisphere pairs are excluded (no geodesic path)."""
    edges = np.arange(0, maxd + binw, binw); cen = edges[:-1] + binw/2
    csum = np.zeros(len(cen)); ccnt = np.zeros(len(cen)); var_acc = []
    G = _geo_meshes() if geodesic else None
    for cp, z in subs:
        if len(z) < 10: continue
        z = z - z.mean(); var_acc.append(np.mean(z**2))
        n = len(z); iu, ju = np.triu_indices(n, 1)
        if geodesic:
            cp = np.asarray(cp, float)
            hid = (G["rh"][2].query(cp)[0] < G["lh"][2].query(cp)[0]).astype(int)   # closer hemi: 0 lh, 1 rh
            esph = np.empty_like(cp)
            for hi, h in enumerate(("lh", "rh")):
                mm = hid == hi
                if mm.any(): esph[mm] = G[h][1][G[h][2].query(cp[mm])[1]]
            D = geodesic_mm(esph, esph); d = D[iu, ju].copy()
            d[hid[iu] != hid[ju]] = np.inf                                          # exclude cross-hemi pairs
        else:
            d = pdist(cp)
        zz = z[iu]*z[ju]
        idx = np.digitize(d, edges) - 1; m = (idx >= 0) & (idx < len(cen))
        np.add.at(csum, idx[m], zz[m]); np.add.at(ccnt, idx[m], 1.0)
    var0 = float(np.mean(var_acc))
    return cen, csum/np.maximum(ccnt, 1)/var0, ccnt, var0


def load_fs5(hemi="lh"):
    vpos, faces = read_fs_surface(f"{FS5}/{hemi}.pial")        # geometry (electrode space)
    ipos, _     = read_fs_surface(f"{FS5}/{hemi}.inflated")    # for display
    return vpos, faces, ipos


def build_gaussian_W(vpos, sigma):
    """Row-normalized sparse Gaussian smoothing operator over vertices (euclidean 3D)."""
    tree = cKDTree(vpos); nb = tree.query_ball_point(vpos, 3*sigma)
    rows, cols, vals = [], [], []
    for i, n in enumerate(nb):
        n = np.asarray(n); dd = np.linalg.norm(vpos[n] - vpos[i], axis=1)
        w = np.exp(-dd**2/(2*sigma**2))
        rows.append(np.full(len(n), i)); cols.append(n); vals.append(w)
    W = csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                   shape=(len(vpos), len(vpos)))
    return W.multiply(1.0/np.asarray(W.sum(1)))


def plot_surf_row(fields, titles, ipos, faces, outpath, clim=None, cmap="RdBu_r",
                  suptitle="", cbar_label=""):
    """Render a row of scalar fields on the inflated lh surface (lateral, cheap depth shading)."""
    tris = ipos[faces]
    tn = np.cross(tris[:, 1]-tris[:, 0], tris[:, 2]-tris[:, 0])
    tn /= (np.linalg.norm(tn, axis=1, keepdims=True) + 1e-9)
    sh = np.clip(0.55 + 0.45*(-tn[:, 0]), 0.35, 1.0)[:, None]
    cm = plt.get_cmap(cmap); k = len(fields)
    fig = plt.figure(figsize=(3.8*k, 4.4))
    for j, (f, t) in enumerate(zip(fields, titles)):
        ax = fig.add_subplot(1, k, j+1, projection="3d")
        lo, hi = clim if clim is not None else np.nanpercentile(f, [2, 98])
        fc = cm(plt.Normalize(lo, hi)(f[faces].mean(1))); fc = fc.copy(); fc[:, :3] *= sh
        ax.add_collection3d(Poly3DCollection(tris, facecolors=fc, linewidths=0, shade=False))
        mn, mx = ipos.min(0), ipos.max(0)
        ax.set_xlim(mn[0], mx[0]); ax.set_ylim(mn[1], mx[1]); ax.set_zlim(mn[2], mx[2])
        ax.set_box_aspect((mx-mn)); ax.view_init(elev=0, azim=180); ax.axis("off")
        ax.set_title(t, fontsize=10)
    if cbar_label:
        sm = plt.cm.ScalarMappable(cmap=cm, norm=plt.Normalize(*(clim if clim else (lo, hi))))
        fig.colorbar(sm, ax=fig.axes, shrink=0.5, pad=0.01, label=cbar_label)
    fig.suptitle(suptitle, fontsize=12); fig.savefig(outpath, dpi=140, bbox_inches="tight"); plt.close(fig)
