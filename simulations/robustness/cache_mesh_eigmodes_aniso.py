#!/usr/bin/env python3
"""Anisotropic fs5 Laplacian eigenmodes for the anisotropic 1/f noise model. The metric is warped by the
area-preserving B=diag(rho^-1/2, rho^1/2, 1) (same as the SE aniso model), so we build a Gaussian RADIUS
graph in the warped mm space (sigma_g mm; a plain mesh-edge Laplacian is ~combinatorial and blind to the
warp). rho=1 is the isotropic baseline of this family. K modes -> <config.OUT_DIR>/mesh_eig_aniso_<h>_rho<rho>.npz.
  python cache_mesh_eigmodes_aniso.py [K=1500] [sigma_g_mm=8] [rho1 rho2 ...]
"""
import os as _os, sys as _sys, time
_HERE=_os.path.dirname(_os.path.abspath(__file__)); _sys.path.insert(0,_HERE)
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from scipy.spatial import cKDTree
import sensitivity_core as SC, config as C

argv=_sys.argv[1:]
K = int(argv[0]) if len(argv)>0 else 1500
SG = float(argv[1]) if len(argv)>1 else 8.0
RHOS = [float(x) for x in argv[2:]] if len(argv)>2 else [1.0,2.0,3.0,5.0,10.0]
MM = 100.0/SC.MM2SPHERE                                   # unit-sphere chordal -> mm
base=C.FS_SURF; OUT=C.OUT_DIR; _os.makedirs(OUT, exist_ok=True)

for h in ("lh","rh"):
    sph,_ = SC.read_fs_surface(f"{base}/{h}.sphere.reg"); nV=len(sph)
    U = sph/np.linalg.norm(sph,axis=1,keepdims=True)
    for rho in RHOS:
        dst=_os.path.join(OUT, f"mesh_eig_aniso_{h}_rho{rho:g}.npz")
        if _os.path.exists(dst) and int(np.load(dst)["Q"].shape[1])>=K: print(f"{h} rho{rho:g}: cached -> skip"); continue
        X = (U*np.array([rho**-0.5, rho**0.5, 1.0])[None,:]) * MM        # warped coords in mm
        tree=cKDTree(X); R=3*SG; pairs=tree.query_pairs(R, output_type="ndarray")
        i,j=pairs[:,0],pairs[:,1]; dd=np.linalg.norm(X[i]-X[j],axis=1); w=np.exp(-dd**2/(2*SG**2))
        W=sp.csr_matrix((np.r_[w,w],(np.r_[i,j],np.r_[j,i])),shape=(nV,nV))
        dg=np.asarray(W.sum(1)).ravel(); dg[dg==0]=1; Dm=sp.diags(1/np.sqrt(dg))
        Lsym=sp.eye(nV)-Dm@W@Dm; Lsym=((Lsym+Lsym.T)/2).tocsc()
        t=time.time(); lam,Q=eigsh(Lsym,k=K,sigma=-1e-3,which="LM")
        o=np.argsort(lam); lam=np.maximum(lam[o],0.0); Q=Q[:,o]
        np.savez(dst, lam=lam, Q=Q.astype(np.float32), rho=rho, sigma_g=SG)
        print(f"{h} rho{rho:g}: {nV} vtx, {K} modes, {len(i)} edges in {time.time()-t:.0f}s, lam<= {lam.max():.3f} -> {_os.path.basename(dst)}", flush=True)
print("DONE")
