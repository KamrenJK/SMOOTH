#!/usr/bin/env python3
"""Cache the fsaverage5 registration-sphere graph-Laplacian eigenmodes (per hemi) for the isotropic and
nonstationary 1/f noise models (mesh-edge graph with Gaussian arc weights, sigma = 50 rad, i.e. ~the
combinatorial mesh Laplacian). Kee_1f = Q_e diag((lam+c)^-beta) Q_e^T uses these.
  python cache_mesh_eigmodes.py [K=1500]      # K modes -> <config.OUT_DIR>/mesh_eig_<hemi>.npz"""
import os as _os, sys as _sys, time
_HERE=_os.path.dirname(_os.path.abspath(__file__)); _sys.path.insert(0,_HERE)
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import sensitivity_core as SC, config as C

K = int(_sys.argv[1]) if len(_sys.argv)>1 else 1500
base = C.FS_SURF
OUT = _os.path.join(C.OUT_DIR); _os.makedirs(OUT, exist_ok=True)
for h in ("lh","rh"):
    dst = _os.path.join(OUT, f"mesh_eig_{h}.npz")
    if _os.path.exists(dst):
        d=np.load(dst);
        if int(d["Q"].shape[1])>=K: print(f"{h}: cached ({d['Q'].shape[1]} modes) -> skip"); continue
    sph,faces = SC.read_fs_surface(f"{base}/{h}.sphere.reg"); nV=len(sph)
    U = sph/np.linalg.norm(sph,axis=1,keepdims=True)
    E = np.sort(np.vstack([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1); E=np.unique(E,axis=0)
    i,j = E[:,0],E[:,1]
    th = np.arccos(np.clip((U[i]*U[j]).sum(1),-1,1)); w=np.exp(-th**2/(2*50.0**2))
    W = sp.csr_matrix((np.r_[w,w],(np.r_[i,j],np.r_[j,i])),shape=(nV,nV))
    dg = np.asarray(W.sum(1)).ravel(); dg[dg==0]=1; Dm=sp.diags(1/np.sqrt(dg))
    Lsym = sp.eye(nV)-Dm@W@Dm; Lsym=((Lsym+Lsym.T)/2).tocsc()
    t=time.time(); lam,Q = eigsh(Lsym,k=K,sigma=-1e-3,which="LM")
    o=np.argsort(lam); lam=np.maximum(lam[o],0.0); Q=Q[:,o]
    np.savez(dst, lam=lam, Q=Q.astype(np.float32))
    print(f"{h}: {nV} vtx, {K} modes in {time.time()-t:.0f}s, lam<= {lam.max():.3f} -> {dst}", flush=True)
print("DONE")
