"""Setup: build the fs5 graph-Laplacian eigenmode caches used by the 1/f GRF noise, into $PROJ/out/eigmodes.
  isotropic : mesh_eig<K>_<h>.npz            (Gaussian arc weights ~ combinatorial mesh Laplacian)
  warped    : mesh_eig<K>_aniso_<h>_rho<r>.npz  (B=diag(rho^-1/2,rho^1/2,1)-warped Gaussian radius graph)
Dense subset eigh (fast for a large fraction of the spectrum). K from SENS_MESH_EIG_TAG (default 5000).
  python build_eigmodes.py iso                  # both hemis, isotropic
  python build_eigmodes.py aniso 1 2 3 ... 10   # both hemis, given rho (default 1..10)
  python build_eigmodes.py all                  # iso + aniso rho=1..10
Idempotent (skips existing). Run once as a setup job; caches are dataset-independent (shared)."""
import os as _os, sys as _sys, time
_sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..")))
import numpy as np, scipy.sparse as sp, scipy.linalg as sla
from scipy.spatial import cKDTree
import sensitivity_core as SC, config as C

K = int(_os.environ.get("SENS_MESH_EIG_TAG", "5000")); SG = 8.0; MM = 100.0/SC.MM2SPHERE
base = C.FS_SURF; _os.makedirs(C.EIG_DIR, exist_ok=True)

def _laplacian_eig(nV, i, j, w):
    W = sp.csr_matrix((np.r_[w, w], (np.r_[i, j], np.r_[j, i])), shape=(nV, nV))
    dg = np.asarray(W.sum(1)).ravel(); dg[dg == 0] = 1; Dm = sp.diags(1/np.sqrt(dg))
    L = (sp.eye(nV) - Dm@W@Dm); L = ((L + L.T)/2).toarray()
    lam, Q = sla.eigh(L, subset_by_index=[0, K-1]); return np.maximum(lam, 0.0), Q.astype(np.float32)

def build_iso(h):
    dst = _os.path.join(C.EIG_DIR, f"mesh_eig{K}_{h}.npz")
    if _os.path.exists(dst): print(f"iso {h}: skip"); return
    sph, faces = SC.read_fs_surface(f"{base}/{h}.sphere.reg"); nV = len(sph)
    U = sph/np.linalg.norm(sph, axis=1, keepdims=True)
    E = np.unique(np.sort(np.vstack([faces[:, [0,1]], faces[:, [1,2]], faces[:, [2,0]]]), axis=1), axis=0)
    i, j = E[:, 0], E[:, 1]; th = np.arccos(np.clip((U[i]*U[j]).sum(1), -1, 1)); w = np.exp(-th**2/(2*50.0**2))
    t = time.time(); lam, Q = _laplacian_eig(nV, i, j, w)
    np.savez(dst, lam=lam, Q=Q); print(f"iso {h}: {K} modes {time.time()-t:.0f}s -> {_os.path.basename(dst)}", flush=True)

def build_aniso(h, rho):
    dst = _os.path.join(C.EIG_DIR, f"mesh_eig{K}_aniso_{h}_rho{rho:g}.npz")
    if _os.path.exists(dst): print(f"aniso {h} rho{rho:g}: skip"); return
    sph, _ = SC.read_fs_surface(f"{base}/{h}.sphere.reg"); nV = len(sph)
    U = sph/np.linalg.norm(sph, axis=1, keepdims=True)
    X = (U*np.array([rho**-0.5, rho**0.5, 1.0])[None, :])*MM
    tree = cKDTree(X)
    try:
        pr = tree.query_pairs(3*SG, output_type="ndarray")                 # scipy >= 1.9
    except TypeError:                                                      # older scipy: no output_type kwarg
        s = tree.query_pairs(3*SG); pr = np.array(sorted(s), dtype=np.int64) if s else np.empty((0, 2), np.int64)
    i, j = pr[:, 0], pr[:, 1]
    w = np.exp(-np.linalg.norm(X[i]-X[j], axis=1)**2/(2*SG**2))
    t = time.time(); lam, Q = _laplacian_eig(nV, i, j, w)
    np.savez(dst, lam=lam, Q=Q); print(f"aniso {h} rho{rho:g}: {K} modes {time.time()-t:.0f}s -> {_os.path.basename(dst)}", flush=True)

if __name__ == "__main__":
    mode = _sys.argv[1] if len(_sys.argv) > 1 else "all"
    rhos = [float(x) for x in _sys.argv[2:]] if len(_sys.argv) > 2 else list(range(1, 11))
    for h in ("lh", "rh"):
        if mode in ("iso", "all"): build_iso(h)
        if mode in ("aniso", "all"):
            for r in rhos: build_aniso(h, r)
    print("DONE")
