"""Build paper/data/anisononstat_fig_data_v3.mat for paper/code/genFigS3.m (Supplementary Figure 3:
false-positive control under anisotropic and nonstationary noise).

FWER curves: read from the cached sweep files in paper/data/anisononstat/, produced by
simulations/robustness/run_sweeps.sh (fwer_sweep_aniso1f_nonug.npz: beta = 8 on the 1500-mode anisotropic
bases, rho = 1, 2, 3, 5, 10; fwer_sweep_nonstat1f_nonug.npz: beta = 10 on the 1500-mode isotropic mesh
basis, q = 1, 2, 4, 6, 8; 1000 null datasets each). The
stored per-simulation p-values are one-tailed positive (only positive clusters are formed), so FWER is
mean(p <= THRESH) at the one-tailed 0.05 threshold of the main FWER figure, with a 95% binomial CI.

Exemplar fields: one illustrative draw per rho / q at a fixed spectral slope per row, so the panels differ
only in anisotropy / nonstationarity. They are not the realizations behind the FWER numbers. beta is not
comparable across eigenbases -- the spectral weight is (lam + median(lam))^-beta, and median(lam) is 0.996
on the anisotropic basis, 0.247 on the 1500-mode mesh basis and 0.712 on the 5000-mode mesh basis -- so
each row uses its own basis and slope, checked against the sweeps' recorded `observed_SA` (electrode-level
Moran's I, 15 mm Gaussian weights):

  ANISOTROPY row: beta = 8 on the 1500-mode anisotropic (warped radius-graph) bases, the sweep's own
  basis and slope. It reproduces the
  simulated fields' SA at every rho, including its gradual decrease as rho grows:
        rho          1      2      3      5     10
        exemplar  0.483  0.457  0.446  0.420  0.369
        simulated 0.482  0.473  0.459  0.432  0.392
  (A per-rho beta would hold the exemplar SA constant and so track the simulated fields less closely.)

  NONSTATIONARITY row: beta = 29 on the 5000-mode isotropic basis (the main simulations' basis and
  slope), used only for these illustrative maps. The sweep itself ran at beta = 10 on the 1500-mode
  isotropic basis; beta = 29 on 5000 modes was chosen because it reproduces the sweep's recorded
  observed_SA to ~0.01 across q:
        q               1      2      4      8
        recorded     0.452  0.464  0.469  0.470
        5000/b=29    0.461  0.474  0.482  0.487

Nugget condition: a = 1.0 (no nugget, the adversarial pure-1/f case).
Inputs: paper/data/anisononstat/*.npz (sweeps + 1500-mode lh aniso bases), paper/data/methodscomparison_pvals.mat
(lh fsaverage5 mesh + curvature), and the 5000-mode isotropic basis from autocorr/build_eigmodes.py (config.EIG_DIR).

Usage: python export_anisononstat_v3.py      # -> <repo>/paper/data/anisononstat_fig_data_v3.mat (override: SENS_AN_OUT)
"""
import os, sys
import numpy as np
import scipy.io as sio
from scipy.sparse import csr_matrix
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")))
import config as C

THRESH = 0.05          # one-tailed positive; matches the main FWER figure
NUGGET = "a=1.0 (no nugget)"

VF   = C.PAPER_DATA
SRC  = os.path.join(VF, "anisononstat")
EIG  = C.EIG_DIR
DST  = os.environ.get("SENS_AN_OUT") or os.path.join(VF, "anisononstat_fig_data_v3.mat")

# Exemplar slopes, per row, each on its own basis (see the module docstring).
BETA_ANISO,   NMODE_ANISO   = 8.0,  1500     # aniso radius-graph bases, median(lam) ~ 0.996
BETA_NONSTAT, NMODE_NONSTAT = 29.0, 5000     # isotropic mesh basis, median(lam) = 0.712 -- the main
                                             # simulations' condition; matches the sweep's recorded SA

# ---------------------------------------------------------------- surface (lh)
mc   = sio.loadmat(os.path.join(VF, "methodscomparison_pvals.mat"), squeeze_me=True)
nL   = int(mc["nL"])
pos  = np.asarray(mc["fs5_pos"], float)[:nL]
curv = np.asarray(mc["fs5_curv"], float)[:nL]
tri  = np.asarray(mc["fs5_tri"], float).astype(int) - 1
tri  = tri[(tri < nL).all(1)]

def adjacency(tri, n):
    i = tri[:, [0, 1, 2]].ravel(); j = tri[:, [1, 2, 0]].ravel()
    A = csr_matrix((np.ones(len(i) * 2), (np.r_[i, j], np.r_[j, i])), shape=(n, n))
    return (A > 0).astype(float)

ADJ = adjacency(tri, nL)

def moran(x, A=ADJ):
    """Graph Moran's I on the mesh adjacency. Printed for reference only -- it is NOT the metric the
    sweeps recorded (that is electrode-level with 15 mm Gaussian weights) and the two are not on the
    same scale."""
    x = x - x.mean(); w = A.sum()
    return (len(x) / w) * float(x @ (A @ x)) / float(x @ x)

# ---------------------------------------------------------------- field synthesis
def onef(lam, beta):
    P = 1.0 / np.power(lam + np.median(lam), beta)
    P[0] = 0.0
    return P

def draw(Q, lam, beta, rng):
    a = np.sqrt(onef(lam, beta)) * rng.standard_normal(len(lam))
    y = Q @ a
    return (y - y.mean()) / (y.std() + 1e-12)

# ---------------------------------------------------------------- FWER sweeps
def sweep(fname):
    d = np.load(os.path.join(SRC, fname), allow_pickle=True)
    x = np.asarray(d["vals"], float); n = int(d["nsim"])
    out = {}
    for arm, key in [("smooth", "smooth_p"), ("naive", "naive_p")]:
        p = (np.asarray(d[key]) <= THRESH).mean(1)
        se = np.sqrt(p * (1 - p) / n)
        out[arm] = p; out[arm + "_lo"] = np.clip(p - 1.96 * se, 0, 1); out[arm + "_hi"] = np.clip(p + 1.96 * se, 0, 1)
    out["vals"] = x; out["n"] = n; out["SA"] = np.asarray(d["observed_SA"], float).mean(1)
    print(f"  {fname}: vals={x.astype(int).tolist()}  SMOOTH={np.round(out['smooth'],3).tolist()}  "
          f"naive={np.round(out['naive'],3).tolist()}")
    return out

print(f"FWER at one-tailed p <= {THRESH}  [{NUGGET}]")
A = sweep("fwer_sweep_aniso1f_nonug.npz")
N = sweep("fwer_sweep_nonstat1f_nonug.npz")

# ---------------------------------------------------------------- exemplars
RHOS = [1, 2, 3, 5, 10]
QS   = [1, 2, 4, 6, 8]

print(f"\nexemplar fields")
print(f"  aniso   beta = {BETA_ANISO:.1f} on the {NMODE_ANISO}-mode aniso bases")
aniso_fields, aniso_SA_ex = [], []
for r in RHOS:
    e = np.load(os.path.join(SRC, f"mesh_eig_aniso_lh_rho{r}.npz"))
    f = draw(e["Q"].astype(float), e["lam"].astype(float), BETA_ANISO, np.random.default_rng(11))
    aniso_fields.append(f); aniso_SA_ex.append(moran(f))
print(f"    mesh-adjacency I: {np.round(aniso_SA_ex,3).tolist()}")

print(f"  nonstat beta = {BETA_NONSTAT:.1f} on the {NMODE_NONSTAT}-mode isotropic basis")
em   = np.load(os.path.join(EIG, f"mesh_eig{NMODE_NONSTAT}_lh.npz"))
lamM = em["lam"].astype(float); QM = em["Q"].astype(float)
assert QM.shape[0] == nL, f"basis has {QM.shape[0]} vertices, surface has {nL}"

ap = pos[:, 1]
lo, hi = np.percentile(ap, [5, 95])
w = np.clip((hi - ap) / (hi - lo), 0, 1)          # 1 posterior -> 0 anterior

nonstat_fields, nonstat_SA_ex = [], []
for q in QS:
    sq = np.sqrt(q); rng = np.random.default_rng(23)
    rough   = draw(QM, lamM, BETA_NONSTAT / sq, rng)          # posterior, rougher
    smoothf = draw(QM, lamM, BETA_NONSTAT * sq, rng)          # anterior, smoother
    f = np.sqrt(w) * rough + np.sqrt(1 - w) * smoothf
    f = (f - f.mean()) / (f.std() + 1e-12)
    nonstat_fields.append(f); nonstat_SA_ex.append(moran(f))
print(f"    mesh-adjacency I: {np.round(nonstat_SA_ex,3).tolist()}")

# ---------------------------------------------------------------- write
sio.savemat(DST, dict(
    thresh=THRESH, thresh_label="one-tailed p <= 0.05", nugget=NUGGET, nsim=A["n"],
    lh_pos=pos, lh_tri=(tri + 1).astype(float), lh_curv=curv,
    aniso_rho=A["vals"], aniso_smooth=A["smooth"], aniso_smooth_lo=A["smooth_lo"], aniso_smooth_hi=A["smooth_hi"],
    aniso_naive=A["naive"], aniso_naive_lo=A["naive_lo"], aniso_naive_hi=A["naive_hi"], aniso_SA=A["SA"],
    nonstat_q=N["vals"], nonstat_smooth=N["smooth"], nonstat_smooth_lo=N["smooth_lo"], nonstat_smooth_hi=N["smooth_hi"],
    nonstat_naive=N["naive"], nonstat_naive_lo=N["naive_lo"], nonstat_naive_hi=N["naive_hi"], nonstat_SA=N["SA"],
    aniso_fields=np.column_stack(aniso_fields), aniso_field_rho=np.array(RHOS, float),
    aniso_beta=BETA_ANISO, aniso_nmodes=float(NMODE_ANISO),
    nonstat_fields=np.column_stack(nonstat_fields), nonstat_field_q=np.array(QS, float),
    nonstat_beta=BETA_NONSTAT, nonstat_nmodes=float(NMODE_NONSTAT),
    ap_weight=w,
))
print(f"\nwrote {DST}")
