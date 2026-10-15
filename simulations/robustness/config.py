#!/usr/bin/env python3
"""Configuration for the anisotropy / nonstationarity FWER sweeps (Supplementary Figure 3). Paths are
env-overridable:
    SENS_DATASET   - output label recorded in the sweep files (default socialcog, the primary dataset)
    SENS_FS_SURF   - fsaverage5 surf/ directory holding {lh,rh}.{pial,sphere.reg} (default:
                     $SUBJECTS_DIR/fsaverage5/surf if present, else $FREESURFER_HOME/subjects/fsaverage5/surf)
    SENS_ELEC      - .mat with per-subject electrode coords (default <repo>/demo/source.mat)
    SENS_ATLAS_DIR - directory holding <atlas>_fs5.npz parcel lookups (default <repo>/paper/data; used
                     only by the per-hemisphere ROI arm of HemiSim, not by the FWER sweeps)
    SENS_OUT       - output directory for eigenmode caches and sweep results (default robustness/out)
    SENS_MESH, SENS_TAIL, SENS_CLUSTERALPHA, SENS_A0, SENS_NITER, SENS_NPERM - optional parameter overrides"""
import os
DIR  = os.path.dirname(os.path.abspath(__file__))                       # simulations/robustness/
REPO = os.environ.get("SMOOTH_REPO") or os.path.abspath(os.path.join(DIR, "..", ".."))

# ---- dataset + paths (override via env) ----
DATASET   = os.environ.get("SENS_DATASET", "socialcog")   # output label
MESH      = os.environ.get("SENS_MESH", "fsaverage5")
def _fs_surf():
    if os.environ.get("SENS_FS_SURF"): return os.environ["SENS_FS_SURF"]
    sd = os.environ.get("SUBJECTS_DIR")
    if sd and os.path.isdir(os.path.join(sd, MESH, "surf")): return os.path.join(sd, MESH, "surf")
    return os.path.join(os.environ.get("FREESURFER_HOME", "/usr/local/freesurfer"), "subjects", MESH, "surf")
FS_SURF   = _fs_surf()                                    # <FS_SURF>/{lh,rh}.{pial,sphere.reg}
ATLAS_DIR = os.environ.get("SENS_ATLAS_DIR", os.path.join(REPO, "paper", "data"))
ELEC_FILE = os.environ.get("SENS_ELEC", os.path.join(REPO, "demo", "source.mat"))
OUT_DIR   = os.environ.get("SENS_OUT",  os.path.join(DIR, "out"))

# ---- mesh + SMOOTH parameters (validated against MATLAB SMOOTHstat) ----
ATLAS       = os.environ.get("SENS_ATLAS", "aparc")   # ROI-arm parcellation: aparc (Desikan-Killiany) | schaefer300
KERNELWIDTH = 12.0     # mm FWHM, electrode->surface Gaussian smoothing
GRAPHSIGMA  = 8.0      # mm, electrode-graph Laplacian sigma for the eigenmode null
MINNBSUB    = 3        # vertices with < this many subjects are masked
CLUSTERALPHA = float(os.environ.get("SENS_CLUSTERALPHA", "0.01"))   # cluster-forming z = norminv(1-clusteralpha); run_sweeps.sh sets 0.05 (the SMOOTHstat default)
ALPHA       = 0.05     # cluster-level significance
TAIL        = os.environ.get("SENS_TAIL", "positive")   # 'both' tests clusters at ALPHA/2; run_sweeps.sh sets 'both'
FT_OFFSET   = (1.0, 1.0, 1.0)   # offset (mm) added to the fsaverage5 pial and sphere.reg vertex coordinates

# ---- squared-exponential noise parameters (SE noise models; the 1/f models ignore them, SA_L is
#      recorded in the sweep file as baseL) ----
SA_A = 0.521           # correlated variance fraction (1-SA_A = nugget); geodesic fit
SA_L = 13.1            # mm, Gaussian correlation length (geodesic / registration sphere)

# ---- ground-truth event + effect size (sensitivity mode only) ----
BUMP_FWHM = 20.0       # mm, simulated functional-region size
A0        = float(os.environ.get("SENS_A0", "4.0"))   # bump amplitude / noise SD (SNR); operating point
R_HIT     = 12.0       # mm, a significant cluster within this of the seed counts as a true positive

# ---- simulation size (sensitivity mode only; the FWER sweep takes --nsim / --nperm) ----
NITER     = int(os.environ.get("SENS_NITER", "100"))  # iterations per seed vertex
NPERM     = int(os.environ.get("SENS_NPERM", "100"))  # SMOOTH permutations per iteration (exact null)
BASE_SEED = 42                                         # RNG base; per-seed streams are deterministic
