"""Simulation-pipeline configuration. Paths derive from $PROJ (this simulations/ folder) and the SMOOTH
repository root that contains it; fsaverage5 surfaces come from a FreeSurfer installation. Analysis
parameters are locked to the manuscript's matched specificity+sensitivity design (1/f GRF, beta=29,
ca=0.05, whole-brain); each is env-overridable for sweeps.

  PROJ            simulations/ folder (default: parent of code/); outputs go to $PROJ/out
  SMOOTH_REPO     repository root (default: parent of $PROJ); electrode data + atlas lookups live here
  SENS_FS_SURF    fsaverage5 surf/ directory (default: $SUBJECTS_DIR/fsaverage5/surf if present,
                  else $FREESURFER_HOME/subjects/fsaverage5/surf)
  SENS_ELEC       electrode file (default: <repo>/demo/source.mat)
  SENS_DATASET    output label / subdirectory (default socialcog, the primary dataset)
  SENS_A          correlated fraction a: 0.53 (with-nugget, primary) | 1.0 (no-nugget)
  SENS_BETA       1/f slope (default 29, the rho(d) fit at 5000 modes)
  SENS_MESH_EIG_TAG  eigenbasis mode count tag (default 5000 -> eigmodes/mesh_eig5000_<h>.npz)
  SENS_CA         cluster-forming alpha (LOCKED 0.05)
  SENS_ATLAS      aparc (Desikan-Killiany) | schaefer300   (ROI comparator pass)
  SENS_NSIM/NITER/NPERM   sim sizes (default 1000/100/1000)
  SLURM_CPUS_PER_TASK / NCORES   worker count (flexible; defaults to all cores)
"""
import os
PROJ = os.environ.get("PROJ") or os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
REPO = os.environ.get("SMOOTH_REPO") or os.path.abspath(os.path.join(PROJ, ".."))
PAPER_DATA = os.path.join(REPO, "paper", "data")              # cached figure data (and the export destination)

# ---- noise: 1/f GRF (locked to the rho(d) fit) ----
os.environ.setdefault("SENS_MESH_EIG_TAG", "5000")            # 5000-mode eigenbasis; read by sensitivity_core.load_mesh_eig
BETA = float(os.environ.get("SENS_BETA", "29.0"))             # empirical 1/f slope (5000-mode rho(d) fit)
SA_A = float(os.environ.get("SENS_A", "0.53"))                # correlated fraction (with-nugget primary); set 1.0 for no-nugget
SA_L = 13.1                                                   # SE fallback only (unused for 1/f)

# ---- mesh + surfaces ----
MESH = os.environ.get("SENS_MESH", "fsaverage5")
def _fs_surf():
    if os.environ.get("SENS_FS_SURF"): return os.environ["SENS_FS_SURF"]
    sd = os.environ.get("SUBJECTS_DIR")
    if sd and os.path.isdir(os.path.join(sd, MESH, "surf")): return os.path.join(sd, MESH, "surf")
    return os.path.join(os.environ.get("FREESURFER_HOME", "/usr/local/freesurfer"), "subjects", MESH, "surf")
FS_SURF   = _fs_surf()                                        # <FS_SURF>/{lh,rh}.{pial,sphere.reg,inflated,curv}
ATLAS_DIR = os.environ.get("SENS_ATLAS_DIR", PAPER_DATA)       # <ATLAS_DIR>/<atlas>_fs5.npz (per-vertex parcel labels)

# ---- dataset + paths ----
DATASET   = os.environ.get("SENS_DATASET", "socialcog")       # output label (primary dataset)
ELEC_FILE = os.environ.get("SENS_ELEC", os.path.join(REPO, "demo", "source.mat"))
OUT_DIR   = os.environ.get("SENS_OUT",  os.path.join(PROJ, "out", DATASET))         # per-dataset results
EIG_DIR   = os.environ.get("SENS_EIG",  os.path.join(PROJ, "out", "eigmodes"))      # SHARED mesh-Laplacian eigenmode caches

# ---- SMOOTH parameters (locked; validated vs MATLAB SMOOTHstat) ----
ATLAS       = os.environ.get("SENS_ATLAS", "aparc")           # ROI comparator: aparc | schaefer300
KERNELWIDTH = 12.0
GRAPHSIGMA  = 8.0
MINNBSUB    = 3
CLUSTERALPHA = float(os.environ.get("SENS_CA", "0.05"))       # LOCKED 0.05 (specificity+sensitivity matched)
ALPHA       = 0.05
TAIL        = os.environ.get("SENS_TAIL", "positive")
FT_OFFSET   = (1.0, 1.0, 1.0)   # offset (mm) added to the fsaverage5 pial and sphere.reg vertex coordinates

# ---- ground-truth event (sensitivity) ----
BUMP_FWHM = 20.0
A0        = float(os.environ.get("SENS_A0", "4.0"))           # peak SNR
R_HIT     = 12.0

# ---- simulation sizes ----
NSIM  = int(os.environ.get("SENS_NSIM",  "1000"))             # FWER null datasets / regime
NITER = int(os.environ.get("SENS_NITER", "100"))             # sensitivity iterations / seed
NPERM = int(os.environ.get("SENS_NPERM", "1000"))            # SMOOTH permutations / test
BASE_SEED = 42

# ---- resources (flexible: scales to whatever the node grants) ----
NCORES = int(os.environ.get("SLURM_CPUS_PER_TASK", os.environ.get("NCORES", str(os.cpu_count() or 8))))
