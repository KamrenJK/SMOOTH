# SMOOTH

**S**urrogate **M**ap **O**ptimized **O**bservation and **T**opographic **H**ypothesis testing —
a MATLAB toolbox for **group-level statistical inference on intracranial EEG (ECoG/sEEG)**.

Electrode coverage in iEEG is sparse, anatomically uneven, and differs between every
participant. That breaks the assumptions behind the permutation schemes normally used for
group inference. SMOOTH addresses this by building null distributions from **spatially
constrained surrogate maps**: each subject's cortical map is decomposed into graph Laplacian
eigenmodes, the modal coefficients are randomly sign-flipped, and the map is reconstructed.
This preserves the intrinsic spatial smoothness of the data while destroying its topographic
alignment across subjects — so cluster-based tests ask whether observed spatial convergence
exceeds chance, without relying on a predefined atlas.

---

## Requirements

| | |
|---|---|
| MATLAB | R2022a or newer (developed and tested on R2024b), with the Statistics and Machine Learning Toolbox |
| [FieldTrip](https://www.fieldtriptoolbox.org) | any recent release |
| [FreeSurfer](https://surfer.nmr.mgh.harvard.edu) | 7.x, with the `fsaverage` subject installed |

Figures 3-4 and most of the supplementary-figure simulations additionally use `fsaverage5`.
FreeSurfer is required for its surface files only — no FreeSurfer processing is run.

Two third-party MATLAB packages ship inside `external/` and need no separate install:
[slanCM](https://github.com/slandarer) (colormaps) and
[`fdr_bh`](https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh)
(Groppe's Benjamini-Hochberg FDR, used by Supplementary Figure 1's channel-level correction).
Both are BSD-licensed; see [NOTICE](NOTICE).

## Installation

```matlab
% 1. clone the repository
%    git clone https://github.com/KamrenJK/SMOOTH.git

% 2. open smooth_setup.m and edit two lines:
%       paths.fieldtrip  = '/path/to/fieldtrip';
%       paths.freesurfer = '/path/to/freesurfer';   % e.g. /Applications/freesurfer/7.3.2
%    (or set the FIELDTRIP_HOME / FREESURFER_HOME environment variables instead)

% 3. from the repository root:
paths = smooth_setup();
```

`smooth_setup` derives the rest of the layout from its own location, validates that
FieldTrip and `fsaverage` are actually present, and puts everything on the MATLAB path.
Nothing else needs editing — the demo, figure scripts, and tests all call it.

## Quick start

```matlab
paths = smooth_setup();
load(fullfile(paths.demo,'source.mat'));   % 38 subjects, 3950 channels

cfg                  = [];
cfg.fshome           = paths.freesurfer;   % required
cfg.numrandomization = 1000;
cfg.randomseed       = 42;                 % omit and results are NOT reproducible
stat = SMOOTHstat(cfg, source{:});

sig = find([stat.posclusters.prob] <= 0.05);
```

Or just run the annotated walkthrough:

```matlab
run(fullfile(paths.demo,'SMOOTHdemo.m'))
```

**Runtime**: ~3 minutes for 1000 permutations on the bundled 38-subject dataset
(Apple M-series, MATLAB R2024b). Enabling surrogate retention for diagnostics adds
memory pressure — see `cfg.keepsubsurrogates` below.

---

## Preparing your own data

`SMOOTHstat` takes one struct per subject:

```matlab
source{i}.stat            % nchan x 1, one statistic per channel (e.g. a t or z value)
source{i}.elec.chanpos    % nchan x 3, electrode coordinates in fsaverage surface space
source{i}.label           % nchan x 1 cell of channel labels (optional)
```

`source{i}.stat` is whatever per-channel effect you want to test — SMOOTH is agnostic to how
it was computed. **Coordinates must already be in fsaverage space**; hemisphere is assigned by
the sign of x. Inputs are validated up front, so mismatched sizes, NaNs, or coordinates that
do not look like fsaverage produce an explicit error or warning rather than silently wrong maps.

## Configuration

Only `cfg.fshome` is required.

### Core
| option | default | meaning |
|---|---|---|
| `cfg.fshome` | — | FreeSurfer home (must contain `subjects/fsaverage`) |
| `cfg.parameter` | `'stat'` | field of the source struct holding the statistic |
| `cfg.numrandomization` | `1000` | number of surrogates |
| `cfg.randomseed` | `[]` | integer for reproducible results. If unset, a warning is raised and the RNG state is stored in `stat.cfg.randomseed` |

### Spatial model
| option | default | meaning |
|---|---|---|
| `cfg.smooth` | `'sphere'` | `'sphere'` measures distance on `sphere.reg` (fast; correctly separates opposite banks of a sulcus). `'pial'` uses true geodesic distance in mm (5–10× slower) |
| `cfg.kernel` | `'gaussian'` | `'gaussian'` or `'boxcar'` |
| `cfg.kernelwidth` | `12` | smoothing FWHM in mm, matched to the empirical HFB spatial-autocorrelation scale |
| `cfg.graphsigma` | `8` | Gaussian σ (mm) for the electrode graph. Larger → smoother surrogates |
| `cfg.randmethod` | `'signflip'` | `'signflip'`, or `'bandrotate'` for band-wise orthogonal rotation |
| `cfg.rankrescale` | `'no'` | match surrogate marginals to the empirical map: `'exact'`, `'quantile'`, `'zscore'` |

### Inference
| option | default | meaning |
|---|---|---|
| `cfg.tail` | `'both'` | which tail `stat.mask` reflects. One-tailed uses `cfg.alpha`; `'both'` uses `cfg.alpha/2` per tail (Bonferroni across tails) |
| `cfg.clusteralpha` | `0.05` | cluster-forming threshold per tail |
| `cfg.minnbsub` | `3` | minimum subjects contributing to a vertex |
| `cfg.equalize_n` | `'no'` | downweight high-coverage vertices in the t-statistic |
| `cfg.normalize` | `'no'` | per-subject z-scoring before the group test |

> **On `cfg.tail`.** The default is two-tailed because iEEG contrasts often show widespread
> suppression as well as activation — in the bundled dataset 69% of channel-level statistics
> are negative. A one-tailed default would silently hide half the effects. Set
> `cfg.tail='positive'` when your hypothesis is directional.

### Output control
| option | default | meaning |
|---|---|---|
| `cfg.keepmaps` | `'no'` | return subject-level cortical maps |
| `cfg.keepsubsurrogates` | `'no'` | return subject-level surrogates — **memory expensive** |
| `cfg.keepgroupsurrogates` | `'no'` | return group-level surrogate t-maps |
| `cfg.subsamplesurr` | `'yes'` | keep only `nsurrsamples` surrogates instead of all |
| `cfg.nsurrsamples` | `100` | how many to keep |
| `cfg.feedback` | `'text'` | `'no'` to silence progress |

Subject-level surrogates are stored as `nsubjects × nvertices × nsurrsamples`. For 38 subjects
that is ~9 GB at `nsurrsamples=100`. SMOOTHstat reports the allocation and warns above 4 GB.

## Outputs

```
stat.stat            group-level t-map over fsaverage vertices
stat.mask            significance mask, FWE-controlled according to cfg.tail
stat.prob            vertex-wise cluster p-value, ONE-TAILED and not corrected across
                     tails -- prefer stat.mask; see cfg.tail before thresholding this
stat.posclusters     positive clusters, each with .clusterstat and .prob
stat.negclusters     negative clusters, likewise
stat.coverage        subjects contributing to each vertex
stat.cfg             configuration actually used, including the resolved random seed
```

## Other functions

| | |
|---|---|
| `SMOOTHstat` | group-level inference — the main entry point |
| `SMOOTHsub` | single-subject spectral decomposition and surrogates, for illustration/QC. Performs **no** inference |
| `SMOOTHdiag` | diagnostics: checks surrogate spatial autocorrelation matches the empirical maps |
| `SMOOTHdummy` | random sign flipping (RSF), for comparison against SMOOTH |

Surrogates are only valid if their spatial autocorrelation matches the empirical data. Use
`SMOOTHdiag` to check; if smoothness is mismatched, adjust `cfg.kernelwidth` or `cfg.graphsigma`.

## Python

`python/smoothstat.py` is a self-contained Python port of `SMOOTHstat` (only numpy and scipy required
— no FieldTrip, no MATLAB). Useful if your pipeline is already Python-based, or to independently
reproduce a result without a MATLAB license.

```python
import numpy as np
from scipy.io import loadmat
import smoothstat as ss

src  = loadmat("demo/source.mat", squeeze_me=True, struct_as_record=False)["source"]
subs = [(np.asarray(s.elec.chanpos, float), np.asarray(s.stat, float).ravel()) for s in np.atleast_1d(src)]

surf, nL, nR, adj = ss.load_mesh("fsaverage")           # needs $SUBJECTS_DIR or $FREESURFER_HOME
res = ss.run_smoothstat(subs, surf, nL, nR, adj, numrand=1000, tail="positive", seed=42)

for mass, p, members in res["posclusters"]:
    print(f"cluster mass={mass:.1f}  p={p:.4f}  n_vertices={len(members)}")
```

For a runnable introduction instead of a snippet, see [`demo/smoothdemo.py`](demo/smoothdemo.py) —
reproduces the primary dataset's reported clusters (pMTG p=0.001, rTPJ p=0.009) end to end
(`python demo/smoothdemo.py`; add `--quick` for a ~15 s fsaverage5 sanity check instead of the full
~5 min fsaverage run).

---

## Reproducing the manuscript

Every `genFig*` script loads a cached `.mat`/`.npz` from `paper/data/` and renders directly from
it — none require re-running an analysis to reproduce a figure. Some of those cached files were
themselves produced by a one-off `export_*` script that depends on private, unpublished
infrastructure (a separate simulation/analysis codebase, raw clinical recordings, or a compute
cluster); those exporters are documented in each `genFig*` script's header for provenance, but are
not required, and will not run standalone outside that infrastructure. Run with a display unless a
script's header says otherwise.

```matlab
paths = smooth_setup();
cd(fullfile(paths.root,'paper','code'));

% Main figures
genFig1;          % Fig 1  graphical overview (real rTPJ workflow example)
genFig1_inputs;   % Fig 1  inset: toy dependent-variable icons
genFig2;            % Fig 2  spatial autocorrelation of iEEG signals -- panel b
genFig2_seedcorr;   % Fig 2  panel a (seed-correlation berry map)
genFig2_surrogates; % Fig 2  panels c-f (surrogate generation / semivariogram / Moran's I)
genFig3;            % Fig 3  false-positive control under spatial autocorrelation
genFig4;          % Fig 4  detection sensitivity vs. atlas-based analyses
genFig5;          % Fig 5  recovers effects across three iEEG datasets (ALSO writes SF7, see below)

% Supplementary figures
genFigS1;  % SF1  communication task and channel-level effect estimation
genFigS2;  % SF2  participant-specific spatial basis (spectral decomposition)
genFigS3;  % SF3  false-positive control under anisotropic/nonstationary noise
genFigS4;  % SF4  ground-truth sensitivity-simulation schematic
genFigS5;  % SF5  detection sensitivity, FDR-corrected atlas analyses
genFigS6;  % SF6  comparison of analysis approaches across datasets
%          % SF7  RSF reconstruction -- produced by genFig5 (above), not a separate script:
%          % it reuses the same SMOOTH/RSF run Figure 5 computes, so splitting it out would mean
%          % either duplicating that run or shipping the multi-GB surrogate ensembles it depends on.
```

Set `savefigs = true` near the top of a script (where present) to write PDFs to `paper/figs/<n>/`
or `paper/figs/S<n>/`.

Figure 3's false-positive-control simulations and most supplementary-figure simulations were run on
a computing cluster with a separate simulation pipeline; `paper/data/` ships their cached results
so the figures regenerate without it.

## Data

`demo/source.mat` is the full dataset from the manuscript: 38 participants, 3,950 bipolar
channels, one HFB statistic per channel, with electrode coordinates in fsaverage space.
Recordings were collected under IRB-approved protocols with written informed consent, and
the distributed file is de-identified.

`paper/data/corrsource.mat` carries the same cohort's native-space electrode positions and
per-channel LFP/low-frequency/HFB correlation matrices (Figure 2's spatial-autocorrelation
analysis). `paper/data/natcortex_P19.mat` is one participant's native pial surface reconstruction,
used only for Figure 2's seed-correlation berry map (`genFig2_seedcorr.m`) — a cortical surface
mesh carries no facial or scalp geometry. Both are de-identified consistent with the IRB-approved
protocols above.

## Citation

If you use SMOOTH, please cite the manuscript and the software (see `CITATION.cff`).

## Third-party code

`external/slanCM` — colormaps by Zhaoxu Liu, BSD-3-Clause (see
`external/slanCM/license.txt`), incorporating colormaps derived from
[matplotlib](https://matplotlib.org) and [scicomap](https://github.com/ThomasBury/scicomap)
(MIT). Full notices in [NOTICE](NOTICE).

## License

See [LICENSE](LICENSE).

## Contact

Questions, feedback, and bug reports: open an issue, or email
kamren[dot]khan[at]pennmedicine[dot]upenn[dot]edu

Mutual Understanding Lab — https://www.mutualunderstanding.nl/
