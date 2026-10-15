# SMOOTH — Python

A self-contained Python port of the SMOOTH group-level inference algorithm (`../smooth/SMOOTHstat.m`).

## Install

```bash
pip install -r requirements.txt          # numpy, scipy (+ matplotlib for the demo)
```

FreeSurfer surfaces are required (the `fsaverage` subject; `fsaverage5` for a fast run). Point the code at
them with either environment variable:

```bash
export SUBJECTS_DIR=/path/to/freesurfer/subjects      # or FREESURFER_HOME
```

## Quick start

```python
import numpy as np
from scipy.io import loadmat
import smoothstat as ss

# per-subject dependent variable (channel statistic) + electrode positions in fsaverage coordinates
src  = loadmat("../demo/source.mat", squeeze_me=True, struct_as_record=False)["source"]
subs = [(np.asarray(s.elec.chanpos, float), np.asarray(s.stat, float).ravel()) for s in np.atleast_1d(src)]

surf, nL, nR, adj = ss.load_mesh("fsaverage")           # "fsaverage5" is faster
res = ss.run_smoothstat(subs, surf, nL, nR, adj, numrand=1000, tail="positive", seed=42)

for mass, p, members in res["posclusters"]:
    print(f"cluster mass={mass:.1f}  p={p:.4f}  n_vertices={len(members)}")
```

`res` contains: `stat` (group t-map, length nL+nR), `coverage`, `mask`, `posclusters` / `negclusters`
(`[(mass, p, member_vertices), ...]`, strongest first), `sig` (boolean significance mask), the null cluster-mass
distributions `posdist`/`negdist`, and the per-subject operators (`subjects`) for the surrogate spatial-
autocorrelation diagnostic (`ss.moran_sa`).

## Notes

- **Coordinates**: electrode `chanpos` must be in fsaverage surface space (the frame the MATLAB pipeline
  uses). Surfaces are used exactly as stored, as in MATLAB. On `fsaverage`, a warning is issued if a
  subject's electrodes sit on the pial surface shifted by +[1,1,1] mm, the signature of coordinates made
  with a modified `read_surf.m` (see the main README); `fsaverage5` is too coarse for this check.
- **Agreement with MATLAB**: maps, t-statistics, coverage and observed cluster masses match `SMOOTHstat.m`
  to floating-point precision. Null distributions, and hence p-values, agree statistically rather than
  bit for bit: the random number generators differ, and so can the eigenvector basis that numpy and MATLAB
  pick when an electrode graph splits into disconnected parts.
- **Resolution**: `fsaverage` reproduces the manuscript; `fsaverage5` is a fast lower-resolution option
  (the simulations in Figures 3-4 use it). Run with `numrand=1000` for reported p-values.
- **Tail**: positive and negative clusters are each tested against their own null (one-tailed p-values).
  `sig` applies `alpha` for `'positive'` / `'negative'`, and `alpha/2` per tail for `'both'`. Note that
  `run_smoothstat` defaults to `tail="positive"` (the manuscript's directional test), whereas `SMOOTHstat.m`
  defaults to `'both'`; pass `tail` explicitly.
- See [`../demo/smoothdemo.py`](../demo/smoothdemo.py) for a runnable walkthrough that reproduces the
  manuscript's primary-dataset clusters end to end (`python ../demo/smoothdemo.py`, or `--quick` for
  a fast fsaverage5 sanity check); [`../demo/SMOOTHdemo.ipynb`](../demo/SMOOTHdemo.ipynb) for a
  worked example with visualization, and `../demo/SMOOTHdemo.m` for the MATLAB equivalent.
