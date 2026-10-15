# Changelog

## 0.2.0

Results from this version differ from 0.1.0 for most users; see "Changed defaults" and "Fixed".

### Fixed
- **Surface geometry no longer depends on the MATLAB path.** SMOOTH now reads the fsaverage
  surfaces with its own FreeSurfer reader (`smooth/utils/smooth_read_surf.m`). In 0.1.0 it called
  whichever `read_surf.m` came first on the path. Some toolboxes ship a modified copy (iELVis's
  `read_surf.m` adds 1 to every vertex), which shifted the pial and registration-sphere surfaces
  by +1 mm on each axis relative to the electrodes. `SMOOTHstat`, `SMOOTHsub`, `SMOOTHdummy`,
  `SMOOTHdiag`, `load_fsaverage_mesh` and `fsavg_coord_to_dk` all use the new reader.
- **Python port:** removed the +[1, 1, 1] mm surface offset that mirrored the shifted MATLAB
  geometry, and uses the exact FWHM-to-sigma constant, 2·sqrt(2·ln 2).

### Changed defaults
- `cfg.rankrescale` now defaults to `'exact'` in `SMOOTHstat` and `SMOOTHsub` (was `'no'`). This is
  the setting used throughout the paper and the one under which surrogates match the empirical maps'
  distribution and spatial smoothness. Calls that set `cfg.rankrescale` explicitly are unaffected.
  `SMOOTHdummy`, the random-sign-flip comparison, still defaults to `'no'`.

### Added
- A coordinate-frame check: `SMOOTHstat` warns (`SMOOTHstat:shiftedFrame`), and the Python port
  issues a matching warning, when a participant's electrodes sit on the fsaverage pial surface
  shifted by +[1, 1, 1] mm, the signature of coordinates made with a modified `read_surf.m`.
- Input validation: an invalid `cfg.tail` now raises `SMOOTHstat:badTail` / `SMOOTHdummy:badTail`.
- `simulations/`: the simulation pipeline behind Figures 3 and 4 and Supplementary Figures 3–5.
- Paper scripts: `genFig5_nullhist` (redraws the Figure 5 null-distribution panels from the saved
  results), `genFigS7_tencomponents`, `export_fig5_sa` and `export_dist2surf.py`.

### Incompatible
- Caches built with 0.1.0 (`cfg.returncache`) are rejected with `SMOOTHstat:badCache`, because they
  hold the old geometry. Rebuild them with `cfg.returncache = 'yes'`.

### Paper
- All empirical results and figures were regenerated with the corrected geometry; the cached
  results in `paper/data/` are the regenerated ones.
- Figure 2c–f is now seeded (seed 42); a right-hemisphere indexing error in the exemplar panels was
  fixed, and the Figure 2e lag axis is now in mm.

## 0.1.0

Initial public release.
