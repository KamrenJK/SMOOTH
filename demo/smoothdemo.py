#!/usr/bin/env python3
"""
SMOOTH demo (Python) -- reproduces the primary ("social cognition") result from the manuscript
================================================================================================
A minimal end-to-end walkthrough of the Python port, for anyone whose pipeline is already
Python-based:
  1) load the bundled example dataset (demo/source.mat: 38 subjects, 3,950 bipolar channels)
  2) run SMOOTHstat (group-level cluster-based permutation inference)
  3) check the result against the manuscript's reported clusters
  4) a null-distribution plot for the top cluster, and the surrogate spatial-autocorrelation
     robustness check (same two diagnostics demo/SMOOTHdemo.m ends with)

This uses cfg.tail="positive" to match the manuscript's directional hypothesis for this dataset
(genFig5.m, the script that actually generates the published figure) -- NOT the two-tailed default
in demo/SMOOTHdemo.m, which is written for general/exploratory use. With a directional hypothesis,
one-tailed is the correct and more powerful choice; with no prior hypothesis about direction, prefer
SMOOTHdemo.m's two-tailed default.

Expect two significant clusters, matching Figure 5a/d of the manuscript:
    pMTG   p = 0.001
    rTPJ   p = 0.009  (permutation-based, so this can read ~0.008-0.011 depending on the random draw)

Requirements: numpy, scipy (see ../python/requirements.txt) and a FreeSurfer installation (for the
fsaverage surfaces -- point SUBJECTS_DIR or FREESURFER_HOME at it; no FreeSurfer processing is run,
only its shipped surface files are read).

Runtime: a few minutes for the full 1,000-permutation run below (full fsaverage, 38 subjects --
about 5 minutes on a recent laptop). For a quick sanity check instead of a publication-matching run,
pass --quick (fsaverage5, 100 permutations, ~15 seconds) -- cluster masses and p-values will not
match the paper at that resolution, only the qualitative pattern (two significant clusters).
"""
import argparse, os, sys, time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
import numpy as np
from scipy.io import loadmat
import smoothstat as ss


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true",
                    help="fsaverage5 + 100 permutations (~15s): a fast sanity check, not a reproduction")
    args = ap.parse_args()

    mesh = "fsaverage5" if args.quick else "fsaverage"
    numrand = 100 if args.quick else 1000

    demo_source = os.path.join(os.path.dirname(__file__), "source.mat")
    src = loadmat(demo_source, squeeze_me=True, struct_as_record=False)["source"]
    subs = [(np.asarray(s.elec.chanpos, float), np.asarray(s.stat, float).ravel()) for s in np.atleast_1d(src)]
    nchan = sum(len(y) for _, y in subs)
    print(f"Loaded demo dataset: {len(subs)} subjects, {nchan} total electrodes.")

    print(f"Loading {mesh} surfaces...")
    surf, nL, nR, adj = ss.load_mesh(mesh)

    print(f"Running SMOOTHstat ({numrand} permutations, seed=42)...")
    t0 = time.time()
    res = ss.run_smoothstat(subs, surf, nL, nR, adj, numrand=numrand, clusteralpha=0.05, alpha=0.05,
                            minnbsub=3, tail="positive", seed=42, kernelwidth=12.0, graphsigma=8.0,
                            verbose=False)
    print(f"SMOOTHstat completed in {time.time() - t0:.1f}s.")

    sig = [(m, p, mem) for m, p, mem in res["posclusters"] if p <= 0.05]
    print(f"\n{len(sig)} significant positive cluster(s):")
    for m, p, mem in sig:
        print(f"   p = {p:.4f}   mass = {m:.1f}   n_vertices = {len(mem)}")
    if not args.quick:
        print("\nCompare with the manuscript (Figure 5a,d): pMTG p=0.001, rTPJ p=0.009.")

    # ---------------------------------------------------------------- null-distribution plot
    try:
        import matplotlib.pyplot as plt
        top_mass = res["posclusters"][0][0]
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(res["posdist"], bins=40, color="0.7", edgecolor="none")
        ax.axvline(top_mass, color="C3", lw=2, label=f"observed top cluster (mass={top_mass:.0f})")
        ax.set_xlabel("max cluster mass (null)"); ax.set_ylabel("count")
        ax.set_title("SMOOTH null distribution vs. observed top cluster"); ax.legend()
        fig.tight_layout(); fig.savefig(os.path.join(os.path.dirname(__file__), "smoothdemo_nulldist.png"), dpi=150)
        print("\nWrote demo/smoothdemo_nulldist.png")
    except ImportError:
        print("\n(matplotlib not installed -- skipping the null-distribution plot)")

    # ---------------------------------------------------------------- surrogate SA robustness check
    print("\nChecking surrogate spatial autocorrelation against the empirical maps (SMOOTHdiag equivalent)...")
    emp, surr_mu, surr_sd = ss.moran_sa(res, adj, nsurr=10, seed=7)
    ok = np.isfinite(emp) & np.isfinite(surr_mu)
    from scipy.stats import ttest_rel
    tstat, p = ttest_rel(emp[ok], surr_mu[ok])
    print(f"  paired t-test, empirical vs. surrogate Moran's I: t = {tstat:.2f}, p = {p:.3g}")
    print("  (p > 0.05 here means surrogate smoothness is indistinguishable from the empirical data --")
    print("   the expected result. A significant difference would mean cfg.kernelwidth/graphsigma need adjusting.)")


if __name__ == "__main__":
    main()
