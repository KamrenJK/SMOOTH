#!/usr/bin/env python3
"""Whole-brain FWER (specificity) sweep: the no-bump twin of the sensitivity simulation. For each swept
value (correlation length L, anisotropy ratio rho, nonstationarity ratio q, or 1/f exponent beta), draw
`nsim` INDEPENDENT null datasets (per-subject spatially autocorrelated noise, NO effect), and on BOTH
hemispheres (one pooled test) score the observed max positive-cluster mass under TWO nulls: SMOOTH
(eigenmode-rotation, preserves spatial autocorrelation) and NAIVE (DV<->electrode shuffle, destroys it ->
anti-conservative). Records the whole-brain FWER for each arm plus per-trial quantities (p-values,
null-mass 95th pctiles, observed max-mass, observed electrode SA / Moran's I).

Supplementary Figure 3 uses the 1/f modes (see run_sweeps.sh):
  python fwer_sweep.py --noise aniso1f --rhos 1 2 3 5 10 --a 1.0 --beta 8 --outtag aniso1f_nonug --nsim 1000 --nperm 1000
  python fwer_sweep.py --noise nonstat1f --qs 1 2 4 6 8 --a 1.0 --beta 10 --outtag nonstat1f_nonug --nsim 1000 --nperm 1000

Writes <config.OUT_DIR>/fwer_sweep_<outtag>.npz: fwer_smooth / fwer_naive vs the swept value, and
per-(value, trial) arrays smooth_p, naive_p, obs_max, smooth_p95, naive_p95, observed_SA, smooth_sig, naive_sig.
"""
import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _sys.path.insert(0, _HERE)
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # 1 BLAS thread/worker so N workers use N cores
import argparse, time, numpy as np
from multiprocessing import Pool
import config as C, sensitivity_core as SC

_S = {}
def _init(noise, val, baseL, a, beta, nperm):
    if noise == "aniso":
        _S["sim"] = SC.BrainSim(sa_a=a, sa_l=baseL, noise="aniso", rho=val)          # SE, sweep anisotropy ratio rho
    elif noise == "nonstat":
        _S["sim"] = SC.BrainSim(sa_a=a, sa_l=baseL, noise="nonstat", q=val)          # SE, sweep nonstationarity ratio q
    elif noise == "aniso1f":
        _S["sim"] = SC.BrainSim(sa_a=a, noise="aniso1f", rho=val, beta=beta)         # 1/f on warped Laplacian, sweep rho
    elif noise == "nonstat1f":
        _S["sim"] = SC.BrainSim(sa_a=a, noise="nonstat1f", q=val, beta=beta)         # A-P mixture of two 1/f spectra, sweep q
    elif noise == "1f":
        _S["sim"] = SC.BrainSim(sa_a=a, noise="1f", beta=val)                        # isotropic 1/f, sweep beta
    else:
        _S["sim"] = SC.BrainSim(sa_a=a, sa_l=val)                                    # sweep correlation length L (isotropic SE)
    _S["nperm"] = nperm

def _trial(key):
    rng = np.random.default_rng(list(key))
    return _S["sim"].run_null(rng, _S["nperm"])                             # both hemispheres; SMOOTH + naive arms

# per-trial quantities saved for each swept value
FIELDS = ["smooth_p", "naive_p", "obs_max", "smooth_p95", "naive_p95", "observed_SA", "smooth_sig", "naive_sig"]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Ls", type=float, nargs="+", default=[1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
                    help="correlation lengths (mm) to sweep")
    ap.add_argument("--a", type=float, default=None, help="correlated variance fraction (default 1.0 for se, SA_A for aniso/nonstat)")
    ap.add_argument("--noise", choices=["se", "aniso", "nonstat", "1f", "aniso1f", "nonstat1f"], default="se",
                    help="null noise model: se|aniso|nonstat (metric-space SE) | 1f|aniso1f|nonstat1f (graph-Laplacian power-law)")
    ap.add_argument("--rhos", type=float, nargs="+", default=[1, 2, 3, 5, 10], help="anisotropy ratios (noise=aniso/aniso1f)")
    ap.add_argument("--qs", type=float, nargs="+", default=[1, 2, 4, 6, 8], help="nonstationarity ratios (noise=nonstat/nonstat1f)")
    ap.add_argument("--betas", type=float, nargs="+", default=[10.0], help="1/f exponents to sweep (noise=1f)")
    ap.add_argument("--beta", type=float, default=10.0, help="fixed 1/f exponent for aniso1f/nonstat1f")
    ap.add_argument("--baseL", type=float, default=C.SA_L, help="fixed base correlation length (mm) for aniso/nonstat")
    ap.add_argument("--outtag", default="", help="output suffix -> fwer_sweep_<outtag>.npz (e.g. aniso1f_nonug); overrides the noise-based default")
    ap.add_argument("--nsim", type=int, default=1000, help="null datasets per swept value")
    ap.add_argument("--nperm", type=int, default=1000, help="permutations per dataset (SMOOTH + naive nulls)")
    ap.add_argument("--ncores", type=int, default=1)
    args = ap.parse_args()
    noise = args.noise
    if noise in ("aniso", "aniso1f"):     vals, label = args.rhos, "rho"
    elif noise in ("nonstat", "nonstat1f"): vals, label = args.qs, "q"
    elif noise == "1f":                   vals, label = args.betas, "beta"
    else:                                 vals, label = args.Ls, "L"
    aval = args.a if args.a is not None else (1.0 if noise == "se" else C.SA_A)     # default: empirical nugget for the non-SE models
    os.makedirs(C.OUT_DIR, exist_ok=True)
    print(f"FWER sweep (SMOOTH + naive-shuffle) | dataset={C.DATASET} | noise={noise} sweep {label}={vals} | "
          f"a={aval} baseL={args.baseL} | nsim={args.nsim} nperm={args.nperm} | "
          f"kw={C.KERNELWIDTH} gs={C.GRAPHSIGMA} alpha={C.ALPHA} tail={C.TAIL} ca={C.CLUSTERALPHA}", flush=True)
    tag = ("_"+args.outtag) if args.outtag else ("" if noise == "se" else f"_{noise}")
    out = os.path.join(C.OUT_DIR, f"fwer_sweep{tag}.npz")
    xs = []; fwer_s = []; fwer_n = []; acc = {f: [] for f in FIELDS}; t0 = time.time()
    for val in vals:
        tv = time.time(); vid = int(round(val*1000))
        with Pool(args.ncores, initializer=_init, initargs=(noise, val, args.baseL, aval, args.beta, args.nperm)) as pool:
            res = pool.map(_trial, [(vid, i) for i in range(args.nsim)])
        xs.append(float(val))
        for f in FIELDS: acc[f].append([r[f] for r in res])
        fwer_s.append(float(np.mean([r["smooth_sig"] for r in res])))
        fwer_n.append(float(np.mean([r["naive_sig"] for r in res])))
        print(f"  {label}={val:6.2f}  FWER smooth={fwer_s[-1]:.3f} naive={fwer_n[-1]:.3f}  "
              f"| median obs-SA={np.median(acc['observed_SA'][-1]):.3f}  [{(time.time()-tv)/60:.1f} min]", flush=True)
        np.savez(out, vals=np.array(xs), param=label, noise=noise, Ls=np.array(xs),   # Ls: alias of vals (swept axis)
                 fwer_smooth=np.array(fwer_s), fwer_naive=np.array(fwer_n), fwer=np.array(fwer_s),
                 **{f: np.array(acc[f], float) for f in FIELDS},                        # each (nvals, nsim)
                 a=aval, baseL=args.baseL, nsim=args.nsim, nperm=args.nperm, alpha=C.ALPHA,
                 clusteralpha=C.CLUSTERALPHA, kernelwidth=C.KERNELWIDTH, graphsigma=C.GRAPHSIGMA,
                 tail=C.TAIL, dataset=C.DATASET)
    print(f"\nDONE ({(time.time()-t0)/3600:.2f} h)  ->  {out}", flush=True)

if __name__ == "__main__":
    main()
