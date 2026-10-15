#!/usr/bin/env python3
"""Whole-brain FWER (false-positive control) under 1/f-GRF noise, four arms per null dataset:
SMOOTH (eigenmode sign-flip cluster), NAIVE (electrode shuffle), and two ROI comparators
(Desikan-Killiany 'aparc', Schaefer300). One driver, three sweep MODES:

  beta     noise=1f,       sweep 20 log-spaced integer beta in [1,50]   (Figure 3: FP control vs smoothness)
  aniso    noise=aniso1f,  per-rho SA-matched beta, sweep rho           (robustness: anisotropy)
  nonstat  noise=nonstat1f, symmetric mix about beta=29, sweep q        (robustness: nonstationarity)

Each mode is run at the nugget conditions in SENS_NUGS (default BOTH: a=0.53 with-nugget, matched to
sensitivity; and a=1.0 no-nugget, the adversarial pure-1/f case). SENS_NUGS='1.0' or '0.53' runs one
condition; output files are keyed by a, so staged runs never collide and a rerun skips what is already
done; merge_fwer.py aggregates whatever partials are present.
ca=0.05, one-tailed alpha=0.05 (config-locked).

Cluster layout -- a SLURM job array over a flat task grid  (cell = one noise config x nugget):
    task = cell*CHUNKS + chunk        cell in [0, NCELL),  chunk in [0, CHUNKS)
each task runs its strided subset of NSIM null datasets on all NCORES and writes a partial-count npz;
merge_fwer.py sums counts across chunks -> FWER + Wilson CI -> figure. CHUNKS defaults to 1 (one task
per cell); raise SENS_FWER_CHUNKS to spread a cell across more array tasks when nodes are plentiful.

  python fwer.py --mode beta --task $SLURM_ARRAY_TASK_ID           # one array task
  python fwer.py --mode beta --task all                           # whole sweep locally (all cells x chunks)
"""
import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _sys.path.insert(0, _os.path.join(_HERE, ".."))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    _os.environ.setdefault(_v, "1")
import argparse, time, numpy as np
from multiprocessing import get_context
import sensitivity_core as SC, config as C

ATL = ["aparc", "schaefer300"]
def _nugs():
    """Nugget conditions (correlated fraction a) to sweep; default both. SENS_NUGS='1.0' runs only
    no-nugget, '0.53' only with-nugget. Staged runs coexist (files keyed by a)."""
    return [float(x) for x in _os.environ.get("SENS_NUGS", "0.53,1.0").split(",") if x.strip()]
def _betas():
    """20 strictly-increasing integer betas in [1,50], log-spaced (rounding collisions bumped up), with the
    empirical slope C.BETA snapped in so FWER is scored exactly at the fitted value."""
    b = []
    for x in np.geomspace(1, 50, 20):
        v = int(round(x))
        if b and v <= b[-1]: v = b[-1] + 1
        b.append(v)
    emp = int(round(C.BETA))
    if emp not in b: b[int(np.argmin(np.abs(np.array(b) - emp)))] = emp
    return sorted(set(b))
BETAS = _betas()                                                     # 20 log-spaced integer beta in [1,50], incl. empirical
CHUNKS = int(_os.environ.get("SENS_FWER_CHUNKS", "1"))

def _aniso_beta_match():
    """rho -> SA-matched 1/f slope from calibrate_aniso_beta.py. The aniso field lives on a different
    (radius-graph) basis than the beta sweep, so C.BETA=29 is far too smooth there; each rho gets its own
    beta fit so its Moran's I matches the empirical DVs (holds smoothness fixed, varies only the shape)."""
    p = _os.path.join(C.OUT_DIR, "aniso_beta_match.npz")
    if not _os.path.exists(p):
        raise SystemExit(f"aniso needs SA-matched betas -- run calibrate_aniso_beta.py first (missing {p})")
    d = np.load(p); return {int(round(r)): float(b) for r, b in zip(d["rho"], d["beta"])}

def grid(mode, nugs=None):
    """List of noise cells for a mode: one per (config value x nugget). kwargs feed BrainSim."""
    nugs = _nugs() if nugs is None else nugs
    if mode == "beta":   vals = [("beta", b, dict(noise="1f", beta=float(b))) for b in BETAS]
    elif mode == "aniso":
        bm = _aniso_beta_match()                                       # {rho: SA-matched beta}, calibrated for rho=1..10
        rvals = [int(x) for x in _os.environ.get("SENS_ANISO_RHO", "1,4,7,10").split(",") if x.strip()]
        miss = [r for r in rvals if r not in bm]
        if miss: raise SystemExit(f"aniso rho {miss} not in aniso_beta_match.npz -- run calibrate_aniso_beta.py for them (or pick calibrated rho)")
        vals = [("rho", r, dict(noise="aniso1f", beta=bm[r], rho=float(r))) for r in rvals]
    elif mode == "nonstat":     # SA-preserving symmetric mix around SENS_NONSTAT_BETA; q capped so beta*sqrt(q) stays non-degenerate
        nb = float(_os.environ.get("SENS_NONSTAT_BETA", C.BETA))                 # anchor (default 29 -> structure SA; e.g. 9.09 -> raw-empirical SA)
        qvals = [float(x) for x in _os.environ.get("SENS_NONSTAT_Q", "1,2,3,4").split(",") if x.strip()]
        vals = [("q", q, dict(noise="nonstat1f", beta=nb, q=q)) for q in qvals]
    else: raise SystemExit(f"bad mode {mode}")
    cells = []
    for a in nugs:
        for pname, pval, kw in vals:
            cells.append(dict(pname=pname, pval=pval, a=a, kw={**kw, "sa_a": a, "roi_atlases": ATL}))
    return cells

_S = {}
# per-simulation fields kept for the manuscript QQ / null-max-cluster-mass distribution panels (Figure 3),
# not just the counts.
PERSIM = ["smooth_p", "naive_p", "obs_max", "smooth_p95", "naive_p95", "observed_SA"]
def _trial(i):
    r = _S["sim"].run_null(np.random.default_rng([C.BASE_SEED, 909, i]), _S["nperm"])
    return tuple(r[k] for k in PERSIM) + (int(r["roi_aparc_sig"]), int(r["roi_schaefer300_sig"]))

def run_cell(cell, chunk, nsim, nperm, ncores):
    idx = list(range(chunk, nsim, CHUNKS))
    _S["sim"] = SC.BrainSim(**cell["kw"]); _S["nperm"] = nperm   # built ONCE in parent; fork workers share the 5000-mode caches copy-on-write
    thr = _S["sim"].sig_thresh                                    # one-tailed cutoff used to count false positives
    t0 = time.time()
    with get_context("fork").Pool(ncores) as pool:                # fork (not spawn) so workers inherit the shared BrainSim
        rows = pool.map(_trial, idx)
    cols = list(zip(*rows)); A = {k: np.asarray(cols[j], float) for j, k in enumerate(PERSIM)}
    roi_ap = np.asarray(cols[len(PERSIM)], float); roi_s3 = np.asarray(cols[len(PERSIM)+1], float)
    return dict(pname=cell["pname"], pval=cell["pval"], a=cell["a"], n=len(idx),
                k_smooth=int((A["smooth_p"] <= thr).sum()), k_naive=int((A["naive_p"] <= thr).sum()),
                k_aparc=int(roi_ap.sum()), k_schaefer300=int(roi_s3.sum()),
                SA=float(np.median(A["observed_SA"])), secs=time.time()-t0,
                **A)                                              # full per-sim arrays (smooth_p, naive_p, obs_max, p95s, observed_SA)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["beta", "aniso", "nonstat"])
    ap.add_argument("--task", default="all", help="SLURM array task id, or 'all' to run every cell x chunk locally")
    ap.add_argument("--nsim", type=int, default=C.NSIM)
    ap.add_argument("--nperm", type=int, default=C.NPERM)
    ap.add_argument("--ncores", type=int, default=C.NCORES)
    ap.add_argument("--nugs", default=None, help="comma-separated nugget a-values (default env SENS_NUGS or 0.53,1.0)")
    ap.add_argument("--count", action="store_true", help="print the array-task count (ncells*chunks) and exit")
    a = ap.parse_args()
    nugs = [float(x) for x in a.nugs.split(",") if x.strip()] if a.nugs else None
    cells = grid(a.mode, nugs); ncell = len(cells); ntask = ncell*CHUNKS
    if a.count: print(ntask); return
    od = _os.path.join(C.OUT_DIR, f"fwer_{a.mode}")
    if a.mode == "nonstat":                                          # keep different anchors in separate dirs so a re-anchor never overwrites a prior run
        _nb = float(_os.environ.get("SENS_NONSTAT_BETA", C.BETA))
        if abs(_nb - C.BETA) > 1e-9: od += f"_b{int(round(_nb))}"    # e.g. beta=9.09 -> out/<ds>/fwer_nonstat_b9/ ; default 29 stays fwer_nonstat/
    _os.makedirs(od, exist_ok=True)
    tasks = range(ntask) if a.task == "all" else [int(a.task)]
    print(f"FWER mode={a.mode} | nugs a={sorted(set(c['a'] for c in cells))} "
          f"| {ncell} cells x {CHUNKS} chunks = {ntask} tasks | nsim={a.nsim} nperm={a.nperm} "
          f"ncores={a.ncores} | ca={C.CLUSTERALPHA} tail={C.TAIL} thr={C.ALPHA}", flush=True)
    for task in tasks:
        cell = cells[task // CHUNKS]; chunk = task % CHUNKS
        out = _os.path.join(od, f"a{cell['a']:g}_{cell['pname']}{cell['pval']:g}_chunk{chunk:02d}.npz")
        if _os.path.exists(out): print(f"  task {task}: exists, skip"); continue
        r = run_cell(cell, chunk, a.nsim, a.nperm, a.ncores)
        np.savez(out, **{k: v for k, v in r.items()})
        print(f"  task {task:>3} [{r['pname']}={r['pval']:g} a={r['a']:g} chunk {chunk}] n={r['n']} "
              f"SA={r['SA']:.3f} smooth={r['k_smooth']}/{r['n']} naive={r['k_naive']}/{r['n']} "
              f"aparc={r['k_aparc']} s300={r['k_schaefer300']} | {r['secs']:.0f}s", flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
