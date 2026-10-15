#!/usr/bin/env python3
"""Whole-brain detection sensitivity under 1/f-GRF noise (beta=29, with-nugget a=0.53, ca=0.05).
For every testable fs5 vertex (seed) a geodesic bump is planted in ONE hemisphere and the whole brain
goes through ONE SMOOTH cluster test with a single pooled permutation null; per seed we record SMOOTH
sensitivity / FPR / localization error plus the two ROI comparators (Desikan-Killiany, Schaefer300).
Electrode coverage comes from config.ELEC_FILE (the primary dataset).

Cluster-ready: a SLURM job array strides the seed list (--chunk/--nchunks); each task runs on all
NCORES and is resumable (skips seeds already in its chunk file). merge_sensitivity.py stitches chunks
into whole-brain per-vertex maps (per_vertex_maps.npz).

  python run_sensitivity.py --chunk $SLURM_ARRAY_TASK_ID --nchunks 200 --ncores $SLURM_CPUS_PER_TASK
  python run_sensitivity.py --scan            # A0 effect-size calibration (writes nothing)
"""
import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _sys.path.insert(0, _os.path.join(_HERE, ".."))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    _os.environ.setdefault(_v, "1")
import argparse, time, numpy as np
from multiprocessing import get_context
import sensitivity_core as SC, config as C

ATL = ["aparc", "schaefer300"]
_S = {}
def _work(arg):
    seed, niter, nperm = arg
    return seed, _S["sim"].run_seed(seed, niter, nperm)

def _save(outfile, res):
    ks = sorted(res); keys = sorted(res[ks[0]]) if ks else []
    d = {"seed": np.array(ks, int)}
    for m in keys: d[m] = np.array([res[k][m] for k in ks], float)
    np.savez(outfile, **d)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunk", type=int, default=0); ap.add_argument("--nchunks", type=int, default=1)
    ap.add_argument("--ncores", type=int, default=C.NCORES)
    ap.add_argument("--niter", type=int, default=C.NITER); ap.add_argument("--nperm", type=int, default=C.NPERM)
    ap.add_argument("--beta", type=float, default=C.BETA); ap.add_argument("--sa_a", type=float, default=C.SA_A)
    ap.add_argument("--scan", action="store_true")
    a = ap.parse_args()
    _os.makedirs(C.OUT_DIR, exist_ok=True)
    kw = dict(noise="1f", beta=a.beta, sa_a=a.sa_a, roi_atlases=ATL)
    sim = SC.BrainSim(**kw); _S["sim"] = sim                     # built ONCE in parent; fork workers share it copy-on-write
    print(f"[{C.DATASET}] whole-brain 1/f | beta={a.beta} a={a.sa_a} | {len(sim.subs)} subjects | {sim.nV} vertices "
          f"| {len(sim.seeds)} testable seeds | A0={C.A0} bumpFWHM={C.BUMP_FWHM} ca={C.CLUSTERALPHA} tail={C.TAIL}", flush=True)

    if a.scan:
        order = np.argsort(sim.coverage[sim.seeds]); sub = sim.seeds[order][::max(1, len(sim.seeds)//24)][:24]
        print(f"A-scan on {len(sub)} seeds spanning coverage (niter=15, nperm={a.nperm}):", flush=True)
        for A in [2, 3, 4, 6, 8]:
            C.A0 = float(A); s = [sim.run_seed(int(v), 15, a.nperm)["sens"] for v in sub]
            print(f"  A0={A}: median sens={np.median(s):.2f}  mean={np.mean(s):.2f}", flush=True)
        return

    my = sim.seeds[a.chunk::a.nchunks]
    outfile = _os.path.join(C.OUT_DIR, f"sens_wb_c{a.chunk}of{a.nchunks}.npz")
    res = {}
    if _os.path.exists(outfile):
        d = np.load(outfile); ms = [k for k in d.files if k != "seed"]
        res = {int(s): {m: float(d[m][i]) for m in ms} for i, s in enumerate(d["seed"])}
    todo = [int(s) for s in my if int(s) not in res]
    print(f"[chunk {a.chunk}/{a.nchunks}] {len(my)} seeds | {len(todo)} to run, {len(res)} cached "
          f"| ncores={a.ncores} niter={a.niter} nperm={a.nperm}", flush=True)
    if not todo: _save(outfile, res); print("nothing to do; wrote", outfile); return
    t0 = time.time(); n = 0
    with get_context("fork").Pool(a.ncores) as pool:               # fork so workers inherit the shared BrainSim
        for seed, metrics in pool.imap_unordered(_work, [(v, a.niter, a.nperm) for v in todo], chunksize=1):
            res[seed] = metrics; n += 1
            if n % 20 == 0 or n == len(todo):
                _save(outfile, res); el = time.time()-t0; rate = n/el
                print(f"  {n}/{len(todo)} | {rate*3600:.0f} seeds/h | ETA {(len(todo)-n)/max(rate,1e-9)/3600:.2f}h", flush=True)
    _save(outfile, res)
    print(f"[chunk {a.chunk}] DONE {len(res)} seeds -> {outfile}", flush=True)

if __name__ == "__main__":
    main()
