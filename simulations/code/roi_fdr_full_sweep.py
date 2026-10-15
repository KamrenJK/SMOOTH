"""Full 14,075-testable-vertex FDR-only sweep for DK (aparc) and Schaefer-300, niter=1000, no
permutations, no SMOOTH recomputation (Supplementary Figure 5). SMOOTH's own sensitivity is reused as-is
from the cached per_vertex_maps.npz; this script only produces the roi_{atlas}_fdr_{sens,fpr} arrays,
index-aligned to that file's `seeds` array.

Parallelized across local cores via multiprocessing (spawn, macOS default). Each worker builds one
ROIOnlySim once and processes its chunk of seeds sequentially (round-robin split, so any systematic
cost variation across the seed list doesn't unbalance one worker), checkpointing to its own .npz
every CHECKPOINT_EVERY seeds so a mid-run interruption only loses a few seeds of work, not the chunk.

Usage: python roi_fdr_full_sweep.py [nworkers]        # -> $PROJ/out/roi_fdr_sweep/chunk_*.npz
Then viz/export_benchmark_fdr_mat.py merges the chunk files into paper/data/benchmark_fig_data_fdr.mat.
"""
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

PROJ = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
os.environ.setdefault("PROJ", PROJ)
os.environ.setdefault("SENS_MESH_EIG_TAG", "5000")
CODE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, CODE_DIR)

import numpy as np
import multiprocessing as mp

NITER = 1000
ATLASES = ("aparc", "schaefer300")
OUT_DIR = os.path.join(PROJ, "out", "roi_fdr_sweep")
CHECKPOINT_EVERY = 25
PER_VERTEX_MAPS = os.environ.get("SENS_PER_VERTEX_MAPS") or os.path.join(
    os.environ.get("SMOOTH_REPO") or os.path.dirname(PROJ), "paper", "data", "per_vertex_maps.npz")


def worker(args):
    worker_id, seed_chunk = args
    sys.path.insert(0, CODE_DIR)
    import config as C
    from roi_fdr_lean import ROIOnlySim
    sim = ROIOnlySim(noise="1f", beta=C.BETA, sa_a=C.SA_A, roi_atlases=list(ATLASES))
    n = len(seed_chunk)
    res = {f"roi_{atl}_fdr_sens": np.full(n, np.nan) for atl in ATLASES}
    res.update({f"roi_{atl}_fdr_fpr": np.full(n, np.nan) for atl in ATLASES})
    out_path = os.path.join(OUT_DIR, f"chunk_{worker_id:03d}.npz")
    t0 = time.time()
    for i, seed in enumerate(seed_chunk):
        seed = int(seed)
        rng = np.random.default_rng([C.BASE_SEED, 99, seed])
        hit = sim._hit_mask(seed)
        acc = {f"{atl}_tp": [] for atl in ATLASES}
        acc.update({f"{atl}_fp": [] for atl in ATLASES})
        for _ in range(NITER):
            ys = sim._draw_ys(rng, seed)
            for atl in ATLASES:
                tp, fp = sim._roi_detect_fdr(atl, ys, hit)
                acc[f"{atl}_tp"].append(tp); acc[f"{atl}_fp"].append(fp)
        for atl in ATLASES:
            res[f"roi_{atl}_fdr_sens"][i] = float(np.mean(acc[f"{atl}_tp"]))
            res[f"roi_{atl}_fdr_fpr"][i] = float(np.mean(acc[f"{atl}_fp"]))
        if (i + 1) % CHECKPOINT_EVERY == 0 or (i + 1) == n:
            np.savez(out_path, seeds=np.asarray(seed_chunk), **res)
            print(f"[worker {worker_id:03d}] {i + 1}/{n} seeds done, {time.time() - t0:.0f}s elapsed", flush=True)
    return worker_id


if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)
    nworkers = int(sys.argv[1]) if len(sys.argv) > 1 else 11
    pv = np.load(PER_VERTEX_MAPS)
    seeds = pv["seeds"]
    chunks = [(w, seeds[w::nworkers]) for w in range(nworkers)]
    print(f"{len(seeds)} testable seeds, {nworkers} workers, ~{len(seeds) // nworkers} seeds/worker, niter={NITER}")
    t0 = time.time()
    with mp.get_context("spawn").Pool(nworkers) as pool:
        for wid in pool.imap_unordered(worker, chunks):
            print(f"worker {wid:03d} finished, {time.time() - t0:.0f}s since start", flush=True)
    print(f"ALL DONE in {time.time() - t0:.0f}s")
