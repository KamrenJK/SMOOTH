"""ROI-only (Desikan-Killiany / Schaefer-300) detection for the FDR atlas benchmark (Supplementary
Figure 5). ROIOnlySim skips SMOOTH's whole-cortex smoothing + cluster machinery entirely (the P@y
smoothing-operator application, and its nperm-times-per-iteration reconstruction in _null_dist, are the
dominant per-iteration cost and are never touched by the ROI/atlas arms -- they only need the raw
per-electrode planted signal, parcel-averaged). Each atlas is tested two ways on the same draw: the
sign-flip max-statistic permutation FWER of the main benchmark, and a per-parcel one-tailed t-test with
Benjamini-Hochberg FDR across parcels.

SMOOTH's own sensitivity is NOT recomputed here -- the SMOOTH column comes from the cached
paper/data/per_vertex_maps.npz. This module draws its own independent (same-model) 1/f realizations per
seed, shared by both atlas tests in each iteration but not identical to the draws behind the cached
SMOOTH numbers. roi_fdr_full_sweep.py runs it over every testable vertex.

Usage (single seed): python roi_fdr_lean.py <seed> [niter] [nperm]
"""
import os, sys, time

PROJ = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
os.environ.setdefault("PROJ", PROJ)
os.environ.setdefault("SENS_MESH_EIG_TAG", "5000")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from scipy.stats import t as tdist
import config as C
import sensitivity_core as SC

ATLASES = ["aparc", "schaefer300"]


class ROIOnlySim(SC.BrainSim):
    def _draw_ys(self, rng, seed):
        seed_h = "lh" if seed < self.nL else "rh"
        sp_sph = self.sph[seed]
        ys = {}
        for si, sub in enumerate(self.subs):
            ys[si] = {}
            for h, _ in self._slices():
                o = sub[h]
                if o is None: continue
                y = self._draw_noise(o, rng)
                if h == seed_h:
                    y = y + C.A0 * np.exp(-np.sum((o["esph"] - sp_sph) ** 2, 1) / (2 * self.sig_b_sph ** 2))
                ys[si][h] = y
        return ys

    def _hit_mask(self, seed):
        sp_sph = self.sph[seed]
        dvert = np.full(self.nV, np.inf)
        sl = slice(0, self.nL) if seed < self.nL else slice(self.nL, self.nV)
        dvert[sl] = np.linalg.norm(self.sph[sl] - sp_sph, axis=1)
        return dvert < SC.MM2SPHERE * C.R_HIT

    def _roi_detect_fdr(self, atl, ys, hit_mask):
        P = self._parc[atl]; ncol = P["ncol"]; nsub = len(self.subs)
        PM = np.full((nsub, ncol), np.nan)
        for si in range(nsub):
            ssum = np.zeros(ncol); scnt = np.zeros(ncol)
            for h in ("lh", "rh"):
                ec = P["ecols"][si][h]
                if ec is None: continue
                yv = ys[si][h]
                for c in np.unique(ec[ec >= 0]):
                    mm = ec == c; ssum[c] += yv[mm].sum(); scnt[c] += mm.sum()
            g = scnt > 0; PM[si, g] = ssum[g] / scnt[g]
        n = np.sum(~np.isnan(PM), 0)
        with np.errstate(invalid="ignore", divide="ignore"):
            mu = np.nanmean(PM, 0); sd = np.sqrt(np.nanvar(PM, 0, ddof=1))
            t = mu / (sd / np.sqrt(n))
        t[(n < C.MINNBSUB) | ~np.isfinite(t)] = np.nan
        fin = np.where(np.isfinite(t))[0]
        if len(fin) == 0: return 0.0, 0.0
        df = n[fin] - 1
        pv = tdist.sf(t[fin], df)
        order = np.argsort(pv); m = len(pv)
        ranked = pv[order]
        thresh = (np.arange(1, m + 1) / m) * self.sig_thresh
        below = ranked <= thresh
        if not below.any(): return 0.0, 0.0
        cutoff = ranked[np.max(np.where(below)[0])]
        sig_cols = fin[pv <= cutoff]
        on = set(int(c) for c in np.unique(P["vcol"][hit_mask & (P["vcol"] >= 0)]))
        is_tp = is_fp = False
        for i in sig_cols:
            if int(i) in on: is_tp = True
            else: is_fp = True
        return float(is_tp), float(is_fp)

    def run_seed_roi(self, seed, niter, nperm):
        rng = np.random.default_rng([C.BASE_SEED, 99, int(seed)])   # independent stream tag (not matched to cached SMOOTH draws)
        hit = self._hit_mask(seed)
        keys = []
        for atl in self.roi_atlases:
            keys += [f"roi_{atl}_tp", f"roi_{atl}_fp", f"roi_{atl}_fdr_tp", f"roi_{atl}_fdr_fp"]
        acc = {k: [] for k in keys}
        for _ in range(niter):
            ys = self._draw_ys(rng, seed)
            for atl in self.roi_atlases:
                rtp, rfp = self._roi_detect(atl, ys, hit, rng, nperm)
                acc[f"roi_{atl}_tp"].append(rtp); acc[f"roi_{atl}_fp"].append(rfp)
                ftp, ffp = self._roi_detect_fdr(atl, ys, hit)
                acc[f"roi_{atl}_fdr_tp"].append(ftp); acc[f"roi_{atl}_fdr_fp"].append(ffp)
        out = {}
        for atl in self.roi_atlases:
            out[f"roi_{atl}_sens"] = float(np.mean(acc[f"roi_{atl}_tp"]))
            out[f"roi_{atl}_fpr"] = float(np.mean(acc[f"roi_{atl}_fp"]))
            out[f"roi_{atl}_fdr_sens"] = float(np.mean(acc[f"roi_{atl}_fdr_tp"]))
            out[f"roi_{atl}_fdr_fpr"] = float(np.mean(acc[f"roi_{atl}_fdr_fp"]))
        return out


if __name__ == "__main__":
    seed = int(sys.argv[1])
    niter = int(sys.argv[2]) if len(sys.argv) > 2 else C.NITER
    nperm = int(sys.argv[3]) if len(sys.argv) > 3 else C.NPERM
    sim = ROIOnlySim(noise="1f", beta=C.BETA, sa_a=C.SA_A, roi_atlases=ATLASES)
    assert seed in set(sim.seeds.tolist())
    t0 = time.time()
    res = sim.run_seed_roi(seed, niter, nperm)
    dt = time.time() - t0
    print(f"seed={seed} niter={niter} nperm={nperm}  wall={dt:.2f}s")
    for k in sorted(res): print(f"  {k:24s} {res[k]:.4f}")
