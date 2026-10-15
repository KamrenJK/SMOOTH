#!/usr/bin/env python3
"""Aggregate the fwer.py partial-count npz files for a MODE (beta|aniso|nonstat), compute FWER +
Wilson 95% CI per (config value, nugget, arm), write a summary npz, and render the figure: two panels
(with-nugget a=0.53 | no-nugget a=1.0), four arms each, alpha line at 0.05. In beta mode the empirical
slope beta=29 is annotated; x is log-beta. In aniso/nonstat modes x is rho / q.
  python merge_fwer.py --mode beta"""
import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _sys.path.insert(0, _os.path.join(_HERE, ".."))
import argparse, glob, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import config as C

ARMS = [("smooth", "SMOOTH", "#1f6f2e"), ("naive", "naive shuffle", "#8a2417"),
        ("aparc", "ROI: Desikan-Killiany", "#8a6d1f"), ("schaefer300", "ROI: Schaefer300", "#3a4a8a")]
XLAB = {"beta": "1/f slope β  (rough → smooth)", "aniso": "anisotropy ρ", "nonstat": "nonstationarity q"}

def wilson(k, n, z=1.96):
    if n == 0: return np.nan, np.nan, np.nan
    p = k/n; d = 1 + z*z/n
    lo = (p + z*z/(2*n) - z*np.sqrt(p*(1-p)/n + z*z/(4*n*n)))/d
    hi = (p + z*z/(2*n) + z*np.sqrt(p*(1-p)/n + z*z/(4*n*n)))/d
    return p, max(0, lo), min(1, hi)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--mode", required=True, choices=["beta", "aniso", "nonstat"]); a = ap.parse_args()
    od = _os.path.join(C.OUT_DIR, f"fwer_{a.mode}")
    files = sorted(glob.glob(_os.path.join(od, "*_chunk*.npz")))     # a<a>_<param><val>_chunk<c>.npz (any nugget present)
    if not files: raise SystemExit(f"no partial files in {od}")
    agg = {}                                                             # (pval,a) -> dict of summed counts + n + SA-list
    for f in files:
        d = np.load(f); key = (float(d["pval"]), float(d["a"]))
        e = agg.setdefault(key, dict(n=0, SA=[], **{f"k_{arm}": 0 for arm, _, _ in ARMS}))
        e["n"] += int(d["n"]); e["SA"].append(float(d["SA"]))
        for arm, _, _ in ARMS: e[f"k_{arm}"] += int(d[f"k_{arm}"])
    keys = sorted(agg)                                                   # sorted by (pval, a)
    avals = sorted(set(k[1] for k in keys), reverse=True)               # [0.53, 1.0] -> with-nug left
    fig, axes = plt.subplots(1, len(avals), figsize=(6.6*len(avals), 5.0), squeeze=False)
    rows = []                                                            # for summary npz
    for pan, aval in enumerate(avals):
        ax = axes[0][pan]; xs = sorted(k[0] for k in keys if k[1] == aval)
        for arm, lbl, col in ARMS:
            P, LO, HI = [], [], []
            for x in xs:
                e = agg[(x, aval)]; p, lo, hi = wilson(e[f"k_{arm}"], e["n"]); P.append(p); LO.append(lo); HI.append(hi)
            P, LO, HI = map(np.array, (P, LO, HI))
            ax.errorbar(xs, P, [P-LO, HI-P], fmt="-o", color=col, capsize=3, lw=1.8, ms=6, label=lbl)
            for x, p, lo, hi in zip(xs, P, LO, HI):
                rows.append((a.mode, x, aval, arm, agg[(x, aval)]["n"], agg[(x, aval)][f"k_{arm}"], p, lo, hi))
        ax.axhline(C.ALPHA, ls=":", color="0.4", lw=1.3, label=f"α = {C.ALPHA:g}")
        if a.mode == "beta":
            ax.set_xscale("log"); ax.set_xticks(xs); ax.set_xticklabels([f"{x:g}" for x in xs], fontsize=7, rotation=45)
            ax.axvline(C.BETA, ls="--", color="0.55", lw=1.2); ax.text(C.BETA, 0.98, " empirical β", transform=ax.get_xaxis_transform(), va="top", fontsize=8, color="0.4")
        else:
            ax.set_xticks(xs)
        ax.set_xlabel(XLAB[a.mode]); ax.set_ylabel(f"FWER (one-tailed {C.ALPHA:g})"); ax.set_ylim(0, 1)
        ax.set_title(f"{'with nugget (a=0.53)' if aval < 1 else 'no nugget (a=1.0)'}")
        ax.spines[["top", "right"]].set_visible(False); ax.grid(alpha=0.25)
        if pan == 0: ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.suptitle(f"False-positive control — {a.mode} sweep (1/f GRF, whole-brain, ca={C.CLUSTERALPHA:g})", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    figdir = _os.path.join(C.OUT_DIR, "figs"); _os.makedirs(figdir, exist_ok=True)
    for ext in ("png", "pdf"): fig.savefig(_os.path.join(figdir, f"fwer_{a.mode}.{ext}"), dpi=150)
    dt = np.dtype([("mode", "U8"), ("x", float), ("a", float), ("arm", "U12"), ("n", int), ("k", int), ("fwer", float), ("lo", float), ("hi", float)])
    np.savez(_os.path.join(od, f"fwer_{a.mode}_summary.npz"), table=np.array(rows, dtype=dt))
    print(f"wrote {figdir}/fwer_{a.mode}.png  and  {od}/fwer_{a.mode}_summary.npz  ({len(files)} partials, {len(keys)} cells)")

if __name__ == "__main__":
    main()
