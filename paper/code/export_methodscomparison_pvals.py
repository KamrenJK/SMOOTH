#!/usr/bin/env python3
"""Export FDR-corrected channel-level + ROI (Desikan-Killiany and Schaefer-300) results for the
methods-comparison supplementary figure (genFigS6.m), for the three datasets. Everything
is ONE-TAILED POSITIVE (matching the manuscript's cfg.tail='positive'); FDR is Benjamini-Hochberg (q<0.05).

  channel-level : per-electrode p from the dataset's DV null; BH-FDR across ALL electrodes -> chan_sig.
  ROI (DK/Sch)  : each electrode -> its fs5 parcel; per parcel, one-sample t across subjects on the
                  parcel-mean DV; one-tailed positive; BH-FDR across parcels (>=3 subjects) -> sig parcels.

Per-channel p by dataset (from the same DVs SMOOTH is run on):
  socialcog    : DV = permutation z-score      -> p = 1 - Phi(z)
  berezutskaya : DV = Fisher-z; stored pearsonr pvalue + sign(r)  -> one-tailed positive (n=nblocks)
  cogitate     : DV = onset_r = point-biserial R; t = R*sqrt((N-2)/(1-R^2)), N = 2*n_trials EXACT
                 (n_trials = non-target stimulus-onset trials, read per subject from events.tsv; no rejection)

Renders on fsaverage5 (both atlases exist as fs5 npz). Emits per-fs5-vertex significant-parcel t-maps so
the MATLAB figure only loads + renders. Writes paper/data/methodscomparison_pvals.mat.

Requires the Berezutskaya and COGITATE per-subject DV derivatives and fsaverage5 atlas lookups, which
are not distributed; only the cached output is needed to regenerate the figure.
"""
import os, glob, numpy as np, pandas as pd
from scipy.io import loadmat, savemat
from scipy.stats import norm, t as tdist
from scipy.spatial import cKDTree
import nibabel as nib

FS   = os.environ.get("SUBJECTS_DIR") or (os.path.join(os.environ["FREESURFER_HOME"], "subjects")
       if os.environ.get("FREESURFER_HOME") else "/Applications/freesurfer/7.3.2/subjects")
_HERE = os.path.abspath(__file__)
PUB  = os.path.dirname(os.path.dirname(os.path.dirname(_HERE)))                    # repository root
DATA = os.path.join(PUB, "paper", "data")
if not os.environ.get("SMOOTH_DERIVATIVES"):
    raise SystemExit("Set SMOOTH_DERIVATIVES to the per-dataset DV derivatives root (not distributed).")
ROOT = os.environ["SMOOTH_DERIVATIVES"]
FT   = np.array([1.0, 1.0, 1.0])
DVLIM = {"socialcog": 10.0, "berezutskaya": 2.0, "cogitate": 1.0}                  # per-dataset DV colour limit
def _bad(n): n = n.lower(); return ("unknown" in n or "medial" in n or "background" in n or "corpuscallosum" in n or n == "")

def bh_fdr(p):
    p = np.asarray(p, float); adj = np.full(p.shape, np.nan); ok = np.isfinite(p)
    q = p[ok]; n = q.size
    if n == 0: return adj
    o = np.argsort(q); r = np.empty(n); r[o] = np.minimum.accumulate((q[o]*n/np.arange(1, n+1))[::-1])[::-1]
    adj[ok] = np.clip(r, 0, 1); return adj

# ---------------------------------------------------------------- fs5 geometry + atlases
def load_fs5():
    P = {}; trees = {}
    for h in ("lh", "rh"):
        pos, tri = nib.freesurfer.read_geometry(f"{FS}/fsaverage5/surf/{h}.pial")
        curv = nib.freesurfer.read_morph_data(f"{FS}/fsaverage5/surf/{h}.curv")
        P[h] = (pos + FT, tri, curv); trees[h] = cKDTree(pos + FT)
    nL = len(P["lh"][0])
    pos = np.vstack([P["lh"][0], P["rh"][0]]); tri = np.vstack([P["lh"][1], P["rh"][1] + nL])
    curv = np.concatenate([P["lh"][2], P["rh"][2]])
    return pos, tri, curv, nL, trees

def load_atlas_key(atlas, nL):
    """Per-fs5-vertex parcel key 'h:name' (''=medial-wall/unknown), whole-brain length 2*nL."""
    d = np.load(os.path.join(DATA, f"{atlas}_fs5.npz"), allow_pickle=True)
    key = np.empty(2*nL, dtype=object); key[:] = ""
    for h, off in (("lh", 0), ("rh", nL)):
        lab = d[f"{h}_labels"].astype(int); nm = [str(x) for x in d[f"{h}_names"]]
        for v, l in enumerate(lab):
            if l > 0 and not _bad(nm[l]): key[off+v] = f"{h}:{nm[l]}"
    return key

def elec_keys(chanpos, trees, key, nL):
    out = np.empty(len(chanpos), dtype=object); out[:] = ""
    for h, off in (("lh", 0), ("rh", nL)):
        m = (chanpos[:, 0] < 0) if h == "lh" else (chanpos[:, 0] >= 0)
        if not m.any(): continue
        out[np.where(m)[0]] = key[off + trees[h].query(chanpos[m])[1]]
    return out

def roi_sig_tmap(chanpos, dv, sub, trees, key, nL):
    """Per-parcel one-sample t across subjects (one-tailed +), BH-FDR; return a per-fs5-vertex t-map
    (t where the vertex's parcel is FDR-significant, NaN elsewhere) + the parcel table."""
    ek = elec_keys(chanpos, trees, key, nL)
    parcel_t = {}; names = []; ts = []; ps = []; ns = []
    for k in sorted(set(ek[ek != ""])):
        m = ek == k; subs = np.unique(sub[m])
        if len(subs) < 3: continue
        vals = np.array([dv[m & (sub == s)].mean() for s in subs])
        mu, sd = vals.mean(), vals.std(ddof=1)
        tv = mu/(sd/np.sqrt(len(subs))) if sd > 0 else np.nan
        parcel_t[k] = tv; names.append(k); ts.append(tv); ps.append(tdist.sf(tv, len(subs)-1)); ns.append(len(subs))
    padj = bh_fdr(ps); sig = padj <= 0.05
    sigkeys = {names[i]: ts[i] for i in range(len(names)) if sig[i]}
    vtx = np.full(2*nL, np.nan)
    for v, k in enumerate(key):
        if k in sigkeys: vtx[v] = sigkeys[k]
    tab = dict(roi=np.array(names, dtype=object), t=np.array(ts), p=np.array(ps), p_fdr=padj,
               sig=sig.astype(float), nsub=np.array(ns))
    return vtx, tab

# ---------------------------------------------------------------- per-dataset DV + channel p (one-tailed +)
def get_socialcog():
    S = np.atleast_1d(loadmat(os.path.join(PUB, "demo", "source.mat"), squeeze_me=True, struct_as_record=False)["source"])
    cp, dv, p, sub = [], [], [], []
    for i, s in enumerate(S):
        z = np.asarray(s.stat, float).ravel(); cp.append(np.asarray(s.elec.chanpos, float))
        dv.append(z); p.append(1 - norm.cdf(z)); sub.append(np.full(len(z), i))
    return list(map(np.concatenate, (cp, dv, p, sub)))

def get_berezutskaya():
    cp, dv, p, sub = [], [], [], []
    for i, f in enumerate(sorted(glob.glob(os.path.join(ROOT, "ds003688-pilot", "out", "subjects", "sub-*.npz")))):
        d = np.load(f); r = d["r"]; pv = d["pvalue"]
        cp.append(np.asarray(d["chanpos"], float)); dv.append(d["fisher_z"])
        p.append(np.where(r >= 0, pv/2.0, 1 - pv/2.0)); sub.append(np.full(len(r), i))
    return list(map(np.concatenate, (cp, dv, p, sub)))

def get_cogitate():
    cp, dv, p, sub = [], [], [], []
    for i, f in enumerate(sorted(glob.glob(os.path.join(ROOT, "cogitate-smooth", "out", "dv", "sub-*_dv.mat")))):
        d = loadmat(f, squeeze_me=True); R = np.asarray(d["onset_r"], float).ravel()
        if R.size < 3 or "chanpos" not in d: continue
        sid = os.path.basename(f).replace("_dv.mat", "")
        ntr = 0
        for evf in glob.glob(os.path.join(ROOT, "cogitate-smooth", "data", sid, "**", "*_events.tsv"), recursive=True):
            tt = pd.read_csv(evf, sep="\t")["trial_type"].astype(str)
            ntr += int((tt.str.contains("stimulus onset") &
                        tt.apply(lambda t: any(k in t for k in ("face", "object", "letter", "false"))) &
                        ~tt.str.contains("Relevant target", case=False)).sum())
        df = 2*ntr - 2                                                            # onset point-biserial: 2*n_trials samples
        tval = R*np.sqrt(df/np.clip(1 - R**2, 1e-12, None))
        cp.append(np.asarray(d["chanpos"], float)); dv.append(R); p.append(tdist.sf(tval, df)); sub.append(np.full(len(R), i))
    return list(map(np.concatenate, (cp, dv, p, sub)))

# ---------------------------------------------------------------- run
pos, tri, curv, nL, trees = load_fs5()
KEY = {"dk": load_atlas_key("aparc", nL), "sch": load_atlas_key("schaefer300", nL)}
out = dict(fs5_pos=pos, fs5_tri=(tri+1).astype(float), fs5_curv=curv, nL=float(nL))
for name, getter in (("socialcog", get_socialcog), ("berezutskaya", get_berezutskaya), ("cogitate", get_cogitate)):
    print(f"[{name}]", flush=True)
    cp, dv, p, sub = getter()
    chan_fdr = bh_fdr(p); chan_sig = chan_fdr <= 0.05
    dk_vtx, dk_tab = roi_sig_tmap(cp, dv, sub, trees, KEY["dk"], nL)
    sch_vtx, sch_tab = roi_sig_tmap(cp, dv, sub, trees, KEY["sch"], nL)
    print(f"  channels {len(p)} | chan sig {int(np.nansum(chan_sig))} ({100*np.nanmean(chan_sig):.0f}%) | "
          f"DK sig {int(dk_tab['sig'].sum())}/{dk_tab['roi'].size} | Sch sig {int(sch_tab['sig'].sum())}/{sch_tab['roi'].size}", flush=True)
    out.update({
        f"{name}_chanpos": cp, f"{name}_dv": dv, f"{name}_sub": sub.astype(float),
        f"{name}_chan_p": p, f"{name}_chan_p_fdr": chan_fdr, f"{name}_chan_sig": chan_sig.astype(float),
        f"{name}_dvlim": DVLIM[name], f"{name}_dk_sigt": dk_vtx, f"{name}_sch_sigt": sch_vtx,
        f"{name}_dk_roi": dk_tab["roi"], f"{name}_dk_sig": dk_tab["sig"],
        f"{name}_sch_roi": sch_tab["roi"], f"{name}_sch_sig": sch_tab["sig"],
    })
dst = os.path.join(DATA, "methodscomparison_pvals.mat"); savemat(dst, out); print("wrote", dst)
