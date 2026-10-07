#!/usr/bin/env python3
"""Export paper/data/fig1_data.mat for genFig1.m (Figure 1: graphical overview).

Runs the SMOOTH group test on the RIGHT hemisphere over all socialcog (primary dataset) subjects
using each subject's REAL per-electrode z-stat as the DV (no simulated bump), finds the significant
cluster on right TPJ, and picks 3 exemplar subjects covering that cluster with positive effects.
Exports per-subject coverage/smoothed maps, 3 sign-flip surrogates each, and the graph + eigenmodes +
spectra for each exemplar -- all the pieces for the Figure 1 workflow schematic. Reuses
sensitivity_core from a separate simulation pipeline (not distributed), and needs the fsaverage5
Schaefer-300 (17-network) right-hemisphere annot, so -- like every other one-off exporter in this codebase -- this script cannot
be run from the published repository alone; only its cached output (fig1_data.mat) is needed to
regenerate the figure.

  SMOOTH_SCHAEFER_ANNOT=<rh.Schaefer2018_300Parcels_17Networks_order.annot> SMOOTH_PRIVATECODE=<simulation-pipeline code dir>
  SENS_MESH_EIG_TAG=5000 CS_SURF_FS5=<.../fsaverage5> python export_fig1.py
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
if not os.environ.get("SMOOTH_SCHAEFER_ANNOT"): sys.exit("Set SMOOTH_SCHAEFER_ANNOT to the fsaverage5 rh Schaefer-300 (17-network) annot file.")
if not os.environ.get("SMOOTH_PRIVATECODE"): sys.exit("Set SMOOTH_PRIVATECODE to the simulation-pipeline code directory (not distributed).")
sys.path.insert(0, os.environ["SMOOTH_PRIVATECODE"])
import numpy as np, nibabel as nib
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.stats import norm
from scipy.io import savemat
import config as C, cs_common as CC, sensitivity_core as SC

DATA = os.path.abspath(os.path.join(HERE, "..", "data"))
FS = os.environ.get("SENS_FS_SURF", "/Applications/freesurfer/7.3.2/subjects/fsaverage5/surf")
HEMI = "rh"; NPERM = int(os.environ.get("OV_NPERM", "300")); NSURR = 3; FLIP_SEED = 11
rng = np.random.default_rng(7)
base = os.path.join(C.SURF_DIR, C.MESH); off = np.array(C.FT_OFFSET)

# ---------------- geometry ----------------
vpos, faces = SC.read_fs_surface(f"{base}/{HEMI}.pial"); vpos = vpos + off
sph, _ = SC.read_fs_surface(f"{base}/{HEMI}.sphere.reg"); sph = sph + off
curv = nib.freesurfer.read_morph_data(f"{FS}/{HEMI}.curv").astype(float)
nV = len(vpos)
Ei = faces[:, [0, 1, 2]].ravel(); Ej = faces[:, [1, 2, 0]].ravel()
adj = (csr_matrix((np.ones(len(Ei)*2), (np.r_[Ei, Ej], np.r_[Ej, Ei])), shape=(nV, nV)) > 0).astype(np.int8)
zpos = norm.ppf(1 - C.CLUSTERALPHA)
lam_m, Qm = SC.load_mesh_eig(HEMI); Pw = 1.0/np.power(lam_m + np.median(lam_m), C.BETA); Pw[0] = 0.0
vt = cKDTree(vpos)
ap = np.load(os.path.join(DATA, "aparc_fs5.npz"), allow_pickle=True)
rh_names = [str(x) for x in ap["rh_names"]]; rh_lab = ap["rh_labels"].astype(int)
tpj_parcels = [i for i, nm in enumerate(rh_names) if nm in ("supramarginal", "inferiorparietal")]
tpj_vtx = np.isin(rh_lab, tpj_parcels)                                          # DK rTPJ vertices (reference overlap)

# ---------------- all rh subjects: real DV (socialcog z-stat) -> operators + smoothed map ----------------
subs = CC.load_socialcog(None)
ops = []
for sid, (cp, st) in enumerate(subs):
    cp = np.asarray(cp, float); st = np.asarray(st, float).ravel()
    m = cp[:, 0] >= 0; elec = cp[m]; dv = st[m]
    if len(elec) < 3: continue
    P, esph, sup = SC.hemi_smoother(elec, vpos, sph, C.KERNELWIDTH)
    Q, lam = SC.graph_eigenbasis(esph, C.GRAPHSIGMA)
    emp = np.asarray(P @ dv).ravel(); emp[~sup] = np.nan
    ops.append(dict(id=sid+1, elec=elec, esph=esph, P=P, Q=Q, lam=lam, sup=sup, dv=dv,
                    mc=Q.T @ dv, emp=emp, tgt=np.sort(emp[sup]), vidx=vt.query(elec)[1]))
print(f"{len(ops)} rh subjects | real socialcog z-stat DV", flush=True)

# ---------------- real SMOOTH group test over ALL subjects ----------------
allmaps = np.column_stack([o["emp"] for o in ops]); coverage = np.sum(~np.isnan(allmaps), 1); mask = coverage < C.MINNBSUB
tmap = SC.group_t(allmaps, mask, C.MINNBSUB); obs = SC.pos_clusters(tmap, adj, zpos)

def surrogate_map(o):
    sf = np.sign(rng.standard_normal(len(o["mc"]))); sf[sf == 0] = 1
    m = np.asarray(o["P"] @ (o["Q"] @ (o["mc"]*sf))).ravel(); m[~o["sup"]] = np.nan
    fin = o["sup"]; order = np.argsort(m[fin]); tmp = m[fin].copy(); tmp[order] = o["tgt"]; m[fin] = tmp
    return m
posdist = np.empty(NPERM); surr_tmaps = []
for p in range(NPERM):
    pt = SC.group_t(np.column_stack([surrogate_map(o) for o in ops]), mask, C.MINNBSUB)
    posdist[p] = max((mm for mm, _ in SC.pos_clusters(pt, adj, zpos)), default=0.0)
    if p < NSURR: surr_tmaps.append(pt)
obs.sort(key=lambda x: -x[0])
sig = np.zeros(nV, bool); clusters = []
for mass, mem in obs:
    pv = (np.sum(posdist >= mass)+1)/(NPERM+1)
    if pv <= C.ALPHA: sig[mem] = True; clusters.append((mass, mem, pv))
# rTPJ cluster = the significant cluster with the most overlap with the DK rTPJ parcels
tpj_cluster = max(clusters, key=lambda cm: np.sum(tpj_vtx[cm[1]]), default=None)
tpj_mem = tpj_cluster[1]; tpj_p = tpj_cluster[2]; tpj_mass = tpj_cluster[0]
sig_t = np.where(sig, tmap, np.nan)
in_tpj = np.zeros(nV, bool); in_tpj[tpj_mem] = True
print(f"  {len(clusters)} significant clusters | rTPJ cluster: mass {tpj_mass:.0f}, p {tpj_p:.3f}, "
      f"{len(tpj_mem)} vtx ({100*np.mean(tpj_vtx[tpj_mem]):.0f}% in DK supramarginal/IP)", flush=True)

# ---------------- pick 3 exemplars: most electrodes inside the rTPJ cluster with positive DV ----------------
for o in ops:
    inc = in_tpj[o["vidx"]]; o["n_tpj_pos"] = int(np.sum(inc & (o["dv"] > 0))); o["mean_tpj"] = float(np.mean(o["dv"][inc])) if inc.any() else -np.inf
shown = sorted(ops, key=lambda o: -(o["n_tpj_pos"] + 0.1*(o["mean_tpj"] if np.isfinite(o["mean_tpj"]) else 0)))[:3]
print(f"  exemplars: {[(o['id'], o['n_tpj_pos'], round(o['mean_tpj'],1)) for o in shown]}", flush=True)

# ---------------- per-exemplar: surrogates (3), graph edges, eigenmodes, spectra ----------------
def graph_edges(o):
    Un = o["esph"]/np.linalg.norm(o["esph"], axis=1, keepdims=True)
    D = 100.0*np.arccos(np.clip(Un @ Un.T, -1, 1)); gs = SC.MM2SPHERE*C.GRAPHSIGMA
    W = np.exp(-D**2/(2*gs**2)); np.fill_diagonal(W, 0.0); ed = []
    for i in range(len(o["elec"])):
        for j in np.argsort(W[i])[::-1][:3]:
            a, b = (int(i), int(j)) if i < j else (int(j), int(i)); ed.append((a, b, float(W[i, j])))
    g = np.unique(np.array(ed), axis=0); g[:, :2] += 1; return g
MODES = [1, 3, 6, 12, 22]
ex_elec = np.array([o["elec"] for o in shown], dtype=object)
ex_dv   = np.array([o["dv"] for o in shown], dtype=object)
ex_emp  = np.column_stack([o["emp"] for o in shown])
ex_surr = np.array([np.column_stack([surrogate_map(o) for _ in range(NSURR)]) for o in shown])   # 3 x nV x NSURR
ex_graph = np.array([graph_edges(o) for o in shown], dtype=object)
ex_eig  = np.array([np.column_stack([o["Q"][:, k]/(np.max(np.abs(o["Q"][:, k]))+1e-12) for k in MODES]) for o in shown], dtype=object)
ex_mc   = np.array([o["mc"] for o in shown], dtype=object)
ex_mcsurr = np.array([o["mc"]*(lambda sf: np.where(sf == 0, 1, sf))(np.sign(np.random.default_rng(FLIP_SEED+i).standard_normal(len(o["mc"]))))
                     for i, o in enumerate(shown)], dtype=object)

# ---------------- all electrodes: subject index + DK/Schaefer parcel colours (electrodes coloured by parcel) ----------------
allE = np.vstack([o["elec"] for o in ops]); vidx_all = vt.query(allE)[1]
all_sub = np.concatenate([np.full(len(o["elec"]), si) for si, o in enumerate(ops)]).astype(float)
FSLAB = os.path.join(os.path.dirname(FS), "label")
dk_lab, dk_ctab, _ = nib.freesurfer.read_annot(os.path.join(FSLAB, "rh.aparc.annot"))
sch_annot = os.environ["SMOOTH_SCHAEFER_ANNOT"]
sch_lab, sch_ctab, _ = nib.freesurfer.read_annot(sch_annot)
def parcel_rgb(lab, ctab, vi):
    l = lab[vi].copy(); l[l < 0] = 0; return ctab[l, :3].astype(float)/255.0
all_dk_rgb = parcel_rgb(dk_lab, dk_ctab, vidx_all); all_schaefer_rgb = parcel_rgb(sch_lab, sch_ctab, vidx_all)
print(f"  all electrodes: {len(allE)} | {int(all_sub.max())+1} subjects | DK {len(dk_ctab)} parcels, Schaefer {len(sch_ctab)} parcels", flush=True)

savemat(os.path.join(DATA, "fig1_data.mat"), dict(
    pial_pos=vpos, pial_tri=(faces+1).astype(float), pial_curv=curv, hemi=HEMI, n_all=float(len(ops)),
    tmap=tmap, sig_t=sig_t, tpj_cluster=in_tpj.astype(float), tpj_p=float(tpj_p), tpj_mass=float(tpj_mass),
    coverage=coverage.astype(float), tpj_vtx=tpj_vtx.astype(float),
    surr_tmap=np.column_stack(surr_tmaps), posdist=posdist, obs_mass=float(obs[0][0]),
    all_elec=allE, all_dv=np.concatenate([o["dv"] for o in ops]), all_sub=all_sub,
    all_dk_rgb=all_dk_rgb, all_schaefer_rgb=all_schaefer_rgb,
    sub_ids=np.array([o["id"] for o in shown], float), n_sub=float(len(shown)),
    sub_elec=ex_elec, sub_dv=ex_dv, sub_emp=ex_emp, sub_surr=ex_surr,
    graph_edges=ex_graph, eig=ex_eig, mode_idx=np.array(MODES, float), mc=ex_mc, mc_surr=ex_mcsurr,
    params=dict(kernelwidth=C.KERNELWIDTH, graphsigma=C.GRAPHSIGMA, beta=C.BETA,
                minnbsub=C.MINNBSUB, clusteralpha=C.CLUSTERALPHA, alpha=C.ALPHA),
))
print(f"wrote {os.path.join(DATA, 'fig1_data.mat')}", flush=True)
