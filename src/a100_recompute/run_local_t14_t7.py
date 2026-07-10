"""
Complete local reproduction of T14 (H3 activation control) and T7 (graph ablation).

T14: Uses gc_nk_subset.h5ad for gene-level NK activation signature.
T7:  Graph-structure ablation -- measures impact of metabolic_crosstalk/SST edges
     on the spectral gene embedding itself (no TCGA classifier needed -- directly
     tests the embedding quality that the manuscript claims as the graph's value).

Run: cd G:/cc/GC-NKGraph-Atlas && python src/a100_recompute/run_local_t14_t7.py
"""
import sys, os, time, io
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import linregress
from scipy.sparse.linalg import svds
from pathlib import Path

# Fix Windows GBK console encoding
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
T = "results/tables/"
GRAPH_DIR = "data/processed/graph"

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

print(f"T14 T7 local reproduction, {time.strftime('%Y-%m-%d %H:%M:%S')}")

# ============================================================================
# T14 -- H3 activation control with REAL gene-level NK expression from h5ad
# ============================================================================
log("=" * 60)
log("T14: H3 activation control (REAL gene expression from h5ad)")
log("=" * 60)

import anndata
ad = anndata.read_h5ad("data/processed/scrna/gc_nk_subset.h5ad")
log(f"Loaded h5ad: {ad.n_obs} NK cells x {ad.n_vars} genes")

NK_ACT = ["CD69","TNF","XCL1","XCL2","CCL3","CCL4","CCL5",
          "CSF2","IL2RA","ICOS","TNFSF10","FASLG","CD38",
          "HLA-DRA","HLA-DRB1","MKI67"]
act_present = [g for g in NK_ACT if g in ad.var_names]
log(f"Activation genes: {len(act_present)}/{len(NK_ACT)}")

PROTRUSION_GENES = ["EZR","MSN","RDX","ACTR2","ACTR3","ARPC2","ARPC3",
                     "WAS","WASL","CDC42","RAC1","RHOA","DIAPH1","FMNL1"]
CYTOTOXICITY_GENES = ["NKG7","GNLY","GZMB","PRF1","IFNG",
                       "LCP2","LAT","VAV1","TLN1","ITGAL","ITGB2"]
prot_genes = [g for g in PROTRUSION_GENES if g in ad.var_names]
cyto_genes = [g for g in CYTOTOXICITY_GENES if g in ad.var_names]
log(f"Protrusion genes: {len(prot_genes)}  Cytotoxicity genes: {len(cyto_genes)}")

# Get raw expression
expr = ad[:, prot_genes + cyto_genes + act_present].X
if hasattr(expr, 'toarray'):
    expr = expr.toarray()
elif hasattr(expr, 'todense'):
    expr = np.array(expr.todense())

# Z-score per gene
from scipy.stats import zscore
expr_z = zscore(expr, axis=0, nan_policy='omit')
expr_z = np.nan_to_num(expr_z, 0)

prot_score = expr_z[:, :len(prot_genes)].mean(axis=1)
cyto_score = expr_z[:, len(prot_genes):len(prot_genes)+len(cyto_genes)].mean(axis=1)
act_score  = expr_z[:, len(prot_genes)+len(cyto_genes):].mean(axis=1)

# Raw H3
r_raw, p_raw = stats.pearsonr(prot_score, cyto_score)
log(f"Raw H3:    r={r_raw:.4f}  P={p_raw:.2e}  r2={r_raw**2:.4f}  n={len(prot_score)}")

r_ap, _ = stats.pearsonr(act_score, prot_score)
r_ac, _ = stats.pearsonr(act_score, cyto_score)
log(f"Activation~protrusion: r={r_ap:.4f}  Activation~cytotoxicity: r={r_ac:.4f}")

# Partial: protrusion ~ cytotoxicity | activation
slope_xa, intercept_xa, _, _, _ = linregress(act_score, prot_score)
prot_resid = prot_score - (slope_xa * act_score + intercept_xa)
slope_ya, intercept_ya, _, _, _ = linregress(act_score, cyto_score)
cyto_resid = cyto_score - (slope_ya * act_score + intercept_ya)
r_partial, p_partial = stats.pearsonr(prot_resid, cyto_resid)

r2_raw = r_raw**2
r2_partial = r_partial**2
r2_act = r2_raw - r2_partial

log(f"Partial H3: r={r_partial:.4f}  P={p_partial:.2e}  r2={r2_partial:.4f}")
log(f"Activation component: r2={r2_act:.4f} ({100*r2_act/max(r2_raw,1e-10):.1f}% of raw r2)")

# Bulk check
log("Bulk liver check:")
try:
    liver = pd.read_csv(T + "sst_axis_scores_liver_bulk.tsv", sep="\t", index_col=0)
    p_b, c_b = liver["nk_protrusion_machinery_score"], liver["nk_synapse_cytotoxicity_outcome_score"]
    r_b, p_b_val = stats.pearsonr(p_b, c_b)
    log(f"  Bulk raw H3: r={r_b:.4f} P={p_b_val:.2e}")
    if "checkpoint_link_score" in liver.columns:
        act_b = -liver["checkpoint_link_score"].values
        v = ~(np.isnan(p_b) | np.isnan(c_b) | np.isnan(act_b))
        pb, cb, ab = p_b[v], c_b[v], act_b[v]
        sxa, _, _, _, _ = linregress(ab, pb)
        sya, _, _, _, _ = linregress(ab, cb)
        r_bp, p_bp = stats.pearsonr(pb - sxa*ab, cb - sya*ab)
        log(f"  Bulk partial H3: r={r_bp:.4f} P={p_bp:.2e}")
except Exception as e:
    log(f"  Skipped: {e}")

# Verdict
print()
if r_partial > 0.15 and p_partial < 0.001:
    print("VERDICT: H3 ROBUST to activation control")
    print("  The 'effector arm recovers' claim SURVIVES.")
    print(f"  Only {100*r2_act/max(r2_raw,1e-10):.0f}% of r2 shared with activation.")
    print(f"  Residual protrusion~cytotoxicity r={r_partial:.3f} (r2={r2_partial:.3f}) is a real, independent effect.")
    rec = "robust"
elif r_partial > 0.05 and p_partial < 0.05:
    print("VERDICT: H3 PARTIALLY survives activation control")
    rec = "partial"
else:
    print("VERDICT: H3 DOES NOT SURVIVE activation control")
    rec = "not_robust"

pd.DataFrame([{
    "resolution":"single-cell NK","n":len(prot_score),
    "r_raw":round(r_raw,4),"r2_raw":round(r2_raw,5),
    "r_partial_activation":round(r_partial,4),"r2_partial":round(r2_partial,5),
    "r2_activation_component":round(r2_act,5),
    "p_raw":p_raw,"p_partial":p_partial,
    "activation_genes_n":len(act_present),
    "activation_genes":",".join(act_present),
    "verdict":rec,
}]).to_csv(T + "h3_activation_control.tsv", sep="\t", index=False)
log(f"T14 PASS -- verdict: {rec}")

# ============================================================================
# T7 -- Graph-structure ablation (spectral embedding quality)
# ============================================================================
log("\n" + "=" * 60)
log("T7: Graph ablation -- spectral embedding structure")
log("=" * 60)

nodes = pd.read_csv(f"{GRAPH_DIR}/nodes.tsv", sep="\t")
edges = pd.read_csv(f"{GRAPH_DIR}/edges.tsv", sep="\t")
log(f"Graph: {len(nodes)} nodes, {len(edges)} edges")
for et, n in edges["edge_type"].value_counts().items():
    log(f"  {et}: {n}")

TUMOR_SERINE = ["PHGDH","PSAT1","PSPH","SHMT1","SHMT2","MTHFD1","MTHFD2","MTHFD1L","SLC1A4","SLC1A5"]
NK_SM = ["SGMS1","SGMS2","SMPD1","SMPD2","SMPD3","SMPD4"]
NK_PROT = ["EZR","MSN","RDX","ACTR2","ACTR3","ARPC2","ARPC3","WAS","WASL","CDC42","RAC1","RHOA","DIAPH1","FMNL1"]
NK_CYTO = ["NKG7","GNLY","GZMB","PRF1","IFNG","LCP2","LAT","VAV1","TLN1","ITGAL","ITGB2"]

def build_adj(edges_df, nodes_df):
    node_ids = nodes_df["node_id"].tolist()
    n2i = {str(n): i for i, n in enumerate(node_ids)}
    n = len(node_ids)
    adj = np.zeros((n,n), dtype=np.float64)
    for _, r in edges_df.iterrows():
        s, d = str(r["src"]), str(r["dst"])
        if s in n2i and d in n2i:
            i, j = n2i[s], n2i[d]
            w = float(r.get("weight", 1.0))
            adj[i,j] = w; adj[j,i] = w
    deg = adj.sum(axis=1); deg = np.where(deg > 0, deg, 1.0)
    deg_is = np.diag(1.0 / np.sqrt(deg))
    return deg_is @ adj @ deg_is, n2i

def spectral_embed(adj_norm, dim=128):
    k = min(dim, adj_norm.shape[0] - 2)
    u, s, vt = svds(adj_norm.astype(np.float64), k=k)
    idx = np.argsort(s)[::-1]
    emb = u[:, idx] @ np.diag(np.sqrt(np.maximum(s[idx], 0)))
    if emb.shape[1] < dim:
        pad = np.zeros((emb.shape[0], dim - emb.shape[1])); emb = np.hstack([emb, pad])
    else:
        emb = emb[:, :dim]
    return emb

def module_coupling(emb, node_to_idx, genes_a, genes_b):
    idx_a = [node_to_idx[g] for g in genes_a if g in node_to_idx]
    idx_b = [node_to_idx[g] for g in genes_b if g in node_to_idx]
    if len(idx_a) < 2 or len(idx_b) < 2:
        return np.nan
    cross = 0.0; npairs = 0
    for ia in idx_a:
        for ib in idx_b:
            cross += np.dot(emb[ia], emb[ib]) / (np.linalg.norm(emb[ia]) * np.linalg.norm(emb[ib]) + 1e-10)
            npairs += 1
    return cross / npairs

VARIANTS = {
    "FULL":  set(),
    "-MC":   {"metabolic_crosstalk"},
    "-SST":  {"metabolic_crosstalk", "sm_topology_axis"},
}
LABELS = {
    "FULL": "All edges",
    "-MC":  "Without metabolic_crosstalk",
    "-SST": "Without all SST edges",
}

results = []
for variant_key, drop_set in VARIANTS.items():
    log(f"\n{variant_key}: {LABELS[variant_key]}")
    masked = edges[~edges["edge_type"].isin(drop_set)]
    log(f"  Edges: {len(masked)} (dropped {len(edges)-len(masked)})")

    adj_norm, n2i = build_adj(masked, nodes)
    emb = spectral_embed(adj_norm, dim=128)

    h2 = module_coupling(emb, n2i, NK_SM, NK_PROT)
    h3 = module_coupling(emb, n2i, NK_PROT, NK_CYTO)
    h1 = module_coupling(emb, n2i, TUMOR_SERINE, NK_SM)

    # Modularity: intra vs inter module similarity
    modules = {"serine": TUMOR_SERINE, "sm": NK_SM, "prot": NK_PROT, "cyto": NK_CYTO}
    mod_genes = {m: [g for g in gs if g in n2i] for m, gs in modules.items()}
    mnames = list(mod_genes.keys())
    intra_sim, inter_sim = 0.0, 0.0
    ni, nj = 0, 0
    for mi in range(len(mnames)):
        for mj in range(mi, len(mnames)):
            c = module_coupling(emb, n2i, mod_genes[mnames[mi]], mod_genes[mnames[mj]])
            if np.isnan(c): continue
            if mi == mj: intra_sim += c; ni += 1
            else: inter_sim += c; nj += 1
    intra_sim /= max(ni, 1); inter_sim /= max(nj, 1)
    modularity = intra_sim - inter_sim

    log(f"  H2(SM<->prot)={h2:.4f}  H3(prot<->cyto)={h3:.4f}  H1(serine<->SM)={h1:.4f}")
    log(f"  Intra={intra_sim:.4f}  Inter={inter_sim:.4f}  Modularity={modularity:.4f}")

    results.append({
        "variant": variant_key, "n_edges": len(masked),
        "H2_SM_prot_coupling": round(h2,5),
        "H3_prot_cyto_coupling": round(h3,5),
        "H1_serine_SM_coupling": round(h1,5),
        "intra_module_sim": round(intra_sim,5),
        "inter_module_sim": round(inter_sim,5),
        "modularity": round(modularity,5),
    })

df = pd.DataFrame(results)
df.to_csv(T + "ablation_results.tsv", sep="\t", index=False)

print(f"\n{'='*70}")
print("ABLATION RESULTS -- Spectral Embedding Quality")
print(f"{'='*70}")
for _, r in df.iterrows():
    print(f"  {r['variant']:6s} edges={r['n_edges']:5d}  H2={r['H2_SM_prot_coupling']:+.4f}  H3={r['H3_prot_cyto_coupling']:+.4f}  H1={r['H1_serine_SM_coupling']:+.4f}  mod={r['modularity']:+.4f}")

full = df[df["variant"] == "FULL"]
minus_mc = df[df["variant"] == "-MC"]
if len(full) and len(minus_mc):
    print(f"\n--- Delta (FULL - -MC) ---")
    for col in ["H2_SM_prot_coupling","H3_prot_cyto_coupling","H1_serine_SM_coupling","modularity"]:
        d = full[col].values[0] - minus_mc[col].values[0]
        print(f"  {col}: {d:+.5f}")
    mod_d = full["modularity"].values[0] - minus_mc["modularity"].values[0]
    h2_d = full["H2_SM_prot_coupling"].values[0] - minus_mc["H2_SM_prot_coupling"].values[0]
    print()
    if mod_d > 0.01 or h2_d > 0.005:
        print(f"VERDICT: metabolic_crosstalk edge ADDS structure (modularity +{mod_d:.4f}, H2 +{h2_d:.4f})")
        print("  -> The edge measurably shapes the gene embedding.")
        print("  -> The graph design is NOT cosmetic -- it has a confirmable structural effect.")
    elif abs(mod_d) < 0.002 and abs(h2_d) < 0.001:
        print(f"VERDICT: metabolic_crosstalk edge has NEGLIGIBLE impact (modularity {mod_d:+.4f}, H2 {h2_d:+.4f})")
        print("  -> Edge weight (0.5) may be too low relative to PPI (0.7-1.0).")
        print("  -> Consider: increase weight OR argue that topological existence matters more than weight magnitude.")
    else:
        print(f"VERDICT: metabolic_crosstalk edge has MODEST impact (modularity {mod_d:+.4f}, H2 {h2_d:+.4f})")

log("T7 PASS")
print(f"\n{'='*70}")
print("T14 + T7 COMPLETE")
print(f"Outputs: {T}h3_activation_control.tsv  +  {T}ablation_results.tsv")
print(f"{'='*70}")
