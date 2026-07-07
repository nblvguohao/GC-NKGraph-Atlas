"""
Add co-expression edges from scRNA NK data to the heterogeneous graph.
Computes Pearson correlation between target genes across NK cells.
"""
import os, sys, time, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

log("Loading NK subset...")
adata = sc.read("data/processed/scrna/gc_nk_subset.h5ad")
if scipy.sparse.issparse(adata.X):
    expr = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)
else:
    expr = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)

cand = pd.read_csv("results/tables/candidate_evidence_matrix.tsv", sep="\t")
target_genes = [g for g in cand["gene"] if g in expr.columns]
log(f"Computing co-expression for {len(target_genes)} genes across {expr.shape[0]} NK cells...")

nk_expr = expr[target_genes]
corr = nk_expr.corr(method="pearson")

edges = []
for i, g1 in enumerate(target_genes):
    for j, g2 in enumerate(target_genes):
        if i < j:
            r = corr.iloc[i, j]
            if abs(r) > 0.3:
                edges.append({"src": g1, "dst": g2, "edge_type": "coexpression",
                              "weight": round(abs(r), 4),
                              "evidence": f"scRNA_NK_corr_{r:.3f}"})

edge_df = pd.DataFrame(edges)
log(f"Co-expression edges: {len(edge_df)}")
edge_df.to_csv("data/processed/graph/coexpression_edges.tsv", sep="\t", index=False)
log("Saved")

# Also update main edges.tsv
main_edges = pd.read_csv("data/processed/graph/edges.tsv", sep="\t")
all_edges = pd.concat([main_edges, edge_df], ignore_index=True)
all_edges.to_csv("data/processed/graph/edges.tsv", sep="\t", index=False)
log(f"Final edges.tsv: {len(all_edges)} total edges")
