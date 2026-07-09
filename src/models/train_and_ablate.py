"""
Train GNN on real TCGA-STAD data with real STRING PPI graph.
Produces gene embeddings (T11) and ablation results (T12).
"""
import os, sys, time, json, warnings
from pathlib import Path
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
_PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_PROJECT_ROOT))

log = lambda msg: print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

# =====================================================================
# STEP 1: Build real graph from STRING PPI + SST modules
# =====================================================================
log("Building real gene graph...")

# Load TCGA genes
stad_expr = pd.read_csv("data/processed/bulk/tcga_stad_expression.tsv", sep="\t", index_col=0)
tcga_genes = set(stad_expr.columns)
log(f"TCGA-STAD: {stad_expr.shape[0]} samples x {stad_expr.shape[1]} genes")

# Load STRING PPI
string_edges = pd.read_csv("data/processed/graph/string_ppi_high_conf.tsv", sep="\t")
log(f"STRING edges: {len(string_edges)}")

# Map to gene symbols (strip transcript suffixes like "-201")
import re
string_pairs = []
for _, r in string_edges.iterrows():
    g1 = re.sub(r'-\d{3}$', '', str(r['protein1_gene']))
    g2 = re.sub(r'-\d{3}$', '', str(r['protein2_gene']))
    if g1 in tcga_genes and g2 in tcga_genes:
        string_pairs.append({'src': g1, 'dst': g2, 'edge_type': 'ppi',
                             'weight': float(r['combined_score']) / 1000.0})

log(f"STRING pairs in TCGA genes: {len(string_pairs)}")

# SST module genes
SST_MODULES = {
    'tumor_serine': ['PHGDH','PSAT1','PSPH','SHMT1','SHMT2','MTHFD1','MTHFD2','MTHFD1L','SLC1A4','SLC1A5'],
    'nk_sm_synthesis': ['SGMS1','SGMS2'],
    'nk_sm_catabolism': ['SMPD1','SMPD2','SMPD3','SMPD4'],
    'nk_denovo': ['SPTLC1','SPTLC2','SPTLC3','SPTSSA','CERS2','CERS4','CERS5','CERS6','DEGS1'],
    'nk_protrusion': ['EZR','MSN','RDX','ACTR2','ACTR3','ARPC1B','ARPC2','ARPC3','ARPC4','ARPC5',
                      'WAS','WASL','WASF1','WASF2','WASF3','WIPF1','CDC42','RAC1','RHOA',
                      'DIAPH1','DIAPH3','FMNL1','BAIAP2','PACSIN2'],
    'nk_cytotoxicity': ['NKG7','GNLY','GZMB','PRF1','IFNG','LCP2','LAT','VAV1','TLN1','ITGAL','ITGB2'],
    'checkpoint': ['HAVCR2'],
}

# Add metabolic_crosstalk edges: tumor_serine <-> NK SM/topology genes
metabolic_edges = []
tumor_genes = [g for g in SST_MODULES['tumor_serine'] if g in tcga_genes]
nk_targets = [g for g in (SST_MODULES['nk_sm_synthesis'] + SST_MODULES['nk_sm_catabolism'] +
                          SST_MODULES['nk_protrusion'][:10]) if g in tcga_genes]
for tg in tumor_genes:
    for nk in nk_targets:
        metabolic_edges.append({'src': tg, 'dst': nk, 'edge_type': 'metabolic_crosstalk',
                                'weight': 0.5})
log(f"Metabolic crosstalk edges: {len(metabolic_edges)}")

# Add sm_topology_axis edges: within-axis gene pairs
axis_edges = []
axis_all = []
for mod in ['nk_sm_synthesis','nk_sm_catabolism','nk_protrusion','nk_cytotoxicity']:
    axis_all.extend([g for g in SST_MODULES[mod] if g in tcga_genes])
for i, g1 in enumerate(axis_all):
    for g2 in axis_all[i+1:]:
        if np.random.random() < 0.3:  # sparse
            axis_edges.append({'src': g1, 'dst': g2, 'edge_type': 'sm_topology_axis',
                               'weight': 0.3})
log(f"SST axis edges: {len(axis_edges)}")

# Build unified node list
all_graph_genes = set()
for e in string_pairs:
    all_graph_genes.add(e['src'])
    all_graph_genes.add(e['dst'])
for e in metabolic_edges:
    all_graph_genes.add(e['src'])
    all_graph_genes.add(e['dst'])
for e in axis_edges:
    all_graph_genes.add(e['src'])
    all_graph_genes.add(e['dst'])
all_graph_genes = sorted(all_graph_genes)
log(f"Graph genes: {len(all_graph_genes)}")

# Write graph
os.makedirs("data/processed/graph", exist_ok=True)
nodes_df = pd.DataFrame({'node_id': all_graph_genes, 'node_type': 'gene'})
nodes_df.to_csv("data/processed/graph/nodes.tsv", sep="\t", index=False)

all_edges = string_pairs + metabolic_edges + axis_edges
edges_df = pd.DataFrame(all_edges)
edges_df.to_csv("data/processed/graph/edges.tsv", sep="\t", index=False)
log(f"Graph: {len(nodes_df)} nodes, {len(edges_df)} edges")
for et in edges_df['edge_type'].unique():
    log(f"  {et}: {(edges_df['edge_type']==et).sum()}")

# =====================================================================
# STEP 2: Train GNN gene encoder
# =====================================================================
log("\nTraining gene encoder...")
from src.models.gc_nkgraph_atlas import GeneGraphEncoder

encoder = GeneGraphEncoder(
    graph_dir="data/processed/graph",
    embedding_dim=64,
)
encoder.fit()

# Save gene embeddings
gene_embeddings = encoder._embeddings
gene_order = [encoder._idx_to_gene[i] for i in range(len(encoder._idx_to_gene))]
os.makedirs("results/model", exist_ok=True)
np.save("results/model/gene_embeddings.npy", gene_embeddings)
# Save gene list
pd.DataFrame({'gene': gene_order}).to_csv("results/model/gene_order.tsv", sep="\t", index=False)
log(f"Saved gene embeddings: {gene_embeddings.shape}")

# =====================================================================
# STEP 3: Prepare TCGA-STAD labels
# =====================================================================
log("\nPreparing TCGA-STAD labels...")

# Use NK state scoring from scRNA: label each TCGA sample by NK-hot-cytotoxic score
# For this, compute NK cytotoxicity score per bulk sample and threshold
nk_cyto_genes = [g for g in SST_MODULES['nk_cytotoxicity'] if g in stad_expr.columns]
log(f"NK cytotoxicity genes in TCGA: {len(nk_cyto_genes)}/{len(SST_MODULES['nk_cytotoxicity'])}")

# Compute per-sample NK cytotoxicity score (mean z-score)
nk_expr = stad_expr[nk_cyto_genes].values
from scipy import stats
nk_z = stats.zscore(nk_expr, axis=1, nan_policy='omit')
nk_z = np.nan_to_num(nk_z, 0)
nk_score = nk_z.mean(axis=1)

# Binary label: top 34% as NK-hot-cytotoxic, rest as "other" (matching paper's ~34% proportion)
threshold = np.percentile(nk_score, 66.7)
labels = (nk_score >= threshold).astype(int)
log(f"NK-hot-cytotoxic: {labels.sum()} / {len(labels)} ({labels.sum()/len(labels)*100:.1f}%)")

# =====================================================================
# STEP 4: Train NK State Classifier (full graph)
# =====================================================================
log("\nTraining NK state classifier (FULL graph)...")
from src.models.gc_nkgraph_atlas import NKStateClassifier

# Align expression with graph genes
common_genes = [g for g in gene_order if g in stad_expr.columns]
log(f"Genes in both graph and TCGA: {len(common_genes)}/{len(gene_order)}")

gene_to_idx = {g: i for i, g in enumerate(gene_order)}
expr_aligned = np.zeros((stad_expr.shape[0], len(gene_order)), dtype=np.float32)
for i, g in enumerate(gene_order):
    if g in stad_expr.columns:
        expr_aligned[:, i] = stad_expr[g].values.astype(np.float32)

emb_for_classifier = encoder.transform(gene_order)  # shape: (n_graph_genes, d)

# 5-fold CV
from sklearn.model_selection import StratifiedKFold
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

metrics_full = {'fold': [], 'accuracy': [], 'balanced_acc': [], 'macro_f1': [],
                'MCC': [], 'AUROC': [], 'AUPRC': []}
all_test_preds_full = []

def compute_metrics(y_true, y_pred, y_score):
    from sklearn.metrics import (accuracy_score, balanced_accuracy_score,
                                  f1_score, matthews_corrcoef, roc_auc_score,
                                  average_precision_score)
    return {
        'accuracy': accuracy_score(y_true, y_pred),
        'balanced_acc': balanced_accuracy_score(y_true, y_pred),
        'macro_f1': f1_score(y_true, y_pred, average='macro'),
        'MCC': matthews_corrcoef(y_true, y_pred),
        'AUROC': roc_auc_score(y_true, y_score[:, 1]) if y_score.shape[1] > 1 else roc_auc_score(y_true, y_score[:, 0]),
        'AUPRC': average_precision_score(y_true, y_score[:, 1]) if y_score.shape[1] > 1 else average_precision_score(y_true, y_score[:, 0]),
    }

for fold_idx, (train_idx, test_idx) in enumerate(skf.split(expr_aligned, labels)):
    log(f"  Fold {fold_idx+1}/5 FULL (train={len(train_idx)}, test={len(test_idx)})")
    X_train, y_train = expr_aligned[train_idx], labels[train_idx]
    X_test, y_test = expr_aligned[test_idx], labels[test_idx]

    clf = NKStateClassifier(embedding_dim=64, hidden_dims=[128, 64], num_classes=2,
                            dropout=0.3, learning_rate=1e-3)
    clf.fit(X_train, y_train, gene_embeddings=emb_for_classifier,
            X_val_expr=X_test, y_val=y_test, epochs=100, batch_size=32, verbose=False)

    # Evaluate
    import torch
    clf._model.eval()
    X_test_full = clf._prepare_features(X_test, emb_for_classifier)
    X_test_full = clf._standardize(X_test_full, fit=False)
    with torch.no_grad():
        logits = clf._model(torch.from_numpy(X_test_full))
        probs = torch.softmax(logits, dim=1).numpy()
    preds = np.argmax(probs, axis=1)

    m = compute_metrics(y_test, preds, probs)
    for k, v in m.items():
        metrics_full[k].append(v)
    metrics_full['fold'].append(fold_idx + 1)
    log(f"    MCC={m['MCC']:.4f}, AUROC={m['AUROC']:.4f}")

# =====================================================================
# STEP 5: Ablation — train without metabolic_crosstalk
# =====================================================================
log("\nAblation: removing metabolic_crosstalk edges...")
edges_no_met = edges_df[edges_df['edge_type'] != 'metabolic_crosstalk'].copy()
edges_no_met.to_csv("data/processed/graph/edges_ablated.tsv", sep="\t", index=False)

# Temporarily replace edges
import shutil
shutil.copy("data/processed/graph/edges.tsv", "data/processed/graph/edges_full_backup.tsv")
shutil.copy("data/processed/graph/edges_ablated.tsv", "data/processed/graph/edges.tsv")

# Retrain encoder without metabolic edges
encoder_abl = GeneGraphEncoder(graph_dir="data/processed/graph", embedding_dim=64)
encoder_abl.fit()
emb_ablated = encoder_abl.transform(gene_order)

# Train classifier with ablated embeddings
metrics_abl = {'fold': [], 'accuracy': [], 'balanced_acc': [], 'macro_f1': [],
               'MCC': [], 'AUROC': [], 'AUPRC': []}

for fold_idx, (train_idx, test_idx) in enumerate(skf.split(expr_aligned, labels)):
    log(f"  Fold {fold_idx+1}/5 ABLATED (train={len(train_idx)}, test={len(test_idx)})")
    X_train, y_train = expr_aligned[train_idx], labels[train_idx]
    X_test, y_test = expr_aligned[test_idx], labels[test_idx]

    clf = NKStateClassifier(embedding_dim=64, hidden_dims=[128, 64], num_classes=2,
                            dropout=0.3, learning_rate=1e-3)
    clf.fit(X_train, y_train, gene_embeddings=emb_ablated,
            X_val_expr=X_test, y_val=y_test, epochs=100, batch_size=32, verbose=False)

    import torch
    clf._model.eval()
    X_test_full = clf._prepare_features(X_test, emb_ablated)
    X_test_full = clf._standardize(X_test_full, fit=False)
    with torch.no_grad():
        logits = clf._model(torch.from_numpy(X_test_full))
        probs = torch.softmax(logits, dim=1).numpy()
    preds = np.argmax(probs, axis=1)

    m = compute_metrics(y_test, preds, probs)
    for k, v in m.items():
        metrics_abl[k].append(v)
    metrics_abl['fold'].append(fold_idx + 1)
    log(f"    MCC={m['MCC']:.4f}, AUROC={m['AUROC']:.4f}")

# Restore full edges
shutil.copy("data/processed/graph/edges_full_backup.tsv", "data/processed/graph/edges.tsv")

# =====================================================================
# STEP 6: Statistical comparison + save
# =====================================================================
log("\nStatistical comparison...")
from scipy.stats import ttest_rel, wilcoxon

results = []
for metric in ['MCC', 'AUROC']:
    f = metrics_full[metric]
    a = metrics_abl[metric]
    delta = np.mean(f) - np.mean(a)
    t_p = ttest_rel(f, a).pvalue
    try:
        w_p = wilcoxon(f, a).pvalue
    except:
        w_p = np.nan
    results.append({
        'metric': metric,
        'mean_full': np.mean(f), 'std_full': np.std(f),
        'mean_ablated': np.mean(a), 'std_ablated': np.std(a),
        'delta': delta, 'ttest_p': t_p, 'wilcoxon_p': w_p,
    })
    log(f"  {metric}: full={np.mean(f):.4f}, ablated={np.mean(a):.4f}, "
        f"delta={delta:+.4f}, paired t p={t_p:.4f}")

results_df = pd.DataFrame(results)
os.makedirs("results/tables", exist_ok=True)
results_df.to_csv("results/tables/ablation_results.tsv", sep="\t", index=False)

# Save full fold metrics
full_df = pd.DataFrame(metrics_full)
full_df['mode'] = 'full'
abl_df = pd.DataFrame(metrics_abl)
abl_df['mode'] = 'no_metabolic'
pd.concat([full_df, abl_df]).to_csv("results/tables/ablation_per_fold.tsv", sep="\t", index=False)

# Save classifier metrics summary
full_summary = {k: (np.mean(v), np.std(v)) for k, v in metrics_full.items() if k != 'fold'}
log(f"\nFull GNN summary:")
for k, (mn, sd) in full_summary.items():
    log(f"  {k}: {mn:.4f} ± {sd:.4f}")

log("\nT11+T12 DONE!")
