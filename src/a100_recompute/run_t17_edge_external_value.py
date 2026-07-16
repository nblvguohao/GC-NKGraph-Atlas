"""
T17 — External-sample value of the metabolic_crosstalk edge (held-out cohort test).

Risk (R2): The current ablation (§3.7) is near-tautological — removing
mechanism-defined edges naturally removes mechanism-defined structure. It proves
"the edge is in the graph," not "the edge encodes independently verifiable biology."

REVISION NOTE (post-hoc): the originally reported result (FULL MCC=0.416 vs
-MC MCC=0.509, bootstrap p=0.003, "measurably degrades") was produced by (a) a
single training seed, and (b) a last-write-wins bug in build_adj_matrix where,
for node pairs connected by more than one edge type, whichever edge type's row
came later in edges.tsv silently overwrote the weight of an earlier one instead
of the two being combined -- so the FULL/-MC adjacency matrices did not purely
differ by the presence/absence of metabolic_crosstalk edges. Both issues are
fixed here: build_adj_matrix now takes max(weight) across colliding edge types,
and this script runs SEEDS (10 seeds) instead of one, using a paired test across
seeds rather than a single bootstrap over the fixed LIHC eval set (which only
captures resampling uncertainty, not the seed-to-seed training variance that
turned out to dominate -- see per-seed spread below).

Test design (held-out cohort prediction):
  1. Build two graph variants: FULL (all edges incl. metabolic_crosstalk) and -MC
     (minus metabolic_crosstalk edges only).
  2. For each variant, for each of SEEDS:
     a. Compute gene embeddings via SVD on the heterogeneous graph (deterministic,
        shared across seeds).
     b. Train NK-state classifier on TCGA-STAD (gastric).
     c. Evaluate on TCGA-LIHC (held-out, liver — entirely different tissue/organ).
  3. Paired t-test and Wilcoxon signed-rank test across seeds on the MCC
     difference (FULL - -MC).

Criterion:
  - If FULL > -MC in cross-cohort MCC with p < 0.05 (either test) → the
    metabolic_crosstalk edge provides externally measurable predictive value.
  - If FULL < -MC with p < 0.05 → the edge measurably degrades cross-cohort
    prediction.
  - Otherwise → inconclusive; §3.7 and §4 should describe the edge as shaping
    the embedding without a demonstrated external predictive effect either way.

Output: results/tables/t17_edge_external_value.tsv (aggregate),
        results/tables/t17_edge_external_value_by_seed.tsv (per-seed)

Run:  conda activate gc-nkgraph && python src/a100_recompute/run_t17_edge_external_value.py
"""
import sys, os, copy, time, json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
from sklearn.metrics import matthews_corrcoef, roc_auc_score
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.io_utils import ensure_dir, load_table, load_config
from src.common.seed import set_seed
from src.common.sst_config import load_sst_modules

set_seed(42)

GRAPH_DIR = "data/processed/graph"
RESULTS = "results/tables"
OUT_PATH = f"{RESULTS}/t17_edge_external_value.tsv"
OUT_PATH_SEEDS = f"{RESULTS}/t17_edge_external_value_by_seed.tsv"

def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# =============================================================================
# 1. Build two graph variants in-memory
# =============================================================================

def build_adj_matrix(edges_df, nodes_df, edge_types_to_keep=None, edge_types_to_drop=None):
    """Build combined adjacency matrix from edge list, filtering by edge type.

    Args:
        edges_df: DataFrame with src, dst, edge_type, weight columns.
        nodes_df: DataFrame with node_id column.
        edge_types_to_keep: If set, only include these edge types.
        edge_types_to_drop: If set, exclude these edge types.

    Returns: (adj_combined, node_to_idx, idx_to_node)

    Note: when multiple edge types connect the same node pair, the pair's
    weight is the max across those edge types (not last-write-wins by row
    order) so that e.g. a low-weight generic-prior edge never silently
    overwrites a higher-weight mechanism/PPI edge on the same pair.
    """
    node_ids = nodes_df["node_id"].astype(str).tolist()
    n = len(node_ids)
    node_to_idx = {nid: i for i, nid in enumerate(node_ids)}
    idx_to_node = {i: nid for i, nid in enumerate(node_ids)}

    adj_combined = np.zeros((n, n), dtype=np.float32)

    for _, row in edges_df.iterrows():
        etype = str(row["edge_type"])

        # Apply filters
        if edge_types_to_keep is not None and etype not in edge_types_to_keep:
            continue
        if edge_types_to_drop is not None and etype in edge_types_to_drop:
            continue

        src, dst = str(row["src"]), str(row["dst"])
        if src in node_to_idx and dst in node_to_idx:
            i, j = node_to_idx[src], node_to_idx[dst]
            w = float(row.get("weight", 1.0))
            if w > adj_combined[i, j]:
                adj_combined[i, j] = w
                adj_combined[j, i] = w  # symmetric

    return adj_combined, node_to_idx, idx_to_node


def compute_svd_embeddings(adj, embedding_dim=64):
    """Compute SVD embeddings from adjacency (same method as GeneGraphEncoder)."""
    # Symmetric normalized Laplacian
    deg = adj.sum(axis=1)
    deg = np.where(deg > 0, deg, 1.0)
    deg_inv_sqrt = np.diag(1.0 / np.sqrt(deg))
    adj_norm = deg_inv_sqrt @ adj @ deg_inv_sqrt

    from scipy.sparse.linalg import svds
    k = min(embedding_dim, adj.shape[0] - 2)
    if k < 2:
        # Fallback: random
        rng = np.random.RandomState(42)
        return rng.randn(adj.shape[0], embedding_dim).astype(np.float32)

    u, s, vt = svds(adj_norm.astype(np.float64), k=k)
    idx = np.argsort(s)[::-1]
    embeddings = u[:, idx] @ np.diag(np.sqrt(np.maximum(s[idx], 0)))

    # Pad/truncate
    if embeddings.shape[1] < embedding_dim:
        pad = np.zeros((embeddings.shape[0], embedding_dim - embeddings.shape[1]))
        embeddings = np.hstack([embeddings, pad])
    else:
        embeddings = embeddings[:, :embedding_dim]
    return embeddings.astype(np.float32)


# =============================================================================
# 2. Load data
# =============================================================================

def load_cohort_expression(labels_df, dataset_name):
    """Load expression matrix for a specific dataset cohort."""
    # Find the expression file
    expr_dir = "data/processed/bulk"

    # Map dataset names in labels to expression file prefixes
    cohort_samples = labels_df[labels_df["dataset"] == dataset_name]
    if len(cohort_samples) == 0:
        return None, None

    # Determine expression file
    file_map = {
        "TCGA-LIHC": "tcga_lihc_expression.tsv",
        "TCGA-STAD": "tcga_stad_expression.tsv",
        "GSE62254": "gse62254_expression.tsv",
        "GSE84437": "gse84437_expression.tsv",
    }

    fname = file_map.get(dataset_name)
    if fname is None:
        log(f"  No expression file mapping for {dataset_name}")
        return None, None

    fpath = os.path.join(expr_dir, fname)
    if not os.path.exists(fpath):
        log(f"  Expression file not found: {fpath}")
        return None, None

    log(f"  Loading {fpath}...")
    expr = load_table(fpath)

    # Align samples
    common = expr.index.intersection(cohort_samples.index)
    expr = expr.loc[common]
    labels_sub = cohort_samples.loc[common]

    # Binary label: NK-hot-cytotoxic = 1, else = 0
    y = (labels_sub["nk_immune_state"] == "NK-hot-cytotoxic").astype(int).values

    log(f"    {len(common)} samples, {y.sum()} NK-hot-cytotoxic")
    return expr, y


def prepare_features(expr_df, gene_embeddings, node_to_idx):
    """Project gene expression through graph embeddings (same as NKStateClassifier)."""
    # Get overlapping genes
    genes_in_graph = [g for g in expr_df.columns if g in node_to_idx]
    if len(genes_in_graph) == 0:
        # Fallback: use all genes that intersect
        genes_in_graph = list(expr_df.columns)

    log(f"    {len(genes_in_graph)} genes in graph overlap")

    X_expr = expr_df[genes_in_graph].values.astype(np.float32)

    # Get embeddings for these genes
    gene_indices = [node_to_idx[g] for g in genes_in_graph]
    E = gene_embeddings[gene_indices]  # (n_genes, embedding_dim)

    # Graph projection
    graph_proj = X_expr @ E
    graph_proj = graph_proj / max(np.std(graph_proj), 1e-8)

    # Concatenate raw expression + graph projection
    X = np.hstack([X_expr, graph_proj]).astype(np.float32)

    # Standardize
    X_mean = X.mean(axis=0)
    X_std = X.std(axis=0) + 1e-8
    X = (X - X_mean) / X_std

    return X, X_mean, X_std


def standardize_with_params(X, mean, std):
    """Apply pre-computed standardization."""
    return ((X - mean) / (std + 1e-8)).astype(np.float32)


# =============================================================================
# 3. Classifier (lightweight — mirrors NKStateClassifier core)
# =============================================================================

def train_and_evaluate_cross_cohort(
    X_train, y_train, X_test, y_test,
    hidden_dims=(128, 64), epochs=200, batch_size=32, lr=1e-3, weight_decay=1e-5
):
    """Train MLP on train cohort, evaluate on held-out test cohort.

    Returns: dict with train_mcc, train_auroc, test_mcc, test_auroc, test_preds, test_probs
    """
    import torch
    import torch.nn as nn

    device = torch.device("cpu")

    input_dim = X_train.shape[1]
    n_classes = 2

    # Build model
    dims = [input_dim] + list(hidden_dims) + [n_classes]
    layers = []
    for i in range(len(dims) - 1):
        layers.append(nn.Linear(dims[i], dims[i + 1]))
        if i < len(dims) - 2:
            layers.append(nn.BatchNorm1d(dims[i + 1]))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(0.5))
    model = nn.Sequential(*layers).to(device)

    # Loss with mild class weighting
    n_pos = max((y_train == 1).sum(), 1)
    n_neg = max((y_train == 0).sum(), 1)
    total = n_pos + n_neg
    class_weights = torch.tensor([total / (2.0 * n_neg), total / (2.0 * n_pos)], dtype=torch.float32)
    criterion = nn.CrossEntropyLoss(weight=class_weights)

    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)

    X_train_t = torch.from_numpy(X_train).float()
    y_train_t = torch.from_numpy(y_train.astype(np.int64))

    # Train
    model.train()
    n_samples = X_train.shape[0]

    for epoch in range(epochs):
        # Mini-batch SGD
        perm = torch.randperm(n_samples)
        total_loss = 0.0
        n_batches = 0

        for start in range(0, n_samples, batch_size):
            idx = perm[start:start + batch_size]
            xb = X_train_t[idx].to(device)
            yb = y_train_t[idx].to(device)

            optimizer.zero_grad()
            logits = model(xb)
            loss = criterion(logits, yb)
            loss.backward()
            optimizer.step()

            total_loss += loss.item()
            n_batches += 1

        if (epoch + 1) % 50 == 0:
            log(f"    epoch {epoch+1}/{epochs}  loss={total_loss/max(n_batches,1):.4f}")

    # Evaluate
    model.eval()
    with torch.no_grad():
        # Train metrics
        train_logits = model(X_train_t.to(device)).numpy()
        train_probs = torch.softmax(torch.from_numpy(train_logits), dim=1).numpy()[:, 1]
        train_preds = train_logits.argmax(axis=1)
        train_mcc = matthews_corrcoef(y_train, train_preds)
        train_auroc = roc_auc_score(y_train, train_probs)

        # Test metrics
        X_test_t = torch.from_numpy(X_test).float().to(device)
        test_logits = model(X_test_t).numpy()
        test_probs = torch.softmax(torch.from_numpy(test_logits), dim=1).numpy()[:, 1]
        test_preds = test_logits.argmax(axis=1)
        test_mcc = matthews_corrcoef(y_test, test_preds)
        test_auroc = roc_auc_score(y_test, test_probs)

    return {
        "train_mcc": train_mcc, "train_auroc": train_auroc,
        "test_mcc": test_mcc, "test_auroc": test_auroc,
        "test_preds": test_preds, "test_probs": test_probs,
    }


# =============================================================================
# 4. Bootstrap test
# =============================================================================

def bootstrap_mcc_difference(y_true, preds_full, probs_full, preds_mc, probs_mc,
                              n_bootstrap=1000, seed=42):
    """Bootstrap the difference in MCC (FULL - -MC) on held-out predictions.

    Returns: observed_diff, bootstrap_diffs, p_value, ci_95
    """
    rng = np.random.RandomState(seed)
    n = len(y_true)

    mcc_full_obs = matthews_corrcoef(y_true, preds_full)
    mcc_mc_obs = matthews_corrcoef(y_true, preds_mc)
    obs_diff = mcc_full_obs - mcc_mc_obs

    bootstrap_diffs = []
    for _ in range(n_bootstrap):
        idx = rng.choice(n, size=n, replace=True)
        yb = y_true[idx]

        # Predictions are deterministic given the model — bootstrap the evaluation
        mcc_f = matthews_corrcoef(yb, preds_full[idx])
        mcc_m = matthews_corrcoef(yb, preds_mc[idx])
        bootstrap_diffs.append(mcc_f - mcc_m)

    bootstrap_diffs = np.array(bootstrap_diffs)

    # Two-sided p-value
    p_value = np.mean(np.abs(bootstrap_diffs - obs_diff) >= np.abs(obs_diff))
    # One-sided (FULL > -MC)
    p_one_sided = np.mean(bootstrap_diffs <= 0)  # H0: diff <= 0

    ci_95 = np.percentile(bootstrap_diffs, [2.5, 97.5])

    return {
        "mcc_full": mcc_full_obs,
        "mcc_minus_mc": mcc_mc_obs,
        "mcc_diff": obs_diff,
        "bootstrap_mean_diff": float(bootstrap_diffs.mean()),
        "bootstrap_std_diff": float(bootstrap_diffs.std()),
        "p_two_sided": float(p_value),
        "p_one_sided_full_gt": float(p_one_sided),
        "ci_95_low": float(ci_95[0]),
        "ci_95_high": float(ci_95[1]),
    }


# =============================================================================
# 5. Main
# =============================================================================

SEEDS = [1234, 2345, 3456, 4567, 5678, 6789, 7890, 8901, 9012, 1357]


def main():
    log("=" * 60)
    log("T17 — metabolic_crosstalk EDGE EXTERNAL-SAMPLE VALUE TEST")
    log(f"(multi-seed, n={len(SEEDS)}; revised after fixing a last-write-wins")
    log(" collision bug in build_adj_matrix that a single-seed run had masked)")
    log("=" * 60)

    # ---- Load graph ----
    log("\n[1] Loading graph...")
    nodes = pd.read_csv(f"{GRAPH_DIR}/nodes.tsv", sep="\t")
    edges = pd.read_csv(f"{GRAPH_DIR}/edges.tsv", sep="\t")
    log(f"  Nodes: {len(nodes)}, Edges: {len(edges)}")
    for etype, n in edges["edge_type"].value_counts().items():
        log(f"    {etype}: {n}")

    # ---- Load labels ----
    log("\n[2] Loading NK state labels...")
    labels = pd.read_csv("results/tables/nk_state_labels.tsv", sep="\t", index_col=0, comment="#")
    # Skip comment/provenance rows that may have been ingested
    labels = labels[labels.index.notnull()]
    log(f"  Datasets: {labels['dataset'].value_counts().to_dict()}")

    # ---- Build two graph variants ----
    log("\n[3] Building graph variants...")

    # FULL: all edges
    adj_full, node_to_idx, idx_to_node = build_adj_matrix(
        edges, nodes, edge_types_to_drop=set()
    )
    log(f"  FULL adj: {adj_full.shape}, {int((adj_full > 0).sum()/2)} edges")

    # -MC: drop metabolic_crosstalk only
    adj_mc, _, _ = build_adj_matrix(
        edges, nodes, edge_types_to_drop={"metabolic_crosstalk"}
    )
    log(f"  -MC adj:  {adj_mc.shape}, {int((adj_mc > 0).sum()/2)} edges")
    log(f"  edges removed: {int((adj_full > 0).sum()/2) - int((adj_mc > 0).sum()/2)}")

    # ---- Compute embeddings (deterministic given adjacency; shared across seeds) ----
    log("\n[4] Computing SVD gene embeddings...")
    emb_dim = 64

    emb_full = compute_svd_embeddings(adj_full, emb_dim)
    log(f"  FULL embeddings: {emb_full.shape}")

    emb_mc = compute_svd_embeddings(adj_mc, emb_dim)
    log(f"  -MC embeddings:  {emb_mc.shape}")

    # ---- Prepare train (STAD — gastric, train_primary) and test (LIHC — liver, positive control) ----
    # This is cross-cancer transfer: gastric → liver. The SST axis was originally
    # discovered in liver cancer, so predicting liver NK states from a gastric-trained
    # model with the axis-grounded edges is a genuine external test.
    log("\n[5] Preparing train (STAD → gastric) and test (LIHC → liver) data...")

    expr_train, y_train_data = load_cohort_expression(labels, "TCGA-STAD")
    expr_test, y_test_data = load_cohort_expression(labels, "TCGA-LIHC")

    if expr_train is None or expr_test is None:
        log("ERROR: Missing cohort expression data")
        return

    log(f"  STAD train: {expr_train.shape[0]} samples, {y_train_data.sum()} pos")
    log(f"  LIHC test:  {expr_test.shape[0]} samples, {y_test_data.sum()} pos")

    # ---- Train & evaluate for each variant, across seeds ----
    log(f"\n[6] Training and evaluating across {len(SEEDS)} seeds...")

    per_seed_rows = []
    results_by_seed = {"FULL": [], "minus_MC": []}

    for seed in SEEDS:
      for variant_name, adj, emb in [
        ("FULL", adj_full, emb_full),
        ("minus_MC", adj_mc, emb_mc),
      ]:
        set_seed(seed)
        log(f"\n  --- seed={seed}  {variant_name} ---")

        # Prepare features
        X_train_full, mean_tr, std_tr = prepare_features(expr_train, emb, node_to_idx)
        X_test_full = prepare_features_raw(expr_test, emb, node_to_idx)
        X_test_full_s = standardize_with_params(X_test_full, mean_tr, std_tr)

        # Train and evaluate (cross-cohort)
        result = train_and_evaluate_cross_cohort(
            X_train_full, y_train_data, X_test_full_s, y_test_data,
            hidden_dims=(128, 64), epochs=200, batch_size=32, lr=1e-3
        )

        log(f"    Train MCC={result['train_mcc']:.4f}  AUROC={result['train_auroc']:.4f}")
        log(f"    Test  MCC={result['test_mcc']:.4f}  AUROC={result['test_auroc']:.4f}")

        results_by_seed[variant_name].append(result)
        per_seed_rows.append({
            "seed": seed, "variant": variant_name,
            "train_mcc": round(result["train_mcc"], 4),
            "test_mcc": round(result["test_mcc"], 4),
            "test_auroc": round(result["test_auroc"], 4),
        })

    ensure_dir(RESULTS)
    pd.DataFrame(per_seed_rows).to_csv(OUT_PATH_SEEDS, sep="\t", index=False)
    log(f"\n[7] Per-seed results written to {OUT_PATH_SEEDS}")

    # ---- Paired test across seeds (same seed -> same train/test split & init) ----
    log(f"\n[8] Paired comparison across {len(SEEDS)} seeds (FULL vs -MC)...")

    full_mccs = np.array([r["test_mcc"] for r in results_by_seed["FULL"]])
    mc_mccs = np.array([r["test_mcc"] for r in results_by_seed["minus_MC"]])
    diff = full_mccs - mc_mccs

    from scipy.stats import wilcoxon, ttest_rel
    try:
        _, p_wilcoxon = wilcoxon(full_mccs, mc_mccs)
    except ValueError:
        p_wilcoxon = float("nan")
    _, p_ttest = ttest_rel(full_mccs, mc_mccs)

    log(f"  FULL  test_mcc: {full_mccs.mean():.4f} +/- {full_mccs.std():.4f}  (per-seed: {np.round(full_mccs, 4).tolist()})")
    log(f"  -MC   test_mcc: {mc_mccs.mean():.4f} +/- {mc_mccs.std():.4f}  (per-seed: {np.round(mc_mccs, 4).tolist()})")
    log(f"  mean diff (FULL - -MC): {diff.mean():.4f} +/- {diff.std():.4f}")
    log(f"  paired t-test p:  {p_ttest:.4f}")
    log(f"  Wilcoxon p:       {p_wilcoxon:.4f}")

    # ---- Verdict (n=10 seeds; treat as indicative per this repo's own small-n convention) ----
    significant = (p_ttest == p_ttest and p_ttest < 0.05) or (p_wilcoxon == p_wilcoxon and p_wilcoxon < 0.05)
    if significant and diff.mean() > 0:
        verdict = "PASS — edge has external predictive value"
    elif significant and diff.mean() < 0:
        verdict = "FAIL — edge measurably degrades cross-cohort prediction"
    else:
        verdict = "INCONCLUSIVE — no significant difference across seeds"
    log(f"\n  VERDICT: {verdict}")

    # ---- Write aggregate results ----
    log(f"\n[9] Writing aggregate results to {OUT_PATH}...")

    out_rows = [
        {"metric": "variant", "FULL": "all_edges", "minus_MC": "no_metabolic_crosstalk", "delta": "", "p_value": "", "verdict": ""},
        {"metric": "n_seeds", "FULL": len(SEEDS), "minus_MC": len(SEEDS), "delta": "", "p_value": "", "verdict": ""},
        {"metric": "n_nodes", "FULL": len(nodes), "minus_MC": len(nodes), "delta": "", "p_value": "", "verdict": ""},
        {"metric": "n_edges", "FULL": int((adj_full > 0).sum() / 2), "minus_MC": int((adj_mc > 0).sum() / 2), "delta": int((adj_full > 0).sum() / 2) - int((adj_mc > 0).sum() / 2), "p_value": "", "verdict": ""},
        {"metric": "stad_n_samples", "FULL": expr_train.shape[0], "minus_MC": expr_train.shape[0], "delta": "", "p_value": "", "verdict": ""},
        {"metric": "lihc_n_samples", "FULL": expr_test.shape[0], "minus_MC": expr_test.shape[0], "delta": "", "p_value": "", "verdict": ""},
        {"metric": "test_mcc_mean", "FULL": round(full_mccs.mean(), 4), "minus_MC": round(mc_mccs.mean(), 4), "delta": round(diff.mean(), 4), "p_value": round(p_ttest, 4), "verdict": verdict},
        {"metric": "test_mcc_std", "FULL": round(full_mccs.std(), 4), "minus_MC": round(mc_mccs.std(), 4), "delta": "", "p_value": "", "verdict": ""},
        {"metric": "test_mcc_wilcoxon_p", "FULL": "", "minus_MC": "", "delta": "", "p_value": round(p_wilcoxon, 4) if p_wilcoxon == p_wilcoxon else "NA", "verdict": ""},
    ]

    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(OUT_PATH, sep="\t", index=False)
    log(f"  Written {len(out_df)} rows to {OUT_PATH}")

    log("\n" + "=" * 60)
    log("T17 COMPLETE (revised: multi-seed, collision-bug-fixed adjacency)")
    log("=" * 60)

    return verdict


def prepare_features_raw(expr_df, gene_embeddings, node_to_idx):
    """Prepare features without standardization (for test set)."""
    genes_in_graph = [g for g in expr_df.columns if g in node_to_idx]
    if len(genes_in_graph) == 0:
        genes_in_graph = list(expr_df.columns)

    X_expr = expr_df[genes_in_graph].values.astype(np.float32)
    gene_indices = [node_to_idx[g] for g in genes_in_graph]
    E = gene_embeddings[gene_indices]

    graph_proj = X_expr @ E
    graph_proj = graph_proj / max(np.std(graph_proj), 1e-8)

    X = np.hstack([X_expr, graph_proj]).astype(np.float32)
    return X


if __name__ == "__main__":
    main()
