"""Contracts for the lightweight TREE/GRAFT-inspired multiview analysis."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import torch


def _nodes() -> pd.DataFrame:
    return pd.DataFrame({"node_id": ["PHGDH", "SGMS1", "EZR", "NKG7", "KLRK1"]})


def _edges() -> pd.DataFrame:
    return pd.DataFrame(
        [
            ("PHGDH", "SGMS1", "ppi", 0.8),
            ("EZR", "NKG7", "coexpression", 0.4),
            ("NKG7", "KLRK1", "ligand_receptor", 0.9),
            ("PHGDH", "SGMS1", "metabolic_crosstalk", 0.5),
            ("EZR", "NKG7", "sm_topology_axis", 0.3),
            ("PHGDH", "EZR", "go_prior", 0.2),
            ("SGMS1", "KLRK1", "msigdb_prior", 0.2),
        ],
        columns=["src", "dst", "edge_type", "weight"],
    )


def test_default_views_partition_each_supported_edge_once():
    from src.models.multiview_fusion import DEFAULT_VIEW_SPEC, build_graph_views

    views = build_graph_views(_edges(), _nodes(), DEFAULT_VIEW_SPEC)

    assert list(views) == [
        "generic_interaction",
        "ligand_receptor",
        "metabolic_crosstalk",
        "sm_topology_axis",
        "go_prior",
        "msigdb_prior",
    ]
    assert views["generic_interaction"].n_edges == 2
    assert views["metabolic_crosstalk"].n_edges == 1
    assert views["sm_topology_axis"].n_edges == 1
    assert sum(view.n_edges for view in views.values()) == len(_edges())


def test_missing_go_or_msigdb_view_fails_loudly():
    from src.models.multiview_fusion import DEFAULT_VIEW_SPEC, build_graph_views

    edges = _edges().query("edge_type != 'go_prior'")

    with pytest.raises(ValueError, match="go_prior"):
        build_graph_views(edges, _nodes(), DEFAULT_VIEW_SPEC, require_nonempty=True)


def test_unknown_edge_type_is_rejected_in_strict_mode():
    from src.models.multiview_fusion import DEFAULT_VIEW_SPEC, build_graph_views

    edges = pd.concat(
        [_edges(), pd.DataFrame([("PHGDH", "NKG7", "unknown", 1.0)], columns=_edges().columns)],
        ignore_index=True,
    )

    with pytest.raises(ValueError, match="Unassigned edge types"):
        build_graph_views(edges, _nodes(), DEFAULT_VIEW_SPEC, strict=True)


def test_label_defining_genes_is_the_union_of_all_four_signatures():
    from src.models.multiview_fusion import label_defining_genes

    masked = label_defining_genes()

    assert {"NKG7", "GNLY", "TIGIT", "NT5E", "COL1A1", "FAP", "NCR1"} <= masked
    assert "PHGDH" not in masked


def test_masked_expression_removes_label_genes_but_full_preserves_them():
    from src.models.multiview_fusion import filter_expression_features

    expression = pd.DataFrame(
        {"NKG7": [1.0, 2.0], "TIGIT": [2.0, 1.0], "PHGDH": [3.0, 4.0]},
        index=["s1", "s2"],
    )

    masked, removed = filter_expression_features(expression, mode="masked")
    full, full_removed = filter_expression_features(expression, mode="full")

    assert list(masked.columns) == ["PHGDH"]
    assert removed == {"NKG7", "TIGIT"}
    pd.testing.assert_frame_equal(full, expression)
    assert full_removed == set()


def test_spectral_embedding_is_deterministic_and_has_requested_shape():
    from src.models.multiview_fusion import DEFAULT_VIEW_SPEC, build_graph_views, compute_spectral_embedding

    view = build_graph_views(_edges(), _nodes(), DEFAULT_VIEW_SPEC)["generic_interaction"]
    first = compute_spectral_embedding(view, embedding_dim=4)
    second = compute_spectral_embedding(view, embedding_dim=4)

    assert first.shape == (5, 4)
    assert first.dtype == np.float32
    np.testing.assert_allclose(first, second)


def test_uniform_fusion_is_fixed_and_learned_fusion_normalizes_weights():
    from src.models.multiview_fusion import MultiViewFusionClassifier

    projections = torch.tensor(
        [
            [[1.0, 0.0], [0.0, 1.0]],
            [[3.0, 2.0], [2.0, 3.0]],
        ]
    )
    uniform = MultiViewFusionClassifier(expr_dim=1, emb_dim=2, n_views=2, mode="uniform")
    learned = MultiViewFusionClassifier(expr_dim=1, emb_dim=2, n_views=2, mode="learned")

    np.testing.assert_allclose(uniform.view_weights().detach().numpy(), [0.5, 0.5])
    np.testing.assert_allclose(uniform.fused_projection(projections).detach().numpy(), [[2, 1], [1, 2]])
    assert uniform.view_logits.requires_grad is False
    assert learned.view_logits.requires_grad is True
    assert float(learned.view_weights().sum().detach()) == pytest.approx(1.0)


def test_train_standardizer_never_uses_external_values():
    from src.models.multiview_fusion import fit_standardizer, transform_standardized

    train = np.array([[0.0], [2.0]], dtype=np.float32)
    external = np.array([[1000.0]], dtype=np.float32)
    mean, scale = fit_standardizer(train)

    assert mean.tolist() == [1.0]
    assert scale.tolist() == [1.0]
    assert transform_standardized(external, mean, scale).item() == pytest.approx(999.0)


def test_paired_stratified_bootstrap_detects_positive_and_null_differences():
    from src.baselines.run_multiview_fusion_benchmark import paired_stratified_bootstrap

    y = np.array([0] * 20 + [1] * 20)
    strong = np.r_[np.linspace(0.01, 0.20, 20), np.linspace(0.80, 0.99, 20)]
    weak = np.r_[np.linspace(0.30, 0.65, 20), np.linspace(0.35, 0.70, 20)]
    positive = paired_stratified_bootstrap(y, strong, weak, metric="auroc", n_bootstrap=300, seed=7)
    null = paired_stratified_bootstrap(y, strong, strong, metric="auroc", n_bootstrap=100, seed=7)

    assert positive["observed_delta"] > 0
    assert positive["ci_low"] > 0
    assert null["observed_delta"] == pytest.approx(0.0)
    assert null["ci_low"] == pytest.approx(0.0)
    assert null["ci_high"] == pytest.approx(0.0)


def test_stable_gain_verdict_requires_both_metrics_against_all_comparators():
    from src.baselines.run_multiview_fusion_benchmark import determine_external_gain_verdict

    comparisons = pd.DataFrame(
        {
            "comparator": ["no_graph", "merged_svd", "uniform_multiview"] * 2,
            "metric": ["AUROC"] * 3 + ["AUPRC"] * 3,
            "observed_delta": [0.02] * 6,
            "ci_low": [0.01] * 6,
            "ci_high": [0.03] * 6,
        }
    )

    assert determine_external_gain_verdict(comparisons) == "stable_external_gain"
    comparisons.loc[5, "ci_low"] = -0.001
    assert determine_external_gain_verdict(comparisons) == "no_stable_external_gain"


def test_model_matrix_contains_all_pre_registered_variants():
    from src.baselines.run_multiview_fusion_benchmark import model_variant_specs

    view_names = list(
        [
            "generic_interaction",
            "ligand_receptor",
            "metabolic_crosstalk",
            "sm_topology_axis",
            "go_prior",
            "msigdb_prior",
        ]
    )
    specs = model_variant_specs(view_names)

    assert set(specs) >= {"no_graph", "merged_svd", "uniform_multiview", "learned_multiview"}
    assert {f"single__{name}" for name in view_names} <= set(specs)
    assert {f"leave_out__{name}" for name in view_names} <= set(specs)
    assert specs["uniform_multiview"].fusion_mode == "uniform"
    assert specs["learned_multiview"].fusion_mode == "learned"
    assert specs["leave_out__metabolic_crosstalk"].views == tuple(
        name for name in view_names if name != "metabolic_crosstalk"
    )


def test_result_coverage_rejects_missing_mode_variant_or_seed():
    from src.baselines.run_multiview_fusion_benchmark import validate_result_coverage

    rows = pd.DataFrame(
        [
            {"feature_mode": mode, "variant": variant, "seed": seed, "AUROC": 0.7}
            for mode in ("masked", "full")
            for variant in ("no_graph", "learned_multiview")
            for seed in (1, 2)
        ]
    )
    validate_result_coverage(
        rows,
        feature_modes=("masked", "full"),
        variants=("no_graph", "learned_multiview"),
        seeds=(1, 2),
    )

    with pytest.raises(ValueError, match="Missing benchmark rows"):
        validate_result_coverage(
            rows.iloc[:-1],
            feature_modes=("masked", "full"),
            variants=("no_graph", "learned_multiview"),
            seeds=(1, 2),
        )


def test_projection_uses_only_expression_columns_that_survive_masking():
    from src.models.multiview_fusion import raw_graph_projection

    expression = pd.DataFrame({"PHGDH": [2.0], "NKG7": [100.0]})
    embedding = np.array([[3.0, 0.0], [0.0, 7.0]], dtype=np.float32)
    node_to_idx = {"PHGDH": 0, "NKG7": 1}

    projected, used_genes = raw_graph_projection(expression[["PHGDH"]], embedding, node_to_idx)

    assert used_genes == ("PHGDH",)
    np.testing.assert_allclose(projected, [[6.0, 0.0]])


def test_degree_sequence_randomization_preserves_degrees_and_weights():
    from src.models.multiview_fusion import GraphView, permute_view_node_labels

    adjacency = np.array(
        [
            [0.0, 0.2, 0.0, 0.7],
            [0.2, 0.0, 0.4, 0.0],
            [0.0, 0.4, 0.0, 0.9],
            [0.7, 0.0, 0.9, 0.0],
        ],
        dtype=np.float32,
    )
    view = GraphView("mechanism", ("mechanism",), adjacency, ("A", "B", "C", "D"), 4)
    randomized = permute_view_node_labels(view, np.random.RandomState(3))

    np.testing.assert_array_equal(
        np.sort((view.adjacency > 0).sum(axis=1)),
        np.sort((randomized.adjacency > 0).sum(axis=1)),
    )
    np.testing.assert_allclose(
        np.sort(view.adjacency[np.triu_indices(4, 1)]),
        np.sort(randomized.adjacency[np.triu_indices(4, 1)]),
    )
    assert randomized.node_ids == view.node_ids


def test_module_coupling_calibration_exports_pre_registered_null_contract():
    from src.models.multiview_fusion import GraphView, calibrate_module_coupling

    adjacency = np.array(
        [
            [0.0, 1.0, 0.0, 0.0],
            [1.0, 0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 0.0],
        ],
        dtype=np.float32,
    )
    view = GraphView("mechanism", ("mechanism",), adjacency, ("A", "B", "C", "D"), 3)
    result, null_values = calibrate_module_coupling(
        view,
        source_genes={"A", "B"},
        target_genes={"C", "D"},
        embedding_dim=3,
        n_randomizations=20,
        seed=9,
    )

    assert result["view"] == "mechanism"
    assert result["n_randomizations"] == 20
    assert 0 < result["empirical_p"] <= 1
    assert len(null_values) == 20


def test_contextual_calibration_randomizes_only_the_target_view():
    from src.models.multiview_fusion import GraphView, calibrate_view_in_context

    nodes = ("A", "B", "C", "D")
    generic = GraphView(
        "generic",
        ("ppi",),
        np.array(
            [[0, 1, 0, 0], [1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 1, 0]],
            dtype=np.float32,
        ),
        nodes,
        3,
    )
    mechanism = GraphView(
        "mechanism",
        ("mechanism",),
        np.array(
            [[0, 0, 1, 1], [0, 0, 1, 1], [1, 1, 0, 0], [1, 1, 0, 0]],
            dtype=np.float32,
        ),
        nodes,
        4,
    )

    result, null_values = calibrate_view_in_context(
        {"generic": generic, "mechanism": mechanism},
        randomized_view_name="mechanism",
        source_genes={"A", "B"},
        target_genes={"C", "D"},
        embedding_dim=3,
        n_randomizations=15,
        seed=2,
    )

    assert result["view"] == "mechanism"
    assert result["null_type"] == "target_view_node_labels_permuted_other_views_fixed"
    assert len(null_values) == 15


def test_synthetic_external_benchmark_emits_complete_result_contract(tmp_path):
    from src.baselines.run_multiview_fusion_benchmark import (
        BenchmarkConfig,
        model_variant_specs,
        run_external_benchmark,
    )
    from src.models.multiview_fusion import DEFAULT_VIEW_SPEC

    graph_dir = tmp_path / "graph"
    expression_dir = tmp_path / "expression"
    output_dir = tmp_path / "output"
    graph_dir.mkdir()
    expression_dir.mkdir()
    _nodes().to_csv(graph_dir / "nodes.tsv", sep="\t", index=False)
    _edges().to_csv(graph_dir / "edges.tsv", sep="\t", index=False)

    rng = np.random.RandomState(4)
    stad_ids = [f"STAD_{i}" for i in range(16)]
    lihc_ids = [f"LIHC_{i}" for i in range(12)]
    columns = list(_nodes()["node_id"])
    stad = pd.DataFrame(rng.normal(size=(16, len(columns))), index=stad_ids, columns=columns)
    lihc = pd.DataFrame(rng.normal(size=(12, len(columns))), index=lihc_ids, columns=columns)
    stad.to_csv(expression_dir / "tcga_stad_expression.tsv", sep="\t")
    lihc.to_csv(expression_dir / "tcga_lihc_expression.tsv", sep="\t")
    labels = pd.DataFrame(
        {
            "dataset": ["TCGA-STAD"] * 16 + ["TCGA-LIHC"] * 12,
            "nk_immune_state": (["NK-hot-cytotoxic", "NK-hot-dysfunctional"] * 8)
            + (["NK-hot-cytotoxic", "NK-hot-dysfunctional"] * 6),
        },
        index=stad_ids + lihc_ids,
    )
    labels_path = tmp_path / "labels.tsv"
    labels.to_csv(labels_path, sep="\t")

    config = BenchmarkConfig(
        graph_dir=graph_dir,
        labels_path=labels_path,
        expression_dir=expression_dir,
        output_dir=output_dir,
        seeds=(17,),
        feature_modes=("masked",),
        embedding_dim=4,
        max_epochs=2,
        patience=1,
        batch_size=8,
        n_bootstrap=20,
    )
    result = run_external_benchmark(config)

    expected_variants = model_variant_specs(list(DEFAULT_VIEW_SPEC))
    assert len(result.per_seed) == len(expected_variants)
    assert set(result.per_seed["variant"]) == set(expected_variants)
    assert result.per_seed[["AUROC", "AUPRC", "MCC", "BalancedAccuracy"]].notna().all().all()
    assert result.verdict in {"stable_external_gain", "no_stable_external_gain"}
    assert (output_dir / "multiview_benchmark_per_seed.tsv").exists()
    assert (output_dir / "multiview_external_predictions.tsv").exists()
    ensemble = pd.read_csv(output_dir / "multiview_external_predictions_ensemble.tsv", sep="\t")
    assert ensemble.groupby(["feature_mode", "variant", "sample_id"]).size().eq(1).all()
    assert (output_dir / "multiview_leave_one_out_bootstrap.tsv").exists()
    assert (output_dir / "multiview_weight_stability.tsv").exists()


def test_t19_compatibility_wrapper_delegates_to_canonical_entrypoint():
    from src.a100_recompute import run_t19_multiview_fusion as compatibility
    from src.baselines import run_multiview_fusion_benchmark as canonical

    assert compatibility.main is canonical.main


def test_calibration_module_pairs_keep_upstream_and_effector_layers_separate():
    from src.interpretation.run_multiview_calibration import calibration_module_pairs

    modules = {
        "tumor_serine_capacity": {"genes": ["PHGDH"]},
        "nk_sm_synthesis": {"genes": ["SGMS1"]},
        "nk_sm_catabolism": {"genes": ["SMPD1"]},
        "nk_protrusion_machinery": {"genes": ["EZR"]},
        "nk_synapse_cytotoxicity_outcome": {"genes": ["NKG7"]},
    }

    pairs = calibration_module_pairs(modules)

    assert pairs["metabolic_crosstalk"] == ({"PHGDH"}, {"SGMS1", "SMPD1", "EZR"})
    assert pairs["sm_topology_axis"] == ({"SGMS1", "SMPD1", "EZR"}, {"NKG7"})


def test_supplementary_figure_writes_png_and_pdf(tmp_path):
    from src.figures.make_multiview_supplement import make_multiview_supplement

    summary = pd.DataFrame(
        {
            "feature_mode": ["masked"] * 4,
            "variant": ["no_graph", "merged_svd", "uniform_multiview", "learned_multiview"],
            "AUROC_mean": [0.70, 0.72, 0.73, 0.74],
            "AUROC_std": [0.02] * 4,
            "AUPRC_mean": [0.60, 0.61, 0.62, 0.63],
            "AUPRC_std": [0.03] * 4,
        }
    )
    weights = pd.DataFrame(
        {
            "feature_mode": ["masked"] * 4,
            "variant": ["learned_multiview"] * 4,
            "seed": [1, 1, 2, 2],
            "view": ["generic_interaction", "metabolic_crosstalk"] * 2,
            "weight": [0.6, 0.4, 0.55, 0.45],
        }
    )

    png_path, pdf_path = make_multiview_supplement(summary, weights, tmp_path / "figS1_multiview")

    assert png_path.exists() and png_path.stat().st_size > 0
    assert pdf_path.exists() and pdf_path.stat().st_size > 0


def test_weight_stability_requires_same_top_view_in_eight_of_ten_seeds():
    from src.baselines.run_multiview_fusion_benchmark import summarize_weight_stability

    rows = []
    for seed in range(10):
        top = "ppi" if seed < 8 else "mechanism"
        rows.extend(
            [
                {
                    "feature_mode": "masked",
                    "variant": "learned_multiview",
                    "seed": seed,
                    "view": "ppi",
                    "weight": 0.6 if top == "ppi" else 0.4,
                },
                {
                    "feature_mode": "masked",
                    "variant": "learned_multiview",
                    "seed": seed,
                    "view": "mechanism",
                    "weight": 0.4 if top == "ppi" else 0.6,
                },
            ]
        )

    summary = summarize_weight_stability(pd.DataFrame(rows))
    ppi = summary.set_index("view").loc["ppi"]
    mechanism = summary.set_index("view").loc["mechanism"]

    assert ppi["top_rank_count"] == 8
    assert bool(ppi["stable_preference"]) is True
    assert mechanism["top_rank_count"] == 2
    assert bool(mechanism["stable_preference"]) is False
