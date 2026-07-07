"""
Tests for src/models/gc_nkgraph_atlas.py — GeneGraphEncoder and NKStateClassifier.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest


class TestGeneGraphEncoder:
    """Test the graph encoder that learns gene embeddings from TSV graph files."""

    def test_fit_returns_self(self, sample_graph_dir):
        from src.models.gc_nkgraph_atlas import GeneGraphEncoder
        encoder = GeneGraphEncoder(
            graph_dir=str(sample_graph_dir),
            embedding_dim=8,
        )
        result = encoder.fit()
        assert result is encoder

    def test_fit_populates_embeddings(self, sample_graph_dir):
        from src.models.gc_nkgraph_atlas import GeneGraphEncoder
        encoder = GeneGraphEncoder(
            graph_dir=str(sample_graph_dir),
            embedding_dim=8,
        )
        encoder.fit()
        assert encoder._embeddings is not None
        assert encoder._embeddings.shape[0] == 7  # 7 nodes
        assert encoder._embeddings.shape[1] == 8  # embedding_dim

    def test_transform_returns_correct_shape(self, sample_graph_dir):
        from src.models.gc_nkgraph_atlas import GeneGraphEncoder
        encoder = GeneGraphEncoder(
            graph_dir=str(sample_graph_dir),
            embedding_dim=8,
        )
        encoder.fit()
        genes = ["PHGDH", "SGMS1", "NKG7"]
        result = encoder.transform(genes)
        assert result.shape == (3, 8)
        assert result.dtype == np.float32

    def test_transform_unknown_gene_returns_zeros(self, sample_graph_dir):
        from src.models.gc_nkgraph_atlas import GeneGraphEncoder
        encoder = GeneGraphEncoder(
            graph_dir=str(sample_graph_dir),
            embedding_dim=8,
        )
        encoder.fit()
        result = encoder.transform(["NOT_IN_GRAPH"])
        assert (result == 0.0).all()

    def test_get_gene_embedding(self, sample_graph_dir):
        from src.models.gc_nkgraph_atlas import GeneGraphEncoder
        encoder = GeneGraphEncoder(
            graph_dir=str(sample_graph_dir),
            embedding_dim=8,
        )
        encoder.fit()
        emb = encoder.get_gene_embedding("PHGDH")
        assert emb.shape == (8,)
        assert not np.allclose(emb, 0)  # Should be non-zero for a gene in the graph

    def test_call_transform_before_fit_raises(self, sample_graph_dir):
        from src.models.gc_nkgraph_atlas import GeneGraphEncoder
        encoder = GeneGraphEncoder(
            graph_dir=str(sample_graph_dir),
            embedding_dim=8,
        )
        with pytest.raises(RuntimeError):
            encoder.transform(["PHGDH"])

    def test_known_gene_has_nonzero_embedding(self, sample_graph_dir):
        """PHGDH is connected in the graph → should have non-trivial embedding."""
        from src.models.gc_nkgraph_atlas import GeneGraphEncoder
        encoder = GeneGraphEncoder(
            graph_dir=str(sample_graph_dir),
            embedding_dim=8,
        )
        encoder.fit()
        emb_phgdh = encoder.get_gene_embedding("PHGDH")
        # At least some components should be meaningful
        assert np.any(np.abs(emb_phgdh) > 1e-6)


class TestNKStateClassifier:
    """Test the MLP classifier that uses graph-informed features."""

    @pytest.fixture
    def gene_embeddings(self):
        """8 random gene embeddings for 5 genes."""
        rng = np.random.RandomState(42)
        return rng.randn(5, 8).astype(np.float32)

    @pytest.fixture
    def X_expr(self):
        """50 samples × 5 genes."""
        rng = np.random.RandomState(42)
        return rng.randn(50, 5).astype(np.float32)

    @pytest.fixture
    def y(self):
        """Binary labels."""
        rng = np.random.RandomState(42)
        return (rng.rand(50) > 0.6).astype(np.int64)

    def test_fit_returns_self(self, X_expr, y, gene_embeddings):
        from src.models.gc_nkgraph_atlas import NKStateClassifier
        clf = NKStateClassifier(embedding_dim=8, hidden_dims=[16, 8], num_classes=2)
        result = clf.fit(X_expr, y, gene_embeddings, epochs=10, verbose=False)
        assert result is clf

    def test_predict_returns_correct_shape(self, X_expr, y, gene_embeddings):
        from src.models.gc_nkgraph_atlas import NKStateClassifier
        clf = NKStateClassifier(embedding_dim=8, hidden_dims=[16, 8], num_classes=2)
        clf.fit(X_expr, y, gene_embeddings, epochs=10, verbose=False)
        pred = clf.predict(X_expr, gene_embeddings)
        assert pred.shape == (50,)
        assert set(np.unique(pred)).issubset({0, 1})

    def test_predict_proba_sums_to_one(self, X_expr, y, gene_embeddings):
        from src.models.gc_nkgraph_atlas import NKStateClassifier
        clf = NKStateClassifier(embedding_dim=8, hidden_dims=[16, 8], num_classes=2)
        clf.fit(X_expr, y, gene_embeddings, epochs=10, verbose=False)
        proba = clf.predict_proba(X_expr, gene_embeddings)
        assert proba.shape == (50, 2)
        assert np.allclose(proba.sum(axis=1), 1.0)

    def test_fit_with_validation(self, X_expr, y, gene_embeddings):
        from src.models.gc_nkgraph_atlas import NKStateClassifier
        clf = NKStateClassifier(embedding_dim=8, hidden_dims=[16, 8], num_classes=2)
        clf.fit(
            X_expr[:30], y[:30], gene_embeddings,
            X_val_expr=X_expr[30:], y_val=y[30:],
            epochs=50, verbose=False,
        )
        pred = clf.predict(X_expr[30:], gene_embeddings)
        assert pred.shape == (20,)

    def test_training_improves_accuracy(self, X_expr, gene_embeddings):
        """On separable data, the classifier should achieve near-perfect accuracy."""
        from src.models.gc_nkgraph_atlas import NKStateClassifier

        # Create linearly separable data in the combined feature space
        rng = np.random.RandomState(42)
        n_samples = 80
        n_genes = 5
        emb_dim = 8
        emb = rng.randn(n_genes, emb_dim).astype(np.float32)

        # Generate separable data: class 0 has negative means, class 1 positive
        X0 = rng.randn(n_samples // 2, n_genes).astype(np.float32) - 2.0
        X1 = rng.randn(n_samples // 2, n_genes).astype(np.float32) + 2.0
        X = np.vstack([X0, X1])
        y = np.array([0] * (n_samples // 2) + [1] * (n_samples // 2), dtype=np.int64)

        # Split
        idx = np.arange(n_samples)
        rng.shuffle(idx)
        train_idx = idx[:50]
        test_idx = idx[50:]

        clf = NKStateClassifier(
            embedding_dim=emb_dim, hidden_dims=[32, 16],
            num_classes=2, dropout=0.1, learning_rate=1e-2,
        )
        clf.fit(
            X[train_idx], y[train_idx], emb,
            X_val_expr=X[test_idx], y_val=y[test_idx],
            epochs=100, batch_size=16, verbose=False,
        )
        pred = clf.predict(X[test_idx], emb)
        acc = (pred == y[test_idx]).mean()
        assert acc > 0.7, f"Expected accuracy > 0.7, got {acc:.3f}"
