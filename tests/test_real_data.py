from pathlib import Path

import pytest

from src.common.real_data import RealDataAsset, assert_real_asset, load_real_data_manifest


def test_manifest_requires_accession_url_hash_and_human_species(tmp_path):
    path = tmp_path / "manifest.yaml"
    path.write_text("assets:\n  x: {accession: GSE1}\n", encoding="utf-8")

    with pytest.raises(ValueError, match="source_url"):
        load_real_data_manifest(path)


def test_guard_rejects_synthetic_mock_and_demo_paths():
    asset = RealDataAsset(
        accession="GSE1",
        source_url="https://example.org/x",
        modality="RNA",
        species="Homo sapiens",
        sample_count=2,
        sha256="a" * 64,
        local_path="x.tsv",
    )

    with pytest.raises(ValueError, match="non-real"):
        assert_real_asset(asset, Path("synthetic.tsv"))


def test_manifest_allows_pending_download_hash_only_before_retrieval(tmp_path):
    path = tmp_path / "manifest.yaml"
    path.write_text(
        """assets:
  x:
    accession: GSE1
    source_url: https://example.org/x
    modality: RNA
    species: Homo sapiens
    sample_count: 2
    local_path: data/external/recoverability/GSE1/x.tsv
    status: pending_download
""",
        encoding="utf-8",
    )

    asset = load_real_data_manifest(path)["x"]
    assert asset.status == "pending_download"
    assert asset.sha256 is None


def test_manifest_allows_failed_integrity_metadata_but_guard_rejects_it(tmp_path):
    path = tmp_path / "manifest.yaml"
    path.write_text(
        """assets:
  failed_outer:
    accession: GSE251950
    source_url: https://example.org/GSE251950_RAW.tar
    modality: visium_spatial_transcriptomics
    species: Homo sapiens
    sample_count: 10
    sha256: aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
    local_path: data/external/recoverability/GSE251950/GSE251950_RAW.tar
    status: failed_integrity
    integrity_reason: nested archives are damaged
""",
        encoding="utf-8",
    )
    asset = load_real_data_manifest(path)["failed_outer"]
    assert asset.status == "failed_integrity"

    with pytest.raises(ValueError, match="available"):
        assert_real_asset(asset, path)
