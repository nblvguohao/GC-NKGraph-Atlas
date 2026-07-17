from pathlib import Path

import pytest
import pandas as pd

from src.interpretation.recoverability_modalities import (
    REQUIRED_GSE251950_GSMS,
    not_measured,
    run_visium_spot_module_adjacency,
)


def test_missing_protein_matrix_is_not_measured_not_negative():
    row = not_measured("protein", "MICA abundance", "public matrix unavailable")
    assert row["status"] == "not_measured"
    assert row["not_measured_reason"] == "public matrix unavailable"


def test_real_gse251950_archives_produce_four_spot_adjacency_rows():
    """Only verified GEO per-GSM archives are eligible for spatial evidence."""
    archive_dir = Path("data/external/recoverability/GSE251950/per_gsm")
    archives = sorted(archive_dir.glob("GSM*.tar.gz"))
    if len(archives) != 4:
        pytest.skip("four verified GSE251950 per-GSM archives are not available locally")

    result = run_visium_spot_module_adjacency(archives, n_permutations=10, seed=17)

    assert result["gsm"].tolist() == sorted(REQUIRED_GSE251950_GSMS)
    assert set(result["status"]) == {"measured"}
    assert set(result["scope"]) == {"four_verified_per_gsm_subset"}
    assert (result["coordinate_label_permutation_pvalue"].between(0, 1)).all()
    assert (result["caf_ecm_feature_coverage"] > 0).all()
    assert (result["nk_cytolytic_feature_coverage"] > 0).all()
    assert (result["spot_grid_edge_count"] > 0).all()


def test_real_spatial_analysis_refuses_an_incomplete_verified_subset():
    archive_dir = Path("data/external/recoverability/GSE251950/per_gsm")
    archives = sorted(archive_dir.glob("GSM*.tar.gz"))
    if len(archives) != 4:
        pytest.skip("four verified GSE251950 per-GSM archives are not available locally")

    with pytest.raises(ValueError, match="exactly the four verified"):
        run_visium_spot_module_adjacency(archives[:3], n_permutations=2, seed=17)


def test_direct_modality_writer_preserves_four_real_spatial_rows(tmp_path):
    from src.interpretation.run_recoverability_modalities import write_direct_modality_table

    archive_dir = Path("data/external/recoverability/GSE251950/per_gsm")
    if len(list(archive_dir.glob("GSM*.tar.gz"))) != 4:
        pytest.skip("four verified GSE251950 per-GSM archives are not available locally")
    output = tmp_path / "direct_modality.tsv"
    write_direct_modality_table(output, archive_dir=archive_dir, n_permutations=2, seed=17)

    table = pd.read_csv(output, sep="\t")
    spatial = table.loc[table["modality"] == "spatial_transcriptomics"]
    assert len(spatial) == 4
    assert spatial["gsm"].notna().all()
    assert set(spatial["scope"]) == {"four_verified_per_gsm_subset"}
    assert not table.isna().any().any()
