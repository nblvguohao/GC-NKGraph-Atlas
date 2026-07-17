from pathlib import Path

import pytest

from src.data.download_recoverability_sources import record_download, resolve_download_url


def test_only_allowlisted_accessions_have_urls():
    assert "GSE122401" in resolve_download_url("GSE122401")
    with pytest.raises(ValueError, match="not allowlisted"):
        resolve_download_url("GSE999999")


def test_record_requires_nonempty_download(tmp_path):
    with pytest.raises(FileNotFoundError):
        record_download("GSE122401", tmp_path / "missing.tar")


def test_record_returns_sha256_for_real_nonempty_file(tmp_path):
    archive = tmp_path / "source.tar"
    archive.write_bytes(b"real public archive bytes")

    record = record_download("GSE122401", archive)
    assert record.accession == "GSE122401"
    assert len(record.sha256) == 64
    assert record.size_bytes == archive.stat().st_size
