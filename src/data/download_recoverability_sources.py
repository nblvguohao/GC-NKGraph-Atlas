"""Download and verify allowlisted public inputs for the recoverability atlas."""

from __future__ import annotations

import argparse
import shutil
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.real_data import assert_real_asset, file_sha256, load_real_data_manifest


ROOT = Path(__file__).resolve().parents[2]
MANIFEST_PATH = ROOT / "configs" / "recoverability_atlas" / "real_data_manifest.yaml"
ALLOWLIST = {"GSE122401", "MTBLS3303", "GSE251950"}


@dataclass(frozen=True)
class DownloadRecord:
    accession: str
    path: Path
    sha256: str
    size_bytes: int
    retrieved_at: str


def _asset_for_accession(accession: str):
    if accession not in ALLOWLIST:
        raise ValueError(f"accession is not allowlisted: {accession}")
    assets = load_real_data_manifest(MANIFEST_PATH)
    for asset in assets.values():
        if asset.accession == accession:
            return asset
    raise ValueError(f"allowlisted accession missing from manifest: {accession}")


def resolve_download_url(accession: str) -> str:
    """Return the pre-reviewed official URL for one public accession."""
    return _asset_for_accession(accession).source_url


def record_download(accession: str, path: str | Path) -> DownloadRecord:
    """Create a provenance record for an existing, nonempty allowlisted file."""
    _asset_for_accession(accession)
    archive = Path(path)
    if not archive.is_file() or archive.stat().st_size == 0:
        raise FileNotFoundError(f"download is absent or empty: {archive}")
    return DownloadRecord(
        accession=accession,
        path=archive,
        sha256=file_sha256(archive),
        size_bytes=archive.stat().st_size,
        retrieved_at=datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
    )


def _write_record_to_manifest(record: DownloadRecord) -> None:
    with MANIFEST_PATH.open(encoding="utf-8") as handle:
        manifest = yaml.safe_load(handle)
    for payload in manifest["assets"].values():
        if payload.get("accession") == record.accession:
            payload.update({
                "sha256": record.sha256,
                "local_path": str(record.path).replace("\\", "/"),
                "status": "available",
                "retrieved_at": record.retrieved_at,
            })
            break
    else:
        raise ValueError(f"accession missing from manifest: {record.accession}")
    with MANIFEST_PATH.open("w", encoding="utf-8", newline="\n") as handle:
        yaml.safe_dump(manifest, handle, sort_keys=False, allow_unicode=True)


def download_asset(accession: str, output_root: str | Path) -> DownloadRecord:
    """Stream one reviewed public archive and update its manifest provenance."""
    asset = _asset_for_accession(accession)
    destination = Path(output_root) / accession / Path(asset.local_path).name
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = Request(asset.source_url, headers={"User-Agent": "recoverability-atlas/1.0"})
    with urlopen(request, timeout=60) as response, destination.open("wb") as handle:
        if getattr(response, "status", 200) >= 400:
            raise RuntimeError(f"HTTP error {response.status} for {accession}")
        shutil.copyfileobj(response, handle, length=1024 * 1024)
    record = record_download(accession, destination)
    _write_record_to_manifest(record)
    return record


def verify_all() -> list[DownloadRecord]:
    """Verify every available asset, failing rather than inventing replacement data."""
    records: list[DownloadRecord] = []
    for asset in load_real_data_manifest(MANIFEST_PATH).values():
        if asset.status == "pending_download":
            continue
        assert_real_asset(asset, asset.local_path)
        path = Path(asset.local_path)
        records.append(DownloadRecord(
            accession=asset.accession,
            path=path,
            sha256=file_sha256(path),
            size_bytes=path.stat().st_size,
            retrieved_at=asset.retrieved_at or "pre-existing verified input",
        ))
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--asset", choices=sorted(ALLOWLIST))
    parser.add_argument("--output-root", default="data/external/recoverability")
    parser.add_argument("--metadata-only", action="store_true")
    parser.add_argument("--verify-all", action="store_true")
    args = parser.parse_args()
    if args.verify_all:
        for record in verify_all():
            print(f"{record.accession}\t{record.size_bytes}\t{record.sha256}\t{record.path}")
        return
    if not args.asset:
        parser.error("--asset is required unless --verify-all is used")
    if args.metadata_only:
        print(resolve_download_url(args.asset))
        return
    record = download_asset(args.asset, args.output_root)
    print(f"{record.accession}\t{record.size_bytes}\t{record.sha256}\t{record.path}")


if __name__ == "__main__":
    main()
