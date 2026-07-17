"""Direct-modality evidence helpers; missing inputs are never imputed.

Spatial results in this module are *spot-module* adjacency measurements on the
Visium array grid.  They neither resolve individual cells nor infer cell-cell
contact distances.
"""

from __future__ import annotations

import gzip
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread

from src.common.real_data import assert_real_asset, load_real_data_manifest


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_REAL_DATA_MANIFEST = ROOT / "configs/recoverability_atlas/real_data_manifest.yaml"


REQUIRED_GSE251950_GSMS = (
    "GSM7990474", "GSM7990476", "GSM7990479", "GSM7990481",
)
CAF_ECM_PROGRAM = (
    "COL1A1", "COL1A2", "COL3A1", "FN1", "FAP", "ACTA2", "TGFBR1",
    "TGFBR2", "LOX", "LOXL2", "POSTN", "THBS2", "SPARC", "VCAN",
)
NK_CYTOLYTIC_MACHINERY = (
    "NKG7", "GNLY", "GZMB", "PRF1", "IFNG", "LCP2", "LAT", "VAV1",
    "TLN1", "ITGAL", "ITGB2",
)

def not_measured(modality: str, endpoint: str, reason: str, accession: str = "") -> dict:
    return {"modality": modality, "direct_endpoint": endpoint, "status": "not_measured",
            "not_measured_reason": reason, "accession": accession}


def _member_by_suffix(archive: tarfile.TarFile, suffix: str) -> tarfile.TarInfo:
    matches = [member for member in archive.getmembers() if member.name.endswith(suffix)]
    if len(matches) != 1:
        raise ValueError(f"archive must contain exactly one {suffix}; found {len(matches)}")
    return matches[0]


def _open_nested_gzip(archive: tarfile.TarFile, suffix: str):
    handle = archive.extractfile(_member_by_suffix(archive, suffix))
    if handle is None:
        raise ValueError(f"could not read {suffix}")
    return gzip.GzipFile(fileobj=handle)


def _read_visium_archive(archive_path: Path) -> tuple[np.ndarray, list[str], list[str], pd.DataFrame]:
    """Read the genuine supplied MTX, barcodes, features, and array coordinates."""
    with tarfile.open(archive_path, "r:gz") as archive:
        matrix = mmread(_open_nested_gzip(archive, "_matrix.mtx.gz")).tocsr()
        barcodes = pd.read_csv(_open_nested_gzip(archive, "_barcodes.tsv.gz"), header=None, sep="\t").iloc[:, 0].astype(str).tolist()
        features = pd.read_csv(_open_nested_gzip(archive, "_features.tsv.gz"), header=None, sep="\t")
        if features.shape[1] < 2:
            raise ValueError("Visium feature table lacks gene-name column")
        genes = features.iloc[:, 1].astype(str).str.upper().tolist()

        coordinate_handle = archive.extractfile(_member_by_suffix(archive, "_tissue_positions_list.csv"))
        if coordinate_handle is None:
            raise ValueError("could not read supplied Visium coordinate file")
        coordinates = pd.read_csv(
            coordinate_handle, header=None,
            names=["barcode", "in_tissue", "array_row", "array_col", "pixel_row", "pixel_col"],
        )

    if matrix.shape != (len(genes), len(barcodes)):
        raise ValueError("Visium matrix dimensions do not match supplied features/barcodes")
    coordinates["barcode"] = coordinates["barcode"].astype(str)
    for column in ("in_tissue", "array_row", "array_col"):
        coordinates[column] = pd.to_numeric(coordinates[column], errors="coerce")
    coordinates = coordinates.dropna(subset=["barcode", "in_tissue", "array_row", "array_col"])
    coordinates = coordinates.loc[coordinates["in_tissue"].astype(int) == 1].copy()
    if coordinates.empty:
        raise ValueError("supplied Visium coordinates could not be parsed into in-tissue spots")
    validate_visium_identifiers(barcodes, coordinates)
    return matrix, genes, barcodes, coordinates


def validate_visium_identifiers(barcodes: list[str], coordinates: pd.DataFrame) -> None:
    """Reject ambiguous expression or supplied in-tissue spot-grid identifiers."""
    barcode_index = pd.Index(barcodes)
    if barcode_index.has_duplicates:
        raise ValueError("duplicate expression barcodes are not permitted")
    if coordinates["barcode"].duplicated().any():
        raise ValueError("duplicate in-tissue coordinate barcodes are not permitted")
    if coordinates.duplicated(["array_row", "array_col"]).any():
        raise ValueError("duplicate in-tissue spatial grid coordinates are not permitted")


def _module_scores(matrix, genes: list[str], barcodes: list[str], coordinates: pd.DataFrame, module: tuple[str, ...]):
    gene_to_indices: dict[str, list[int]] = {}
    for index, gene in enumerate(genes):
        gene_to_indices.setdefault(gene, []).append(index)
    indices = [index for gene in module for index in gene_to_indices.get(gene, [])]
    coverage = len(set(genes).intersection(module)) / len(module)
    if not indices:
        raise ValueError("required module has zero feature coverage in this Visium archive")

    barcode_to_index = {barcode: index for index, barcode in enumerate(barcodes)}
    coordinates = coordinates.loc[coordinates["barcode"].isin(barcode_to_index)].drop_duplicates("barcode").copy()
    if coordinates.empty:
        raise ValueError("supplied Visium coordinates have no intersecting expression barcodes")
    expression_indices = np.array([barcode_to_index[barcode] for barcode in coordinates["barcode"]])
    library_sizes = np.asarray(matrix[:, expression_indices].sum(axis=0)).ravel()
    library_sizes[library_sizes == 0] = 1.0
    module_counts = np.asarray(matrix[indices, :][:, expression_indices].sum(axis=0)).ravel()
    # Per-spot library scaling makes scores comparable within the supplied array.
    score = np.log1p(10_000.0 * module_counts / library_sizes) / len(indices)
    return score, coverage, coordinates.reset_index(drop=True)


def _strict_first_order_hex_grid_edges(coordinates: pd.DataFrame) -> np.ndarray:
    """Construct strict first-order Visium hex-grid adjacency from legal offsets.

    Valid neighbours are only `(0, +/-2)` and `(+/-1, +/-1)` in the supplied
    `array_row`/`array_col` convention. Missing spots remain missing: this
    intentionally does not use k-nearest-neighbour edges across spatial holes.
    These are Visium *spot-grid* relations, not inferred cell contacts.
    """
    array_grid = coordinates[["array_row", "array_col"]].astype(int).to_numpy()
    grid_to_index = {tuple(grid): index for index, grid in enumerate(array_grid)}
    legal_offsets = ((0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1))
    edges = {
        tuple(sorted((index, neighbour_index)))
        for index, (row, column) in enumerate(array_grid)
        for row_offset, column_offset in legal_offsets
        if (neighbour_index := grid_to_index.get((row + row_offset, column + column_offset))) is not None
    }
    if not edges:
        raise ValueError("no strict first-order Visium hex-grid edges were found")
    return np.asarray(sorted(edges), dtype=int)


def _cross_module_adjacency(caf_scores: np.ndarray, nk_scores: np.ndarray, edges: np.ndarray) -> float:
    left, right = edges[:, 0], edges[:, 1]
    return float(np.mean((caf_scores[left] * nk_scores[right] + caf_scores[right] * nk_scores[left]) / 2.0))


def _gsm_from_archive_path(path: Path) -> str:
    match = re.search(r"(GSM\d+)", path.name)
    if match is None:
        raise ValueError(f"cannot identify GSM accession from archive name: {path.name}")
    return match.group(1)


def _validate_verified_visium_archives(paths: list[Path], manifest_path: Path) -> dict[str, Path]:
    """Bind formal spatial reads to the four manifest-tracked real GEO assets."""
    gsm_to_path = {_gsm_from_archive_path(path): path for path in paths}
    if len(gsm_to_path) != len(paths) or set(gsm_to_path) != set(REQUIRED_GSE251950_GSMS):
        raise ValueError("analysis requires exactly the four verified GSE251950 per-GSM archives")
    assets = load_real_data_manifest(manifest_path)
    for gsm, path in gsm_to_path.items():
        asset_key = f"gse251950_{gsm.lower()}"
        if asset_key not in assets:
            raise ValueError(f"manifest is missing verified asset {asset_key}")
        asset = assets[asset_key]
        manifest_file = Path(asset.local_path)
        if not manifest_file.is_absolute():
            manifest_file = ROOT / manifest_file
        if path.resolve() != manifest_file.resolve():
            raise ValueError(f"archive path does not match manifest for {gsm}")
        assert_real_asset(asset, path)
    return gsm_to_path


def _benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """Return monotone BH-FDR adjusted p values in the original sample order."""
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    adjusted_ranked = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = np.minimum(adjusted_ranked, 1.0)
    return adjusted


def run_visium_spot_module_adjacency(
    archive_paths: list[Path], *, n_permutations: int = 1000, seed: int = 20260717,
    manifest_path: Path = DEFAULT_REAL_DATA_MANIFEST,
) -> pd.DataFrame:
    """Measure CAF--NK module adjacency for exactly four verified GSE251950 GSMs.

    The p value is a one-sided coordinate-label permutation test: NK spot scores
    are shuffled across the fixed Visium spot grid while CAF scores and grid edges
    remain fixed.  This is spatial calibration, not a cell-contact measurement.
    """
    if n_permutations < 1:
        raise ValueError("n_permutations must be positive")
    paths = [Path(path) for path in archive_paths]
    gsm_to_path = _validate_verified_visium_archives(paths, Path(manifest_path))

    rows: list[dict] = []
    for sample_offset, gsm in enumerate(sorted(gsm_to_path)):
        matrix, genes, barcodes, coordinates = _read_visium_archive(gsm_to_path[gsm])
        caf_scores, caf_coverage, coordinates = _module_scores(
            matrix, genes, barcodes, coordinates, CAF_ECM_PROGRAM,
        )
        nk_scores, nk_coverage, nk_coordinates = _module_scores(
            matrix, genes, barcodes, coordinates, NK_CYTOLYTIC_MACHINERY,
        )
        if not coordinates["barcode"].equals(nk_coordinates["barcode"]):
            raise ValueError("module scoring produced inconsistent Visium barcode intersections")
        edges = _strict_first_order_hex_grid_edges(coordinates)
        observed = _cross_module_adjacency(caf_scores, nk_scores, edges)
        rng = np.random.default_rng(seed + sample_offset)
        null = np.fromiter(
            (_cross_module_adjacency(caf_scores, rng.permutation(nk_scores), edges) for _ in range(n_permutations)),
            dtype=float, count=n_permutations,
        )
        pvalue = (1.0 + np.count_nonzero(null >= observed)) / (n_permutations + 1.0)
        rows.append({
            "card_id": "tgfb_nk_exclusion",
            "modality": "spatial_transcriptomics",
            "direct_endpoint": "CAF--NK Visium spot-module adjacency",
            "status": "measured",
            "accession": "GSE251950",
            "gsm": gsm,
            "scope": "exploratory_four_verified_per_gsm_subset",
            "n_in_tissue_spots": len(coordinates),
            "spot_grid_edge_count": len(edges),
            "observed_caf_nk_spot_module_adjacency": observed,
            "coordinate_label_permutation_pvalue": pvalue,
            "coordinate_label_permutation_null_mean": float(null.mean()),
            "coordinate_label_permutation_null_sd": float(null.std(ddof=1)),
            "n_coordinate_label_permutations": n_permutations,
            "caf_ecm_feature_coverage": caf_coverage,
            "nk_cytolytic_feature_coverage": nk_coverage,
            "adjacency_definition": "strict first-order Visium hex-grid legal offsets; no cross-hole kNN edges",
            "permutation_definition": "shuffle NK module labels across fixed supplied spot grid",
            "not_measured_reason": "",
        })
    result = pd.DataFrame(rows)
    result["coordinate_label_permutation_bh_fdr"] = _benjamini_hochberg(
        result["coordinate_label_permutation_pvalue"].to_numpy(dtype=float)
    )
    return result
