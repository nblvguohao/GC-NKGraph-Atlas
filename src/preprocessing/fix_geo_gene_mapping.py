"""
Fix GEO external-cohort probe->gene-symbol mapping (Blocker C).

PROBLEM
-------
`run_bulk_preprocessing.load_geo_expression()` keeps microarray *probe IDs*
(e.g. Affymetrix `1007_s_at`, Illumina `ILMN_1343291`) as the matrix index.
`standardize_gene_symbols()` only remaps ~20 curated aliases, so downstream
gene-symbol lookups (NK markers, SST-axis genes) find nothing. LOG.md shows
`NK_MARKERS: 0` for GSE62254/GSE84437 and every sample collapses to
`NK-intermediate`. The two external-validation cohorts are unusable.

FIX
---
Map probe IDs -> HGNC symbols using the GEO *platform* (GPL) annotation, then
collapse probes to gene level (keep highest-mean-expression probe per gene,
matching the project's `duplicate_strategy: keep_highest_mean_expression`).

Writes corrected, gene-indexed, sample-by-gene matrices to
`data/processed/bulk/<gse>_expression.tsv`, overwriting the broken ones.

USAGE (run on the compute server, raw data present)
---------------------------------------------------
    python src/preprocessing/fix_geo_gene_mapping.py --gse GSE62254 GSE84437

Requires GEOparse (already a project dependency). If GEOparse cannot fetch the
platform, drop the GPL annotation SOFT file into data/raw/bulk/ and pass
--gpl-file; the script will parse the ID/Symbol columns directly.

NOTE: untested on the author's workstation (no raw data locally). Verify the
first run's LOG output shows non-zero NK gene counts before trusting results.
"""

import argparse
import gzip
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.logging import Logger

RAW_DIR = "data/raw/bulk"
OUT_DIR = "data/processed/bulk"

# Known platform per cohort (from configs + probe-count sanity check in LOG.md:
# GSE62254 = 54675 features = GPL570; GSE84437 = 49576 features ~ Illumina HT-12).
DEFAULT_PLATFORM = {
    "GSE62254": "GPL570",     # Affymetrix HG-U133 Plus 2.0
    "GSE84437": "GPL6947",    # Illumina HumanHT-12 V3 (verify: 49576 features)
}

# Candidate column names that hold the gene symbol in a GPL annotation table,
# in priority order. GPL SOFT tables are inconsistent across platforms.
SYMBOL_COLUMN_CANDIDATES = [
    "Gene Symbol", "Gene symbol", "GENE_SYMBOL", "Symbol", "gene_symbol",
    "ILMN_Gene", "GENE", "Gene_Symbol", "geneSymbol",
]


def _open_maybe_gzip(path, mode="rt"):
    if str(path).endswith((".gz", ".gzip")):
        return gzip.open(path, mode)
    with open(path, "rb") as f:
        if f.read(2) == b"\x1f\x8b":
            return gzip.open(path, mode)
    return open(path, mode)


def load_geo_matrix(gse, logger):
    """Load the probe-level series matrix (samples-by-probe not yet transposed)."""
    candidates = [
        os.path.join(RAW_DIR, f"{gse}_series_matrix.txt.gz"),
        os.path.join(RAW_DIR, f"{gse}_series_matrix.txt"),
    ]
    path = next((p for p in candidates if os.path.exists(p)), None)
    if path is None:
        raise FileNotFoundError(f"{gse}: no series matrix in {RAW_DIR}")

    with _open_maybe_gzip(path, "rt") as f:
        lines = f.readlines()
    start = end = None
    for i, ln in enumerate(lines):
        if "!series_matrix_table_begin" in ln:
            start = i + 1
        elif "!series_matrix_table_end" in ln:
            end = i
    if start is None or end is None:
        raise ValueError(f"{gse}: cannot locate expression table")

    header = lines[start].rstrip("\n").split("\t")
    rows = [ln.rstrip("\n").split("\t") for ln in lines[start + 1:end] if "\t" in ln]
    df = pd.DataFrame(rows, columns=header).set_index(header[0])
    df.index = [p.strip().strip('"') for p in df.index]
    df.columns = [c.strip().strip('"') for c in df.columns]
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(how="all")
    logger.ok("PREPROCESSING", f"{gse}: loaded probe matrix {df.shape}", script=__file__)
    return df  # index = probe IDs, columns = samples


def build_probe_to_symbol(gse, platform, gpl_file, logger):
    """Return a dict probe_id -> gene_symbol from the GPL annotation."""
    # 1) explicit SOFT file, if provided
    if gpl_file and os.path.exists(gpl_file):
        return _parse_gpl_softfile(gpl_file, logger)

    # 2) GEOparse (fetches/caches the platform annotation)
    try:
        import GEOparse
        gpl = GEOparse.get_GEO(geo=platform, destdir=RAW_DIR, silent=True)
        tbl = gpl.table
        id_col = "ID" if "ID" in tbl.columns else tbl.columns[0]
        sym_col = next((c for c in SYMBOL_COLUMN_CANDIDATES if c in tbl.columns), None)
        if sym_col is None:
            raise ValueError(f"{platform}: no symbol column in GPL table "
                             f"(columns: {list(tbl.columns)})")
        mapping = _clean_mapping(tbl.set_index(id_col)[sym_col])
        logger.ok("PREPROCESSING",
                   f"{gse}: {platform} probe->symbol map = {len(mapping)} probes "
                   f"(symbol col '{sym_col}')", script=__file__)
        return mapping
    except Exception as e:
        raise RuntimeError(
            f"{gse}: could not build probe->symbol map via GEOparse ({e}). "
            f"Provide the GPL SOFT annotation with --gpl-file."
        )


def _parse_gpl_softfile(path, logger):
    """Parse a GPLxxxx SOFT annotation file's platform table directly."""
    with _open_maybe_gzip(path, "rt") as f:
        lines = f.readlines()
    start = None
    for i, ln in enumerate(lines):
        if ln.startswith("!platform_table_begin"):
            start = i + 1
            break
    if start is None:
        # some annot files are plain TSV with a header row containing "ID"
        tbl = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    else:
        header = lines[start].rstrip("\n").split("\t")
        rows = [ln.rstrip("\n").split("\t") for ln in lines[start + 1:]
                if not ln.startswith("!platform_table_end") and "\t" in ln]
        tbl = pd.DataFrame(rows, columns=header)
    id_col = "ID" if "ID" in tbl.columns else tbl.columns[0]
    sym_col = next((c for c in SYMBOL_COLUMN_CANDIDATES if c in tbl.columns), None)
    if sym_col is None:
        raise ValueError(f"GPL file {path}: no recognizable symbol column")
    return _clean_mapping(tbl.set_index(id_col)[sym_col])


def _clean_mapping(series):
    """Normalize a probe->symbol Series into a clean dict (first symbol wins)."""
    out = {}
    for probe, sym in series.items():
        if not isinstance(sym, str):
            continue
        sym = sym.strip().strip('"')
        if not sym or sym in {"---", "NA", "NaN"}:
            continue
        # multi-mapping probes: "GENEA /// GENEB" -> take the first symbol
        sym = sym.split("///")[0].split("//")[0].strip()
        if sym:
            out[str(probe).strip().strip('"')] = sym
    return out


def collapse_to_genes(expr_probe, probe2sym, logger, gse):
    """Map probe rows -> symbols and collapse to gene level (highest mean)."""
    symbols = pd.Series(expr_probe.index, index=expr_probe.index).map(probe2sym)
    mapped = expr_probe[symbols.notna().values].copy()
    mapped.index = symbols.dropna().values
    n_mapped = mapped.shape[0]

    # collapse duplicate symbols: keep the probe with the highest mean expression
    means = mapped.mean(axis=1)
    keep = means.groupby(mapped.index).idxmax()
    # idxmax returns positional-safe labels only if index unique; use groupby-transform
    gene_level = mapped.loc[~mapped.index.duplicated(keep=False)]
    dup_genes = mapped.index[mapped.index.duplicated(keep=False)].unique()
    dup_rows = []
    for g in dup_genes:
        sub = mapped.loc[g]
        if isinstance(sub, pd.DataFrame):
            best = sub.iloc[sub.mean(axis=1).values.argmax()]
        else:
            best = sub
        best.name = g
        dup_rows.append(best)
    if dup_rows:
        gene_level = pd.concat([gene_level, pd.DataFrame(dup_rows)])
    gene_level = gene_level[~gene_level.index.duplicated(keep="first")].sort_index()

    logger.ok("PREPROCESSING",
              f"{gse}: {expr_probe.shape[0]} probes -> {n_mapped} mapped -> "
              f"{gene_level.shape[0]} unique genes", script=__file__)
    return gene_level  # index = gene symbols, columns = samples


def sanity_check(gene_level, gse, logger):
    """Warn loudly if canonical NK markers are still absent."""
    nk_markers = ["NKG7", "GNLY", "GZMB", "PRF1", "KLRD1", "NCAM1", "KLRK1"]
    found = [g for g in nk_markers if g in gene_level.index]
    msg = f"{gse}: NK marker check {len(found)}/{len(nk_markers)} present ({found})"
    if len(found) >= 4:
        logger.ok("PREPROCESSING", msg, script=__file__)
    else:
        logger.fail("PREPROCESSING", msg + " -- mapping likely wrong platform",
                    script=__file__)


def main():
    ap = argparse.ArgumentParser(description="Fix GEO probe->gene mapping (Blocker C)")
    ap.add_argument("--gse", nargs="+", required=True, help="GEO accessions, e.g. GSE62254 GSE84437")
    ap.add_argument("--platform", nargs="*", default=None,
                    help="GPL id per --gse (default: built-in guess)")
    ap.add_argument("--gpl-file", default=None,
                    help="Path to a single GPL SOFT annotation file (overrides GEOparse)")
    args = ap.parse_args()

    logger = Logger()
    os.makedirs(OUT_DIR, exist_ok=True)

    for i, gse in enumerate(args.gse):
        platform = (args.platform[i] if args.platform and i < len(args.platform)
                    else DEFAULT_PLATFORM.get(gse))
        if platform is None and not args.gpl_file:
            logger.fail("PREPROCESSING",
                        f"{gse}: no platform known; pass --platform or --gpl-file",
                        script=__file__)
            continue
        try:
            expr_probe = load_geo_matrix(gse, logger)
            probe2sym = build_probe_to_symbol(gse, platform, args.gpl_file, logger)
            gene_level = collapse_to_genes(expr_probe, probe2sym, logger, gse)
            sanity_check(gene_level, gse, logger)

            # transpose to samples-by-genes (project convention) and save
            out = gene_level.T
            out.index.name = "sample_id"
            out_path = os.path.join(OUT_DIR, f"{gse.lower()}_expression.tsv")
            out.to_csv(out_path, sep="\t")
            logger.ok("PREPROCESSING",
                      f"{gse}: saved corrected matrix {out.shape} -> {out_path}",
                      script=__file__)
        except Exception as e:
            logger.fail("PREPROCESSING", f"{gse}: FAILED - {e}", script=__file__)


if __name__ == "__main__":
    main()
