"""
One-shot data processing: GSE246662 scRNA + GEO series matrices + STRING PPI.
Produces all intermediate files needed by T9/T11/T12/T14.
"""
import gzip, pandas as pd, numpy as np, os, re
from scipy import stats

# =============================================================================
# STEP 1: GSE246662 scRNA-seq — merge, QC, NK identification, SST scores
# =============================================================================
print("=" * 60)
print("STEP 1: GSE246662 scRNA-seq processing")
print("=" * 60)

SAMPLES = {
    'HL1': 'GSM7874169_HL1.csv.gz', 'HL2': 'GSM7874170_HL2.csv.gz', 'HL3': 'GSM7874171_HL3.csv.gz',
    'GC1': 'GSM7874172_GC1.csv.gz', 'GC2': 'GSM7874173_GC2.csv.gz', 'GC3': 'GSM7874174_GC3.csv.gz',
    'LM1': 'GSM7874175_LM1.csv.gz', 'LM2': 'GSM7874176_LM2.csv.gz', 'LM3': 'GSM7874177_LM3.csv.gz',
}
TISSUE_MAP = {'HL': 'healthy_liver', 'GC': 'gastric_cancer', 'LM': 'liver_metastasis'}

# First pass: find common genes (efficient: read only the gene name column)
# NOTE: Some samples are genes×cells (genes=rows), others are cells×genes (genes=cols).
# We need to detect orientation and unify.

all_gene_sets = []
ORIENTATION = {}  # 'genes_x_cells' or 'cells_x_genes'

for sample_name, filename in SAMPLES.items():
    # Read header to detect orientation
    with gzip.open(f"data/raw/scrna/{filename}", "rt") as fh:
        header_line = fh.readline().strip()

    header_cols = header_line.replace('"', '').split(',')
    # Check if second column looks like a gene ID (Ensembl: starts with A, contains ".")
    second_col = header_cols[1] if len(header_cols) > 1 else ""
    is_gene_header = (second_col.startswith('A') and '.' in second_col) or second_col in ('NKG7', 'GNLY')

    # Also check first data row
    with gzip.open(f"data/raw/scrna/{filename}", "rt") as fh:
        fh.readline()  # header
        first_line = fh.readline().strip()
    first_val = first_line.replace('"', '').split(',')[0]
    is_gene_first_col = (first_val.startswith('A') and '.' in first_val) or first_val in ('NKG7', 'GNLY', 'GZMB')

    if is_gene_first_col and not is_gene_header:
        # genes × cells: genes in first column, cells in header
        orientation = 'genes_x_cells'
        # Read only the gene name column
        df_genes = pd.read_csv(f"data/raw/scrna/{filename}", compression="gzip",
                               usecols=[0], skiprows=1, header=None, names=['gene'])
        gene_ids = [str(g) for g in df_genes['gene'].values]
    elif is_gene_header and not is_gene_first_col:
        # cells × genes: genes in header, cells in first column — need to transpose
        orientation = 'cells_x_genes'
        df_tmp = pd.read_csv(f"data/raw/scrna/{filename}", compression="gzip",
                             index_col=0, nrows=0)
        gene_ids = [str(c) for c in df_tmp.columns]
    else:
        # Ambiguous: try reading both ways and pick the one with more genes
        # Check column headers first
        df_cols = pd.read_csv(f"data/raw/scrna/{filename}", compression="gzip", nrows=0)
        col_headers = [str(c) for c in df_cols.columns[1:]]
        # Check first column
        df_rows = pd.read_csv(f"data/raw/scrna/{filename}", compression="gzip",
                              usecols=[0], skiprows=1, header=None, names=['gene'])
        row_labels = [str(g) for g in df_rows['gene'].values]

        # Heuristic: the list with gene-like entries (short, no dashes, includes known genes) is genes
        def gene_likeness(names):
            score = 0
            for n in names[:100]:
                if n in ('NKG7', 'GNLY', 'GZMB', 'PHGDH', 'SGMS1', 'ACTB'):
                    score += 100
                if '.' in n and not '-' in n and len(n) < 15:
                    score += 1
                if '-' in n and len(n) > 15:  # cell barcode
                    score -= 1
            return score

        col_score = gene_likeness(col_headers)
        row_score = gene_likeness(row_labels)

        if row_score > col_score:
            orientation = 'genes_x_cells'
            gene_ids = row_labels
        else:
            orientation = 'cells_x_genes'
            gene_ids = col_headers

    ORIENTATION[sample_name] = orientation
    gene_set = set(gene_ids)
    all_gene_sets.append(gene_set)
    print(f"  {sample_name}: {len(gene_set)} genes, orientation={orientation}")

common_genes = sorted(set.intersection(*all_gene_sets))
print(f"Common genes across all samples: {len(common_genes)}")

# Build gene index
gene_to_idx = {g: i for i, g in enumerate(common_genes)}

# Second pass: load all cells
cell_meta = []
all_expr_chunks = []

for sample_name, filename in SAMPLES.items():
    tissue_code = sample_name[:2]
    tissue = TISSUE_MAP[tissue_code]
    orientation = ORIENTATION[sample_name]

    df = pd.read_csv(f"data/raw/scrna/{filename}", compression="gzip", index_col=0)

    if orientation == 'cells_x_genes':
        # Transpose: cells × genes → genes × cells
        df = df.T
        # Now df is genes × cells with gene names as index

    # Subset to common genes
    df = df.loc[df.index.isin(common_genes)]
    df = df.reindex(common_genes)
    n_cells = df.shape[1]
    print(f"  {sample_name} ({tissue}): {n_cells} cells")

    for col in df.columns:
        cell_meta.append({
            'cell_id': f"{sample_name}_{col}",
            'sample': sample_name,
            'tissue': tissue,
        })

    all_expr_chunks.append(df.values.T.astype(np.float64))  # cells × genes

X_all = np.concatenate(all_expr_chunks, axis=0)
n_total = X_all.shape[0]
print(f"\nTotal cells: {n_total}")

# QC filtering
genes_per_cell = (X_all > 0).sum(axis=1)
mito_genes = [g for g in common_genes if g.startswith('MT-')]
mito_idx = [gene_to_idx[g] for g in mito_genes]
pct_mito = np.sum(X_all[:, mito_idx], axis=1) / (np.sum(X_all, axis=1) + 1e-10) * 100

qc_mask = (genes_per_cell >= 200) & (genes_per_cell <= 6000) & (pct_mito <= 20)
print(f"Cells passing QC: {qc_mask.sum()} / {n_total} ({qc_mask.sum()/n_total*100:.1f}%)")

X = X_all[qc_mask]
meta = [cell_meta[i] for i in range(n_total) if qc_mask[i]]
n_qc = X.shape[0]

# Normalize
libsize = X.sum(axis=1)
X_norm = X / (libsize[:, None] + 1e-10) * 10000
X_log = np.log1p(X_norm)

# Marker-based cell type identification
NK_MARKERS = ['NKG7', 'GNLY', 'KLRD1', 'KLRF1', 'NCR1', 'NCAM1', 'FCGR3A', 'KLRK1', 'EOMES']
T_MARKERS = ['CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A']
MONO_MARKERS = ['CD14', 'CD68', 'CSF1R']
B_MARKERS = ['MS4A1', 'CD19', 'CD79A']

def mean_expr(X, markers):
    idxs = [gene_to_idx.get(g) for g in markers if g in gene_to_idx]
    if not idxs: return np.zeros(X.shape[0])
    return X[:, idxs].mean(axis=1)

nk_score = mean_expr(X_log, NK_MARKERS)
t_score = mean_expr(X_log, T_MARKERS)
mono_score = mean_expr(X_log, MONO_MARKERS)
b_score = mean_expr(X_log, B_MARKERS)

is_nk = (nk_score > 0.5) & (nk_score > np.maximum(np.maximum(t_score, mono_score), b_score))
nk_count = is_nk.sum()
print(f"NK cells: {nk_count} ({nk_count/n_qc*100:.1f}%)")

nk_X = X_log[is_nk]
nk_meta_list = [meta[i] for i in range(n_qc) if is_nk[i]]

# SST module scores
def module_zscore(expr, markers):
    idxs = [gene_to_idx.get(g) for g in markers if g in gene_to_idx]
    if len(idxs) < 1:
        return np.full(expr.shape[0], np.nan)
    sub = expr[:, idxs].astype(np.float64)
    z = stats.zscore(sub, axis=0, nan_policy='omit')
    z = np.nan_to_num(z, 0)
    return np.nanmean(z, axis=1)

SST_MODULES = {
    'tumor_serine_capacity': ['PHGDH','PSAT1','PSPH','SHMT1','SHMT2','MTHFD1','MTHFD2','MTHFD1L','SLC1A4','SLC1A5'],
    'nk_sm_synthesis': ['SGMS1','SGMS2'],
    'nk_sm_catabolism': ['SMPD1','SMPD2','SMPD3','SMPD4'],
    'nk_protrusion_machinery': ['EZR','MSN','RDX','ACTR2','ACTR3','ARPC1B','ARPC2','ARPC3','ARPC4','ARPC5',
                                'WAS','WASL','WASF1','WASF2','WASF3','WIPF1','CDC42','RAC1','RHOA',
                                'DIAPH1','DIAPH3','FMNL1','BAIAP2','PACSIN2'],
    'nk_cytotoxicity_outcome': ['NKG7','GNLY','GZMB','PRF1','IFNG','LCP2','LAT','VAV1','TLN1','ITGAL','ITGB2'],
    'checkpoint_link': ['HAVCR2'],
}

module_scores = {}
for mod_name, genes in SST_MODULES.items():
    s = module_zscore(nk_X, genes)
    module_scores[mod_name] = s
    present = sum(1 for g in genes if g in gene_to_idx)
    print(f"  {mod_name}: {present}/{len(genes)} genes, mean={np.nanmean(s):.4f}")

# Build output
df_out = pd.DataFrame({
    'nk_sm_balance_score': module_scores['nk_sm_synthesis'] - module_scores['nk_sm_catabolism'],
    'nk_protrusion_machinery_score': module_scores['nk_protrusion_machinery'],
    'nk_cytotoxicity_outcome_score': module_scores['nk_cytotoxicity_outcome'],
    'nk_sm_synthesis_score': module_scores['nk_sm_synthesis'],
    'nk_sm_catabolism_score': module_scores['nk_sm_catabolism'],
    'tumor_serine_capacity_score': module_scores['tumor_serine_capacity'],
    'checkpoint_link_score': module_scores['checkpoint_link'],
    'tissue': [m['tissue'] for m in nk_meta_list],
    'sample': [m['sample'] for m in nk_meta_list],
})

os.makedirs("results/tables", exist_ok=True)
df_out.to_csv("results/tables/sst_axis_scores_single_cell.tsv", sep="\t", index=False)
print(f"Wrote {len(df_out)} NK cells to results/tables/sst_axis_scores_single_cell.tsv")

# Correlations
r_h2, p_h2 = stats.pearsonr(df_out['nk_sm_balance_score'], df_out['nk_protrusion_machinery_score'])
r_h3, p_h3 = stats.pearsonr(df_out['nk_protrusion_machinery_score'], df_out['nk_cytotoxicity_outcome_score'])
print(f"\n*** REAL DATA RESULTS ***")
print(f"H2 (sm_balance ~ protrusion): r={r_h2:.5f}, p={p_h2:.2e}, n={len(df_out)}")
print(f"H3 (protrusion ~ cytotoxicity): r={r_h3:.5f}, p={p_h3:.2e}, n={len(df_out)}")

for tissue in df_out['tissue'].unique():
    sub = df_out[df_out['tissue'] == tissue]
    if len(sub) < 10: continue
    r, p = stats.pearsonr(sub['nk_protrusion_machinery_score'], sub['nk_cytotoxicity_outcome_score'])
    print(f"  H3 ({tissue}): r={r:.5f}, p={p:.2e}, n={len(sub)}")

# QC summary
qc_summary = pd.DataFrame([{
    'sample': s,
    'cells_before': sum(1 for m in cell_meta if m['sample'] == s),
    'cells_after_qc': sum(1 for m in meta if m['sample'] == s),
    'nk_cells': sum(1 for m in nk_meta_list if m['sample'] == s),
} for s in SAMPLES.keys()])
qc_summary['retained_frac'] = qc_summary['cells_after_qc'] / qc_summary['cells_before']
qc_summary.to_csv("results/tables/scrna_qc_summary.tsv", sep="\t", index=False)
print(f"\nQC summary:\n{qc_summary.to_string()}")

# =============================================================================
# STEP 2: GEO series matrix processing (GSE62254 + GSE84437)
# =============================================================================
print("\n" + "=" * 60)
print("STEP 2: GEO series matrix processing")
print("=" * 60)

def parse_geo_series_matrix(filepath):
    """Parse GEO series matrix, return (expr_df, probe_col, sample_cols)."""
    lines = []
    in_data = False
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                in_data = True
                continue
            if line.startswith('!series_matrix_table_end'):
                break
            if in_data:
                lines.append(line.strip())

    # Parse header
    header = lines[0].strip('"').split('"\t"')
    header = [h.strip('"') for h in header]
    n_cols = len(header)

    # Parse data — handle rows with mismatched column counts
    data_rows = []
    for line in lines[1:]:
        parts = line.strip('"').split('"\t"')
        parts = [p.strip('"') for p in parts]
        # Pad or truncate to match header length
        if len(parts) < n_cols:
            parts.extend([''] * (n_cols - len(parts)))
        elif len(parts) > n_cols:
            parts = parts[:n_cols]
        data_rows.append(parts)

    df = pd.DataFrame(data_rows, columns=header)
    df = df.set_index('ID_REF')
    # Convert to numeric
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    return df

# Process GSE62254
print("\nGSE62254 (GPL570 - Affymetrix):")
gse62254 = parse_geo_series_matrix("data/raw/bulk/GSE62254_series_matrix.txt.gz")
print(f"  Shape: {gse62254.shape} ({gse62254.shape[0]} probes x {gse62254.shape[1]} samples)")

# Process GSE84437
print("GSE84437 (GPL6947 - Illumina):")
gse84437 = parse_geo_series_matrix("data/raw/bulk/GSE84437_series_matrix.txt.gz")
print(f"  Shape: {gse84437.shape} ({gse84437.shape[0]} probes x {gse84437.shape[1]} samples)")

# Download GPL annotation files and map probes to genes
# GPL570: Affymetrix Human Genome U133 Plus 2.0
# GPL6947: Illumina HumanHT-12 V3.0

def download_gpl_annotation(gpl_id):
    """Download GPL annotation table from GEO."""
    import urllib.request, io
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{gpl_id[:5]}nnn/{gpl_id}/annot/{gpl_id}.annot.gz"
    # Try alternative URL format
    try:
        with urllib.request.urlopen(url) as resp:
            content = resp.read()
        return content
    except:
        # Try the SOFT format
        url2 = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{gpl_id[:5]}nnn/{gpl_id}/soft/{gpl_id}.soft.gz"
        try:
            with urllib.request.urlopen(url2) as resp:
                content = resp.read()
            return content
        except:
            print(f"  WARNING: Could not download GPL annotation for {gpl_id}")
            return None

def parse_gpl_annot(content, gpl_id):
    """Parse GPL annotation, return probe->gene mapping."""
    if content is None:
        return {}
    mapping = {}
    lines = gzip.decompress(content).decode('utf-8', errors='replace').split('\n')
    in_table = False
    probe_col = None
    gene_col = None

    for line in lines:
        if line.startswith('!platform_table_begin'):
            in_table = True
            continue
        if line.startswith('!platform_table_end'):
            break
        if in_table:
            parts = line.strip().split('\t')
            if probe_col is None:
                # Header row
                for i, h in enumerate(parts):
                    h_clean = h.strip('"')
                    if h_clean == 'ID':
                        probe_col = i
                    if h_clean in ('Gene Symbol', 'Symbol', 'GENE_SYMBOL', 'gene_assignment'):
                        gene_col = i
                continue
            if probe_col is not None and gene_col is not None and len(parts) > max(probe_col, gene_col):
                probe = parts[probe_col].strip('"')
                gene_raw = parts[gene_col].strip('"')
                # Clean gene symbol
                if '///' in gene_raw:
                    gene_raw = gene_raw.split('///')[0].strip()
                if '//' in gene_raw:
                    gene_raw = gene_raw.split('//')[0].strip()
                # Some platforms have gene_assignment format
                if gpl_id == 'GPL570' and ' // ' in gene_raw:
                    gene_raw = gene_raw.split(' // ')[1].strip() if ' // ' in gene_raw else gene_raw
                mapping[probe] = gene_raw if gene_raw else probe
    print(f"  Parsed {len(mapping)} probe->gene mappings for {gpl_id}")
    return mapping

# Map probes to genes
print("\nMapping probes to genes...")

# GPL570 for GSE62254
gpl570_content = download_gpl_annotation('GPL570')
gpl570_map = parse_gpl_annot(gpl570_content, 'GPL570')

# GPL6947 for GSE84437
gpl6947_content = download_gpl_annotation('GPL6947')
gpl6947_map = parse_gpl_annot(gpl6947_content, 'GPL6947')

def collapse_to_genes(expr_df, probe_map):
    """Map probes to genes, taking max-mean probe per gene."""
    # Map probe IDs to gene symbols
    gene_names = []
    for probe in expr_df.index:
        gene = probe_map.get(probe, probe)
        gene_names.append(gene)

    df = expr_df.copy()
    df.index = gene_names

    # Remove rows that couldn't be mapped
    df = df[df.index != df.index.name]  # remove unmapped if index changed

    # Collapse: for each gene, keep probe with highest mean expression
    gene_means = df.mean(axis=1)
    df['_mean'] = gene_means
    df = df.sort_values('_mean', ascending=False)
    df = df[~df.index.duplicated(keep='first')]
    df = df.drop(columns=['_mean'])

    return df

print("\nCollapsing probes to genes...")
gse62254_genes = collapse_to_genes(gse62254, gpl570_map)
gse84437_genes = collapse_to_genes(gse84437, gpl6947_map)

print(f"  GSE62254: {gse62254_genes.shape[0]} genes x {gse62254_genes.shape[1]} samples")
print(f"  GSE84437: {gse84437_genes.shape[0]} genes x {gse84437_genes.shape[1]} samples")

# Check NK markers
nk_check = ['NKG7','GNLY','GZMB','PRF1','KLRD1','NCAM1','KLRK1']
for name, df in [('GSE62254', gse62254_genes), ('GSE84437', gse84437_genes)]:
    found = [g for g in nk_check if g in df.index]
    print(f"  {name} NK markers: {len(found)}/{len(nk_check)} — {found}")

# Save
os.makedirs("data/processed/bulk", exist_ok=True)
gse62254_genes.T.to_csv("data/processed/bulk/gse62254_expression.tsv", sep="\t")
gse84437_genes.T.to_csv("data/processed/bulk/gse84437_expression.tsv", sep="\t")
print("Wrote GSE62254 and GSE84437 gene-level expression matrices")

# =============================================================================
# STEP 3: STRING PPI — map ENSEMBL to gene symbols
# =============================================================================
print("\n" + "=" * 60)
print("STEP 3: STRING PPI processing")
print("=" * 60)

# Read STRING PPI
string_df = pd.read_csv(
    "data/raw/prior_networks/9606.protein.physical.links.v12.0.txt.gz",
    sep=" ", compression="gzip"
)
print(f"STRING PPI: {len(string_df)} interactions")

# Map ENSEMBL protein IDs to gene symbols using the STRING alias file
# Download alias file
import urllib.request
alias_url = "https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz"
try:
    with urllib.request.urlopen(alias_url) as resp:
        alias_content = resp.read()

    # Save alias file
    os.makedirs("data/raw/prior_networks", exist_ok=True)
    with open("data/raw/prior_networks/9606.protein.aliases.v12.0.txt.gz", "wb") as f:
        f.write(alias_content)

    alias_df = pd.read_csv(
        "data/raw/prior_networks/9606.protein.aliases.v12.0.txt.gz",
        sep="\t", compression="gzip", header=None,
        names=['protein_id', 'alias', 'source']
    )
    # Filter to gene_symbol source or take first alias
    gene_aliases = alias_df[alias_df['source'] == 'Ensembl_HGNC_entrezgene_symbol(synonym)']
    if len(gene_aliases) == 0:
        gene_aliases = alias_df[alias_df['source'] == 'Ensembl_HGNC_entrezgene_symbol']
    if len(gene_aliases) == 0:
        # Take any alias; prefer shorter ones
        gene_aliases = alias_df.groupby('protein_id').first().reset_index()

    ensembl_to_gene = dict(zip(gene_aliases['protein_id'], gene_aliases['alias']))
    print(f"ENSEMBL->gene mapping: {len(ensembl_to_gene)} proteins")
except Exception as e:
    print(f"  WARNING: Could not download aliases: {e}")
    print("  Will use raw ENSP IDs as gene names")
    ensembl_to_gene = {}

# Map protein IDs
ensp_pattern = re.compile(r'^9606\.(ENSP\d+)')
string_df['protein1'] = string_df['protein1'].apply(
    lambda x: ensembl_to_gene.get(x, x) if x in ensembl_to_gene else (ensp_pattern.sub(r'\1', x) if ensp_pattern.match(x) else x))
string_df['protein2'] = string_df['protein2'].apply(
    lambda x: ensembl_to_gene.get(x, x) if x in ensembl_to_gene else (ensp_pattern.sub(r'\1', x) if ensp_pattern.match(x) else x))

# Filter high-confidence (combined_score >= 700, as in paper Methods §2.5)
string_high = string_df[string_df['combined_score'] >= 700].copy()
print(f"High-confidence (>=700) edges: {len(string_high)} / {len(string_df)}")

# Save processed PPI
string_high.to_csv("data/processed/graph/string_ppi_high_conf.tsv", sep="\t", index=False)
print("Wrote data/processed/graph/string_ppi_high_conf.tsv")

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "=" * 60)
print("ALL DATA PROCESSING COMPLETE")
print("=" * 60)
print(f"  scRNA NK cells: {len(df_out)} (GSE246662)")
print(f"  TCGA-LIHC tumors: 371 samples (pre-processed)")
print(f"  TCGA-STAD tumors: 415 samples (pre-processed)")
print(f"  GSE62254: {gse62254_genes.shape[0]} genes x {gse62254_genes.shape[1]} samples")
print(f"  GSE84437: {gse84437_genes.shape[0]} genes x {gse84437_genes.shape[1]} samples")
print(f"  STRING PPI high-conf: {len(string_high)} edges")
print(f"\nKey output: results/tables/sst_axis_scores_single_cell.tsv")
print(f"  H2 r={r_h2:.5f}, p={p_h2:.2e}")
print(f"  H3 r={r_h3:.5f}, p={p_h3:.2e}")
