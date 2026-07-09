"""Process GEO series matrices + STRING PPI (separate from scRNA)."""
import gzip, pandas as pd, numpy as np, os, urllib.request

print("=" * 60)
print("STEP 2: GEO series matrix processing")
print("=" * 60)

def parse_geo_series_matrix(filepath):
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
    header = lines[0].strip('"').split('"\t"')
    header = [h.strip('"') for h in header]
    n_cols = len(header)
    data_rows = []
    for line in lines[1:]:
        parts = line.strip('"').split('"\t"')
        parts = [p.strip('"') for p in parts]
        if len(parts) < n_cols:
            parts.extend([''] * (n_cols - len(parts)))
        elif len(parts) > n_cols:
            parts = parts[:n_cols]
        data_rows.append(parts)
    df = pd.DataFrame(data_rows, columns=header)
    df = df.set_index('ID_REF')
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    return df

print("GSE62254 (GPL570):")
gse62254 = parse_geo_series_matrix("data/raw/bulk/GSE62254_series_matrix.txt.gz")
print(f"  Shape: {gse62254.shape}")
print("GSE84437 (GPL6947):")
gse84437 = parse_geo_series_matrix("data/raw/bulk/GSE84437_series_matrix.txt.gz")
print(f"  Shape: {gse84437.shape}")

def download_gpl(gpl_id):
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{gpl_id[:5]}nnn/{gpl_id}/annot/{gpl_id}.annot.gz"
    try:
        with urllib.request.urlopen(url) as resp:
            return resp.read()
    except:
        print(f"  WARNING: Could not download {gpl_id}")
        return None

def parse_gpl(content, gpl_id):
    if content is None:
        return {}
    mapping = {}
    lines = gzip.decompress(content).decode('utf-8', 'replace').split('\n')
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
                for i, h in enumerate(parts):
                    hc = h.strip('"')
                    if hc == 'ID':
                        probe_col = i
                    if hc in ('Gene Symbol', 'Symbol', 'GENE_SYMBOL'):
                        gene_col = i
                continue
            if probe_col is not None and gene_col is not None and len(parts) > max(probe_col, gene_col):
                probe = parts[probe_col].strip('"')
                gene = parts[gene_col].strip('"')
                if '///' in gene:
                    gene = gene.split('///')[0].strip()
                mapping[probe] = gene if gene else probe
    print(f"  {gpl_id}: {len(mapping)} mappings")
    return mapping

print("Downloading GPL annotations...")
gpl570_map = parse_gpl(download_gpl('GPL570'), 'GPL570')
gpl6947_map = parse_gpl(download_gpl('GPL6947'), 'GPL6947')

def collapse_to_genes(expr_df, probe_map):
    gene_names = [probe_map.get(p, p) for p in expr_df.index]
    df = expr_df.copy()
    df.index = gene_names
    gene_means = df.mean(axis=1)
    df['_mean'] = gene_means
    df = df.sort_values('_mean', ascending=False)
    df = df[~df.index.duplicated(keep='first')]
    df = df.drop(columns=['_mean'])
    return df

print("Collapsing probes to genes...")
gse62254_genes = collapse_to_genes(gse62254, gpl570_map)
gse84437_genes = collapse_to_genes(gse84437, gpl6947_map)
print(f"  GSE62254: {gse62254_genes.shape[0]} genes x {gse62254_genes.shape[1]} samples")
print(f"  GSE84437: {gse84437_genes.shape[0]} genes x {gse84437_genes.shape[1]} samples")

nk_markers = ['NKG7', 'GNLY', 'GZMB', 'PRF1', 'KLRD1', 'NCAM1', 'KLRK1']
for name, df in [('GSE62254', gse62254_genes), ('GSE84437', gse84437_genes)]:
    found = [g for g in nk_markers if g in df.index]
    print(f"  {name} NK markers: {len(found)}/{len(nk_markers)}: {found}")

os.makedirs("data/processed/bulk", exist_ok=True)
gse62254_genes.T.to_csv("data/processed/bulk/gse62254_expression.tsv", sep="\t")
gse84437_genes.T.to_csv("data/processed/bulk/gse84437_expression.tsv", sep="\t")
print("Saved GEO expression matrices")

# =============================================================================
# STEP 3: STRING PPI
# =============================================================================
print()
print("=" * 60)
print("STEP 3: STRING PPI processing")
print("=" * 60)

string_df = pd.read_csv(
    "data/raw/prior_networks/9606.protein.physical.links.v12.0.txt.gz",
    sep=" ", compression="gzip"
)
print(f"Total STRING interactions: {len(string_df)}")

# Download protein aliases for ENSEMBL -> gene symbol mapping
try:
    alias_url = "https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz"
    with urllib.request.urlopen(alias_url) as resp:
        alias_content = resp.read()
    with open("data/raw/prior_networks/9606.protein.aliases.v12.0.txt.gz", "wb") as f:
        f.write(alias_content)

    alias_df = pd.read_csv(
        "data/raw/prior_networks/9606.protein.aliases.v12.0.txt.gz",
        sep="\t", compression="gzip", header=None,
        names=['protein_id', 'alias', 'source']
    )

    # Show available sources
    src_counts = alias_df['source'].value_counts()
    print(f"Available alias sources: {dict(src_counts.head(10))}")

    # Prefer HGNC symbol, fall back to any
    gene_aliases = alias_df[alias_df['source'].str.contains('symbol', case=False, na=False)]
    if len(gene_aliases) == 0:
        gene_aliases = alias_df[alias_df['source'].str.contains('Ensembl_HGNC', case=False, na=False)]
    if len(gene_aliases) == 0:
        gene_aliases = alias_df.drop_duplicates('protein_id')

    ensembl_to_gene = dict(zip(gene_aliases['protein_id'], gene_aliases['alias']))
    print(f"ENSEMBL->gene mapping: {len(ensembl_to_gene)} proteins")
except Exception as e:
    print(f"WARNING: alias download failed: {e}")
    ensembl_to_gene = {}

# Map and filter
import re
ensp_pat = re.compile(r'^9606\.(ENSP\d+)')
def map_id(pid):
    if pid in ensembl_to_gene:
        return ensembl_to_gene[pid]
    m = ensp_pat.match(str(pid))
    if m:
        return m.group(1)
    return pid

string_df['protein1_gene'] = string_df['protein1'].apply(map_id)
string_df['protein2_gene'] = string_df['protein2'].apply(map_id)
string_high = string_df[string_df['combined_score'] >= 700].copy()
print(f"High-confidence (>=700) edges: {len(string_high)} / {len(string_df)}")

os.makedirs("data/processed/graph", exist_ok=True)
string_high.to_csv("data/processed/graph/string_ppi_high_conf.tsv", sep="\t", index=False)
print("Saved STRING PPI")

print()
print("GEO + STRING processing complete!")
