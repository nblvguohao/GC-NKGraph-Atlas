# GSE251950 per-GSM recovery and validation

**Verification date:** 2026-07-17 (Asia/Shanghai)  
**Accession:** GSE251950  
**Scope:** recovery of only the four damaged nested archives identified in
`gse251950-download-report.md`. The original outer archive was not changed.

## Source and retrieval

Each replacement below was retrieved directly from its official NCBI GEO HTTPS
endpoint into `data/external/recoverability/GSE251950/per_gsm/`. No mirror,
proxy, synthetic data, or RNA-derived spatial substitute was used. A final
`curl.exe -I` request returned HTTP 200 and the listed server content length
for every endpoint.

| GSM | Official GEO URL | HTTP content length (bytes) | Local bytes | SHA-256 | Size match |
|---|---|---:|---:|---|---|
| GSM7990474 | https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7990nnn/GSM7990474/suppl/GSM7990474_21_00731_LI_SING.tar.gz | 99,492,794 | 99,492,794 | `449d36fafed541c48f91e18148201b5a0b40f41f5f739312d2c4e8abcf4294cd` | pass |
| GSM7990479 | https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7990nnn/GSM7990479/suppl/GSM7990479_21_01252_LI_SING.tar.gz | 30,694,361 | 30,694,361 | `cf2df49f33766709a4d5c96995dfe31e3e57fab5a998e7235dafde33a1f51291` | pass |
| GSM7990476 | https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7990nnn/GSM7990476/suppl/GSM7990476_21_00733_LI_SING.tar.gz | 117,303,692 | 117,303,692 | `21f34c1c8ef5f454071c9d243a8b1b05985413b66bb5ac8fd68b1cbfd66303f6` | pass |
| GSM7990481 | https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7990nnn/GSM7990481/suppl/GSM7990481_21_01254_LI_SING.tar.gz | 40,179,023 | 40,179,023 | `7614be8a3fb94d465e5bd0adb6d5cd45ef1efeb170db10a02294f3f7b4ada78b` | pass |

## Archive and Visium-input validation

For every recovered file, `tar -tzf <file>` exited with status 0 and emitted no
warning, error, truncation, checksum, or damaged-header diagnostic. Its
listing contained the following direct Visium files. This verifies archive
integrity and the presence of expression plus supplied spot coordinates; it
does **not** perform spatial analysis or claim a biological result.

| GSM | Filtered matrix H5 | Matrix Market / barcodes / features | Tissue coordinates | Scale factors | TAR validation |
|---|---|---|---|---|---|
| GSM7990474 | `21_00731_LI_SING/21_00731_LI_SING_filtered_feature_bc_matrix.h5` | `21_00731_LI_SING_matrix.mtx.gz`; `21_00731_LI_SING_barcodes.tsv.gz`; `21_00731_LI_SING_features.tsv.gz` | `21_00731_LI_SING_tissue_positions_list.csv` | `21_00731_LI_SING_scalefactors_json.json` | pass |
| GSM7990479 | `21_01252_LI_SING/21_01252_LI_SING_filtered_feature_bc_matrix.h5` | `21_01252_LI_SING_matrix.mtx.gz`; `21_01252_LI_SING_barcodes.tsv.gz`; `21_01252_LI_SING_features.tsv.gz` | `21_01252_LI_SING_tissue_positions_list.csv` | `21_01252_LI_SING_scalefactors_json.json` | pass |
| GSM7990476 | `21_00733_LI_SING/21_00733_LI_SING_filtered_feature_bc_matrix.h5` | `21_00733_LI_SING_matrix.mtx.gz`; `21_00733_LI_SING_barcodes.tsv.gz`; `21_00733_LI_SING_features.tsv.gz` | `21_00733_LI_SING_tissue_positions_list.csv` | `21_00733_LI_SING_scalefactors_json.json` | pass |
| GSM7990481 | `21_01254_LI_SING/21_01254_LI_SING_filtered_feature_bc_matrix.h5` | `21_01254_LI_SING_matrix.mtx.gz`; `21_01254_LI_SING_barcodes.tsv.gz`; `21_01254_LI_SING_features.tsv.gz` | `21_01254_LI_SING_tissue_positions_list.csv` | `21_01254_LI_SING_scalefactors_json.json` | pass |

## Reproducible checks

The verification procedure for each table row was:

```powershell
curl.exe -sSIL <official-url>                 # HTTP 200 and Content-Length
Get-Item <local-file>                         # exact byte count
Get-FileHash -Algorithm SHA256 <local-file>   # recorded digest
tar -tzf <local-file>                         # exit 0; no diagnostics
```

The nested files are independently valid recovery artifacts. They do not by
themselves validate the remaining six samples in the original outer TAR;
full-cohort spatial analysis must retain that limitation until every sample is
independently verified or recovered.

## Status

**Pass for the four specified recovery files.** These files are suitable for
subsequent unpacking and direct coordinate-based Visium processing. This task
did not unpack them, construct adjacency, update any manifest, modify code,
or create analysis results.
