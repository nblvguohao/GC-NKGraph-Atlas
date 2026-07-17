# GSE251950 official-download verification

**Verification date:** 2026-07-17 (Asia/Shanghai)  
**Accession:** GSE251950  
**Source:** <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE251nnn/GSE251950/suppl/GSE251950_RAW.tar>

## Retrieval and integrity

- The official GEO HTTPS endpoint returned `200 OK`, `Accept-Ranges: bytes`,
  and `Content-Length: 668221440`.
- The local archive is
  `data/external/recoverability/GSE251950/GSE251950_RAW.tar`.
- Final local size: **668221440 bytes**, exactly matching the official
  `Content-Length`.
- SHA-256: `d411d5aa3163a5c97a408800516a358113c3ca6b7a2c00849dbc111dd71ccd8e`.
- A prior interrupted resume had appended 3,144,368 bytes past the official
  content length. The local file was safely truncated to the server-declared
  byte length before hash calculation and archive inspection; no data were
  invented, substituted, or proxied.
- `tar -tf` completed with exit status 0 and found ten GEO sample archives.
  This verifies the outer TAR directory only; it is **not** evidence that all
  nested sample archives are intact. See the correction below.

## Archive evidence for direct Visium analysis

The outer archive contains ten sample archives (`GSM7990473` through
`GSM7990482`). Inspection of `GSM7990473_20_00331_LI_SING.tar.gz` found the
direct Visium inputs:

- `20_00331_LI_SING_filtered_feature_bc_matrix.h5`
- `20_00331_LI_SING_matrix.mtx.gz`
- `20_00331_LI_SING_barcodes.tsv.gz`
- `20_00331_LI_SING_features.tsv.gz`
- `20_00331_LI_SING_tissue_positions_list.csv`
- `20_00331_LI_SING_scalefactors_json.json`
- tissue images and aligned-fiducials image files.

These are expression and coordinate artifacts supplied by GEO, sufficient in
principle for a direct spot-level Visium analysis. The present report verifies
download and archive structure only; it does not claim that the data have
already been parsed or that any biological result has been obtained.

## Blockers

The outer archive can be retrieved and has the server-declared byte size, but
the nested-archive integrity correction below blocks full-cohort analysis.
Downstream analysis must still unpack each sample and compute spatial adjacency
from the supplied spot coordinates, never from RNA-only proxy data.

## Correction: nested-archive integrity failure

**Status: the local GSE251950 outer archive is unusable for full-cohort
analysis.** A complete-size outer TAR and a SHA-256 value are necessary but not
sufficient integrity checks when the TAR contains compressed per-sample
archives.

### Reproduction commands and results

The following commands were run against the local archive:

```powershell
tar -tf data/external/recoverability/GSE251950/GSE251950_RAW.tar
tar -xf data/external/recoverability/GSE251950/GSE251950_RAW.tar -C C:\tmp\gse251950_integrity <inner-member>
tar -tzf C:\tmp\gse251950_integrity\<inner-member>
```

Results:

- `GSM7990479_21_01252_LI_SING.tar.gz` is truncated and cannot provide a
  complete sample payload.
- `GSM7990474_21_00731_LI_SING.tar.gz`,
  `GSM7990476_21_00733_LI_SING.tar.gz`, and
  `GSM7990481_21_01675_LI_SING.tar.gz` emit `Damaged tar archive (bad header
  checksum)` warnings when their nested TAR structure is listed.
- `GSM7990473_20_00331_LI_SING.tar.gz` was the intact sample used for the
  structure-only inspection above. Its success does not rescue the incomplete
  cohort.

Accordingly, no spatial result, spot adjacency, or cohort-level claim may be
computed from the local outer archive. The archive remains recorded only as a
failed acquisition artifact, not as analysis input.

### Official per-sample recovery route (not yet retrieved)

The NCBI GEO server exposes independent official sample files, so recovery does
not require a non-official mirror or proxy. `curl.exe -I` returned `200 OK` for
these exact URLs:

- <https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7990nnn/GSM7990474/suppl/GSM7990474_21_00731_LI_SING.tar.gz>
  (`Content-Length: 99492794`)
- <https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7990nnn/GSM7990479/suppl/GSM7990479_21_01252_LI_SING.tar.gz>
  (`Content-Length: 30694361`)
- <https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7990nnn/GSM7990476/suppl/GSM7990476_21_00733_LI_SING.tar.gz>
  (`Content-Length: 117303692`)
- <https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7990nnn/GSM7990481/suppl/GSM7990481_21_01254_LI_SING.tar.gz>
  (`Content-Length: 40179023`)

No replacement file was downloaded in this task. Before any future use, each
independent official download must satisfy: (1) the HEAD `Content-Length`
matches its local byte count, (2) SHA-256 is recorded, (3) `tar -tzf` exits
cleanly without damaged-header or truncation diagnostics, and (4) its matrix,
barcodes, features, and tissue-position files are listed. Only after all ten
samples pass this check may a full-cohort spatial analysis proceed.

## Postscript: independently verified four-sample recovery

The failed outer archive remains unsuitable for full-cohort analysis. After the
above inspection, four of the affected samples were independently recovered
from the official per-GSM URLs and validated as described in
`gse251950-per-gsm-retrieval-report.md`. Those four archives are now recorded
with path, byte size, SHA-256, source URL, and `available` status in
`configs/recoverability_atlas/real_data_manifest.yaml`. They support only the
explicit exploratory four-sample subset analysis; they do not change the
full-cohort blocker documented above.
