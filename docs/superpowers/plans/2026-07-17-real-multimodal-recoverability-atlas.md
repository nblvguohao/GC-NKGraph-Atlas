# Real Multimodal Recoverability Atlas Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a no-synthetic-data pipeline that produces an auditable real-data, cross-mechanism recoverability atlas for four gastric-cancer immune-evasion cards.

**Architecture:** A manifest and reality guard admit only accessioned public human inputs. A card specification converts the four YAML cards into explicit comparisons. Shared transcriptomic and direct-modality adapters produce categorical evidence, which an aggregator turns into a card-by-layer figure and a pre-registered cross-mechanism verdict.

**Tech Stack:** Python 3.10; pandas, numpy, scipy, statsmodels, PyYAML, anndata/scanpy, matplotlib, pytest.

## Global Constraints

- Formal analyses and submission artifacts never use synthetic, mock, or demo data.
- Every included source has accession, URL, modality, species, sample count, timestamp, local file and SHA-256 in a manifest.
- Missing direct evidence is `not_measured`, never substituted with RNA, simulated values, or a literature summary.
- Confirmatory results require expected direction, BH-FDR < 0.05, adequate feature coverage, and the approved cross-cohort rule.
- Visium analyses report spot-module adjacency, never single-cell contact distance.
- Do not alter existing dirty worktree files except those explicitly named below.

---

### Task 1: Enforce a real-data manifest

**Files:**
- Create: `configs/recoverability_atlas/real_data_manifest.yaml`
- Create: `src/common/real_data.py`
- Modify: `src/interpretation/run_multicard_analysis.py`
- Test: `tests/test_real_data.py`

**Interfaces:**
- Produces `RealDataAsset`, `load_real_data_manifest(path) -> dict[str, RealDataAsset]`, and `assert_real_asset(asset, path) -> None`.
- Consumed by every downloader and analysis runner.

- [ ] **Step 1: Write the failing tests**

```python
def test_manifest_requires_accession_url_hash_and_human_species(tmp_path):
    path = tmp_path / "manifest.yaml"
    path.write_text("assets:\n  x: {accession: GSE1}\n", encoding="utf-8")
    with pytest.raises(ValueError, match="source_url"):
        load_real_data_manifest(path)

def test_guard_rejects_synthetic_mock_and_demo_paths():
    asset = RealDataAsset("GSE1", "https://example.org/x", "RNA",
                          "Homo sapiens", 2, "a" * 64, "x.tsv")
    with pytest.raises(ValueError, match="non-real"):
        assert_real_asset(asset, Path("synthetic.tsv"))
```

- [ ] **Step 2: Verify RED**

Run: `python -m pytest tests/test_real_data.py -q`  
Expected: FAIL because `src.common.real_data` does not exist.

- [ ] **Step 3: Implement the minimum interface**

```python
@dataclass(frozen=True)
class RealDataAsset:
    accession: str
    source_url: str
    modality: str
    species: str
    sample_count: int
    sha256: str
    local_path: str

def assert_real_asset(asset: RealDataAsset, path: Path) -> None:
    if any(x in str(path).lower() for x in ("synthetic", "mock", "demo")):
        raise ValueError("non-real input path is forbidden")
    if asset.species != "Homo sapiens" or len(asset.sha256) != 64:
        raise ValueError("asset does not satisfy real-data contract")
```

Remove the production synthetic fallback from `run_multicard_analysis.py`: its default must be `synthetic=False`, and a request for synthetic data must raise outside an explicitly test-only fixture path. Formal card outputs must be derived only through manifest-checked sources.

Seed fixed entries for TCGA-STAD, TCGA-LIHC, GSE62254, GSE84437,
GSE122401, MTBLS3303, and GSE251950. New remote entries begin
`status: pending_download`; their hash is filled only after retrieval.

- [ ] **Step 4: Verify GREEN**

Run: `python -m pytest tests/test_real_data.py -q`  
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add configs/recoverability_atlas/real_data_manifest.yaml src/common/real_data.py src/interpretation/run_multicard_analysis.py tests/test_real_data.py
git commit -m "feat: enforce real-data manifest contract"
```

### Task 2: Retrieve allowlisted public sources without fallback

**Files:**
- Create: `src/data/download_recoverability_sources.py`
- Modify: `configs/recoverability_atlas/real_data_manifest.yaml`
- Test: `tests/test_download_recoverability_sources.py`

**Interfaces:**
- Produces `resolve_download_url(accession) -> str` and `record_download(accession, path) -> DownloadRecord`.
- Supports `--asset ACCESSION --output-root PATH`, `--metadata-only`, and `--verify-all`.

- [ ] **Step 1: Write the failing tests**

```python
def test_only_allowlisted_accessions_have_urls():
    assert "GSE122401" in resolve_download_url("GSE122401")
    with pytest.raises(ValueError, match="not allowlisted"):
        resolve_download_url("GSE999999")

def test_record_requires_nonempty_download(tmp_path):
    with pytest.raises(FileNotFoundError):
        record_download("GSE122401", tmp_path / "missing.tsv")
```

- [ ] **Step 2: Verify RED**

Run: `python -m pytest tests/test_download_recoverability_sources.py -q`  
Expected: FAIL because the downloader module is absent.

- [ ] **Step 3: Implement retrieval**

Allowlist the official processed GSE122401 archive, MTBLS3303 public archive,
and GSE251950 supplementary archive. Stream only HTTPS/FTP content into
`data/external/recoverability/<accession>/`; reject HTTP errors, empty files,
unrecognized archive layouts, missing sample annotation, and unknown accessions.
Write the SHA-256 and retrieval timestamp only after validation. No exception
handler may create replacement data.

- [ ] **Step 4: Verify GREEN and metadata-only**

Run: `python -m pytest tests/test_download_recoverability_sources.py -q`  
Run: `python src/data/download_recoverability_sources.py --asset GSE122401 --metadata-only`  
Expected: tests PASS; metadata command prints the official URL and writes no result table.

- [ ] **Step 5: Commit**

```bash
git add src/data/download_recoverability_sources.py configs/recoverability_atlas/real_data_manifest.yaml tests/test_download_recoverability_sources.py
git commit -m "feat: add reproducible public recovery-data retrieval"
```

### Task 3: Make card-layer comparisons explicit

**Files:**
- Create: `configs/recoverability_atlas/card_layer_spec.yaml`
- Create: `src/interpretation/recoverability_spec.py`
- Test: `tests/test_recoverability_spec.py`

**Interfaces:**
- Produces `Comparison(card_id, comparison_id, layer, left_module, right_module, expected_sign, requires_purity_control)`.
- `load_comparisons(path) -> list[Comparison]` rejects `NEEDS_REVIEW` in any confirmatory row.

- [ ] **Step 1: Write the failing test**

```python
def test_confirmatory_comparison_requires_explicit_sign(tmp_path):
    spec = tmp_path / "spec.yaml"
    spec.write_text("comparisons: [{card_id: x, expected_sign: NEEDS_REVIEW}]", encoding="utf-8")
    with pytest.raises(ValueError, match="explicit expected_sign"):
        load_comparisons(spec)
```

- [ ] **Step 2: Verify RED**

Run: `python -m pytest tests/test_recoverability_spec.py -q`  
Expected: FAIL because the specification loader is absent.

- [ ] **Step 3: Implement specification**

Map each card's explicit downstream comparisons to
`upstream_driver`, `receptor_signaling`, or `downstream_effector`.
Keep cross-cell-type direction-calibration hypotheses as
`calibration_only: true`; they cannot count toward recovery. Validate every
named module against its original card YAML.

- [ ] **Step 4: Verify GREEN**

Run: `python -m pytest tests/test_recoverability_spec.py -q`  
Expected: PASS and all four registered cards load.

- [ ] **Step 5: Commit**

```bash
git add configs/recoverability_atlas/card_layer_spec.yaml src/interpretation/recoverability_spec.py tests/test_recoverability_spec.py
git commit -m "feat: formalize card-layer recovery comparisons"
```

### Task 4: Evaluate the four cards on existing real transcriptomes

**Files:**
- Create: `src/interpretation/run_recoverability_transcriptome.py`
- Modify: `src/interpretation/run_tgfb_card_recovery.py`
- Test: `tests/test_recoverability_transcriptome.py`

**Interfaces:**
- Produces `evaluate_transcriptomic_recovery(...) -> pd.DataFrame`.
- Output columns: `card_id, comparison_id, cohort, r, p, p_fdr, coverage, purity_r, purity_p, direction_ok, status`.

- [ ] **Step 1: Write the failing test**

```python
def test_recovered_needs_direction_fdr_coverage_and_two_cohorts():
    row = pd.Series({"direction_ok": True, "p_fdr": 0.01,
                     "coverage": 1.0, "concordant_cohorts": 2})
    assert assign_recovery_status(row) == "recovered"
    assert assign_recovery_status(row.assign(concordant_cohorts=1)) == "partially_recovered"
```

- [ ] **Step 2: Verify RED**

Run: `python -m pytest tests/test_recoverability_transcriptome.py -q`  
Expected: FAIL because the runner is absent.

- [ ] **Step 3: Implement shared analysis**

Reuse `mean_zscore` and `partial_corr` from
`src/topology/bulk_h3_purity_control.py`; calculate Pearson effects, BH-FDR
over each run, module coverage, and NK-lineage sensitivity for NK endpoints.
Call `assert_real_asset` before every read. Convert the existing TGF-beta
runner into a compatibility wrapper so no card-specific result logic is duplicated.

- [ ] **Step 4: Verify GREEN**

Run: `python -m pytest tests/test_recoverability_transcriptome.py tests/test_sst_config.py -q`  
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/interpretation/run_recoverability_transcriptome.py src/interpretation/run_tgfb_card_recovery.py tests/test_recoverability_transcriptome.py
git commit -m "feat: evaluate card recovery across real cohorts"
```

### Task 5: Add honest direct-modality adapters

**Files:**
- Create: `src/interpretation/recoverability_modalities.py`
- Create: `src/interpretation/run_recoverability_modalities.py`
- Test: `tests/test_recoverability_modalities.py`

**Interfaces:**
- `analyze_matched_rna_protein(rna, protein, pairs, asset) -> pd.DataFrame`
- `analyze_metabolite_contrast(table, feature_map, group_column, asset) -> pd.DataFrame`
- `analyze_visium_adjacency(adata, module_pairs, asset) -> pd.DataFrame`
- Every row includes `modality, direct_endpoint, status, not_measured_reason`.

- [ ] **Step 1: Write the failing test**

```python
def test_missing_protein_matrix_is_not_measured_not_negative():
    row = not_measured("protein", "MICA abundance", "public matrix unavailable")
    assert row["status"] == "not_measured"
    assert row["not_measured_reason"] == "public matrix unavailable"
```

- [ ] **Step 2: Verify RED**

Run: `python -m pytest tests/test_recoverability_modalities.py -q`  
Expected: FAIL because the modality module is absent.

- [ ] **Step 3: Implement direct evidence only**

Use Spearman RNA-protein concordance only for verified sample-matched matrices.
For MTBLS3303, compare identified serine, sphingolipid, and adenosine features
between real tumor and adjacent tissue with effect size and BH-FDR. For
GSE251950, score Visium spots and compute spatial module adjacency from supplied
coordinates. An absent valid matrix writes one `not_measured` row and terminates
that adapter; it cannot be replaced by a proxy.

- [ ] **Step 4: Verify GREEN**

Run: `python -m pytest tests/test_recoverability_modalities.py -q`  
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/interpretation/recoverability_modalities.py src/interpretation/run_recoverability_modalities.py tests/test_recoverability_modalities.py
git commit -m "feat: add direct multimodal recovery evidence"
```

### Task 6: Produce the atlas and gate the global claim

**Files:**
- Create: `src/interpretation/run_recoverability_atlas.py`
- Create: `src/figures/make_recoverability_atlas_figure.py`
- Test: `tests/test_recoverability_atlas.py`

**Interfaces:**
- Produces `recoverability_atlas.tsv`, `recoverability_cross_mechanism_verdict.json`, and `figS2_recoverability_atlas.{png,pdf}`.
- `cross_mechanism_verdict(evidence) -> str` returns only `cross_mechanism_pattern_supported` or `comparative_atlas_only`.

- [ ] **Step 1: Write the failing test**

```python
def test_global_pattern_needs_three_cards_two_cohorts_and_direct_evidence():
    evidence = pd.DataFrame([
      {"card_id": "a", "status": "recovered", "concordant_cohorts": 2, "direct_modality": True},
      {"card_id": "b", "status": "recovered", "concordant_cohorts": 2, "direct_modality": False},
      {"card_id": "c", "status": "recovered", "concordant_cohorts": 2, "direct_modality": False},
    ])
    assert cross_mechanism_verdict(evidence) == "cross_mechanism_pattern_supported"
```

- [ ] **Step 2: Verify RED**

Run: `python -m pytest tests/test_recoverability_atlas.py -q`  
Expected: FAIL because the atlas runner is absent.

- [ ] **Step 3: Implement aggregation and Figure S2**

Aggregate categorical status rather than incompatible effect sizes. Render cards as
rows and fixed biological layers as columns, with a distinct symbol for
`not_measured`, feature coverage, modality, and accession in supplementary
tables. Derive the verdict only from the documented 3-card/2-cohort/direct-modality gate.

- [ ] **Step 4: Verify GREEN**

Run: `python -m pytest tests/test_recoverability_atlas.py -q`  
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/interpretation/run_recoverability_atlas.py src/figures/make_recoverability_atlas_figure.py tests/test_recoverability_atlas.py
git commit -m "feat: summarize real multimodal recoverability atlas"
```

### Task 7: Download verified public inputs and create submission artifacts

**Files:**
- Modify: `configs/recoverability_atlas/real_data_manifest.yaml`
- Create: `submission_bundle_BiB/03_supplementary/tables/recoverability_atlas.tsv`
- Create: `submission_bundle_BiB/03_supplementary/tables/recoverability_source_manifest.tsv`
- Create: `submission_bundle_BiB/03_supplementary/tables/recoverability_transcriptome_per_cohort.tsv`
- Create: `submission_bundle_BiB/03_supplementary/tables/recoverability_direct_modality.tsv`
- Create: `submission_bundle_BiB/03_supplementary/tables/recoverability_cross_mechanism_verdict.json`
- Create: `submission_bundle_BiB/02_figures/figS2_recoverability_atlas.png`
- Create: `submission_bundle_BiB/02_figures/figS2_recoverability_atlas.pdf`

**Interfaces:** Consumes only manifest-verified assets and outputs tables that retain accession and SHA-256.

- [ ] **Step 1: Verify downloads**

Run: `python src/data/download_recoverability_sources.py --verify-all`  
Expected: included files are real human inputs with nonempty content and hashes; inaccessible data are recorded as `not_measured`.

- [ ] **Step 2: Run analyses**

Run: `python src/interpretation/run_recoverability_transcriptome.py --manifest configs/recoverability_atlas/real_data_manifest.yaml`  
Run: `python src/interpretation/run_recoverability_modalities.py --manifest configs/recoverability_atlas/real_data_manifest.yaml`  
Run: `python src/interpretation/run_recoverability_atlas.py --submission-root submission_bundle_BiB`  
Expected: all output rows have accession and one allowed status.

- [ ] **Step 3: Run contracts**

Run: `python -m pytest tests/test_real_data.py tests/test_download_recoverability_sources.py tests/test_recoverability_spec.py tests/test_recoverability_transcriptome.py tests/test_recoverability_modalities.py tests/test_recoverability_atlas.py -q`  
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add configs/recoverability_atlas submission_bundle_BiB/02_figures/figS2_recoverability_atlas.* submission_bundle_BiB/03_supplementary/tables/recoverability_*
git commit -m "data: add real multimodal recoverability results"
```

### Task 8: Gate manuscript language and verify submission

**Files:**
- Modify: `submission_bundle_BiB/01_manuscript/main_manuscript.md`
- Modify: `submission_bundle_BiB/01_manuscript/main.tex`
- Modify: `submission_bundle_BiB/03_supplementary/SUPPLEMENTARY_INDEX.md`
- Create: `tests/test_recoverability_submission_contract.py`

**Interfaces:** Consumes Task 7 verdict JSON. The claim gate is the only allowed source for universal versus comparative wording.

- [ ] **Step 1: Write the failing claim test**

```python
def test_universal_claim_requires_supported_verdict():
    verdict = json.loads(VERDICT.read_text())["verdict"]
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    if verdict != "cross_mechanism_pattern_supported":
        assert "general rule across mechanisms" not in manuscript
        assert "universal recoverability law" not in manuscript
```

- [ ] **Step 2: Verify RED**

Run: `python -m pytest tests/test_recoverability_submission_contract.py -q`  
Expected: FAIL because verdict/artifacts do not exist.

- [ ] **Step 3: Implement verdict-conditioned text**

Add manifest provenance, card-by-layer statuses, direct-modality scope, and Fig. S2
reference. If the global gate fails, report a comparative atlas and name unmeasured
layers; if it passes, use only the approved bounded cross-mechanism conclusion.
Do not add model-performance claims.

- [ ] **Step 4: Verify manuscript and PDF**

Run: `python -m pytest tests/test_recoverability_submission_contract.py tests/test_multiview_submission_contract.py -q`  
Run: `pdflatex -interaction=nonstopmode -halt-on-error main.tex` twice in `submission_bundle_BiB/01_manuscript`  
Expected: PASS; PDF references and Fig. S2 resolve.

- [ ] **Step 5: Full regression and commit**

Run: `python -m pytest -q`  
Run: `python src/interpretation/run_recoverability_atlas.py --verify-submission submission_bundle_BiB`  
Expected: all tests pass; every submission artifact traces to a real manifest entry; no synthetic/mock/demo origin appears.

```bash
git add submission_bundle_BiB/01_manuscript submission_bundle_BiB/03_supplementary/SUPPLEMENTARY_INDEX.md tests/test_recoverability_submission_contract.py
git commit -m "docs: report gated multimodal recovery atlas"
```

## Plan self-review

- Tasks 1--2 implement traceable real inputs and the no-synthetic gate.
- Tasks 3--4 implement pre-specified four-card transcriptomic evidence.
- Task 5 implements direct-modality limits; Task 6 implements the universal-claim gate and figure.
- Tasks 7--8 execute real runs, publish artifacts, gate manuscript wording, compile, and regression-test.
- All interfaces named by later tasks are defined by earlier tasks.
