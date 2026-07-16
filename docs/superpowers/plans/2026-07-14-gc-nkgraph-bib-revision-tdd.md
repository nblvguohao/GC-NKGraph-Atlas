# GC-NKGraph-Atlas BiB Revision and Elevation TDD Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 将 GC-NKGraph-Atlas 从“已有结果的计算研究”升级为一篇统计口径统一、机制边界清楚、可复用性和实验可验证性更强的 Briefings in Bioinformatics Problem Solving Protocol 稿件。

**Architecture:** 先锁定证据边界，再以测试驱动方式重算单细胞相关、图模型消融和候选靶点排名；所有结果写入版本化 TSV/JSON，正文只从这些结果表读取数字。最后做稿件一致性、可复现性和 BiB 投稿格式检查。核心叙事从“重建完整机制”调整为“测量转录组对物理/代谢机制的可观测范围”。

**Tech Stack:** Python 3.8+；pytest；pandas/numpy/scipy/statsmodels；scanpy/scvi-tools（如服务器环境已有）；PyTorch Geometric；Markdown/LaTeX；TSV/JSON；现有 `src/`、`tests/`、`submission_bundle_BiB/` 结构。

## Global Constraints

- 不把转录代理指标写成真实鞘磷脂含量、代谢通量或膜突起物理表型。
- 单细胞相关不得把细胞数当作独立生物学重复；正文必须同时给出样本级或混合模型结果。
- H2 的样本校正结果若仍为 `P>=0.05`，正文只能写“未获得稳健证据/效应极小”。
- 单细胞 H3 只有在表达匹配打分、测序深度校正和模块特异性检验均通过后，才能升级为主要结论；否则只保留为探索性结果。
- GNN 不得宣称优于 LightGBM/XGBoost；默认表述为“预测性能相当，提供机制结构化表示”。
- 37 个基因默认称为“机制交叉候选短名单”，只有经过恶性细胞特异性和独立证据筛选后，才可称为肿瘤内在候选。
- 所有正文数字必须来自 `submission_bundle_BiB/03_supplementary/tables/` 的最终结果表。
- 每一个新增分析遵循 RED → GREEN → REFACTOR：先写会失败的验收测试，确认失败原因正确，再改脚本/结果，最后回归全套测试。
- 投稿前必须解决乱码、ORCID、基金信息、章节编号、PDF 重编译和补充方法缺失问题。

## File Map

- `src/topology/sst_axis_validation.py`: H1–H5 的 bulk、单细胞和样本级验证。
- `src/topology/pseudoreplication_correction.py`: 单细胞到样本级的聚合/随机效应校正。
- `src/topology/h3_scoring_method_diagnostic.py`: H3 的打分方法、测序深度和模块置换诊断。
- `src/topology/h3_leave_one_sample_out.py`: H3 样本敏感性分析。
- `src/baselines/run_model_comparison.py`: GNN 与传统模型的同折比较。
- `src/interpretation/prioritize_targets.py`: 候选靶点综合评分。
- `src/interpretation/split_target_lists.py`: 肿瘤内在候选与 NK 读出面板拆分。
- `src/graph_construction/build_heterograph.py`: 异质图和机制边构建。
- `src/models/gc_nkgraph_atlas.py`: 图模型训练和评估。
- `tests/`: 现有单元测试；新增统计契约、结果表和稿件一致性测试。
- `submission_bundle_BiB/03_supplementary/tables/`: 所有最终 TSV/JSON 结果。
- `submission_bundle_BiB/01_manuscript/main_manuscript.md`: 投稿正文。
- `submission_bundle_BiB/01_manuscript/main.tex`: LaTeX 正文。
- `submission_bundle_BiB/01_manuscript/BiB_submission_checklist.md`: 投稿门槛清单。

---

## Task 1: 建立“结果表是唯一事实源”的统计契约

**Files:**
- Create: `tests/test_result_contracts.py`
- Create: `tests/fixtures/result_contract_fixture.tsv`
- Modify: `submission_bundle_BiB/03_supplementary/tables/README.md`（若不存在则创建）

**Interfaces:**
- Produces test helpers `load_result_table(path)`、`assert_required_columns(df, columns)` 和 `assert_statistical_direction(row)`。

- [ ] **Step 1: Write the failing tests**

```python
def test_h2_uses_sample_corrected_p_value():
    row = load_result_table("submission_bundle_BiB/03_supplementary/tables/sst_axis_pseudoreplication_corrected.tsv")
    h2 = row[(row.hypothesis == "H2") & (row.resolution == "single-cell NK")].iloc[0]
    assert h2.corrected_p < 0.05, "current manuscript claim must fail until corrected evidence is regenerated"

def test_all_headline_results_have_corrected_or_bulk_p_values():
    table = load_result_table("submission_bundle_BiB/03_supplementary/tables/sst_axis_positive_control_recovery.tsv")
    assert {"hypothesis", "resolution", "r", "p", "p_fdr", "method"} <= set(table.columns)
```

- [ ] **Step 2: Run the tests to verify RED**

Run: `pytest tests/test_result_contracts.py -q`

Expected: H2 test fails because the current corrected result is `P=0.195`; this is the intended red state proving the old manuscript claim cannot pass the new contract.

- [ ] **Step 3: Define the minimal contract**

Replace the deliberately failing assertion with the final contract: H2 is allowed to pass only if the manuscript language is non-significant/weak when `corrected_p>=0.05`; the test should inspect a machine-readable `claim_level` field rather than force a positive result.

- [ ] **Step 4: Run GREEN and full regression**

Run: `pytest tests/test_result_contracts.py tests/ -q`

Expected: all existing tests pass and the new contract rejects any result table lacking corrected statistics or claim-level labels.

- [ ] **Step 5: Commit**

```bash
git add tests/test_result_contracts.py tests/fixtures/result_contract_fixture.tsv submission_bundle_BiB/03_supplementary/tables/README.md
git commit -m "test: add result-table statistical contracts"
```

---

## Task 2: 重算 H2，消除单细胞伪重复

**Files:**
- Modify: `src/topology/pseudoreplication_correction.py`
- Modify: `src/topology/sst_axis_validation.py`
- Test: `tests/test_pseudoreplication_correction.py`
- Regenerate: `submission_bundle_BiB/03_supplementary/tables/sst_axis_pseudoreplication_corrected.tsv`
- Regenerate: `submission_bundle_BiB/03_supplementary/tables/sst_axis_positive_control_recovery.tsv`

**Interfaces:**
- `correct_cell_level_association(df, x, y, sample_col, method="sample_meta") -> dict`
- Return keys: `naive_r`, `naive_p`, `corrected_r`, `corrected_p`, `ci_lower`, `ci_upper`, `k_samples`, `method`。

- [ ] **Step 1: Write the failing test**

```python
def test_cell_level_h2_does_not_use_8310_as_effective_n():
    result = correct_cell_level_association(h2_fixture, "sm_balance", "protrusion", "sample_id")
    assert result["k_samples"] == 9
    assert result["corrected_n_effective"] < result["naive_n"]
    assert result["method"] in {"sample_mean_meta", "mixed_effects"}
```

- [ ] **Step 2: Run RED**

Run: `pytest tests/test_pseudoreplication_correction.py::test_cell_level_h2_does_not_use_8310_as_effective_n -q`

Expected: FAIL if the function still reports cell-level P values or omits the effective sample count.

- [ ] **Step 3: Implement the minimal correction**

Aggregate module scores within each `sample_id`, calculate per-sample correlations or a prespecified mixed-effects association, and combine nine sample estimates with DerSimonian–Laird or an explicitly documented alternative. Preserve naive statistics only as transparency columns.

- [ ] **Step 4: Verify GREEN**

Run the server analysis and then:

```bash
pytest tests/test_pseudoreplication_correction.py tests/test_result_contracts.py -q
```

Acceptance:

- H2 corrected `r`, `P`, 95% CI and heterogeneity are present;
- if `P>=0.05`, H2 is labeled `NOT_ROBUST`;
- no abstract or main-text claim uses naive H2 P values.

- [ ] **Step 5: Commit**

```bash
git add src/topology/pseudoreplication_correction.py src/topology/sst_axis_validation.py tests/test_pseudoreplication_correction.py submission_bundle_BiB/03_supplementary/tables/
git commit -m "analysis: correct single-cell metabolic association for sample dependence"
```

---

## Task 3: 重做 H3 的技术混杂和模块特异性验证

**Files:**
- Modify: `src/topology/h3_scoring_method_diagnostic.py`
- Modify: `src/topology/h3_leave_one_sample_out.py`
- Test: `tests/test_h3_diagnostic_contract.py`
- Regenerate: `submission_bundle_BiB/03_supplementary/tables/h3_scoring_method_diagnostic.tsv`
- Regenerate: `submission_bundle_BiB/03_supplementary/tables/h3_scoring_method_diagnostic_summary.md`
- Regenerate: `submission_bundle_BiB/03_supplementary/tables/h3_leave_one_sample_out.tsv`

**Interfaces:**
- `score_modules_expression_matched(adata, module_genes, control_genes, n_bins=25) -> np.ndarray`
- `run_h3_diagnostics(adata, protrusion_genes, cytotoxicity_genes, covariates) -> dict`
- Required outputs: raw, expression-matched, count-depth residualized, sample-corrected, matched-permutation empirical P。

- [ ] **Step 1: Write the failing tests**

```python
def test_h3_requires_specificity_and_technical_control():
    result = run_h3_diagnostics(h3_fixture, PROTRUSION, CYTOTOXICITY, ["total_counts", "n_genes"])
    assert result["expression_matched_empirical_p"] < 0.05
    assert result["residualized_r"] >= 0.20

def test_leave_one_sample_out_reports_all_nine_samples():
    table = run_leave_one_sample_out(h3_fixture, "sample_id")
    assert table["excluded_sample"].nunique() == 9
```

- [ ] **Step 2: Run RED**

Run: `pytest tests/test_h3_diagnostic_contract.py -q`

Expected: the first test fails against the current diagnostic, because the matched permutation result does not exceed the null and residualized `r` is small. This is the intended scientific gate, not a software bug.

- [ ] **Step 3: Implement the minimal analysis**

Use expression-matched controls, regress module scores on `total_counts` and `n_genes`, aggregate at sample level, and run a module-size/expression-detection-matched permutation. Do not tune the gene universe to make the observed correlation pass.

- [ ] **Step 4: Verify GREEN and decide claim level**

Run:

```bash
pytest tests/test_h3_diagnostic_contract.py -q
```

If the gate still fails, retain H3 bulk as the primary result and label scRNA H3 exploratory. If it passes, report the corrected effect size, CI, empirical P and exact scoring method.

- [ ] **Step 5: Commit**

```bash
git add src/topology/h3_scoring_method_diagnostic.py src/topology/h3_leave_one_sample_out.py tests/test_h3_diagnostic_contract.py submission_bundle_BiB/03_supplementary/tables/h3_*
git commit -m "analysis: validate effector coupling against scoring and depth confounding"
```

---

## Task 4: 统一胃癌 scRNA 与外部队列的效应层证据

**Files:**
- Modify: `src/preprocessing/run_geo_external_validation.py`
- Modify: `src/topology/sst_axis_gastric_extension.py`
- Test: `tests/test_external_validation_contract.py`
- Regenerate: `submission_bundle_BiB/03_supplementary/tables/external_validation_results.tsv`
- Regenerate: `submission_bundle_BiB/03_supplementary/tables/sst_axis_scrna_by_tissue.tsv`

**Interfaces:**
- `validate_cohort(expr, platform, cohort_id, resolution="sample") -> dict`
- Required fields: `cohort_id`, `n_samples`, `mapped_gene_count`, `nk_marker_coverage`, `effect`, `ci_lower`, `ci_upper`, `p_adjusted`, `scoring_method`。

- [ ] **Step 1: Write the failing tests**

```python
def test_external_validation_has_platform_mapping_audit():
    table = load_external_results()
    assert {"platform", "mapped_gene_count", "nk_marker_coverage", "scoring_method"} <= set(table.columns)

def test_gastric_sc_h3_uses_same_scoring_rule_as_liver_sc_h3():
    assert gastric_sc_result["scoring_method"] == liver_sc_result["scoring_method"]
```

- [ ] **Step 2: Run RED**

Run: `pytest tests/test_external_validation_contract.py -q`

Expected: FAIL if the current external table lacks mapping/scoring metadata or if gastric and liver single-cell analyses use different scoring rules.

- [ ] **Step 3: Implement minimal harmonization**

Apply the same expression normalization, gene mapping, module scoring, sample-aware statistics and FDR rule to TCGA-LIHC scRNA, gastric scRNA, GSE62254 and GSE84437. Keep probe-to-symbol mapping logs in the table.

- [ ] **Step 4: Verify GREEN**

Acceptance:

- GSE62254 and GSE84437 retain nonzero NK-marker coverage;
- all cohorts have mapped-gene and scoring metadata;
- gastric scRNA is not described as confirmatory unless it passes Task 3’s diagnostic gate;
- external bulk replication is reported at the cohort/sample level.

- [ ] **Step 5: Commit**

```bash
git add src/preprocessing/run_geo_external_validation.py src/topology/sst_axis_gastric_extension.py tests/test_external_validation_contract.py submission_bundle_BiB/03_supplementary/tables/
git commit -m "analysis: harmonize gastric validation and platform mapping audits"
```

---

## Task 5: 证明异质图的结构价值，而不是只证明它改变了图

**Files:**
- Modify: `src/graph_construction/build_heterograph.py`
- Modify: `src/models/gc_nkgraph_atlas.py`
- Modify: `src/baselines/run_model_comparison.py`
- Test: `tests/test_graph_ablation_contract.py`
- Regenerate: `submission_bundle_BiB/03_supplementary/tables/ablation_results.tsv`
- Regenerate: `submission_bundle_BiB/03_supplementary/tables/t17_edge_external_value.tsv`

**Interfaces:**
- `build_graph(edge_types: list[str], seed: int) -> HeteroData`
- `evaluate_variant(variant_name: str, train_cohort: str, test_cohort: str, seed: int) -> dict`
- Required metrics: MCC, AUROC, bootstrap CI, modularity, embedding coupling, edge count。

- [ ] **Step 1: Write the failing tests**

```python
def test_ablation_has_predictive_and_structural_metrics():
    table = load_ablation_results()
    assert {"variant", "mcc", "auroc", "modularity", "embedding_h1", "embedding_h2"} <= set(table.columns)

def test_cross_cohort_ablation_has_uncertainty_interval():
    row = load_transfer_results().iloc[0]
    assert row.ci_lower <= row.delta_mcc <= row.ci_upper
```

- [ ] **Step 2: Run RED**

Run: `pytest tests/test_graph_ablation_contract.py -q`

Expected: FAIL if ablation reports only modularity/embedding changes or lacks held-out uncertainty intervals.

- [ ] **Step 3: Implement minimal complete ablation**

Compare FULL, `-metabolic_crosstalk`, `-SST`, PPI/LR/TF-only and tabular baselines under identical folds, seeds, label definitions and preprocessing. Separate structural effects from downstream predictive effects.

- [ ] **Step 4: Verify GREEN**

The manuscript may claim mechanism-edge value only if it reports both:

- structural change in the embedding; and
- held-out predictive change with CI and paired test.

If predictive delta remains null, write that the edge contributes structured representation but not verified predictive improvement.

- [ ] **Step 5: Commit**

```bash
git add src/graph_construction/build_heterograph.py src/models/gc_nkgraph_atlas.py src/baselines/run_model_comparison.py tests/test_graph_ablation_contract.py submission_bundle_BiB/03_supplementary/tables/ablation_results.tsv submission_bundle_BiB/03_supplementary/tables/t17_edge_external_value.tsv
git commit -m "analysis: separate graph structural and predictive ablation value"
```

---

## Task 6: 重建候选靶点，严格拆分肿瘤内在候选和 NK 读出

**Files:**
- Modify: `src/interpretation/prioritize_targets.py`
- Modify: `src/interpretation/split_target_lists.py`
- Test: `tests/test_target_list_contract.py`
- Regenerate: `submission_bundle_BiB/03_supplementary/tables/tumor_intrinsic_candidates.tsv`
- Regenerate: `submission_bundle_BiB/03_supplementary/tables/axis_confirmation_panel.tsv`
- Regenerate: `submission_bundle_BiB/03_supplementary/tables/trivial_baseline_comparison.tsv`

**Interfaces:**
- `prioritize_targets(evidence, tumor_specificity_threshold=0.10, nk_module_penalty=0.50) -> DataFrame`
- `split_target_lists(candidates) -> tuple[DataFrame, DataFrame]`
- Candidate rows must include `candidate_class`, `tumor_specificity_log2`, `nk_module_overlap`, `depmap_evidence`, `nk_state_de_evidence`, `assay`。

- [ ] **Step 1: Write the failing tests**

```python
def test_tumor_intrinsic_list_excludes_negative_or_zero_specificity():
    candidates, readout = split_target_lists(load_candidates())
    assert (candidates.tumor_specificity_log2 > 0).all()
    assert not set(candidates.gene) & {"NKG7", "PRF1", "GZMB", "GNLY"}

def test_readout_panel_is_not_presented_as_target_list():
    candidates, readout = split_target_lists(load_candidates())
    assert set(readout.candidate_class) == {"axis_readout"}
    assert set(candidates.candidate_class) == {"tumor_intrinsic_candidate"}
```

- [ ] **Step 2: Run RED**

Run: `pytest tests/test_target_list_contract.py -q`

Expected: FAIL if NK-side genes remain in the tumor-intrinsic table or if candidate class labels are absent.

- [ ] **Step 3: Implement minimal de-circularization**

Use signed malignant-cell specificity, penalize NK-side module overlap, require either tumor-serine/SM membership or a graph-neighborhood rationale, and retain DepMap/NK-state-DE evidence as independent dimensions. Do not use axis membership alone as a target criterion.

- [ ] **Step 4: Verify GREEN**

Acceptance:

- top-ranked list is a short, evidence-tiered list rather than an unexplained n=37 list;
- PHGDH, PSAT1, PSPH, SGMS2 and SMPD1/3 have explicit evidence and assay rationale;
- every candidate is classified as target, readout, or low-confidence mechanism intersection;
- trivial-baseline overlap and rank changes are reported.

- [ ] **Step 5: Commit**

```bash
git add src/interpretation/prioritize_targets.py src/interpretation/split_target_lists.py tests/test_target_list_contract.py submission_bundle_BiB/03_supplementary/tables/
git commit -m "analysis: enforce tumor-intrinsic candidate and NK-readout separation"
```

---

## Task 7: 增加第二张机制卡，验证框架可复用性

**Files:**
- Create/modify: `configs/mechanism_cards/adenosine_nk_suppression.yaml`
- Modify: `configs/mechanism_cards/registry.yaml`
- Test: `tests/test_mechanism_card_reuse.py`
- Create: `submission_bundle_BiB/03_supplementary/tables/mechanism_card_second_example.tsv`
- Modify: `submission_bundle_BiB/03_supplementary/SUPPLEMENTARY_INDEX.md`

**Interfaces:**
- `load_mechanism_card(card_path) -> MechanismCard`
- `run_card_smoke_test(card: MechanismCard, synthetic: bool = True) -> dict`
- Required output: module validation, graph-edge declaration, label definition, target-assay mapping and pass/fail status。

- [ ] **Step 1: Write the failing tests**

```python
def test_registry_contains_two_executable_cards():
    cards = load_registry()
    assert len(cards) >= 2
    for card in cards:
        result = run_card_smoke_test(card, synthetic=True)
        assert result["status"] == "PASS"

def test_second_card_does_not_change_core_pipeline_code_path():
    assert run_card_smoke_test(load_card("adenosine_nk_suppression"))["core_pipeline_modified"] is False
```

- [ ] **Step 2: Run RED**

Run: `pytest tests/test_mechanism_card_reuse.py -q`

Expected: FAIL until the second card has a complete schema and smoke-test output.

- [ ] **Step 3: Implement the minimal second-card example**

Use the existing registered adenosine or TGF-beta card. Provide only enough biological modules, cell-type attribution, expected direction, graph integration and assay recommendation for a synthetic end-to-end run. Do not claim biological validation of the second mechanism.

- [ ] **Step 4: Verify GREEN**

Run: `pytest tests/test_mechanism_card_reuse.py tests/ -q` and `python src/pipeline.py --synthetic`.

Acceptance: two cards pass the same loader/schema/synthetic pipeline without modifying core analysis code.

- [ ] **Step 5: Commit**

```bash
git add configs/mechanism_cards tests/test_mechanism_card_reuse.py submission_bundle_BiB/03_supplementary/
git commit -m "feat: demonstrate mechanism-card reuse with a second smoke-test card"
```

---

## Task 8: 建立正文—结果表—图表一致性测试

**Files:**
- Create: `tests/test_manuscript_consistency.py`
- Modify: `manuscript/main_manuscript.md`
- Modify: `submission_bundle_BiB/01_manuscript/main_manuscript.md`
- Modify: `submission_bundle_BiB/01_manuscript/main.tex`
- Modify: `submission_bundle_BiB/03_supplementary/SUPPLEMENTARY_INDEX.md`

**Interfaces:**
- `extract_numeric_claims(markdown_path) -> list[dict]`
- `validate_claim_against_tables(claim, table_dir) -> bool`
- `validate_section_numbering(markdown_path) -> list[str]`

- [ ] **Step 1: Write the failing tests**

```python
def test_no_unresolved_section_numbers():
    assert validate_section_numbering("submission_bundle_BiB/01_manuscript/main_manuscript.md") == []

def test_h2_claim_matches_corrected_result_table():
    claims = extract_numeric_claims("submission_bundle_BiB/01_manuscript/main_manuscript.md")
    h2 = next(c for c in claims if c["id"] == "H2_single_cell")
    assert validate_claim_against_tables(h2, "submission_bundle_BiB/03_supplementary/tables")

def test_manuscript_contains_required_boundary_language():
    text = Path("submission_bundle_BiB/01_manuscript/main_manuscript.md").read_text(encoding="utf-8")
    assert "does not substitute for physical membrane topology" in text
    assert "mechanism-intersecting candidate shortlist" in text
```

- [ ] **Step 2: Run RED**

Run: `pytest tests/test_manuscript_consistency.py -q`

Expected: failures identify the missing `3.6`, stale H2/H3 language, and missing candidate-boundary wording.

- [ ] **Step 3: Rewrite the manuscript minimally**

Update abstract, key points, Results 3.2–3.7, Discussion 4.1–4.4 and Conclusion from final TSV values. Add a short Supplementary Methods section covering label thresholds, module lists, QC, scoring, pseudoreplication correction and paired model tests.

- [ ] **Step 4: Verify GREEN**

Run:

```bash
pytest tests/test_manuscript_consistency.py tests/ -q
```

Acceptance:

- no section-number gaps;
- every headline statistic matches a table within a documented tolerance;
- abstract, key points, discussion and conclusion use the same claim level;
- no “physical topology prediction” wording remains.

- [ ] **Step 5: Commit**

```bash
git add tests/test_manuscript_consistency.py manuscript/ submission_bundle_BiB/01_manuscript/ submission_bundle_BiB/03_supplementary/SUPPLEMENTARY_INDEX.md
git commit -m "docs: synchronize manuscript claims with corrected evidence"
```

---

## Task 9: 修复投稿包、编码和作者/基金信息

**Files:**
- Modify: `作者信息和基金.txt`
- Modify: `submission_bundle_BiB/01_manuscript/main_manuscript.md`
- Modify: `submission_bundle_BiB/01_manuscript/main.tex`
- Modify: `submission_bundle_BiB/01_manuscript/BiB_submission_checklist.md`
- Create: `tests/test_submission_package.py`

**Interfaces:**
- `scan_for_mojibake(paths: list[Path]) -> list[str]`
- `validate_front_matter(path) -> dict`
- Required fields: authors, affiliations, corresponding emails, ORCID status, funding, competing interests, ethics, data availability。

- [ ] **Step 1: Write the failing tests**

```python
def test_submission_package_has_no_mojibake():
    assert scan_for_mojibake(SUBMISSION_TEXT_FILES) == []

def test_author_and_funding_fields_are_complete():
    fields = validate_front_matter("submission_bundle_BiB/01_manuscript/main_manuscript.md")
    assert all(fields[name] for name in ["authors", "affiliations", "corresponding_emails", "funding"])
    assert fields["orcid_placeholders"] == 0
```

- [ ] **Step 2: Run RED**

Run: `pytest tests/test_submission_package.py -q`

Expected: failures expose mojibake, placeholder ORCID values and malformed funding text.

- [ ] **Step 3: Implement the minimal cleanup**

Convert all manuscript/source files to UTF-8, verify Chinese names against author records, add real ORCID values, normalize grant numbers and punctuation, remove temporary editorial notes, and synchronize Markdown/LaTeX front matter.

- [ ] **Step 4: Verify GREEN**

Run:

```bash
pytest tests/test_submission_package.py -q
rg -n "鈥|脳|螖|0000-0000|TODO|TBD|3\.6" manuscript submission_bundle_BiB 作者信息和基金.txt
```

Acceptance: no encoding artifacts or placeholders remain; `3.6` is either restored or intentionally renumbered.

- [ ] **Step 5: Commit**

```bash
git add 作者信息和基金.txt tests/test_submission_package.py manuscript/ submission_bundle_BiB/01_manuscript/
git commit -m "chore: clean submission encoding and front matter"
```

---

## Task 10: 最终验证和 BiB 投前门槛

**Files:**
- Modify: `submission_bundle_BiB/01_manuscript/BiB_submission_checklist.md`
- Modify: `submission_bundle_BiB/00_SUBMISSION_GUIDE.md`
- Verify: `submission_bundle_BiB/01_manuscript/main.pdf`
- Verify: `README.md`, `requirements.txt`, `environment.yml`, `Dockerfile`

- [ ] **Step 1: Write the failing release-gate tests**

```python
def test_release_gate():
    assert manuscript_word_count("submission_bundle_BiB/01_manuscript/main_manuscript.md") <= 5000
    assert all_required_figures_present("submission_bundle_BiB/02_figures")
    assert all_required_tables_present("submission_bundle_BiB/03_supplementary/tables")
    assert synthetic_pipeline_exit_code() == 0
```

- [ ] **Step 2: Run RED**

Run: `pytest tests/test_submission_package.py tests/test_manuscript_consistency.py -q`.

Expected: the word-count or PDF freshness gate may fail until the manuscript is compressed and recompiled.

- [ ] **Step 3: Implement release cleanup**

Compress the main text to the selected BiB article type, recompile LaTeX with two passes, confirm zero undefined references, compare all PDF figures visually with source PNG/PDF, update the supplementary index, and confirm public repository state matches the submitted package.

- [ ] **Step 4: Verify GREEN**

Run:

```bash
pytest tests/ -q
python src/pipeline.py --synthetic
git status --short
```

Release gate:

- all tests pass;
- synthetic pipeline exits 0;
- manuscript and supplementary tables are synchronized;
- PDF is rebuilt from current TeX;
- ORCID/funding/author details are complete;
- no unresolved encoding or section-numbering issue;
- every major claim has one primary table/figure and one explicit limitation.

- [ ] **Step 5: Commit**

```bash
git add submission_bundle_BiB/ README.md tests/
git commit -m "release: finalize BiB revision package"
```

---

## Recommended Execution Order

1. Task 1 — statistical result contracts.
2. Task 2 — H2 pseudoreplication correction.
3. Task 3 — H3 scoring/depth/module-specificity diagnostic.
4. Task 4 — gastric cohort harmonization.
5. Task 5 — graph structural versus predictive ablation.
6. Task 6 — target-list de-circularization.
7. Task 7 — second mechanism-card smoke test.
8. Task 8 — manuscript synchronization.
9. Task 9 — encoding and front matter.
10. Task 10 — release gate.

## Decision Gates

- **Gate A — biology:** If corrected H2 remains non-significant, preserve the negative/weak metabolic result; do not add experiments merely to force significance.
- **Gate B — single-cell H3:** If expression-matched and depth-adjusted H3 remains below the matched null, remove “independent single-cell recovery” from the headline claims.
- **Gate C — graph value:** If the mechanism edge changes embeddings but not held-out prediction, frame this as representation value, not predictive superiority.
- **Gate D — targets:** If tumor-specificity remains near zero, reduce the list to a ranked hypothesis panel and recommend validation rather than therapeutic targeting.
- **Gate E — BiB fit:** If the manuscript remains a 7,000+ word single-dataset analysis without a second card or clear protocol generalization, consider a methods-oriented alternative journal or contact the BiB editorial office before submission.

## Final Self-Review Checklist

- [ ] Every required experiment has a failing test before implementation.
- [ ] H2 uses sample-aware statistics.
- [ ] H3 reports scoring-method and technical-confound diagnostics.
- [ ] Gastric validation uses harmonized scoring and platform mapping.
- [ ] Graph ablation separates structure from prediction.
- [ ] Candidate targets are separated from NK readouts.
- [ ] A second mechanism card runs through the synthetic pipeline.
- [ ] All manuscript numbers match final result tables.
- [ ] Encoding, ORCID, funding, section numbering and PDF are clean.
- [ ] Full test suite and synthetic pipeline pass.
