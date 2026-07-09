# TDD：主编终审 + 同行评审修改清单

> **角色双重性**：本文档上半部分以 *Briefings in Bioinformatics* 主编视角逐项列出格式与完成度门槛；下半部分以同行评审专家视角给出关键学术修改意见。最后合并为可执行的 TDD 任务列表。
>
> **阅读顺序**：主编终审（§A）→ 同行评审（§B）→ TDD 任务（§C）
>
> 日期：2026-07-09

---

# A. 主编终审 — BiB 格式与完成度门槛

## A.1 总体评价

这篇稿件在科学诚实性上优于该期刊收到的大多数投稿——这在计算方法学论文中罕见且值得称赞。两臂设计、预先注册的假设、以及"我们不声称完全恢复；范围图本身即是贡献"的明确声明，构成了一个在智识上诚实的论文骨架。然而，在送审之前，以下项目必须解决。

## A.2 阻断性事项（送审前必须完成）

### A.2.1 §3.7 消融实验占位符

BiB 投稿清单标记为 "⚠️ Ablation (§3.7) — either run a minimal ablation (graph edges on/off) or remove the placeholder heading."

**主编意见**：§3.7 在 BiB_submission_checklist 中被列为待完成项，但在 main_manuscript.md 中我未找到 §3.7 的实际内容。要么执行最小消融（关闭 metabolic_crosstalk 边 vs 全图，报告 ΔMCC），要么在正文中诚实声明"消融实验超出了当前研究范围，原因是……"并完全删除该小节标题。

### A.2.2 手稿中的编辑注释残留

以下注释在送审前必须删除（它们属于内部工作文档，不是学术手稿）：

| 行位置 | 内容 | 操作 |
|--------|------|------|
| 作者列表下方 | `*(ORCID iDs: fill in each author's ORCID...)*` | 删除 |
| 作者传记下方 | `*(Biographies ~30 words each...)*` | 删除 |
| CRediT 声明下方 | `*(CRediT taxonomy; adjust...)*` | 删除 |
| 文档底部 | `> **Document status:** Draft v0.3...` | 删除 |
| 作者列表 | `[0000-0000-0000-0000]`（7 处） | 替换为真实 ORCID |

### A.2.3 标题长度

当前标题（含副标题）约 190 字符。BiB 推荐 ≤ 150 字符。

**建议缩短方案**：
- **选项 A（推荐）**：「Reconstructing the Serine–Sphingomyelin–Membrane-Topology Axis of NK-Cell Immune Evasion from Tumor Transcriptomes: A Single-Cell-Informed Heterogeneous Graph Framework」
- **选项 B**：「A Single-Cell-Informed Heterogeneous Graph Framework Reconstructs NK-Cell Immune-Evasion Mechanisms from Tumor Transcriptomes」
- **选项 C**：「GC-NKGraph-Atlas: Reconstructing NK-Cell Immune Evasion Mechanisms from Tumor Transcriptomes via Heterogeneous Graph Learning」

### A.2.4 补充材料索引

BiB 投稿清单提到 "⚠️ Add a Supplementary Methods note + index before submission." 当前 `submission_package/02_supplementary_tables/` 下有文件，但缺少一个索引文件说明每个补充表的内容及与正文的对应关系。

### A.2.5 图表自足性检查

主编要求所有图注必须自足（不依赖正文即可理解图中的每个缩写、颜色编码和统计量）。当前 Fig 1-5 的图注需要逐一核验。

## A.3 强烈建议（决定是否送审）

### A.3.1 T4 scRNA QC 未在服务器执行

`qc_filter.py` 已编写但未在真实数据上运行。回归测试 (`test_qc_regression.py`) 依赖于 QC 输出文件。这是方法学基本要求——审稿人和编辑都会期望看到标准的 scRNA QC 流程（min_genes, max_genes, pct_mito, doublet removal）。

### A.3.2 预注册校准变更日志未更新

根据 BiB_submission_checklist.md 的 G 项："Pre-registration note: log the calibration change (metabolic edge now anchored on the effector arm, since the topology arm did not recover)." 这项变更必须记录在 `configs/sst_axis_config.yaml` 或单独的预注册日志中——这直接关系到论文的科学诚实性声明。

### A.3.3 干净克隆测试未完成

F1 (`pytest -q`) 和 F2 (`--synthetic` end-to-end) 需要在干净机器上验证。如果代码在编辑的机器上能跑但在干净克隆上崩溃，审稿人会立即报告"无法复现"。

### A.3.4 作者传记内容脱节

周爱莲教授的传记描述为"agricultural informatization and blockchain applications"——这与癌症免疫信息学的论文主题存在明显脱节。审稿人会注意到这一点。建议补充她与论文主题相关的专业背景，或解释跨领域合作的关系。

### A.3.5 通讯作者单位问题

周爱莲的归属单位（中国农业科学院农业信息研究所 / 农业部农业区块链应用重点实验室）同样与论文主题存在领域距离。需确保 ScholarOne 注册时通讯作者邮箱域（@caas.cn）能够通过期刊验证。

---

# B. 同行评审报告

# **Reviewer Report**

**Manuscript:** GC-NKGraph-Atlas: Reconstructing the Serine–Sphingomyelin–Membrane-Topology Axis of NK-Cell Immune Evasion from Tumor Transcriptomes
**Journal:** Briefings in Bioinformatics
**Reviewer expertise:** Computational biology, single-cell transcriptomics, graph neural networks in biomedicine

---

## Overall Assessment

This is an intellectually honest and methodologically rigorous paper. The two-arm design with pre-registered hypotheses, the mechanism-card abstraction, and — most importantly — the willingness to report *partial* recovery and measure the limits of transcriptomic inference are all commendable and, frankly, rare in the computational biology literature. The paper has the bones of an excellent Problem Solving Protocol.

However, the manuscript in its current form has several vulnerabilities that a competent reviewer will identify. My review below is structured to help the authors strengthen these points *before* submission, because I believe the core contribution deserves to be published.

**Recommendation:** Major revision (but structurally sound — all issues are addressable within 1-2 weeks of focused work).

---

## Major Comments

### B.1 The H2 effect size is statistically significant but biologically negligible (CRITICAL)

The paper frames the single-cell H2 result (SM-balance→protrusion, r=+0.030, p=6×10⁻³) as "recovered" — and uses this as evidence that "cell-type resolution is required for the metabolic layer." The p-value is indeed significant, but with n=8,310 cells, a correlation of r=0.030 explains **0.09% of the variance**. Virtually any non-zero correlation would be "significant" at this sample size.

This is the paper's most vulnerable claim. A reviewer will ask:

1. What is the *biologically meaningful* effect size for SM-balance→protrusion coupling? Is r=0.030 within or below that range?
2. The claim "invisible in bulk, significant in single-cell" conflates two different things: (a) the correlation changed (from r=-0.017 to r=+0.030 — a tiny shift) and (b) the p-value crossed the significance threshold (from p=0.72 to p=6×10⁻³ — entirely driven by n).
3. What happens at more realistic NK sample sizes? If you down-sample to n=500 NK cells (a typical scRNA NK yield), is the correlation still detectable?

**Required revision:**
- Qualify the H2 language: "weak but correctly signed" rather than "recovered."
- Report confidence intervals on r (e.g., bootstrap CI).
- Add a down-sampling analysis showing at what NK cell count the signal becomes undetectable.
- Add a brief discussion of biologically meaningful effect sizes for metabolic-enzyme transcript correlations.

### B.2 Alternative explanation for H4/H5 failure is not discussed (MAJOR)

The paper interprets the H4/H5 failures (machinery transcription runs *opposite* to the physical phenotype) as "transcription does not proxy the physical topology phenotype." This is the paper's central scoping claim.

However, an alternative interpretation exists and should be discussed: the NK protrusion machinery module contains 25 genes, many of which are **general actin cytoskeleton regulators** (Arp2/3 complex, Rho GTPases, WASP/WAVE, formins) that are part of a broad NK activation program. Intratumoral NK cells, despite being functionally suppressed, may still exhibit transcriptional activation signatures — a phenomenon well-documented in T-cell exhaustion (where activated/exhausted T cells upregulate effector genes but fail to produce functional protein). The "higher transcription, lower function" pattern for protrusion machinery could therefore reflect **post-transcriptional or protein-level regulation** rather than transcriptional inadequacy per se.

**Required revision:**
- Add a paragraph in §4.1 or §4.3 discussing this alternative interpretation.
- Cite relevant T-cell exhaustion / NK dysfunction literature on transcription-protein discordance.
- Clarify that "transcription does not proxy" could mean either (a) the wrong genes were chosen for the module, (b) post-transcriptional regulation dominates, or (c) the physical phenotype truly has no transcriptional correlate.

### B.3 The metabolic_crosstalk edge calibration is opaque (MAJOR)

The Methods (§2.5) state that the metabolic_crosstalk edge sign is "calibrated on the liver positive-control cohort, not hard-coded." But the Results (§3.2) show that H1 (tumor_serine_capacity ⟂ nk_sm_balance) was null (r=-0.016, p=0.74). If the calibration source failed, how was the edge actually calibrated?

A reviewer will press on this. The pre-registration promise was calibration on liver; liver didn't provide a signal; the downstream analysis nevertheless uses the edge. This needs transparent documentation.

**Required revision:**
- Add a paragraph in Methods §2.5 explicitly stating: "Because H1 was null in the liver control, the metabolic_crosstalk edge weight was set to a hypothesis-driven default derived from the anchor paper's mechanistic direction (tumor serine ↑ → NK SM ↓), flagged as NEEDS_REVIEW in the pre-registration log, and treated as an uncalibrated prior in all downstream analyses."
- Update `configs/sst_axis_config.yaml` to reflect this.
- Consider whether the edge should be downgraded (weight 0.5 → 0.2) to reflect the calibration failure, or whether results with and without the edge should be compared.

### B.4 GNN "interpretability" claim is stated but not demonstrated (MAJOR)

The paper's core defense of the GNN — since it doesn't outperform LightGBM/XGBoost — is that it provides "mechanism-structured gene embeddings" that tabular baselines don't. This is a legitimate argument, but it needs evidence.

Specifically:
1. What downstream analysis *uses* these embeddings that could not be done with XGBoost feature importance or SHAP values?
2. Can you show a concrete example where the gene embedding reveals a relationship (e.g., two genes that are mechanistically related cluster in embedding space) that a correlation matrix or SHAP dependence plot would miss?
3. The axis analyses (§3.2, §3.5) use per-cell module scores (mean z-scores of gene sets), not the learned embeddings. Where exactly do the embeddings enter the biological analysis?

**Required revision:**
- Add one concrete demonstration of the embeddings' utility: e.g., a UMAP of the learned gene embeddings colored by axis membership, showing that mechanistically related genes cluster.
- Or: show that the embedding-space distance between tumor-serine genes and NK-topology genes correlates with the axis coupling strength across tissues.
- Or: if the embeddings are only used for NK-state classification and the axis analyses are independent, say so explicitly and rephrase the "added interpretability" claim to reflect the actual scope.

### B.5 External validation scope is narrow (MODERATE)

The external validation in two independent gastric cohorts (GSE62254, GSE84437) is a genuine strength — but it validates only the effector arm (protrusion~cytotoxicity), which was already the most robust finding. The more interesting claims — cell-type-resolved SM-balance coupling, topology non-recovery — are not externally validated because no external scRNA cohort with NK cells is available.

This is understandable and acceptable, but the limitation should be stated explicitly: "We externally validated the axis layer for which validation data exist (the effector arm in bulk); the cell-type-resolved metabolic coupling and topology non-recovery remain to be validated when independent NK scRNA cohorts become available."

### B.6 Mechanism-card reusability: claimed but not demonstrated (MODERATE)

The mechanism-card abstraction is the paper's most novel conceptual contribution. However, the paper demonstrates only one card (`zheng_nk_sm_topology`). The claim "applying the framework to a new mechanism requires only a new card" is a promise, not an empirical demonstration.

**Suggested revision:**
- Option A (stronger): Create a second minimal mechanism card for a different mechanism (e.g., adenosine-A2AR mediated NK suppression, or TGFβ-driven NK exclusion) — even at the specification level only, with the card YAML included as a supplementary file. This shows the abstraction works for at least one additional mechanism.
- Option B (acceptable): Add an explicit scope limitation: "Demonstrating the mechanism-card abstraction with a second published mechanism is important future work; the current study establishes the formalism and validates it on one well-characterized mechanism."
- Option C (current, risky): Keep as-is but be aware reviewers may ask for Option A.

### B.7 Candidate target rank stability not assessed (MODERATE)

The composite target score uses fixed weights (tumor specificity 0.30, NK dysfunction 0.20, SST-axis membership 0.30, axis core 0.10, literature 0.10). A reviewer will ask: how sensitive are the rankings to these weights? Would PHGDH and SGMS2 still top the list if the weights were ±0.10?

**Required revision:**
- Add a rank stability analysis: perturb each weight by ±0.10 and report whether the top-10 list is stable (e.g., Jaccard index of overlap).
- Or state the weights as "chosen to reflect X, validated by Y" with a sensitivity note in Limitations.

### B.8 Four-state characterization → binary classification gap (MINOR)

Methods §2.4 defines four NK states in detail, then collapses them to binary ("NK-hot-cytotoxic vs rest") for the classification task. The justification is reasonable (biological relevance, underpowered multi-class, embedding evaluation focus). However, a reviewer will ask: at what sample size *would* multi-class classification be adequately powered? A brief power analysis or rule-of-thumb note would help readers assess the framework's scalability.

---

## Minor Comments

### B.9 Methods §2.3 gene count verification

The seven SST modules are listed as 10 + 2 + 4 + 9 + 25 + 11 + 1 = 62 genes. Verified. Consider adding the total gene count explicitly in §2.3 for reader convenience.

### B.10 Table 3 row ordering

In Table 3, methods are ordered by MCC descending. However, the GNN row is bolded — this is reasonable but should be explicitly noted in the table caption as "the proposed method."

### B.11 The term "de-circularized" in §3.5

The term "de-circularized" is used extensively and is a good coinage for this problem class. However, consider defining it briefly on first use (§3.5): "i.e., NK effector markers that are both the axis readout and enriched in tumor samples are separated from genuine tumor-intrinsic targets."

### B.12 GSE62254/GSE84437 probe mapping

The external validation required probe→gene remapping because the processed matrices retained probe IDs. This is a well-executed fix, but the Methods should briefly describe the remapping approach (platform annotation file used, handling of multi-probe mapping).

### B.13 In-silico SM-restoration stratification (§3.6)

The section is present but the `main_claims.md` marks R6 (SM-restoration) as "⬜ Pending." If computed results exist, they should be described; if not, either remove the section or add a note that the readout is defined but stratification awaits a larger cohort.

### B.14 Missing Supplementary Methods

The supplementary tables (`submission_package/02_supplementary_tables/`) need a brief "Supplementary Methods" document explaining the content and structure of each supplementary file.

---

## Summary of Required vs Recommended Changes

| Priority | Item | Effort |
|----------|------|--------|
| **Required** | B.1 — H2 effect size qualification + down-sampling | 2-4 hr |
| **Required** | B.2 — Alternative interpretation of H4/H5 | 1-2 hr writing |
| **Required** | B.3 — metabolic_crosstalk calibration transparency | 1-2 hr + config update |
| **Required** | B.4 — Concrete GNN interpretability demonstration | 4-8 hr |
| **Required** | A.2.2 — Remove editorial notes from manuscript | 30 min |
| **Required** | A.2.1 — Ablation or remove §3.7 | 2-4 hr or 15 min |
| **Required** | A.2.3 — Title length | 15 min |
| **Recommended** | B.5 — External validation scope limitation | 30 min writing |
| **Recommended** | B.6 — Second mechanism card or scope limitation | 1-3 hr or 15 min |
| **Recommended** | B.7 — Rank stability analysis | 2-3 hr |
| **Recommended** | A.3.1 — T4 server execution | 1-2 hr server |
| **Recommended** | A.3.2 — Pre-registration log update | 30 min |
| **Recommended** | A.3.3 — Clean clone test | 1-2 hr |
| **Optional** | B.8 — Multi-class power analysis note | 1 hr |
| **Optional** | A.2.5 — Figure legend audit | 1 hr |
| **Optional** | A.3.4 — Author bio alignment | 30 min |

---

# C. TDD 任务列表（续接 TDD_remaining_work.md）

> 图例：🔴 未开始 · 🟡 进行中 · 🟢 通过
> 任务编号接续现有 T1–T7，新增 T8–T13。

---

## T8 — 手稿格式清洁（主编阻断项）🔴

**RED（验收判据）**
- A1：`grep -nE '\[PLACEHOLDER|\[TBD|To be (expanded|written|designed)' main_manuscript.md` 返回空（已有）✅
- A2：`grep -nE 'ORCID iDs: fill in|Biographies ~30|CRediT taxonomy; adjust|Document status: Draft' main_manuscript.md` 返回空（新增）
- A3：`grep -c '0000-0000-0000-0000' main_manuscript.md` 返回 0
- A4：标题 ≤ 150 字符
- A5：`cover_letter.md` 中 `Internal note (delete before sending)` 已删除

**完成定义**：手稿中不含任何编辑注释、占位符或内部工作标记。

```bash
# 一键验证
grep -nE '\[PLACEHOLDER|\[TBD|To be (expanded|written|designed)|ORCID iDs: fill|Biographies ~30|CRediT taxonomy; adjust|Document status: Draft|Internal note' manuscript/main_manuscript.md manuscript/cover_letter.md
# 预期输出：空
```

---

## T9 — H2 效应量限定 + 下采样分析（审稿人 B.1, CRITICAL）🔴

**RED**
- B1：H2 单细胞相关性汇报 **bootstrap 95% CI**（r=+0.030 的 CI 下界不得跨越零）
- B2：**下采样曲线**：从 n=8310 逐步降至 n=100（步长 ~200），绘制 r ± CI vs n 曲线。标注 r 失去显著性（p≥0.05）时的样本量
- B3：§3.2 和 §4.1 中 H2 措辞从"recovered"改为"weak but correctly signed (r=+0.030 [95% CI: X–Y], p=6×10⁻³, variance explained = 0.09%)"
- B4：§2.9 H2 预期方向处添加注释"效应量可能极小（代谢酶转录是酶活性的噪声代理）"

**GREEN**
```bash
python src/analysis/h2_effect_size.py \
    --input results/tables/sst_axis_scores_single_cell.tsv \
    --output results/tables/h2_bootstrap_ci.tsv \
    --output-fig results/figures/figS_h2_downsample.pdf
python - <<'PY'
import pandas as pd
ci = pd.read_csv("results/tables/h2_bootstrap_ci.tsv", sep="\t")
assert ci["r_lower"].iloc[0] > 0, "H2 CI crosses zero"
assert ci["r"].iloc[0] < 0.05, f"H2 r={ci['r'].iloc[0]} — verify this is the genuinely tiny effect"
print("T9 B1-B2 PASS")
PY
```

**REFACTOR**：将下采样图加入补充材料；在 §4.3 Limitations 添加一条关于"大 n 下 p 值膨胀"的说明。

---

## T10 — H4/H5 替代解释 + metabolic_crosstalk 校准透明化（审稿人 B.2 + B.3）🔴

**RED**
- C1：§4.1 或 §4.3 新增一段（~200 词），讨论 H4/H5 的替代解释：
  - 内肿瘤 NK 的转录激活-功能解耦（类比 T 细胞耗竭中的转录-蛋白不一致）
  - 突出肌动蛋白骨架基因是广谱激活标志而非拓扑状态特异的可能性
  - 引用相关文献（如 T 细胞耗竭中 effector gene 转录上调但蛋白不表达的现象）
- C2：§2.5 新增一段说明 metabolic_crosstalk 边的校准状态：
  - "Because H1 was null, the edge weight defaults to the anchor paper's mechanistic direction, flagged NEEDS_REVIEW."
  - 考虑将权重从 0.5 降为 0.2（或以开关形式做敏感性分析）
- C3：更新 `configs/sst_axis_config.yaml` 中 `tumor_serine_capacity.expected_direction` 为 `CALIBRATION_FAILED_DEFAULTING_TO_ANCHOR_PAPER_DIRECTION`
- C4：在预注册日志 (`configs/sst_axis_prereg_log.md`) 新增条目记录校准变更

**GREEN**
```bash
# 验证 manuscript 中存在替代解释段落
grep -c "alternative" manuscript/main_manuscript.md
# 预期 ≥ 1（§4.1 或 §4.3）
grep -c "NEEDS_REVIEW\|CALIBRATION_FAILED" configs/sst_axis_config.yaml
# 预期 ≥ 1
```

---

## T11 — GNN 可解释性具体演示（审稿人 B.4, MAJOR）🔴

**RED**
- D1：生成**基因嵌入空间 UMAP 图**：用学到的基因嵌入矩阵 E，按 SST 轴模块着色（tumor_serine / nk_sm_synthesis / nk_sm_catabolism / nk_protrusion_machinery / nk_cytotoxicity），展示机械性相关的基因在嵌入空间中聚类
- D2：计算**嵌入空间距离 vs 轴耦合强度**：对每个 SST 基因对，计算嵌入余弦距离；与它们在 scRNA 中的共表达相关性做 Spearman 相关。预期正相关（嵌入距离近 → 共表达高）
- D3：§4.2 中将"added interpretability"替换为具体可引用的发现（例如"在嵌入空间中，PHGDH 与 SGMS2 的余弦距离为 X，而它们之间的代谢耦合是该轴的限速步骤"）

**GREEN**
```bash
python src/interpretation/gene_embedding_analysis.py \
    --embeddings results/model/gene_embeddings.npy \
    --gene-list data/processed/graph/gene_nodes.tsv \
    --sst-config configs/sst_axis_config.yaml \
    --output-dir results/figures/
# 预期输出: results/figures/fig6_gene_embedding_umap.pdf
#           results/tables/embedding_distance_vs_coexpression.tsv
python - <<'PY'
import pandas as pd
df = pd.read_csv("results/tables/embedding_distance_vs_coexpression.tsv", sep="\t")
assert df["spearman_r"].iloc[0] > 0, "embedding distance should positively correlate with coexpression"
print("T11 D2 PASS")
PY
```

**REFACTOR**：更新 §4.2 以引用新的嵌入分析结果；在 Table 3 脚注中添加一条指向嵌入分析的说明。

---

## T12 — 消融实验（主编 A.2.1 + 审稿人 B.3 附带）🔴

**RED**
- E1：运行最小消融实验：比较全图模型 vs 去除 `metabolic_crosstalk` 边的模型（5 折 CV，相同种子）
- E2：报告 ΔMCC 和 ΔAUROC，配对方差检验
- E3：若 ΔMCC 在 ±0.02 范围内 → 报告"去除代谢边不影响分类性能（ΔMCC=X, p=Y），表明该边的价值在嵌入空间而非分类精度"→ 消融结果直接支持 §4.2 的论点

**GREEN**
```bash
python src/baselines/run_ablation.py \
    --full-graph data/processed/graph/edges_full.tsv \
    --ablated-graph data/processed/graph/edges_no_metabolic.tsv \
    --output results/tables/ablation_results.tsv
python - <<'PY'
import pandas as pd
df = pd.read_csv("results/tables/ablation_results.tsv", sep="\t")
full_mcc = df[df.mode=="full"].MCC.mean()
ablated_mcc = df[df.mode=="no_metabolic"].MCC.mean()
delta = full_mcc - ablated_mcc
print(f"ΔMCC = {delta:.4f}")
# 不要求 GNN 在消融中显著更好；诚实地报告结果
assert abs(delta) < 0.1, f"ΔMCC={delta:.4f} unexpectedly large"
print("T12 E1-E3 PASS")
PY
```

---

## T13 — 机制卡片复用性证明（审稿人 B.6）🔴

**RED（选择策略 A：最弱可行方案）**
- F1：新增第二个机制卡片 YAML（如 `adenosine_a2ar_nk_suppression.yaml`），至少包含：
  - 完整的 `origin`、`biology`、`transcriptional_proxy`（不少于 3 个基因模块）
  - `physical_ground_truth`（GATED targets）
  - `validation`（预注册假设）
- F2：将该卡片注册到 `configs/mechanism_cards/registry.yaml`
- F3：在 §4.2 中添加一段（~100 词）："As a preliminary test of the card format's reusability, we authored a second card for [mechanism X]. The card required only specification of [N] gene modules and [M] hypotheses; no pipeline code was modified."
- F4：或选择策略 B：在 §4.3 Limitations 中明确声明"Demonstrating the mechanism-card with a second published mechanism is important future work."

**GREEN**
```bash
python - <<'PY'
import yaml
with open("configs/mechanism_cards/registry.yaml") as f:
    reg = yaml.safe_load(f)
assert len(reg["cards"]) >= 2, f"registry has {len(reg['cards'])} card(s), need ≥ 2"
print([c["id"] for c in reg["cards"]])
print("T13 F1-F2 PASS")
PY
```

**REFACTOR**：更新 Figure 5（mechanism-card 概念图）以展示多卡片注册表。

---

## T14 — 候选靶点排名稳定性分析（审稿人 B.7）🔴

**RED**
- G1：对五个评分维度分别做 ±0.10 的权重扰动（每次仅扰动一个维度，其他维度的权重按比例重归一化）
- G2：计算每次扰动后的 top-10 列表与原始 top-10 的 Jaccard 重叠指数
- G3：报告平均 Jaccard ≥ 0.7（即 top-10 在前 10 个权重配置中基本稳定）
- G4：特别检查 PHGDH、SGMS2 是否在所有扰动下保持在 top-3

**GREEN**
```bash
python src/interpretation/target_rank_stability.py \
    --input results/tables/tumor_intrinsic_candidates.tsv \
    --output results/tables/target_rank_stability.tsv \
    --n-perturbations 10
python - <<'PY'
import pandas as pd
df = pd.read_csv("results/tables/target_rank_stability.tsv", sep="\t")
# 检查平均 Jaccard
avg_j = df[df.metric=="jaccard_top10"]["value"].mean()
assert avg_j >= 0.7, f"average Jaccard = {avg_j:.3f} < 0.7 — top targets unstable to weight perturbation"
# 检查 PHGDH 是否稳定在 top-3
phgdh_top3 = df[df.gene=="PHGDH"]["frac_in_top3"].mean()
assert phgdh_top3 >= 0.8, f"PHGDH in top-3 only {phgdh_top3:.0%} of perturbations"
print(f"T14 PASS: avg Jaccard={avg_j:.3f}, PHGDH top-3 stability={phgdh_top3:.0%}")
PY
```

---

## T15 — 补充材料索引 + 最终润色（主编 A.2.4, A.2.5）🔴

**RED**
- H1：创建 `submission_package/02_supplementary_tables/README.md`，列出每个补充文件的内容、列说明及与正文表格/图片的对应关系
- H2：逐图核验 Fig 1-5 的图注自足性（每个缩写首次出现时给出全称，每个颜色编码标注含义，统计量标注检验方法）
- H3：外部验证方法补充：在 Methods §2.9 或 Supplementary Methods 中添加 GSE62254/GSE84437 探针→基因重映射的方法细节（GPL 平台文件来源、多探针映射策略）
- H4：§3.6（in-silico SM-restoration）若结果已计算则完整描述；若未计算则删除并标注为未来工作

**完成定义**：补充材料索引存在；所有图注自足；外部验证方法透明。

---

## 里程碑与依赖（更新）

```
T4 ──► T9（H2 效应量，依赖 QC'd scRNA 数据）
T8（格式清洁，无依赖，立即执行）
T10（文本修改，无依赖，立即执行）
T11（嵌入分析，需要 GNN checkpoint）
T12（消融实验，需要服务器 GPU）
T13（第二个卡片，纯文本规格）
T14（排名稳定性，依赖 tumor_intrinsic_candidates.tsv）
T15（补充索引 + 润色，最后执行）

关键路径：T12（消融）→ T11（嵌入）→ T15（润色）→ 投稿
可并行：T8 ∥ T10 ∥ T13 ∥ T9（T9 依赖 T4）
```

### 优先级分层

| 层级 | 任务 | 如果时间不够 |
|------|------|------------|
| **P0 必须** | T8（格式清洁）、T9（H2 效应量）、T10（替代解释+校准透明）、T11（GNN 嵌入演示） | — |
| **P1 强烈建议** | T12（消融实验）、T14（排名稳定性） | T12 可用"移除 §3.7 标题"替代 |
| **P2 锦上添花** | T13（第二张卡片）、T15（补充索引） | T13 可退化为 Limitations 中的一句声明 |
| **服务器依赖** | T4（QC）、T9（下采样）、T12（消融） | 需在 A100 上执行 |

---

## 预计新增工作量

| 类别 | 时间 |
|------|------|
| 纯文本修改（T8, T10, T15-H2/H3/H4） | 4-6 hr |
| 数据分析脚本（T9, T14） | 4-6 hr |
| 嵌入分析（T11） | 4-8 hr |
| 消融实验（T12，服务器） | 2-4 hr |
| 第二个卡片规格（T13，纯文本） | 1-3 hr |
| **合计** | **15-27 hr（~2-4 个工作日）** |

---

*本文件由主编视角 + 同行评审视角双审生成，2026-07-09。应与 `TDD_remaining_work.md` 一起阅读——后者覆盖 T1–T7（数据/分析阻断项），本文档覆盖 T8–T15（格式/论证/可解释性）。*
