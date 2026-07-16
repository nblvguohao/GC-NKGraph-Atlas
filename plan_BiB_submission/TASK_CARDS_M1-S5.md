# 修回任务卡 · M1–S5（Nature-reviewer 视角 → 可执行）

> 来源：`/nature-reviewer` 三审 + 交叉综合（2026-07-12）后整理的自查清单。
> 目标稿件：`submission_bundle_BiB/01_manuscript/main_manuscript.md`（v0.5）。
> 与既有 `BiB_ACTION_PLAN.md`（P0–P3）互补：本文件按 Nature 视角的 M（必须）/S（建议）
> 编号，并在每卡标注对应的 P 代码。
> **边界声明**：卡内所列脚本/函数/路径均来自本次对仓库源码的实际抽查（`src/…`），
> 输出文件名部分为脚本内已写死路径、部分为建议新增（标 *新增*）。运行环境以
> `environment.yml` / `.venv` 为准；DepMap、真实 h5ad 等外部数据的可得性需本地确认。

---

## 图例

- **对应顾虑**：三审中的哪条主张 / 综合共识风险。
- **涉及脚本**：需运行或修改的 `.py`。
- **输入**：脚本依赖的数据/表。
- **产出**：预计生成或更新的文件名（`*新增` = 当前不存在，需创建）。
- **DoD**：完成判据（可勾选）。
- **退路**：确实做不到时的合规降级措辞。

---

# 🔴 必须（Must — 阻断项）

## M1 — 用真实数据替换合成 / 文献推导的关键数字
> 对应 R1 顾虑 1 + 综合共识风险 1；稿件局限 3、12。拆为 M1a/M1b/M1c 三张子卡。

### M1a · 真实 DepMap CERES 替换文献推导值 ⭐最高优先
- **对应顾虑**：表 4 Tier 1（SGMS2/SMPD3/SMPD1/NT5E）"体外非必需"判定当前来自
  `_literature_essentiality_fallback()`（DepMap 下载失败）。
- **涉及脚本**：`src/interpretation/validation_v2_nk_depmap.py`
  （函数 `query_depmap()` → 失败时落到 `_literature_essentiality_fallback()`）。
- **输入**：
  - DepMap Public 25Q2 `CRISPRGeneEffect.csv`（脚本内 `depmap_url`，~80 MB；
    自动下载失败则手动下载放本地，改 `query_depmap()` 读本地路径）。
  - `data/processed/bulk/tcga_stad_expression.tsv`、`results/tables/nk_state_labels.tsv`（DE 维度）。
  - `GENES_37`（脚本内已定义，37 基因）。
- **动作**：
  1. 手动下载 `CRISPRGeneEffect.csv` 到 `data/external/depmap/`，在 `query_depmap()`
     增加"本地文件优先"分支，跳过失效的 URL 下载。
  2. 重跑 `python src/interpretation/validation_v2_nk_depmap.py`。
  3. 核对 `data_source` 列不再是 `literature fallback`。
- **产出**：
  - `results/tables/target_validation_depmap.tsv`（更新，`data_source`=DepMap 25Q2）
  - `results/tables/target_validation_nk_state_de.tsv`（更新）
  - `results/tables/target_validation_v2_merged.tsv`（更新）
  - `submission_bundle_BiB/03_supplementary/tables/target_validation_depmap.tsv`（同步）
- **DoD**：
  - [ ] 表 4 CERES 列全部来自真实 DepMap；脚注删除"literature-derived"。
  - [ ] 稿件局限 12 删除或改为"已用真实 DepMap 25Q2 确认"。
- **退路**：跑不通则把表 4"DepMap CERES"整列降级为"literature-annotated, pending
  DepMap confirmation"，并从 Tier 定级依据中移除体外必需性（仅保留 NK-state DE 作定级），
  正文明确此弱化。

### M1b · 真实 h5ad 重跑 count-depth 控制
- **对应顾虑**：局限 3——P0-2 仅合成数据验证。
- **涉及脚本**：`src/topology/count_depth_control.py`（已内置 `--real` 开关；
  `residualize()` 对 `total_counts`/`n_genes` 残差化）。
- **输入**：`data/processed/scrna/gc_integrated.h5ad`（真实，git-ignored，需本地就位）。
- **动作**：`python src/topology/count_depth_control.py --real`。
- **产出**：`results/tables/sst_axis_count_depth_control.tsv`（真实数据版，覆盖合成版）。
- **DoD**：
  - [ ] Methods §2.3 补一句"H2/H3 在计数深度残差化后方向/显著性不变（真实数据）"或如实报告变化。
  - [ ] 局限 3 中删去 count-depth 的"synthetic-only"表述。
- **退路**：真实 h5ad 不可得时，保留合成版但在正文明确标注"method demonstrated on
  synthetic; real-data rerun pending data access"。

### M1c · 真实 h5ad 重跑模块置换检验
- **对应顾虑**：局限 3——P0-3 置换检验仅合成数据；同时支撑 M2。
- **涉及脚本**：`src/topology/module_permutation_test.py`（N=10,000 置换，
  `partial_corr_via_residuals()`）。
- **输入**：`data/processed/scrna/gc_integrated.h5ad`（真实）。
- **动作**：确保真实 h5ad 就位后直接 `python src/topology/module_permutation_test.py`
  （脚本自动优先真实、缺失时回落合成）。
- **产出**：`results/tables/h3_module_permutation_test.tsv`（真实数据版）。
- **DoD**：
  - [ ] §3.2 报告真实数据下 protrusion~cytotoxicity 的经验 P（vs 同尺寸随机模块零分布）。
- **退路**：同 M1b。

---

## M2 — 证明效应臂耦合独立于 NK 通用激活程序
> 对应 R2 顾虑 1 + 综合共识风险 2；扩展 P0-3。现有偏相关 partial r=0.286（33% 共享方差）不足。

- **涉及脚本**：
  - `src/a100_recompute/run_h3_activation_control.py`（现有偏相关，读
    `sst_axis_scores_single_cell.tsv`，产出 `h3_activation_control.tsv`）。
  - `src/topology/module_permutation_test.py`（经验零分布，M1c）。
- **输入**：`results/tables/sst_axis_scores_single_cell.tsv`；真实 scNK h5ad（scVI 潜因子）。
- **动作**（二选一或并用）：
  1. **激活匹配子集**：按 `_activation_score` 分箱，在激活水平匹配的 NK 子集内重算
     protrusion~cytotoxicity，报告组内 r。
  2. **scVI 潜因子残差化**：对 protrusion / cytotoxicity 两模块打分先对 scVI 潜因子
     （或前 k 个 PC）残差化，再算相关。
- **产出**：
  - `results/tables/h3_activation_control.tsv`（扩展：新增激活匹配子集/残差化列）
  - `results/tables/h3_activation_matched_subset.tsv` *新增*
  - `submission_bundle_BiB/03_supplementary/tables/h3_activation_matched_subset.tsv` *新增*
- **DoD**：
  - [ ] §3.2 增补"激活匹配子集 / 潜因子残差化后效应臂仍显著（或减弱到 X）"。
  - [ ] 正文把"independent replication of the functional endpoint"改为与实测残差
        比例一致的克制措辞。
- **退路**：若残差化后显著性明显下降，如实报告"约 X% 可归因共激活"，并将效应臂结论
  从"稳健独立复现"降为"与 NK 激活程序部分耦合的相关"。

---

## M3 — 收紧"可复用引擎/框架"措辞，使之与"仅一卡真实验证"对齐
> 对应 R3 顾虑 1 + 综合共识风险 3；稿件局限 11。等价于 P2-1。

- **涉及脚本**：`src/interpretation/run_multicard_analysis.py`
  （现产出卡间比较/重叠，非端到端真实数据）。
- **输入**：`configs/mechanism_cards/registry.yaml`（adenosine / TGFβ / MICA-B 卡）；
  对应真实表达数据。
- **动作**（(a) 与 (b) 二选一）：
  - **(a) 补真实端到端**：选 adenosine 或 TGFβ 卡，跑通最小端到端（哪怕仅恢复图），
    使"demonstrated on ≥2 cards"成立。
  - **(b) 下调措辞**：摘要、Key Points 第 1 条、§4.2 标题、结论段统一把
    "framework/engine/reusable" 改为 "a mechanism-card formalism demonstrated on
    one mechanism (Zheng 2023), with a schema designed for reuse"。
- **产出**：
  - (a) `results/tables/mechanism_card_<name>_recovery.tsv` *新增* + 补充图
  - (a) 更新 `submission_bundle_BiB/03_supplementary/tables/mechanism_card_comparison.tsv`
  - (b) `submission_bundle_BiB/01_manuscript/main_manuscript.md`（措辞修订）
- **DoD**：
  - [ ] 摘要 + Key Points 不再把"可复用引擎"作为核心卖点，除非有 ≥2 卡真实验证。
- **退路**：(a) 跑不通即执行 (b)——纯改写，零计算成本。

---

## M4 — 澄清 NK 状态标注与下游分析的循环风险
> 对应 R3 顾虑 3；稿件局限 6。P0–P3 未覆盖，属新增。

- **涉及脚本**：
  - `src/immune_scoring/nk_scores.py`（NK-state 标签定义基因集）。
  - `src/topology/sst_axis_scrna.py`、`src/topology/sst_axis.py`（SST-axis 打分基因集）。
  - `src/models/gc_nkgraph_atlas.py`（分类器输入特征）。
- **输入**：三处基因集清单（state-defining vs SST-axis vs 分类输入）。
- **动作**：
  1. 导出三处基因集，计算两两重叠（Jaccard / 交集清单）。
  2. 若 H3 / 分类所用基因与 state-defining 基因重叠，去重后重算并对比。
- **产出**：
  - `results/tables/geneset_separation_audit.tsv` *新增*（三集合成员 + 重叠矩阵）
  - `submission_bundle_BiB/03_supplementary/tables/geneset_separation_audit.tsv` *新增*
- **DoD**：
  - [ ] Methods 或补充新增一段"基因集分离/去重"说明，读者可确认非自证。
  - [ ] 若存在重叠，报告去重后 H3/分类结果与原结果一致（或如实差异）。
- **退路**：无重叠即一句话说明"state-defining 与 axis / 分类基因集互斥"并附审计表即可结题。

---

# 🟠 建议（Should — 显著加分）

## S1 — H3 留一样本敏感性分析（回应 I²=96%）
> 对应 R1 顾虑 2 + 共识风险 4；扩展 P0-1。

- **涉及脚本**：`src/topology/pseudoreplication_correction.py`
  （现产出 `sst_axis_pseudoreplication_corrected.tsv`；per-sample r → Fisher z → DerSimonian–Laird）。
- **输入**：9 样本的 per-sample r（脚本已算）。
- **动作**：增加 leave-one-sample-out 循环，逐次剔除 1 样本重算合并 r 与 95% CI；
  标出驱动效应的样本（per-sample r 跨度 [0.009, 0.560]）。
- **产出**：
  - `results/tables/h3_leave_one_sample_out.tsv` *新增*
  - `submission_bundle_BiB/03_supplementary/tables/h3_leave_one_sample_out.tsv` *新增*
- **DoD**：
  - [ ] §3.2 / §4.1 给出去极端样本后 CI 仍不含 0（或如实说明哪些样本驱动效应）。
  - [ ] "robust"措辞与敏感性结果一致。
- **退路**：若某单一样本主导，改称"effect present but sample-driven; robust to
  leave-one-out except sample X"。

## S2 — 证明 `metabolic_crosstalk` 符号校准无信息泄漏
> 对应 R1 顾虑 3。关联既有 T17 边价值分析。

- **涉及脚本**：
  - `src/graph_construction/build_heterograph.py`（`metabolic_crosstalk` 边、符号校准）。
  - `src/a100_recompute/run_t17_edge_external_value.py`（边的外部价值/迁移测试）。
- **输入**：LIHC（定标）与 STAD（测试）表达数据。
- **动作**：在留出 / 独立队列上定标符号，再在测试队列检验；或证明定标只用方向不接触结局标签。
  一并如实呈现既有 STAD→LIHC 迁移无增益结果（ΔMCC=+0.002, p=0.44）。
- **产出**：
  - `results/tables/mc_edge_sign_calibration_audit.tsv` *新增*
  - 更新引用现有 `03_supplementary/tables/`（边消融相关表）
- **DoD**：
  - [ ] 读者可确认"同队列校准符号"不构成循环（附留出定标结果）。
- **退路**：若无第二队列可定标，明确声明"sign fixed by mechanism direction, not
  fitted to outcome"，并把该边定位为"结构性/可解释性贡献，非下游必需"。

## S3 — NK-state DE 的功效说明与（若可能）外部复现
> 对应 R2 顾虑 + R1 顾虑 4；稿件局限 13（dysfunctional n=20）。

- **涉及脚本**：`src/interpretation/validation_v2_nk_depmap.py`（DE 维度，NK-state 对比）。
- **输入**：TCGA-STAD DE（现有）；若有第二个带 scRNA 的胃癌队列则复现 SGMS2/NT5E 方向。
- **动作**：
  1. Tier 表脚注明确标注"dysfunctional n=20，P 值谨慎解读，log2FC 可靠"。
  2. （可选）在外部胃癌队列复现 SGMS2/NT5E 的 DE 方向。
- **产出**：
  - `results/tables/target_validation_nk_state_de.tsv`（补功效标注列）
  - `results/tables/nk_state_de_external_replication.tsv` *新增，若做外部复现*
- **DoD**：
  - [ ] 表 4 / §3.5 脚注写明功效限制；或外部复现方向一致。
- **退路**：无外部队列时仅做脚注标注即达标。

## S4 — 明确 GNN / 图组件的定位（标题与摘要层面）
> 对应 R3 顾虑 2 + 共识风险 3；关联 P1-1 / P2-3。

- **涉及脚本**：无需运行（已有证据表）；引用现有
  `03_supplementary/tables/model_comparison_stats.tsv`、`ablation_edge_types.tsv`、
  `domain_baselines_summary.tsv`。
- **输入**：已有基线/消融结果（GNN 与 LightGBM/XGBoost、无图 SST-module 逻辑回归无显著差异；
  metabolic_crosstalk 边跨队列无增益）。
- **动作**：标题与摘要不暗示图模型带来分类精度优势；把价值锚定在"机制结构化 embedding +
  主张边界纪律"。正文 §3.4 / §3.7 已诚实，主要改标题/摘要口径。
- **产出**：`submission_bundle_BiB/01_manuscript/main_manuscript.md`（标题 + 摘要修订）。
- **DoD**：
  - [ ] 标题 / 摘要不宣称精度优势；与 §3.4 结论一致。
- **退路**：纯改写，无计算成本。

## S5 — 提供图像本体供审稿核对
> 对应评审边界（本次仅见图注，未见图像）。等价于 P3-3。

- **涉及脚本**：`src/figures/make_figures.py`、`src/figures/make_workflow_figure.py`。
- **输入**：各分析产出表（Fig1↔表 2、Fig2↔外部验证表、Fig3↔靶点表、Fig4↔基线表）。
- **动作**：重生成 fig0–fig4，逐子图核对与正文数值一致（尤其 Fig1D = 表 2 森林图）。
- **产出**：`submission_bundle_BiB/02_figures/fig{0..4}_*.{pdf,png}`（更新）。
- **DoD**：
  - [ ] 每张图每个子图能在正文找到对应数值；面板编号与图注一致。
- **退路**：无（投稿硬性要求）。

---

## 依赖与建议顺序

```
M1a(DepMap) ─┐
M1b/M1c(真实h5ad) ─┼─→ M2(激活独立性,依赖真实h5ad) ─→ S1(留一,依赖伪重复表)
             │
M4(去循环,独立) ── S2(符号校准,独立) ── S3(DE功效,依赖M1a的DE) 
M3(措辞或补卡,独立) ── S4(措辞,依赖已有基线表) ── S5(图,依赖所有表定稿)
```

- **关键路径**：真实数据就位（DepMap + h5ad）→ M1 → M2 → S1。这条链一通，靶点表与
  效应臂两大结论同时升级。
- **零计算成本可先做**：M3(b) 措辞、S4 措辞、M4 若无重叠——可与数据重跑并行推进。
- **最小可投稿集（资源受限）**：M1a + M1b + M2 + M3(b) + S4 + S5。

## 与 BiB_ACTION_PLAN.md 的映射

| 本卡 | 对应 P 代码 | 关系 |
|------|------------|------|
| M1a | （新）§3.5 靶点验证 | 新增，真实 DepMap |
| M1b | P0-2 | 真实数据重跑 |
| M1c | P0-3 | 真实数据重跑 |
| M2  | P0-3 | 更严独立性检验 |
| M3  | P2-1 | 补卡或下调措辞 |
| M4  | （新）循环风险 | 新增 |
| S1  | P0-1（扩展） | 留一敏感性 |
| S2  | T17 边价值 | 符号校准审计 |
| S3  | （新）DE 功效 | 新增 |
| S4  | P2-3 / P1-1 | 定位措辞 |
| S5  | P3-3 | 图像核对 |
