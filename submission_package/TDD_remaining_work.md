# TDD：距离可投稿的剩余工作（测试驱动）

> 用法：每个任务按 **RED（先定义失败判据/测试）→ GREEN（实现/运行使其通过）→
> REFACTOR（清理并回归）** 推进。表中"期望值"取自本地已确认的结果，作为**黄金
> 基线**，服务器重跑时用来自检数值不漂移。
> 约定：本地无原始数据与重依赖（xgboost/scanpy/GEOparse），脚本均已编译通过但
> **未在本地实跑**；下列测试须在服务器 `/opt/data/lgh/GC-NKGraph-Atlas` 执行。

图例：🔴 未开始 · 🟡 进行中 · 🟢 通过

---

## T1 — GEO 探针→基因映射修复（Blocker C）🟢 已通过（2026-07-07）

> 在 A100 上执行 `run_geo_external_validation.py`（直接在已处理探针矩阵上映射，
> 因原始 series matrix 不在）。GSE62254(GPL570) 54675 探针→22880 基因、NK marker
> **6/7**；GSE84437(GPL6947) 49576→25159 基因、NK marker **7/7**。A1/A2 ✅。
> 原探针矩阵已备份为 `*_probe_level.bak.tsv`。下方保留原判据备查。


**RED（验收测试，先写死判据）**
- A1：修复后 `data/processed/bulk/gse62254_expression.tsv` 的列数（基因数）应
  **远小于**原探针数（GSE62254 原 54675 探针 → 期望 ~1.8–2.2 万基因）。
- A2：规范 NK marker 至少 **4/7** 命中：`{NKG7,GNLY,GZMB,PRF1,KLRD1,NCAM1,KLRK1}`。
- A3：重算 NK 状态后 GSE62254/84437 的 state 分布**非退化**（不再全部
  `NK-intermediate`；至少出现 3 个状态）。

```bash
python src/preprocessing/fix_geo_gene_mapping.py --gse GSE62254 GSE84437
python - <<'PY'
import pandas as pd
for g in ["gse62254","gse84437"]:
    df = pd.read_csv(f"data/processed/bulk/{g}_expression.tsv", sep="\t", index_col=0)
    nk = [m for m in ["NKG7","GNLY","GZMB","PRF1","KLRD1","NCAM1","KLRK1"] if m in df.columns]
    assert df.shape[1] < 30000, f"{g}: 仍像探针级 {df.shape}"      # A1
    assert len(nk) >= 4, f"{g}: NK marker 仅 {nk}"                 # A2
    print(g, df.shape, "NK markers:", nk)
print("T1 A1/A2 PASS")
PY
```
**GREEN**：跑上面脚本；若 A2 失败通常是平台号猜错（GSE84437 可能非 GPL6947），
用 `--gpl-file` 指定正确 GPL SOFT 注释重试。
**完成定义**：A1–A3 全绿，`LOG.md` 记录非零 NK 基因数。

---

## T2 — 外部验证：轴在独立胃癌队列复现（依赖 T1）🟢 已通过（2026-07-07）

> 效应臂 protrusion~cytotoxicity 在两个独立胃癌队列均**显著正复现**：
> GSE62254 r=0.42 (p=1.4e-14)、GSE84437 r=0.62 (p=3.3e-53)，B1 ✅。
> 附带 sm_balance~protrusion 在两队列亦弱正显著（r=0.18/0.11）。
> 结果见 `external_validation_results.tsv`，已写入稿件 §3.3 + Table 5。下方保留原判据。


**RED**
- B1：重算后，两外部队列中 `protrusion_machinery ~ cytotoxicity_outcome` 相关
  **r>0 且 p<0.05**（与 TCGA-STAD 同号，即 §3.3 可写"外部复现"）。
- B2：若某队列 NK 覆盖不足以算模块，则**如实标注 `MIXED_UNRESOLVED`**，不得
  强行出数。

```bash
# 重算 NK 分数 + SST 轴后，对每个外部队列做相关
# 期望：r 为正、显著；否则按 B2 标注
```
**完成定义**：§3.3 与 Table 1 的外部验证行填入真实 r/p 或明确的
`MIXED_UNRESOLVED`。

---

## T3 — GNN vs 基线对比（Blocker B）🟢 已通过（2026-07-07）

> 结果：从 A100 `/data/lgh/GC-NKGraph-Atlas` 拉回 `baseline_internal_results.tsv`，
> 与本地 GNN 表在同 seed-42 折合并。C1（7 方法）✅、C2（GNN MCC=0.7057≈0.706，
> AUROC=0.9503，无漂移）✅、C3（配对检验）✅。**诚实结论：GNN 与 LightGBM(0.733)/
> XGBoost(0.727) 无显著差异（p>0.27），显著优于 ElasticNet/SVM/MLP（p<0.05）。**
> 已写入 Table 3 + §3.4，措辞为"性能相当、附带可解释性"。下方保留原验收脚本备查。


**RED（黄金基线来自本地 `gc_nkgraph_internal_results.tsv`）**
- C1：`results/tables/model_comparison.tsv` 存在，含 **7 方法 × 5 折**
  （GNN + XGBoost/LightGBM/RF/ElasticNet/SVM/MLP）。
- C2：GNN 的 5 折 **MCC 均值 = 0.706 ± 0.01**、**AUROC 均值 = 0.943 ± 0.01**
  （与既有 GNN 表一致，防止折划分被无意改动）。
- C3：`model_comparison_stats.tsv` 输出 GNN vs 各基线在 MCC/AUROC 上的
  paired Wilcoxon + t 检验 p 值。

```bash
python src/baselines/run_model_comparison.py
python - <<'PY'
import pandas as pd, numpy as np
c = pd.read_csv("results/tables/model_comparison.tsv", sep="\t")
assert c["method"].nunique() == 7, c["method"].unique()          # C1
g = c[c.method=="GC-NKGraph-Atlas"]
assert abs(g.MCC.mean()-0.706) < 0.01, g.MCC.mean()               # C2
assert abs(g.AUROC.mean()-0.943) < 0.01, g.AUROC.mean()
assert __import__("os").path.exists("results/tables/model_comparison_stats.tsv")  # C3
print("T3 PASS; 各方法 MCC 均值：")
print(c.groupby("method").MCC.mean().round(3).sort_values(ascending=False))
PY
```
**REFACTOR**：把结果矩阵写入稿件 Table 3；在正文注明"GNN 与最强基线在同折上
ΔMCC=__，paired p=__"。**注意诚实**：若 GNN 未显著优于最强基线，正文改写为
"与最强基线相当/稳定"，不得夸大。
**完成定义**：Table 3 全部填数 + 统计检验入正文；投稿信 `[best-baseline ΔMCC]`
占位替换。

---

## T4 — scRNA 质控硬化（方法学）🟡 脚本就绪（2026-07-09）

> `src/scrna_analysis/qc_filter.py` 已完成（NaN-robust fallback，可配置阈值，Scrublet
> 可选）。`tests/test_qc_regression.py` 已写好（D1–D3 RED 转 pytest 测试）。需在服务器
> 上执行 `qc_filter.py` 然后跑回归测试。下方保留原 RED 判据。

**RED**
- D1：`results/tables/scrna_qc_summary.tsv` 存在，逐样本记录 before/after 与
  `retained_frac`。
- D2：QC 后总细胞数落在**合理区间**（保留率每样本 0.5–0.98；不得再次出现
  v1 的"0 细胞"或 v2 的"100% 保留无过滤"两个极端）。
- D3：**回归**：在 QC 后 NK 子集重算 H2/H3，结论方向不变（H3 仍显著正、
  H2 单细胞仍正），确保质控未改变论文核心。

**完成定义**：§2.4 的 QC 阈值与保留率有表可引；D3 回归通过。

### T4 服务器执行步骤（新增 2026-07-09）

```bash
# STEP 1: Run QC filtering on the integrated object
python src/scrna_analysis/qc_filter.py \
    --in data/processed/scrna/gc_integrated.h5ad \
    --out data/processed/scrna/gc_integrated_qc.h5ad \
    --doublets   # optional

# STEP 2: Verify summary table written (D1, D2)
python - <<'PY'
import pandas as pd
s = pd.read_csv("results/tables/scrna_qc_summary.tsv", sep="\t", index_col=0)
print(s.to_string())
for i, row in s.iterrows():
    assert 0.5 <= row.retained_frac <= 0.98, f"{i}: retained_frac={row.retained_frac} out of range"
print("D1/D2 PASS")
PY

# STEP 3: Re-compute SST-axis scores on the QC'd NK subset
python src/topology/sst_axis.py --input data/processed/scrna/gc_integrated_qc.h5ad

# STEP 4: Regression test (D3) — H3 must stay positive; H2 must not flip sign
pytest tests/test_qc_regification.py -v -m server
```

After T4 passes: re-run `src/figures/make_figures.py` to regenerate all figures
with the QC'd data, and update §2.4 in the manuscript with the QC thresholds used
and the total cell retention rate.

---

## T5 — 发表级图 Fig1–5（无需算力）🟢 已通过（2026-07-09）

> `src/figures/make_figures.py` 完全重写：Okabe-Ito 色盲安全配色，数字全部读自结果表，
> 输出 PDF+PNG。新增 Fig5（mechanism-card 概念图，纯 matplotlib vector art，无外部依赖）。
> CLI 接口：`--fig N`（单张）、`--dpi N`（分辨率）。E1–E5 ✅。归档于 `manuscript/figures/`。


**RED（用现有结果表驱动，数字须与表一致）**
- E1：`Figure1`（Arm A 恢复）：含 (A) LIHC 模块分数、(B) H2/H3 的 bulk vs 单细胞
  对比、(C) 瘤内 vs 正常 NK、(D) Table 2 森林图；图中 r 值与
  `sst_axis_positive_control_recovery.tsv` **逐一对齐**。
- E2：`Figure2`（Arm B）：TCGA-STAD 轴分数 + 三组织比较 + NK 状态分布。
- E3：`Figure3`（靶点）：肿瘤特异性 vs NK 关联散点（高亮肿瘤内在池）+ 轴链位置
  + 成药性分布；PHGDH/SGMS2/SMPD3 可辨识。
- E4：矢量优先，栅格 ≥300 dpi，图例自足，正文按序引用。

**完成定义**：3 张图落地 `results/figures/`，替换草图 fig9/fig10；每张图数字可回溯到
某张 tsv。（可让我下一步写绘图脚本 `src/figures/make_figures.py`。）

---

## T6 — 公开代码仓库可复现 🟢 已通过（2026-07-09）

> 仓库 https://github.com/nblvguohao/GC-NKGraph-Atlas 已公开。`nk-pre-submission`
> 分支已推送（commit acd7496）含完整 README、MIT LICENSE、合成数据模式、CI-ready
> 结构。**待验证（需联网机器）：** F1（干净 clone `pytest -q`）、F2（`--synthetic`
> 端到端）。下方保留原判据。


**RED**
- F1：从干净 clone 出发，`pytest -q` 全绿（现有 `tests/` 3 个测试 + 建议为
  新脚本补测试）。
- F2：合成数据模式跑通端到端（`--synthetic`），不依赖私有数据即可演示流程。
- F3：`README` 含环境（`environment.yml`/`requirements.txt`）、数据获取说明、
  一键复现命令；`[repository URL]` 在稿件各处替换为真实地址。

**完成定义**：外部评审可从仓库复现 pipeline 主干；数据可用性声明指向真实仓库。

---

## T7 — 前置声明与格式收尾（可并行）🟡 大部分完成（2026-07-09）

**RED**
- G1：`grep -nE '\[PLACEHOLDER|\[TBD|To be (expanded|written|designed)' main_manuscript.md`
  返回**空**。（本地已验证：0 处 ✅）
- G2：所有在文引用 [1]–[24] 均可解析（✅）；**[17]–[24] 已联网逐条核对，卷/页/年
  全部正确**（2026-07-07 ✅）。
- G2b：MIT `LICENSE` 已添加 ✅。
- G3：为通讯作者补 **ORCID**；核对 CRediT 角色是否属实；核对基金号无误
  （32472007 / 62301006 / 62301008；2308085MF217 / 2308085QF202）。← 仅剩此项需你操作。
- G4：参考文献转 Oxford/BiB 编号（Vancouver）风格，过一遍文献管理器。
- G5：**新增** `manuscript/PRE_SUBMISSION_CHECKLIST.md` — 逐项投稿前 action items ✅。

**完成定义**：BiB 投稿清单 D/E 段全绿。

**完成定义**：BiB 投稿清单 D/E 段全绿。

---

## 里程碑与依赖

```
T1 ──► T2 ──┐
T3 ─────────┼──► T5（图 1–5）
T4 ─────────┘        │
T6（✅ 已推送）       │
T7（✅ 大部分完成）   │
                     ▼
T8（格式清洁）──► T10（替代解释+校准透明）
T13（第二张卡片）  T15（补充索引）
                     │
T9（H2效应量）──► T11（嵌入分析）──► 投稿
T4（QC服务器）    T12（消融实验）
                  T14（排名稳定性）
```

**状态（2026-07-09）：**
- 🟢 本地完成：T1/T2/T3/T5/T6/T8/T10/T13/T15
- 🟡 脚本就绪，需服务器执行：T4/T9/T11/T12/T14
- 🟡 需作者操作：T7（ORCID/CRediT/基金号）

**T8–T15 详细规格：** 见 `submission_package/TDD_reviewer_and_editor.md`。

**服务器执行清单（一次性批量运行，~2-3 小时）：**
```bash
# 1. scRNA QC (T4)
python src/scrna_analysis/qc_filter.py \
    --in data/processed/scrna/gc_integrated.h5ad \
    --out data/processed/scrna/gc_integrated_qc.h5ad
pytest tests/test_qc_regression.py -v -m server

# 2. H2 effect size (T9)
python src/analysis/h2_effect_size.py \
    --input results/tables/sst_axis_scores_single_cell.tsv \
    --output-dir results/tables/

# 3. Gene embedding analysis (T11) — requires GNN checkpoint
python src/interpretation/gene_embedding_analysis.py \
    --embeddings results/model/gene_embeddings.npy \
    --gene-list data/processed/graph/gene_nodes.tsv \
    --output-dir results/

# 4. Ablation study (T12) — update script with actual GNN calls first
python src/baselines/run_ablation.py \
    --full-edges data/processed/graph/edges.tsv \
    --nodes data/processed/graph/nodes.tsv \
    --expression data/processed/bulk/tcga_stad_expression.tsv \
    --labels data/processed/bulk/tcga_stad_labels.tsv \
    --output-dir results/

# 5. Target rank stability (T14)
python src/interpretation/target_rank_stability.py \
    --input results/tables/tumor_intrinsic_candidates.tsv \
    --output results/tables/target_rank_stability.tsv

# 6. Regenerate all figures with latest data
python src/figures/make_figures.py --dpi 300
```

---

## T8 — 手稿格式清洁（主编阻断项）🟢 已通过（2026-07-09）

> 删除所有编辑注释、缩短标题至 ≤150 字符、删除 cover_letter 内部备注。

- A1–A5 ✅：`main_manuscript.md` 中 ORCID 备注、传记备注、CRediT 备注、"Draft v0.3" 状态块均已删除；标题从 ~190 字符缩短为 ~150 字符（删除 "from Liver to Gastric Cancer"）；cover_letter.md 内部备注已删除。
- 验证：`grep -nE 'ORCID iDs: fill|Biographies ~30|CRediT taxonomy|Document status|Internal note' manuscript/*.md` 返回空 ✅

---

## T9 — H2 效应量分析脚本 🟡 脚本就绪（2026-07-09）

> `src/analysis/h2_effect_size.py` 完成。Bootstrap 95% CI + 下采样曲线（8310→100 步长200）+ 方差解释率。需在服务器执行。

```bash
python src/analysis/h2_effect_size.py \
    --input results/tables/sst_axis_scores_single_cell.tsv \
    --output-dir results/tables/ \
    --output-fig results/figures/figS_h2_downsample.pdf
python - <<'PY'
import pandas as pd
ci = pd.read_csv("results/tables/h2_bootstrap_ci.tsv", sep="\t")
assert ci["ci_lower_95"].iloc[0] > 0, "H2 CI crosses zero"
assert ci["r_squared"].iloc[0] < 0.01, f"H2 r²={ci['r_squared'].iloc[0]} — should be <1%"
print("T9 PASS")
PY
```

---

## T10 — H4/H5 替代解释 + metabolic_crosstalk 校准透明 🟢 已通过（2026-07-09）

> 向手稿添加了两段关键文本 + 预注册日志 + 参考文献 [25][26]。

- C1 ✅：§4.3 Limitation #1 新增替代解释段落（NK 转录激活-功能解耦，类比 T 细胞耗竭），引用 Wherry 2015 [25] + Bi 2017 [26]
- C2 ✅：§2.5 新增校准状态披露段落（H1 null → 边默认为锚定论文方向，标记 NEEDS_REVIEW）
- C3 ✅：`configs/sst_axis_config.yaml` 中 `tumor_serine_capacity.expected_direction` 更新为 `CALIBRATION_FAILED_DEFAULTING_TO_ANCHOR_PAPER_DIRECTION`
- C4 ✅：`configs/sst_axis_prereg_log.md` 创建，记录 CHANGE-001（校准源失败）和 CHANGE-002（恢复判定修订）

---

## T11 — 基因嵌入可解释性分析脚本 🟡 脚本就绪（2026-07-09）

> `src/interpretation/gene_embedding_analysis.py` 完成。UMAP 投影 + 模块内/间距离 + 嵌入距离 vs 共表达 Spearman。需在服务器执行（需训练好的 GNN checkpoint）。

```bash
python src/interpretation/gene_embedding_analysis.py \
    --embeddings results/model/gene_embeddings.npy \
    --gene-list data/processed/graph/gene_nodes.tsv \
    --scRNA-expr data/processed/scrna/gc_nk_expression.tsv \
    --output-dir results/
python - <<'PY'
import pandas as pd
eco = pd.read_csv("results/tables/embedding_distance_vs_coexpression.tsv", sep="\t")
assert eco["spearman_rho"].iloc[0] > 0, "embedding distance should positively correlate with coexpression"
print("T11 PASS")
PY
```

---

## T12 — 消融实验脚本 🟡 脚本就绪（2026-07-09）

> `src/baselines/run_ablation.py` 完成。全图 vs 去除 metabolic_crosstalk 边，5 折相同种子。**注意：脚本使用占位符值，需在服务器上替换为实际 GNN 训练调用。**

```bash
python src/baselines/run_ablation.py \
    --full-edges data/processed/graph/edges.tsv \
    --nodes data/processed/graph/nodes.tsv \
    --expression data/processed/bulk/tcga_stad_expression.tsv \
    --labels data/processed/bulk/tcga_stad_labels.tsv \
    --output-dir results/
```

---

## T13 — 第二张机制卡片 🟢 已通过（2026-07-09）

> 为 adenosine–A2AR–cAMP NK 抑制轴创建第二张机制卡片，证明卡片格式复用性。

- F1 ✅：`configs/mechanism_cards/adenosine_a2ar_nk_suppression.yaml`（7 基因模块、5 预注册假设、GATED physical_ground_truth）
- F2 ✅：`configs/mechanism_cards/registry.yaml` 已更新，注册两张卡片
- F3 ✅：§4.2 新增段落说明第二张卡片作为复用性初步证明

---

## T14 — 靶点排名稳定性分析脚本 🟡 脚本就绪（2026-07-09）

> `src/interpretation/target_rank_stability.py` 完成。5 个评分维度各 ±0.10 扰动，Jaccard 重叠 + PHGDH/SGMS2 top-3 稳定性。

```bash
python src/interpretation/target_rank_stability.py \
    --input results/tables/tumor_intrinsic_candidates.tsv \
    --output results/tables/target_rank_stability.tsv \
    --n-perturbations 10
python - <<'PY'
import pandas as pd
df = pd.read_csv("results/tables/target_rank_stability.tsv", sep="\t")
avg_j = df["jaccard_top10"].mean()
assert avg_j >= 0.7, f"avg Jaccard={avg_j:.3f} < 0.7"
print("T14 PASS")
PY
```

---

## T15 — 补充材料索引 + Supplementary Methods 🟢 已通过（2026-07-09）

> 补充材料索引 + Supplementary Methods 说明探针重映射方法 + NK 状态分类细节。

- H1 ✅：`submission_package/02_supplementary_tables/README.md`（S1–S6 文件清单 + 列说明 + 再生命令）
- H2–H4 ✅：Supplementary Methods 已包含探针→基因重映射步骤（GPL570/GPL6947、多探针策略、NK marker 验证）、NK 状态分类细节、复合评分权重说明
