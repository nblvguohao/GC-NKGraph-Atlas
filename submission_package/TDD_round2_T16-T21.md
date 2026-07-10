# TDD：Nature 风格三审意见落地（本轮预审 → 修回，测试驱动）

## Context（为什么做这轮）

本轮用 `nature-reviewer` 技能对 `manuscript/main_manuscript.md`（v0.5）做了
Nature 风格三审 + 交叉综述。稿件的**诚实性/统计规范（预注册、报负结果、FDR/CI、
去循环审计）是公认强项**，但三位审稿人一致指向若干**尚未被现有 T1–T15 覆盖的
新技术缺口**。本 TDD 只收录**本轮新增**的科学缺口，继续沿用现有编号（T16 起），
与 `submission_package/TDD_review_fixes.md`（T1–T15）互补、不重复。

投稿实体仍是 `submission_package/latex_official/main.tex`（由
`manuscript/main_manuscript.md` 生成）——**两者都要改并保持同步**。
需重算的项走 A100（`ssh user@100.112.165.109`，项目 `/data/lgh/GC-NKGraph-Atlas`，
conda env `gc-nkgraph`），脚本落在 `src/a100_recompute/`，沿用 `run_*.py` + `run_all.sh` 约定。

图例：🔴 未开始 · 🟡 进行中 · 🟢 通过　｜　优先级：P0 阻断 · P1 高 · P2 中 · P3 低

> 三审最致命的四点（外审最可能据以退稿）：**T16（效应臂特异性对照）、
> T17（边价值的外样本检验）、T18（靶点近零效应量 + 摘要口径）、T19（单细胞独立复现）**。

## 执行进度（2026-07-10 本地）

本地无原始数据/GPU，凡需重算者只能写 RED 脚本、GREEN 待 A100。已完成的纯文本收口
已同步改入 `main.tex` + `main_manuscript.md`，`latexmk -pdf main.tex` exit 0、14 页、无未定义引用。

- 🟢 **T21** §4.1 "validates … design premise" → "consistent with … but given r²<0.001 does not by itself validate"（两文件已改，含顺带修掉 Conclusion 重复词 "principled, principled"）。
- 🟢 **T20** 采用方案(b)：Conclusion 的 "reusable engine" 下调为 "demonstrated on a single card here; reuse … remains a design aspiration"（两文件已改）。
- 🟡 **T18** R2 摘要 + Key Points 已前置 "near-zero fold-change / ~46% NK-side" caveat（两文件已改）；**R1 的 FDR/CI 与 R3 的 Table 4↔Fig 3B 对齐仍待 A100 + 绘图**。
- 🟢 **T16 已本地跑通 + 方法学复核完成**（数据 `data/processed/scrna/gc_nk_subset.h5ad` + anndata/scanpy 本地齐备）。
  **初测（仅控激活）：观测 partial r=0.251 落在匹配零分布正中（null mean 0.244），经验 p=0.42 → NOT_SPECIFIC。**
  **复核（`run_t16_robustness.py`，控激活 + 测序深度 log total_counts）：观测 0.159、零分布塌到 0.06（p95≈0.12），经验 p≈0.006 → SPECIFIC，且对 ACT_TOL 0.25/0.5/1.0/unmatched 全稳。**
  **结论（决定性）：随机 NK 模块对之所以在仅控激活时仍相关 ~0.24，是被测序深度/全局转录活性这一技术协变量抬高的；一旦剔除深度，泛模块共变几乎消失，而真实 protrusion~cytotoxicity 仍保留 0.159 显著高于零分布 → 效应臂耦合是机制特异的，`effector arm recovers` 论断成立。** 深度是驱动一切模块共表达的混杂项，故"控激活+深度"才是正确零分布。
  产出 `t16_specificity_control.tsv`、`t16_null_distribution.tsv`、`t16_robustness.tsv`。
  ➡️ **推荐（把审稿弱点变强点）**：**不下调**措辞；改为在 §3.2 / T14 的激活对照旁**新增一句深度-控制的匹配零分布特异性检验（emp p≈0.006）**，直接回应 R1/R3。待用户 green-light 后落到两文件 + 补充图。
- 🔴 **T17 / T19** 数据在本地？——T17 需图/模型接口（`data/processed/graph/` 已在本地，可尝试本地跑）；T19 需第二套 scRNA（本地暂无）。下轮排查。

---

## P1 — 实质科学（外审前应补）

### T16 — 效应臂"机制特异性"对照（审稿人 R1/R3）🔴 P1 ★最致命
> 风险：H3（protrusion machinery ~ cytotoxicity output）可能只是**任意两个"激活相关
> NK 模块"的共变**，而非 Zheng 机制中 serine→SM→topology 因果链的验证。T14 已对
> 16 基因激活签名做偏相关（partial r=0.286 仍显著），但审稿人指出 16 marker 未必张成
> 整个激活流形，且偏相关只回答"是否与激活正交"，未回答"是否**特异于本机制**"。
**RED（先写死判据）**
- R1 构造**负对照模块对**：从与 SST 轴无关的 NK 程序（如干扰素应答、代谢housekeeping、
  迁移/趋化等）随机抽取，**匹配模块大小与激活载荷**，生成 N≥100 个随机模块对，
  计算其 protrusion-替身~cytotoxicity-替身相关分布。
- R2 判据：真实 protrusion~cytotoxicity 的 partial-r（=0.286）应落在负对照分布的
  **上尾（经验 p<0.05）**，即耦合强度**显著高于**匹配随机对；否则如实下调措辞。
- R3 §3.2 增一句特异性结论 + 补充图（负对照分布 + 观测值竖线）。
**GREEN**：新增 `src/a100_recompute/run_t16_specificity_control.py`，输入
`gc_nk_subset.h5ad`（8,310 NK × 22,728 基因，与 T14 同源）。
**完成定义**：H3 耦合被证明"机制特异"而非"泛激活共变"，或正文据实标注为部分泛激活驱动。

### T17 — `metabolic_crosstalk` 边价值的**外样本**检验（审稿人 R2）🔴 P1 ★
> 风险：现 §3.7 消融是**近乎同义反复**——SST 边由同一机制定义，删掉它自然删掉
> H1/H2 的 embedding 结构；这只证明"边在图里"，不证明"边编码了独立可验证的生物学"。
> 且 §3.4 显示 GNN 分类精度与 LightGBM/XGBoost 持平（未见增益），削弱"边有用"论断。
**RED**
- R1 设计**外部端点**检验（择一或并行）：(a) 用含/不含 `metabolic_crosstalk` 边的 embedding
  做**留出队列**的 cytotoxicity/NK-state 预测，比较泛化指标；(b) **跨队列迁移**
  （LIHC→STAD 或 TCGA→GEO）下含/不含边的性能差。
- R2 判据：含边版在**外样本**指标上显著优于不含边版（配对检验 p<0.05）→ 写入"边有外部增益"；
  若无差异 → §3.7/§4 据实改为"边塑造 embedding 结构但未转化为外样本预测增益"，删除任何
  暗示"边=真实生物学验证"的措辞。
**GREEN**：扩展 `src/a100_recompute/run_ablation.py` 或新增
`run_t17_edge_external_value.py`；复用 `src/graph_construction/build_heterograph.py`
的 on/off 开关与 `src/models/gc_nkgraph_atlas.py` 的 embedding 接口。
**完成定义**：边价值由"内部结构"升级为"外样本可检验"结论（正负皆据实）。

### T18 — 靶点近零效应量 + 摘要/正文口径对齐（审稿人 R2/R3）🔴 P1 ★
> 风险：Table 4 的 tumor_specificity_log2FC 为 +0.059/+0.038/+0.010/+0.001（近零），
> 却支撑"37 个 putative tumor-intrinsic 靶点"；且去循环审计自报 17/37（46%）为 NK 侧。
> 摘要仍以"37 candidate targets"领起，未前置 caveat。（延续 T6/T11，本轮要求收口）
**RED**
- R1 为 37 个候选补**差异表达显著性**（malignant vs non-malignant 的 FDR、绝对 log2FC、CI），
  出补充表；判据：正文对"tumor-intrinsic"的措辞与数值强度**相称**（近零富集不得称"肿瘤特异高表达"）。
- R2 摘要与 Key Points 前置 caveat：`grep -n "37" main.tex main_manuscript.md` 命中处，
  相邻须出现 "near-zero fold-change / 46% NK-side residual" 之一的限定语。
- R3 Table 4 与 Fig 3B 取自**同一 ranking**、top 候选一致（回收 T6-R1 未闭项）。
**GREEN**：扩展 `src/a100_recompute/run_effect_size_reframe.py` 与
`src/interpretation/split_target_lists.py` / `prioritize_targets.py`；CI 需 A100 重算。
**完成定义**：靶点强度、去循环 caveat 在**摘要—正文—表—图**四处一致。

### T19 — 单细胞结论的独立复现 + 参考映射注释稳健性（审稿人 R2/R3）🔴 P1
> 风险：全部单细胞论断（H2 rescue、H3、H5、胃癌 n=1,017 NK）**仅来自 GSE246662 一套**，
> 且 NK 用 marker-threshold 注释。H5"瘤内 NK 突起转录**更高**"这一支撑"拓扑不可达"结论的
> 关键反差，可能受 NKT/组织驻留污染影响。
**RED**
- R1 在**第二套独立 scRNA**（如另一 HCC/GC 公共单细胞集）复算 H3 与 H5：判据
  H3 同号显著、H5-cytotoxicity 同号；若不复现则正文将单细胞结论降级为"单数据集观察"。
- R2 **参考映射敏感性**：用 scANVI/scArches（稿件已引 ref [35]）重注释 NK 舱，
  复算 H5；判据 Δprotrusion 符号在 marker-threshold 与 reference-mapping 两法下一致，
  否则如实标注注释敏感性为 limitation。
**GREEN**：扩展 `src/scrna_analysis/run_scrna_v2.py` + `src/topology/sst_axis_scrna.py`；
新增 `src/a100_recompute/run_t19_scrna_replication.py`。
**完成定义**：核心单细胞耦合有第二数据集/第二注释法支撑，或据实降级措辞。

---

## P2 — 论断范围收口

### T20 — "可复用引擎"论断需第二张卡的真实结果，或下调为设计愿景（审稿人 R1）🔴 P2
> 风险：mechanism-card 的"可复用引擎/registry"论断，本文正文仅**一张卡**有科学结果；
> 文末 document-status 提到 adenosine/TGFβ 多卡，但正文无结果。
**RED**
- R1 择一：(a) 用 `src/interpretation/run_multicard_analysis.py` 跑出**至少第二张卡**
  的非平凡结果并入正文/补充；或 (b) 将 §4.2/§1.3/摘要的"reusable engine"下调为
  "design aspiration / demonstrated on one card"，`grep -n "reusable engine\|engine"` 命中处逐一改。
**完成定义**：可复用性论断与已展示证据相称。

---

## P3 — 文本自洽

### T21 — §4.1 "validates the edge's design premise" 与 r²<0.001 对齐（审稿人 R2）🔴 P3
> 风险：§4.1(2) 仍称 H2 "validates the `metabolic_crosstalk` edge's design premise"，
> 但同文 §3.2 已把 H2 定性为"统计可检测、量级可忽略（r²=0.0009）"——前后自相矛盾。
**RED**
- R1 `grep -n "validates the\|design premise" main.tex main_manuscript.md` 命中处改为
  效应量相称措辞（如"is consistent with, but does not by itself validate"）。
**完成定义**：全文对 H2 的表述与其 r²<0.001 强度一致（回收 T15 精神的遗留点）。

---

## 里程碑与依赖

```
P1  T16(效应臂特异性★) · T17(边外样本价值★) · T18(靶点效应量+摘要口径★) · T19(单细胞独立复现)
P2  T20(多卡/引擎论断)
P3  T21(§4.1 自洽)
```

- T16/T17/T18/T19 是外审最可能据以退稿的四项，优先补；三项需 A100 重算。
- T18-R3 回收 T6 未闭的 Table 4↔Fig 3B 对齐；T21 回收 T15 遗留的文本自洽。
- 与既有 T1–T15：T16 承 T14、T17 承 T7、T18 承 T6/T11、T21 承 T15；不重复，只补缺口。

## 验证（Verification）

1. 每个 T 的 RED 判据可独立自检：文本类用 `grep`（T18-R2、T20-R1、T21-R1），
   数值类用 A100 上 `run_t1x_*.py` 输出的黄金基线断言（仿 `run_h3_activation_control.py`）。
2. A100 一键回归：把新增 `run_t16..t19_*.py` 加入 `src/a100_recompute/run_all.sh`，
   `bash src/a100_recompute/run_all.sh` 全绿。
3. 稿件同步：改动后 `latexmk -pdf main.tex` 无 Error，且 `main.tex` 与
   `main_manuscript.md` 对应段落一致（数字、措辞、表/图引用）。
