# 投稿前修改意见与执行清单

> 依据 2026-07-07 模拟同行评审整理。目的不是再扩大战果，而是把论文改成“证据强度与主张强度一致”的可投稿版本。

## 总体判断

这篇稿件最稳的定位是一个机制驱动的计算框架和边界图谱：它说明一个已发表的 NK 细胞 serine-SM-topology 免疫逃逸机制中，哪些层能从公开转录组中恢复，哪些层只能在单细胞分辨率下弱检测，哪些层必须依赖物理/代谢实验。投稿前应避免把论文包装成“完整重建物理机制”或“GNN 显著优于强基线”的论文。

## 必须修改的主张边界

| 审稿风险 | 投稿前处理 |
|---|---|
| 标题和摘要中的 “reconstructing/reconstructs” 容易被理解为完整重建机制 | 改为 “mapping the transcriptional reach” 或 “transcriptionally recoverable layers” |
| H2 单细胞相关 r=0.030，统计显著但效应量极小 | 写成 “weak but significant / weakly detectable after cell-type resolution”，不要写成强恢复 |
| H4 和 H5 protrusion 方向错误削弱 topology recovery | 明确作为核心 scoping result：physical topology is not captured by machinery transcription |
| GNN 低于 LightGBM/XGBoost | 主张改为 “on par with strong tree baselines, added mechanism-structured embeddings”，不要声称模型性能优势 |
| 37 个靶点没有实验验证 | 全文使用 “putative tumor-intrinsic candidate targets”，并保持 “recommended assays” 而非 “validated targets” |
| in-silico SM-restoration stratification 未验证 | 只作为 secondary hypothesis，不作为主图或主要贡献 |

## 需要补强的技术说明

1. 在 Methods 或 Supplement 中补充模块定义和 target list 的敏感性分析：至少说明候选靶点排序是否依赖某一组权重或 gene module。
2. 对 scRNA 注释和 QC 给出可审计说明：阈值、marker、批次校正、NK cell 数量、QC 前后是否影响主要结果。
3. 对 `metabolic_crosstalk` edge 加一句限制：它是 mechanism-grounded / hypothesis-weighted edge，不是因果证明。
4. 在模型比较处强调 identical folds 和 paired tests，同时说明图模型价值来自机制化 embedding，而非 raw accuracy gain。
5. 对外部 GEO 队列说明 probe-to-symbol 修正、NK marker coverage 和平台差异，避免 reviewers 质疑微阵列映射。

## 投稿包清理清单

| 项目 | 状态 | 说明 |
|---|---|---|
| 正文标题、摘要、Key Points 降低强主张 | Done | 已从 “reconstructing” 改为 “mapping transcriptional reach” |
| pre-registration recovery definition 改为 partial-recovery scoping | Done | 已说明 full criterion 未满足 |
| cover letter 删除内部备注并同步新标题 | Done | 已将强主张改为 bounded scoping contribution |
| README 从 “submitted/reconstructs” 改为 “prepared/maps” | Done | 更符合投稿前状态 |
| main_claims 中胃癌外部验证状态 | Done | 已从 Pending 改为 Done |
| ORCID、CRediT、基金号 | Author action | 需要作者逐项确认 |
| suggested reviewer 邮箱与利益冲突 | Author action | 需要提交前人工核验 |
| References [17]-[24] 细节 | Author action | 需要逐条核对页码、DOI、期刊格式 |
| PDF/Word 最终排版与字符编码 | Author action | 需要确认无乱码、无内部状态块、图表引用正确 |

## 推荐投稿叙事

一句话版本：

GC-NKGraph-Atlas turns a published wet-lab immune-evasion mechanism into a reusable, single-cell-informed computational workflow that maps which layers of the serine-SM-topology axis are visible from public transcriptomes and which require physical or metabolite-level validation.

投稿时应反复守住三条边界：

- Transcriptomic proxy is not physical membrane topology.
- Weak single-cell metabolic coupling is not a strong mechanistic recovery.
- Candidate targets are hypotheses for wet-lab testing, not validated therapeutic targets.
