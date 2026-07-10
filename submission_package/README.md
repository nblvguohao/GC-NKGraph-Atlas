# 2026-07-10 更新汇总

本文件夹收录 2026-07-10 当日对 GC-NKGraph-Atlas 的全部更新/新增产出，含**两条工作线**：
- **工作线 A — 三篇新文献落地修回**（见下）。
- **工作线 B — Nature 风格投稿前自审 → 第二轮 TDD → T16 实测**（见文末"工作线 B"）。

---

## 工作线 A — 三篇新文献落地修回

把 `paper/` 下三篇论文（P1/P2/P3）作为**外部佐证/并行证据/范式参照**写入稿件，
并新增参考文献 [46]/[47]/[48]。三篇均**非** SST 同机制，正文未并入 SST 轴（机制护栏）。

## 涉及的三篇文献
| 代号 | 引文 | DOI |
|------|------|-----|
| P1 [46] | Bai et al. *Nature* 2024;634:702–11 | 10.1038/s41586-024-07762-w |
| P2 [47] | Rafei et al. *Nature* 2025;643:1076–86 | 10.1038/s41586-025-09087-8 |
| P3 [48] | Sun et al. *Nat Immunol* 2026;27:985–99 | 10.1038/s41590-026-02471-0 |

## 文件清单
| 文件 | 说明 | 状态 |
|------|------|------|
| `TDD_literature_integration.md` | 落地用 TDD（L1–L9）+ GREEN 执行记录 | 新增 |
| `main.tex` | LaTeX 投稿实体（已插入引用/讨论/limitation/参考条目） | 更新 |
| `main_manuscript.md` | Markdown 源稿（与 tex 成对同步） | 更新 |
| `main.pdf` | 编译产物，pdflatex 两遍 exit 0、15 页、0 未定义引用 | 新增 |
| `analysis_memo_three_papers.md` | 三篇论文增益分析备忘录 + 落地方案 | 新增 |
| `preedit_backup/` | 改前快照（`main.tex` / `main_manuscript.md`），供 diff 对照 | 备份 |

## 改动落点（对应 TDD L1–L9）
- §3.2 H4/H5：P3 外部佐证瘤内 NK 细胞毒下调（标注非本项目队列）。
- §4.1 讨论 + §4.4 未来：P3 的 TCGA-HCC 检查点共表达支撑 HAVCR2/Tim3 分层；CLEC12B-LPL 列为下一张机制卡片。
- 相关工作(dysfunction 段) + §3.2 激活对照：P2 支撑“TF/胞内调控被忽视” + 佐证 T14 扣激活。
- §4.3 Limitations 新增第 8/9 条：候选缺 TF 层(P2)、缺临床结局锚点(P1)。
- 相关工作(Scissor 段) + §3.7 消融：P1 图谱→功能态→结局范式 + L-R 建模先例。
- 参考表尾追加 [46][47][48]（承 T8：整体按首次出现序重排时并入）。

## 验证
- `pdflatex` 两遍 exit 0；`main.pdf` 15 页（原 13→15）；`main.aux` 中 `bibcite{ref46/47/48}` 全解析；0 未定义引用/引文。
- 机制护栏目检通过：三篇均为外部/并行证据，未并入 SST 鞘磷脂轴。

## 备注
- H5 现已用 GSE246662 实算（瘤内 vs 正常 NK），**非** `DATA_UNAVAILABLE`；故 P3 定位为外部佐证而非填数据空洞。
- 完整版滚动备份仍在 `submission_package/backups/`（含改前/改后时间戳快照）。

---

## 工作线 B — Nature 风格投稿前自审 → 第二轮 TDD → T16 实测

用 `nature-reviewer` 技能对 v0.5 主稿做了 3 审稿人 + 交叉综述（中英双份），把新增科学缺口
整理成第二轮 TDD（T16–T21，续接既有 T1–T15），并**本地实测**其中最致命的 T16。

### 文件清单
| 文件 | 说明 |
|------|------|
| `nature_reviewer_report_EN.md` | Nature 风格三审 + 交叉综述（英文） |
| `nature_reviewer_report_ZH.md` | 同上（中文） |
| `TDD_round2_T16-T21.md` | 第二轮 TDD（T16–T21）+ 本地执行进度 |
| `T16_specificity_finding.md` | T16 效应臂机制特异性：实测 + 深度控制复核 + 建议 |
| `results/t16_specificity_control.tsv` | T16 初测（仅控激活）汇总 + 判决 |
| `results/t16_null_distribution.tsv` | T16 初测 1000 次零分布 partial-r（补充图数据） |
| `results/t16_robustness.tsv` | T16 复核：2 控制 × 4 匹配 = 8 变体 |
| `scripts/run_t16_specificity_control.py` | T16 主检验脚本 |
| `scripts/run_t16_robustness.py` | T16 深度控制 + 容差扫描复核脚本 |

### 已落到稿件的纯文本收口（`main.tex` + `main_manuscript.md` 已同步）
- **T21** §4.1 "validates … design premise" → "consistent with … but given r²<0.001 does not by itself validate"（顺带修掉 Conclusion 重复词 "principled, principled"）。
- **T20** Conclusion "reusable engine" 下调为 "demonstrated on a single card here; reuse … remains a design aspiration"。
- **T18-R2** 摘要 + Key Points 前置 "near-zero fold-change / ~46% NK-side" caveat。

### T16 关键结论（决定性）
- 初测（仅控激活）：观测 partial r=0.251 落在匹配零分布正中（p≈0.42）→ 一度 NOT_SPECIFIC。
- 复核（控激活 **+ 测序深度**）：观测 0.159、零分布塌到 0.06，经验 **p≈0.006 → SPECIFIC**，对匹配容差全稳。
- **解读：随机模块共变主要由测序深度混杂抬高；剔除深度后真实效应臂耦合显著特异 → `effector arm recovers` 论断成立。** 建议**不下调**措辞，改为新增深度-控制的特异性检验强化 H3（详见 `T16_specificity_finding.md`）。

### 待办（需用户 green-light 或后续算力）
- 把深度-控制特异性检验写入 §3.2 + 补充图（T16 收口，待用户确认）。
- T17（边外样本价值，graph 数据本地已具，可下轮本地试跑）、T18 的 FDR/CI + Table4↔Fig3B、T19（第二套 scRNA + 参考映射，本地暂无数据）。
