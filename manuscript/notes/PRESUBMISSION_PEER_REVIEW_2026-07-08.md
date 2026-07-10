# 投稿前模拟同行评审（Nature 风格审稿人视角）

> 日期：2026-07-08　｜　对象：GC-NKGraph-Atlas（实际投向 *Briefings in Bioinformatics*）
> 说明：以通用大刊标尺做的**审稿人视角**评估，非编辑决定；不断言是否适合 *Nature*。
> 仅基于稿件事实（源 Markdown + 编译 LaTeX + Fig 0–4 + Table 1–5 + 45 篇文献）。
> 结构：3 份审稿意见 + 交叉综述 + 风险清单。中文为主，文末附英文原版。

---

## 评审设置

- **输入范围：** 全文稿、编译 LaTeX、五张图、Table 1–5、45 篇文献。
- **评估边界：** 不可评估代码/仓库可解析性、补充材料、原始数据、数字复现；湿实验验证按作者设计缺失。
- **共享论断：** 量化 Zheng 2023 丝氨酸→鞘磷脂→膜拓扑→细胞毒轴的"转录可及度"；mechanism-card 抽象 + 含 `metabolic_crosstalk` 边的异质图 + GNN；两臂设计给出 scoping 结论（效应臂可恢复、代谢臂仅单细胞层、拓扑不可恢复）；输出 37 个去循环肿瘤内在候选。
- **证据基础：** TCGA-LIHC(423)+8,310 scNK；TCGA-STAD(450)+胃癌 scNK(1,017)+GSE62254(300)/GSE84437(483)；5 折 CV（GNN MCC 0.706/AUROC 0.950，与 LightGBM/XGBoost 相当）；Table 4 靶点。
- **缺失材料：** 无实验验证；无消融；H1 未报；§3.6 无数字；可复用性仅一张 card；膜拓扑从未直接测量；embedding 可解释性只被断言。

---

## 审稿人 1 — 侧重：生物学意义与免疫学

- **总体：** 谨慎、如实的计算再检验，定位克制值得称道；但核心生物学增量偏薄——唯一"稳健恢复"的一层可能近乎自带，代谢臂头条效应量微不足道，且无一结果经实验验证。
- **谁会感兴趣：** NK 与肿瘤免疫代谢研究者；将成像/物理机制外推到转录组规模的团队。专门化。
- **主要优点：** 明确区分转录代理与物理表型；逐假设如实报告（含失败）；效应耦合在两外部胃癌队列独立复现；NK 读出与靶点做了去循环分离。
- **主要担忧：**
  1. **H3 可能是定义性重叠而非验证。** "细胞毒输出"模块（NKG7/GNLY/GZMB/PRF1/IFNG + LAT/VAV1/TLN1/ITGAL/ITGB2/LCP2）与"突起机器"模块都是 NK 激活/突触程序的侧面，正相关可能只是共激活。须剔除通用激活签名后再看。
  2. **被"救回"的代谢耦合可忽略。** 单细胞 H2 r=+0.030，共享方差 <0.1%，p 值是 n=8,310 的产物。"细胞类型分辨率是决定性的"夸大了微弱效应。
  3. **同名生物表型从未被测量。** 膜拓扑/SM 含量按设计排除，本文无法真正谈论其命名机制。
  4. **靶点可信度。** Table 4 肿瘤特异性 log2FC≈0.001–0.06（实质为零），"肿瘤内在"未确立。
- **须先解决：** 对 H3 控制激活签名；报告效应量并重述 H2；补报 H1；解释近零 log2FC 为何算肿瘤内在。
- **对标：** 原创性中；重要性中/有限；跨学科中；技术稳健混合；可读性尚可。
- **姿态：** 专门期刊大修；不足以达广读者显著性门槛。

## 审稿人 2 — 侧重：计算与方法学严谨性

- **总体：** 成熟组件的合理重组 + 一条有生物依据的专设边 + 有价值的纪律；但作者数字动摇其主旨：GNN 无准确率增益，其被声称贡献的"机制结构化图"从未经消融。
- **谁会感兴趣：** 图/多组学方法与单细胞归因开发者，冲着 mechanism-card 范式与去循环协议。
- **主要优点：** 同折基线 + 配对检验 + 诚实"相当、非 SOTA"；预注册；QC/注释透明。
- **主要担忧：**
  1. **无消融 → 核心方法论主张无支撑。** 无 `metabolic_crosstalk`/SST 边 on/off 实验，也无 embedding 下游效用 vs 树基线归因的对比。图的必要性是断言。
  2. **去循环疑未做净。** Fig 3C 显示"肿瘤内在池(n=37)"被 NK 侧模块主导（protrusion machinery 7 / de novo sphingolipid 5 / sm catabolism 4），RAC1 亦在候选，与"分离"表述冲突。
  3. **统计。** ~8–9 项检验无多重校正；"显著"反复由大 n 在可忽略效应量上驱动。
  4. **模型含糊。** 编码器视 torch_geometric 在谱分解/HGT 间二选一，Table 3 由哪个产生不明；边权硬编码无依据无敏感性。
  5. **可复用性仅一张 card。**
- **须先解决：** 图/边消融；定量 embedding 增益；多重校正 + 突出效应量；说明 Table 3 模型；修 Fig 3C；边权敏感性。
- **对标：** 原创性中；图贡献未确立；重要性取决于消融；可读性术语密集。
- **姿态：** 大修；无消融则方法论不成立，现稿不可推荐。

## 审稿人 3 — 侧重：广泛兴趣、意义定位与非专家可读性

- **总体：** 流畅、自省，诚实是最大长处；但贡献偏增量且部分为阴性，可泛化的核心思想论证不足。对通用读者增量不足。
- **谁会感兴趣：** 相对狭窄的计算肿瘤/NK 读者；若证明可泛化，mechanism-card 或吸引更广方法读者。
- **主要优点：** 叙事清晰；"转录可及度地图"是有用重构；图 1/3/4 专业色盲友好；声明完整。
- **主要担忧：**
  1. **意义温和且部分阴性**（近恒真相关 + 微弱代谢效应 + 一个真阴性）。
  2. **新颖性重组式**（复用 PPI/LR/TF + HGT，新元素是一条边 + 流程纪律）；纪律才是真新颖，应前置并跨 ≥2 机制演示。
  3. **可读性**：术语密集、正文含文件/脚本路径、反复"honest"自评显防御性。
  4. **呈现缺陷侵蚀信任**：表/图编号乱、缺模型比较图、符号乱码、工作流阶段数(4)与正文(5)矛盾。
- **须先解决：** 证明通用性；修呈现缺陷；显著性与效应量相称；核验可用性链接。
- **对标：** 原创性中；重要性有限至中；跨学科中；可读性可改进。
- **姿态：** 现稿不支持通用期刊接收；修回后适合专门生信方法期刊。

## 交叉综述

- **一致优点：** 罕见诚实；效应耦合独立外部复现；核心图整洁；声明完整。
- **一致技术风险：**（1）无消融；（2）正向结论被夸大（H3 定义重叠、H2 可忽略效应）；（3）去循环不自洽（Fig 3C）；（4）无实验验证、同名表型未测量；（5）可复用性仅一张 card。
- **侧重差异：** R1 压生物学、R2 压方法、R3 压意义与呈现；三人皆归"大修"。
- **意义读数：** 价值真实但专门化；最具广度的 mechanism-card 恰最欠演示；现稿低于通用大刊广泛兴趣门槛。
- **最需先解决：**（i）图/边消融 + 定量 embedding 增益；（ii）H3 控激活签名、重述 H2；（iii）自洽 Fig 3C；（iv）补 H1、补齐/删 §3.6；（v）card 之外证可复用；（vi）修编译-打包缺陷。

## 风险 / 无支撑或不可评估

- "可复用引擎/可泛化" —— 仅一机制演示，未确立。
- "图 embedding 带来可解释性价值" —— 无消融，未确立。
- "细胞类型分辨率是决定性的" —— 效应量 r=0.03 可忽略，论断强度无支撑。
- "37 个肿瘤内在候选" —— 肿瘤特异性近零 + Fig 3C 冲突，未确立。
- H1 结果 / §3.6 数字 / 代码复现 / 补充方法 —— 不可评估。
- 是否适合 *Nature* —— 不可评估 / 非编辑决定。

---

<details>
<summary>English original (click to expand)</summary>

See conversation transcript 2026-07-08 for the full English version of this 3-reviewer + synthesis assessment. Key grounded critiques (identical substance):
- **R1 (biology):** H3 may be definitional overlap of an NK-activation program, not axis validation; H2 single-cell r=0.030 is a negligible effect inflated by n=8,310; named membrane-topology phenotype never measured; Table 4 tumor-specificity log2FC ≈0.
- **R2 (method):** no ablation for the graph/edges → core methodological claim unsupported; Fig 3C shows NK-side modules dominating the "tumor-intrinsic" pool (de-circularization leak); no multiple-testing control; encoder ambiguity (spectral vs HGT); reusability shown on one card.
- **R3 (broad interest):** significance modest/partly negative; novelty recombinative; jargon + self-labeling ("honest"); compile-package defects erode trust.
- **Synthesis / most important:** (i) edge ablation + quantify embedding value; (ii) control H3 for activation, reframe H2 by effect size; (iii) reconcile Fig 3C; (iv) report H1, fill/remove §3.6; (v) demonstrate reuse beyond one card; (vi) fix package defects.
- **Not assessable:** journal fit (incl. Nature), H1, §3.6 numbers, code reproducibility, supplementary.

</details>
