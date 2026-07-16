# 修回任务卡 · M5 — bulk 效应臂相关的 NK 浸润比例混杂控制

> 来源:`/nature-reviewer` 三审(2026-07-16)交叉综合的**唯一共识技术漏洞**。
> 目标稿件:`submission_bundle_BiB/01_manuscript/main_manuscript.md`。
> 与 `TASK_CARDS_M1-S5.md` 同体例(M=必须/阻断项),编号续接为 M5。
> **边界声明**:卡内脚本/函数/字段名来自对本仓库源码的实际抽查(`src/…`);
> 产出文件名部分为脚本内既有路径、部分为建议新增(标 *新增*)。

---

## 0. 一句话问题

稿件把 **bulk TCGA-LIHC 的 H3 相关(protrusion~cytotoxicity, r=0.55)** 明确定为效应臂的
**头号证据**(单细胞 H3 已因混杂被自行降级),但在 bulk 肿瘤转录组里,
`nk_protrusion_machinery_score` 与 `nk_synapse_cytotoxicity_outcome_score` **两个模块打分都会随样本的 NK 浸润比例同涨同落**。
因此这个正相关可能部分来自 NK 丰度协变,而非 NK 细胞**内**的真实机制耦合。
稿件目前未展示对此已做控制——三位审稿人一致认为这是投稿前应闭合的技术项(BiB 审稿大概率也会问)。

---

## 🔴 M5 — bulk H3(及两个外部队列)对 NK 浸润比例做偏相关/残差化

- **对应顾虑**:nature-reviewer R1 主诉 + 交叉综合"共识技术风险(1)";关联稿件 §3.2/§3.3、Limitations 新增项。
- **涉及脚本**:
  - `src/topology/sst_axis_validation.py`
    (`hypothesis_test(..., "positive_corr", ["nk_protrusion_machinery_score","nk_synapse_cytotoxicity_outcome_score"])`
    → 现为 `scipy.stats.pearsonr` 零阶相关;产出 `sst_axis_positive_control_recovery.tsv`)。
  - `src/preprocessing/run_geo_external_validation.py`
    (外部 GSE62254/GSE84437 的同一 protrusion~cytotoxicity 相关;产出 `external_validation_results.tsv`)。
  - `src/immune_scoring/nk_scores.py`
    (`NK_infiltration_score = mean_zscore(expr, NK_MARKERS)` — 现成的、**无需 R** 的 NK 丰度代理;
    **但见下方⚠️,不可直接使用**)。

> ⚠️ **关键方法学陷阱(务必遵守)**:`NK_MARKERS`(15 基因)里有 **8 个与 cytotoxicity 模块重叠**
> (CCL5, GNLY, GZMB, IFNG, NKG7, PRF1, XCL1, XCL2)。直接用整套 `NK_infiltration_score` 做协变量
> = 把结局变量自己的基因扣掉 → **过度控制**,会人为塌掉 H3、得到虚假阴性。
> 必须改用**纯 lineage/受体代理**:`NKG_LINEAGE = [FCGR3A, KLRD1, KLRF1, KLRK1, NCAM1, NCR1, TYROBP]`
> ——这 7 个都在 `NK_MARKERS` 内、但**不在** cytotoxicity 模块**也不在** protrusion 模块中,
> 因而能代表 NK 丰度而不吞掉两模块的机制信号。路径 B 的去卷积 NK 分数天然满足此要求。
- **输入**:
  - `data/processed/bulk/tcga_lihc_expression.tsv`(LIHC 正对照)。
  - `data/processed/bulk/gse62254_expression.tsv`、`gse84437_expression.tsv`(基因级,已 probe→symbol 重映射)。
  - `NK_MARKERS`(`nk_scores.py` 内已定义,15 基因)。

- **动作**(路径 A 必做、无外部依赖;路径 B 可选、更强):
  1. **路径 A — NK lineage 丰度偏相关(主控制,始终可跑)**
     对每个 bulk 样本算 `nk_lineage_score = mean_zscore(expr, NKG_LINEAGE)`(7 基因,见上⚠️),
     再计算 protrusion~cytotoxicity 的**偏相关**(控制 `nk_lineage_score`);
     同一队列并报零阶 r 与偏相关 r、95% CI、p。三个 bulk 队列(LIHC/GSE62254/GSE84437)各做一遍。
     - 实现建议:在 `sst_axis_validation.py` 内新增 `partial_corr(x, y, covar)`
       (对 x、y 各自对 covar 做 OLS 残差,再对残差取 Pearson),或复用
       `src/topology/count_depth_control.py` 已有的 `residualize()` 逻辑。
     - **稳健性附检**:同时报告以 `NK_infiltration_score`(全 15 基因)为协变量的结果作为**下界**
       (预期会过度衰减),使读者能看到真实 partial r 落在 [lineage 代理, 全浸润代理] 区间内。
  2. **路径 B — 去卷积 NK 分数残差化(可选,若本地有 R)**
     用 quanTIseq 或 CIBERSORTx 估计每样本 NK 细胞分数,以其替代 `NK_infiltration_score`
     重复偏相关;把两种代理的结果并列,证明结论不依赖单一浸润代理的选择。
     - 注:Methods §2.3 已把 CIBERSORTx/quanTIseq 列为 bulk 反卷积回退;§3.4 亦附
       `03_supplementary/CIBERSORTx_quanTIseq_Scissor_roadmap.md`。R 环境不可得时,仅做路径 A 即达标。

- **产出**:
  - `results/tables/bulk_h3_purity_control.tsv` *新增*
    (列:cohort、n、r_zeroorder、p_zeroorder、
    r_partial_lineage、ci_low、ci_high、p_partial_lineage、  ← 主结果(7 基因 lineage 代理)
    r_partial_fullinfil、p_partial_fullinfil、               ← 下界附检(15 基因全浸润代理)
    [r_partial_deconv、p_partial_deconv]、                    ← 路径 B 可选
    covariate_used)。
  - `sst_axis_positive_control_recovery.tsv`(更新:H3 bulk 行新增偏相关列)。
  - `external_validation_results.tsv`(更新:GSE62254/GSE84437 行新增偏相关列)。
  - `submission_bundle_BiB/03_supplementary/tables/bulk_h3_purity_control.tsv`(同步)。

- **DoD**:
  - [ ] §3.2(H3 bulk)与 §3.3(外部复现)各补一句:控制 NK 浸润比例后 protrusion~cytotoxicity 的偏相关 r 与显著性。
  - [ ] 若偏相关后 r 仍显著且效应量可观 → 明确写"稳健于 NK 浸润比例调整(zero-order r=0.55 → partial r=X,p=Y)",效应臂头号证据升级。
  - [ ] Limitations 新增一条(或并入现 item 8/9 邻域):说明 bulk 相关已做浸润代理控制及其代理选择的局限。
  - [ ] 更新投稿终检:把本卡对应的 Risk 项("bulk H3 未展示纯度控制")标记为已闭合。
  - [ ] 三个 bulk 队列结论方向一致(或如实说明哪个队列衰减最大)。

- **退路**(若偏相关后显著衰减):
  如实报告"bulk 效应臂相关部分反映 NK 丰度协变,调整后偏相关 r=X";
  把效应臂措辞从"稳健恢复"下调为"与 NK 丰度部分耦合、调整后仍/不再显著的相关",
  并相应弱化摘要/Key Points/Conclusion 中效应臂的确定性表述。
  这一步是**零新增数据**的诚实降级,不阻断投稿,但须全文口径一致。

---

## 依赖与优先级

```
M5(bulk 纯度混杂) ── 独立,无前置;仅依赖已就位的 bulk 表达矩阵与现成 NK_infiltration_score
                    路径 A 无需 R,数小时可完成;路径 B 视 R 环境可选
```

- **关键性**:这是效应臂**头号证据**的直接稳健性检验;闭合它,稿件最核心的正向结论才算完全确立。
- **与 M1-S5 的关系**:M2/S1 处理的是**单细胞**效应臂的独立性与敏感性;M5 补上的是
  **bulk**层此前未覆盖的浸润混杂——两者合起来,效应臂在单细胞与 bulk 两个分辨率的混杂都被显式控制。

## 与 BiB_ACTION_PLAN.md / TASK_CARDS 的映射

| 本卡 | 对应 P 代码 | 关系 |
|------|------------|------|
| M5 | (新)§3.2/§3.3 bulk 效应臂 | 新增,NK 浸润比例偏相关;补 P0 系列在 bulk 层的空缺 |
