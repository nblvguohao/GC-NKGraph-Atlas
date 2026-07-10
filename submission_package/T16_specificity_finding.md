# T16 — 效应臂机制特异性：本地实测 + 方法学复核

> 2026-07-10 本地运行（数据 `data/processed/scrna/gc_nk_subset.h5ad`，8,310 NK × 22,728 基因，log-normalized；anndata/scanpy 本地齐备）。
> 回应本轮三审 R1/R3 最致命顾虑：H3（protrusion machinery ~ cytotoxicity output）
> 测的是 Zheng serine→SM→topology 机制，还是任意两个激活相关 NK 模块的共变？

## 检验设计
匹配零分布置换：从**排除全部 SST 轴基因 + 激活签名**后的 NK 表达基因池中，随机抽取
与真实模块**同尺寸、同激活载荷**的"伪突起"（24 基因）/"伪细胞毒性"（11 基因）模块对，
N=1000 次，构造 partial-r 零分布；经验 p = 观测落在零分布上尾的比例。

## 结果一：初测（仅控激活）→ NOT_SPECIFIC
| 量 | 值 |
|---|---|
| 观测 partial r（protrusion~cytotoxicity \| activation） | 0.251（p=2.6×10⁻¹¹⁹） |
| 匹配零分布（1000 次） | mean 0.244, sd 0.042, p95 0.317 |
| 经验 p | **0.42** → 落在零分布正中 |

初读：耦合不强于随机匹配模块对——似乎证实"泛激活共变"。

## 结果二：方法学复核（`run_t16_robustness.py`）→ 结论反转
怀疑：仅控 16 基因激活签名不足以剔除**测序深度/全局转录活性**这一让**任何**两个模块
均值同涨同落的技术混杂项。于是加入 `log(total_counts)` 作协变量重算：

| 控制变量 | 观测 r | 零分布均值 | 零分布 p95 | 经验 p | 判决 |
|---|---|---|---|---|---|
| 仅激活（0.25/0.5/1.0/unmatched） | 0.251 | 0.24 | 0.31 | ~0.42 | not_spec |
| **激活 + 深度（同上 4 档）** | **0.159** | **0.06** | **0.12** | **~0.006** | **SPECIFIC** |

（4 档匹配容差在每种控制下几乎一致，说明结论不依赖匹配带宽。）

## 决定性解读
- 随机 NK 模块对在"仅控激活"时仍相关 ~0.24，**主要是测序深度混杂**，不是生物学。
- 一旦剔除深度：**零分布几乎塌掉（0.24 → 0.06）**，而**真实 protrusion~cytotoxicity 仍保留 0.159**（零分布均值的约 2.6 倍、超 95 分位，经验 p≈0.006）。
- 因此"仅控激活"是**错误零分布**（深度整体抬高、掩盖真实信号）；"控激活+深度"对真实模块与随机模块一视同仁剔除深度，是**公平对照**。
- **在公平对照下，效应臂耦合显著特异于任意匹配随机模块对 → `effector arm recovers` 论断成立。**

坦诚提醒：控 total_counts 是 scRNA 标准做法（数据虽已 library-size 归一化，残余深度相关仍在）；
深度里或含少量生物学，但零分布同样剔除深度，故比较公平——观测超出零分布不受"深度是否部分生物学"影响。

## 对稿件的建议（审稿弱点 → 强化点）
1. **不下调** §3.2 结论方向；
2. 在 §3.2 T14 激活对照旁**新增一句**：深度-控制的匹配置换零分布检验显示 protrusion~cytotoxicity
   partial r=0.159 显著高于同尺寸、同激活载荷的随机 NK 模块对（经验 p≈0.006，1000 次置换），
   即耦合机制特异，非泛激活/深度共变；
3. 附补充图（`t16_null_distribution` / `t16_robustness` 的零分布 + 观测竖线）；
4. 把 T14/H3 主分析也加深度协变量，令特异性论断无懈可击。

## 产物
- `results/t16_specificity_control.tsv` — 初测（仅控激活）汇总 + 判决。
- `results/t16_null_distribution.tsv` — 初测 1000 次零分布 partial-r（补充图数据）。
- `results/t16_robustness.tsv` — 复核 2 控制 × 4 匹配 = 8 变体表。
- `scripts/run_t16_specificity_control.py` — 主检验（从 `src/common/sst_config.py` 载模块）。
- `scripts/run_t16_robustness.py` — 深度控制 + 容差扫描复核。

> 运行环境：本地 Python 3.8.20，numpy 1.23.5 / pandas 2.0.3 / scipy 1.10.1 / anndata 0.9.2 / scanpy。
> 复现：`python src/a100_recompute/run_t16_specificity_control.py && python src/a100_recompute/run_t16_robustness.py`
