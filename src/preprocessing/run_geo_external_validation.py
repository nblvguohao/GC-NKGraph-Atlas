"""
外部验证：修正 GEO 探针->基因映射，并检验 SST 轴在独立胃癌队列中的复现。

背景
----
A100 上 `data/raw/bulk/` 为空（原始 series matrix 不在），但
`data/processed/bulk/{gse}_expression.tsv` 是**探针级**矩阵（samples × probes，
列名形如 "1007_s_at" / "ILMN_1343048"，且被 preprocessing 做过 .upper() 和加引号）。
因此本脚本直接在已处理探针矩阵上做映射，不依赖原始 series matrix，
且**不触碰 TCGA-STAD/LIHC 的既有 NK 标签**（保持已完成的 GNN/基线对比可比）。

流程
----
1. 用 GEOparse 下载平台注释，构建 probe(大写)->gene_symbol 映射。
2. 清洗列名（去三重引号/空白/大写），映射到基因，按"均值最高探针"折叠到基因级。
3. 备份原探针矩阵为 *_probe_level.bak.tsv，写回基因级 {gse}_expression.tsv。
4. 计算 SST 轴模块分数（复用 mechanism card 的模块基因），检验：
   - protrusion_machinery ~ cytotoxicity_outcome（期望正，H3 的外部复现）
   - sm_balance ~ protrusion_machinery（期望正，H2 的外部检验）
   并记录 NK marker 命中数。结果写 results/tables/external_validation_results.tsv。

用法（在 A100 项目根，激活 gc-nkgraph 环境）
    python src/preprocessing/run_geo_external_validation.py
"""
import os, sys, shutil
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

# SST 轴模块基因（硬编码自 mechanism card / sst_config，避免依赖 A100 上尚未同步的模块）
MODULES = {
    "nk_protrusion_machinery": ['EZR', 'MSN', 'RDX', 'ACTR2', 'ACTR3', 'ARPC1B',
        'ARPC2', 'ARPC3', 'ARPC4', 'ARPC5', 'WAS', 'WASL', 'WASF1', 'WASF2',
        'WASF3', 'WIPF1', 'CDC42', 'RAC1', 'RHOA', 'DIAPH1', 'DIAPH3', 'FMNL1',
        'BAIAP2', 'PACSIN2'],
    "nk_synapse_cytotoxicity_outcome": ['NKG7', 'GNLY', 'GZMB', 'PRF1', 'IFNG',
        'LCP2', 'LAT', 'VAV1', 'TLN1', 'ITGAL', 'ITGB2'],
    "nk_sm_synthesis": ['SGMS1', 'SGMS2'],
    "nk_sm_catabolism": ['SMPD1', 'SMPD2', 'SMPD3', 'SMPD4'],
}

PROC = "data/processed/bulk"
RAW = "data/raw/bulk"
OUT = "results/tables/external_validation_results.tsv"

# 平台：GSE62254=Affy GPL570；GSE84437=Illumina GPL6947（失败则回退试其它）
PLATFORM = {"gse62254": ["GPL570"], "gse84437": ["GPL6947", "GPL10558", "GPL6884"]}
SYMBOL_COLS = ["Gene Symbol", "Gene symbol", "Symbol", "ILMN_Gene", "GENE_SYMBOL",
               "gene_symbol", "GENE", "geneSymbol", "Gene_Symbol"]
NK_MARKERS = ["NKG7", "GNLY", "GZMB", "PRF1", "KLRD1", "NCAM1", "KLRK1"]


def clean(c):
    return str(c).strip().strip('"').strip().strip('"').strip().upper()


def build_map(gse):
    import GEOparse
    last = None
    for plat in PLATFORM[gse]:
        try:
            gpl = GEOparse.get_GEO(geo=plat, destdir=RAW, silent=True)
            t = gpl.table
            idc = "ID" if "ID" in t.columns else t.columns[0]
            symc = next((c for c in SYMBOL_COLS if c in t.columns), None)
            if symc is None:
                last = f"{plat}: 无符号列 {list(t.columns)[:8]}"; continue
            m = {}
            for pid, sym in zip(t[idc], t[symc]):
                if not isinstance(sym, str):
                    continue
                s = sym.split("///")[0].split("//")[0].strip()
                if s and s not in {"---", "NA"}:
                    m[str(pid).strip().upper()] = s
            if len(m) > 1000:
                print(f"  [{gse}] 平台 {plat}: {len(m)} 探针映射（符号列 '{symc}'）")
                return m
            last = f"{plat}: 映射过少 ({len(m)})"
        except Exception as e:
            last = f"{plat}: {e}"
    raise RuntimeError(f"{gse}: 构建探针->基因映射失败 -> {last}")


def collapse(gse, probe2sym):
    df = pd.read_csv(f"{PROC}/{gse}_expression.tsv", sep="\t", index_col=0)  # samples×probes
    orig_cols = df.shape[1]
    sym = {col: probe2sym.get(clean(col)) for col in df.columns}
    keep = [c for c in df.columns if sym[c] is not None]
    sub = df[keep].copy()
    sub.columns = [sym[c] for c in keep]
    # 每个基因保留"均值最高"的探针列：先按列均值降序排，再去重保第一
    order = sub.mean(axis=0).sort_values(ascending=False).index
    sub = sub.loc[:, order]
    sub = sub.loc[:, ~pd.Index(sub.columns).duplicated(keep="first")]
    sub = sub.sort_index(axis=1)
    print(f"  [{gse}] {orig_cols} 探针 -> {len(keep)} 有映射 -> {sub.shape[1]} 基因")
    return sub  # samples×genes


def module_score(df_genes, genes):
    present = [g for g in genes if g in df_genes.columns]
    if not present:
        return None, 0
    x = df_genes[present]
    z = (x - x.mean(axis=0)) / x.std(axis=0).replace(0, np.nan)
    return z.mean(axis=1), len(present)


def main():
    os.makedirs("results/tables", exist_ok=True)
    prot = MODULES["nk_protrusion_machinery"]
    cyto = MODULES["nk_synapse_cytotoxicity_outcome"]
    smsyn = MODULES["nk_sm_synthesis"]
    smcat = MODULES["nk_sm_catabolism"]

    rows = []
    for gse in ["gse62254", "gse84437"]:
        print(f"\n=== {gse} ===")
        probe2sym = build_map(gse)
        genes_df = collapse(gse, probe2sym)

        # 备份 + 写回基因级
        bak = f"{PROC}/{gse}_expression.probe_level.bak.tsv"
        if not os.path.exists(bak):
            shutil.copy(f"{PROC}/{gse}_expression.tsv", bak)
        out = genes_df.copy(); out.index.name = "sample_id"
        out.to_csv(f"{PROC}/{gse}_expression.tsv", sep="\t")

        nk_hit = [g for g in NK_MARKERS if g in genes_df.columns]
        ps, npr = module_score(genes_df, prot)
        cs, ncy = module_score(genes_df, cyto)
        syn, _ = module_score(genes_df, smsyn)
        cat, _ = module_score(genes_df, smcat)

        def rec(name, a, b, exp):
            if a is None or b is None:
                rows.append(dict(dataset=gse, test=name, r="", p="", n_genes="",
                                 expected=exp, outcome="MIXED_UNRESOLVED"))
                return
            d = pd.DataFrame({"a": a, "b": b}).dropna()
            r, p = stats.pearsonr(d.a, d.b)
            ok = (r > 0 and exp == "+" and p < 0.05)
            rows.append(dict(dataset=gse, test=name, r=round(r, 4), p=f"{p:.2e}",
                             n_genes="", expected=exp,
                             outcome="RECOVERED" if ok else ("weak/ns" if r > 0 else "wrong_sign")))

        rec("protrusion~cytotoxicity", ps, cs, "+")
        if syn is not None and cat is not None:
            smbal = syn - cat
            rec("sm_balance~protrusion", smbal, ps, "+")
        rows.append(dict(dataset=gse, test="NK_marker_coverage",
                         r=f"{len(nk_hit)}/{len(NK_MARKERS)}", p="", n_genes=genes_df.shape[1],
                         expected="", outcome=str(nk_hit)))
        print(f"  NK markers: {len(nk_hit)}/{len(NK_MARKERS)} {nk_hit}; "
              f"protrusion({npr} genes)~cytotox({ncy} genes)")

    res = pd.DataFrame(rows)
    res.to_csv(OUT, sep="\t", index=False)
    print(f"\n写出: {OUT}")
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
