#!/usr/bin/env python3
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "work/high_impact/data"
META_DIR = ROOT / "outputs/研究升级/meta_analysis_v2"
OUT = ROOT / "outputs/研究升级/high_impact_v4/human_bone_localization"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(20260908)
N = 20000


def bh(values):
    p = np.asarray(values, float); order = np.argsort(p); q = np.empty(len(p))
    q[order] = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    return np.minimum(q, 1)


def prepare_meta():
    m = pd.read_csv(META_DIR / "gene_meta_analysis.csv").drop_duplicates("gene")
    m = m[(m.se_hk > 0) & np.isfinite(m.se_hk)].copy()
    m["meta_z"] = m.pooled_effect / m.se_hk
    m["se_bin"] = pd.qcut(m.se_hk.rank(method="first"), 10, labels=False)
    m["stratum"] = m.k.astype(int).astype(str) + "_" + m.se_bin.astype(str)
    return m


def matched_null(meta, genes, value):
    members = meta[meta.gene.isin(genes)]
    candidates = meta[~meta.gene.isin(genes)]
    out = np.zeros(N); expanded = 0
    for s, k in members.stratum.value_counts().items():
        pool = candidates.loc[candidates.stratum == s, value].to_numpy()
        if len(pool) < k:
            kk, sb = map(int, s.split("_"))
            for radius in range(1, 17):
                near = candidates[(candidates.k.sub(kk).abs() + candidates.se_bin.sub(sb).abs()) <= radius]
                pool = near[value].to_numpy()
                if len(pool) >= k:
                    expanded += 1; break
        out += np.vstack([RNG.choice(pool, k, replace=len(pool) < k) for _ in range(N)]).sum(axis=1)
    return out / len(members), len(members), expanded


def test_set(meta, genes, label, value, alternative):
    genes = set(genes) & set(meta.gene)
    obs = float(meta.loc[meta.gene.isin(genes), value].mean())
    null, n, expanded = matched_null(meta, genes, value)
    if alternative == "greater":
        p = (1 + np.sum(null >= obs)) / (N + 1)
    else:
        p = (1 + np.sum(np.abs(null-null.mean()) >= abs(obs-null.mean()))) / (N + 1)
    return {"cluster": label, "endpoint": value, "n_genes": n, "observed": obs,
            "null_mean": float(null.mean()), "effect": float(obs-null.mean()),
            "standardized_effect": float((obs-null.mean())/null.std(ddof=1)),
            "permutation_p": float(p), "expanded_strata": expanded}


def read_program(sheet):
    p = DATA / "NatureGenetics_2026_supp_tables_2_13.xlsx"
    d = pd.read_excel(p, sheet_name=sheet)
    d = d.rename(columns={"GENE": "gene", "CLUSTER NUMBER AND CELL TYPE": "cluster",
                          "AVERAGE LOG2(FC)": "avg_log2fc"})
    d["gene"] = d.gene.astype(str).str.upper()
    return d


meta = prepare_meta(); meta["gene"] = meta.gene.str.upper(); meta["abs_meta_z"] = meta.meta_z.abs()
sc, spatial = read_program("2c"), read_program("2d")

sc_osteo_labels = ["0. MSCS (CXCL12+GALNT17+)", "1. MSCS (CXCL12+STMN2+)",
                   "3. PRE-OSTEOBLASTS", "12. MATURE OSTEOBLASTS"]
sc["cluster_upper"] = sc.cluster.astype(str).str.upper()
sc_union = set(sc.loc[sc.cluster_upper.isin(sc_osteo_labels), "gene"])
primary = test_set(meta, sc_union, "scRNA_osteoblast_lineage_union", "abs_meta_z", "greater")

sc_rows = []
for cluster, d in sc.groupby("cluster"):
    sc_rows += [test_set(meta, d.gene, cluster, "meta_z", "two_sided"),
                test_set(meta, d.gene, cluster, "abs_meta_z", "greater")]
sc_tests = pd.DataFrame(sc_rows)
sc_tests["BH_q"] = sc_tests.groupby("endpoint")["permutation_p"].transform(bh)
sc_tests.to_csv(OUT / "human_scRNA_cluster_tests.csv", index=False)

spatial["cluster_upper"] = spatial.cluster.astype(str).str.upper()
sp_labels = ["0. MSCS", "9. OSTEO-MSCS", "1. OSTEOBLASTS"]
sp_union = set(spatial.loc[spatial.cluster_upper.isin(sp_labels), "gene"])
replication = test_set(meta, sp_union, "spatial_osteoblast_lineage_union", "abs_meta_z", "greater")

sp_rows = []
for cluster, d in spatial.groupby("cluster"):
    sp_rows += [test_set(meta, d.gene, cluster, "meta_z", "two_sided"),
                test_set(meta, d.gene, cluster, "abs_meta_z", "greater")]
sp_tests = pd.DataFrame(sp_rows)
sp_tests["BH_q"] = sp_tests.groupby("endpoint")["permutation_p"].transform(bh)
sp_tests.to_csv(OUT / "human_spatial_cluster_tests.csv", index=False)

corr_rows = []
for source, d in [("scRNA", sc), ("spatial", spatial)]:
    for cluster, g in d.groupby("cluster"):
        z = g.merge(meta[["gene", "meta_z", "abs_meta_z"]], on="gene")
        for endpoint in ["meta_z", "abs_meta_z"]:
            if len(z) >= 10:
                r = spearmanr(z.avg_log2fc, z[endpoint])
                corr_rows.append({"source": source, "cluster": cluster, "endpoint": endpoint,
                                  "n": len(z), "rho": r.statistic, "p": r.pvalue})
corr = pd.DataFrame(corr_rows)
corr["BH_q"] = corr.groupby(["source", "endpoint"])["p"].transform(bh)
corr.to_csv(OUT / "continuous_program_correlations.csv", index=False)

summary = {"primary_scRNA_osteoblast_lineage": primary,
           "prespecified_spatial_replication": replication,
           "scRNA_source_cells": 125063,
           "decision": "support" if primary["effect"] > 0 and primary["permutation_p"] < 0.05
                       else ("trend" if primary["effect"] > 0 else "not_support")}
(OUT / "human_bone_localization_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
