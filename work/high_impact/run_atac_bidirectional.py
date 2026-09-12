#!/usr/bin/env python3
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent/'pydeps'))

import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "work/high_impact/data"
META_DIR = ROOT / "outputs/研究升级/meta_analysis_v2"
OUT = ROOT / "outputs/研究升级/high_impact_v4/atac_validation"
RNG = np.random.default_rng(20260908)
N_PERM = 20000
FILES = {
    "young_1": DATA / "GSM8747711_young_1.bw",
    "young_2": DATA / "GSM8747712_young_2.bw",
    "force_minus_1": DATA / "GSM8747713_Force_-_-1.bw",
    "force_minus_2": DATA / "GSM8747714_Force_-_-2.bw",
}


def extract(frame, name, path):
    bw = pyBigWig.open(str(path))
    chroms = bw.chroms()
    vals = []
    for r in frame.itertuples(index=False):
        a, b = max(0, int(r.tss) - 2000), min(chroms[r.chrom], int(r.tss) + 2001)
        v = bw.stats(r.chrom, a, b, type="mean", exact=True)[0]
        vals.append(0.0 if v is None or not np.isfinite(v) else float(v))
    header = bw.header()
    bw.close()
    frame[name] = vals
    frame[name + "_rank"] = pd.Series(np.log1p(vals)).rank(method="average", pct=True).to_numpy()
    frame[name + "_scaled"] = np.log1p(np.asarray(vals) / max(float(header["sumData"]), 1.0) * 1e9)
    return header


def matched_null(df, members, forbidden, value):
    m = df[df.gene.isin(members)]
    c = df[~df.gene.isin(forbidden)]
    counts = m.stratum.value_counts().to_dict()
    result = np.zeros(N_PERM)
    for s, k in counts.items():
        pool = c.loc[c.stratum == s, value].to_numpy()
        result += np.vstack([RNG.choice(pool, k, replace=len(pool) < k) for _ in range(N_PERM)]).sum(axis=1)
    return result / len(m), len(m)


old = pd.read_csv(OUT / "gene_promoter_atac_results.csv.gz")
cols = ["gene", "gene_id", "chrom", "strand", "start", "end", "tss", "gene_length",
        "pooled_effect", "atac_delta"]
x = old[cols].copy()
headers = {name: extract(x, name, path) for name, path in FILES.items()}
x["young_rank"] = x[["young_1_rank", "young_2_rank"]].mean(axis=1)
x["force_minus_rank"] = x[["force_minus_1_rank", "force_minus_2_rank"]].mean(axis=1)
x["soft_delta"] = x.force_minus_rank - x.young_rank
x["young_scaled"] = x[["young_1_scaled", "young_2_scaled"]].mean(axis=1)
x["force_minus_scaled"] = x[["force_minus_1_scaled", "force_minus_2_scaled"]].mean(axis=1)
x["soft_delta_scaled"] = x.force_minus_scaled - x.young_scaled
x["baseline_bin"] = pd.qcut(x.young_rank.rank(method="first"), 10, labels=False)
x["length_bin"] = pd.qcut(x.gene_length.rank(method="first"), 10, labels=False)
x["stratum"] = x.baseline_bin.astype(str) + "_" + x.length_bin.astype(str)

sig = json.loads((META_DIR / "meta_signature.json").read_text())
pos, neg = set(sig["positive"]) & set(x.gene), set(sig["negative"]) & set(x.gene)
forbidden = pos | neg
pos_obs = float(x.loc[x.gene.isin(pos), "soft_delta"].mean())
neg_obs = float(x.loc[x.gene.isin(neg), "soft_delta"].mean())
d = pos_obs - neg_obs
pos_null, npg = matched_null(x, pos, forbidden, "soft_delta")
neg_null, nng = matched_null(x, neg, forbidden, "soft_delta")
null = pos_null - neg_null
p_expected_negative = float((1 + np.sum(null <= d)) / (N_PERM + 1))

qc = {}
for a, b in [("young_1", "young_2"), ("force_minus_1", "force_minus_2")]:
    qc[a + "__" + b] = float(spearmanr(np.log1p(x[a]), np.log1p(x[b])).statistic)
rho_mech = spearmanr(x.atac_delta, x.soft_delta)
rho_meta = spearmanr(x.pooled_effect, x.soft_delta)
scaled = float(x.loc[x.gene.isin(pos), "soft_delta_scaled"].mean() -
               x.loc[x.gene.isin(neg), "soft_delta_scaled"].mean())

x.to_csv(OUT / "bidirectional_gene_promoter_results.csv.gz", index=False, compression="gzip")
summary = {
    "comparison": "young cells on soft 1.5 kPa substrate (Force-) versus young control",
    "n_genes": int(len(x)),
    "primary_like_directional_score": d,
    "expected_direction": "negative",
    "one_sided_matched_permutation_p": p_expected_negative,
    "n_positive": npg, "n_negative": nng,
    "matched_null_mean": float(null.mean()), "matched_null_sd": float(null.std(ddof=1)),
    "total_signal_scaled_directional_score": scaled,
    "within_group_qc": qc,
    "soft_vs_force_plus_genomewide_spearman": {"rho": float(rho_mech.statistic), "p": float(rho_mech.pvalue)},
    "soft_delta_vs_meta_effect_spearman": {"rho": float(rho_meta.statistic), "p": float(rho_meta.pvalue)},
    "bidirectional_direction_check": bool(d < 0),
    "bidirectional_statistical_support": bool(d < 0 and p_expected_negative < 0.05),
    "bigwig_headers": headers,
}
(OUT / "bidirectional_validation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
