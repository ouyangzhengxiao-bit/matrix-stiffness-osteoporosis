#!/usr/bin/env python3
"""Prespecified chronological leave-future-out transportability analysis."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "outputs/研究升级/meta_analysis_v2/study_level_gene_effects.csv.gz"
OUT = ROOT / "outputs/研究升级/high_impact_v5/temporal_transportability"
OUT.mkdir(parents=True, exist_ok=True)

DATES = {
    "GSE55867": "2014-03-13",
    "GSE181512": "2021-08-05",
    "GSE193021": "2022-01-04",
    "GSE226411": "2023-03-01",
    "GSE255574": "2024-02-12",
    "GSE288678": "2025-02-04",
    "GSE310514": "2025-11-19",
}
ORDER = sorted(DATES, key=lambda s: (DATES[s], s))


def pool(train: pd.DataFrame) -> pd.DataFrame:
    z = train.assign(weight=1 / train.variance)
    g = z.groupby("gene", sort=False)
    out = g.apply(lambda x: pd.Series({
        "training_effect": np.average(x.effect, weights=x.weight),
        "training_k": x.study.nunique(),
        "training_weight": x.weight.sum(),
    }), include_groups=False)
    return out[out.training_k >= 3]


def evaluate(training, heldout, mode):
    tr = pool(df[df.study.isin(training)])
    ho = df[df.study == heldout].set_index("gene")[["effect", "variance"]]
    x = tr.join(ho, how="inner").dropna()
    rho = float(spearmanr(x.training_effect, x.effect).statistic)
    top = x.loc[x.training_effect.abs().nlargest(min(200, len(x))).index].copy()
    agreement = float((np.sign(top.training_effect) == np.sign(top.effect)).mean())
    weighted_mse = float(np.average((x.effect - x.training_effect) ** 2,
                                    weights=1 / x.variance))
    return {
        "mode": mode, "held_out_study": heldout,
        "held_out_submission_date": DATES[heldout],
        "training_studies": ";".join(training), "n_training_studies": len(training),
        "n_shared_genes": len(x), "spearman_rho": rho,
        "top200_directional_agreement": agreement,
        "inverse_variance_weighted_mse": weighted_mse,
    }


df = pd.read_csv(SRC)
df = df[df.study.isin(ORDER)].copy()
rows = []
for i in range(3, len(ORDER)):
    rows.append(evaluate(ORDER[:i], ORDER[i], "sequential_expanding"))
fixed = ORDER[:3]
for held in ORDER[3:]:
    rows.append(evaluate(fixed, held, "fixed_earliest_three"))

res = pd.DataFrame(rows)
res.to_csv(OUT / "chronological_transportability_results.csv", index=False)
pd.DataFrame({"study": ORDER, "geo_submission_date": [DATES[x] for x in ORDER],
              "chronological_rank": range(1, len(ORDER)+1)}).to_csv(
                  OUT / "study_chronology.csv", index=False)

summary = {}
for mode, z in res.groupby("mode"):
    summary[mode] = {
        "tests": len(z),
        "median_spearman_rho": float(z.spearman_rho.median()),
        "rho_range": [float(z.spearman_rho.min()), float(z.spearman_rho.max())],
        "positive_rho_tests": int((z.spearman_rho > 0).sum()),
        "median_top200_directional_agreement": float(z.top200_directional_agreement.median()),
        "agreement_range": [float(z.top200_directional_agreement.min()),
                            float(z.top200_directional_agreement.max())],
    }
(OUT / "temporal_transportability_summary.json").write_text(
    json.dumps(summary, indent=2), encoding="utf-8")
print(res.to_string(index=False))
print(json.dumps(summary, indent=2))
