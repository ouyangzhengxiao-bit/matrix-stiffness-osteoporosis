from pathlib import Path
import gzip
import hashlib
import io
import itertools
import json
import re

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import false_discovery_control

from run_prediction_v1 import prepare_data


DATA = Path("work/upgrade/data")
OUT = Path("outputs/研究升级/cross_experiment_v1")
OUT.mkdir(parents=True, exist_ok=True)


def read_gse22011():
    with gzip.open(DATA / "GSE22011_series_matrix.txt.gz", "rt") as handle:
        text = handle.read()
    table = text.split("!series_matrix_table_begin\n")[1].split("!series_matrix_table_end")[0]
    expression = pd.read_csv(io.StringIO(table), sep="\t", index_col=0)
    expression.index = expression.index.astype(str)
    with gzip.open(DATA / "GPL6244.annot.gz", "rt") as handle:
        header_line = None
        for number, line in enumerate(handle):
            if line.startswith("ID\t"):
                header_line = number
                break
    annotation = pd.read_csv(
        DATA / "GPL6244.annot.gz",
        sep="\t",
        skiprows=header_line,
        dtype=str,
        low_memory=False,
    )[["ID", "Gene symbol"]]
    annotation = annotation.dropna()
    annotation = annotation[
        annotation["Gene symbol"].str.match(r"^[A-Za-z][A-Za-z0-9_.-]*$")
    ].drop_duplicates("ID")
    annotation = annotation.set_index("ID")["Gene symbol"]
    expression["gene"] = annotation.reindex(expression.index)
    expression = expression.dropna().groupby("gene").median()
    if expression.shape[1] != 15:
        raise ValueError("Expected 15 GSE22011 samples")
    return expression


def discover_signature(expression):
    stiffness = np.log2(np.array([100, 400, 1600, 6400, 25600], dtype=float))
    stiffness -= stiffness.mean()
    denominator = np.square(stiffness).sum()
    slopes = {}
    for donor in range(3):
        values = expression.iloc[:, donor * 5 : (donor + 1) * 5].to_numpy()
        slopes[f"donor_{donor + 1}"] = values @ stiffness / denominator
    slopes = pd.DataFrame(slopes, index=expression.index)
    slopes["mean_slope"] = slopes.mean(axis=1)
    slopes["consistent_positive"] = (slopes.iloc[:, :3] > 0).all(axis=1)
    slopes["consistent_negative"] = (slopes.iloc[:, :3] < 0).all(axis=1)
    positive = slopes[slopes.consistent_positive].nlargest(100, "mean_slope").index.tolist()
    negative = slopes[slopes.consistent_negative].nsmallest(100, "mean_slope").index.tolist()
    if len(positive) != 100 or len(negative) != 100 or set(positive) & set(negative):
        raise ValueError("Unable to form the prespecified directional signature")
    slopes.to_csv(OUT / "discovery_gene_slopes.csv")
    signature = {"positive": positive, "negative": negative}
    (OUT / "discovery_signature.json").write_text(json.dumps(signature, indent=2))
    return signature


def gencode_transcript_map():
    mapping = {}
    pattern_transcript = re.compile(r'transcript_id "([^"]+)"')
    pattern_gene = re.compile(r'gene_name "([^"]+)"')
    with gzip.open(DATA / "gencode.v38.annotation.gtf.gz", "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 8)
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            transcript = pattern_transcript.search(fields[8])
            gene = pattern_gene.search(fields[8])
            if transcript and gene:
                mapping[transcript.group(1)] = gene.group(1)
    return pd.Series(mapping, name="gene")


def read_gse255574():
    transcript_map = gencode_transcript_map()
    matrices = []
    metadata = []
    folder = DATA / "GSE255574_DMSO"
    for path in sorted(folder.glob("*.txt.gz")):
        match = re.search(r"_(150|500|t20)-(24|48)-d_(\d)-3", path.name)
        if not match:
            raise ValueError(f"Cannot parse {path.name}")
        stiffness = {"150": 150.0, "500": 500.0, "t20": 2000.0}[match.group(1)]
        time = int(match.group(2))
        replicate = int(match.group(3))
        sample = path.name.split("_")[0]
        table = pd.read_csv(path, sep="\t", usecols=["Name", "NumReads"])
        table["gene"] = transcript_map.reindex(table.Name).to_numpy()
        counts = table.dropna(subset=["gene"]).groupby("gene").NumReads.sum()
        counts.name = sample
        matrices.append(counts)
        metadata.append(
            {"sample": sample, "stiffness_pa": stiffness, "time_hours": time, "replicate": replicate}
        )
    counts = pd.concat(matrices, axis=1).fillna(0)
    metadata = pd.DataFrame(metadata).set_index("sample").loc[counts.columns]
    if counts.shape[1] != 18:
        raise ValueError("Expected 18 GSE255574 samples")
    metadata.to_csv(OUT / "GSE255574_sample_metadata.csv")
    return counts, metadata


def score_samples(matrix, signature):
    ranks = matrix.rank(axis=0, pct=True)
    positive = ranks.index.intersection(signature["positive"])
    negative = ranks.index.intersection(signature["negative"])
    score = ranks.loc[positive].mean() - ranks.loc[negative].mean()
    return score, len(positive), len(negative)


def exact_three_group_slope(values, stiffness):
    centered = stiffness - stiffness.mean()
    observed = np.dot(values - values.mean(), centered) / np.square(centered).sum()
    slopes = []
    positions = set(range(9))
    unique_levels = np.sort(np.unique(stiffness))
    for first in itertools.combinations(range(9), 3):
        remaining = sorted(positions - set(first))
        for second in itertools.combinations(remaining, 3):
            labels = np.full(9, unique_levels[2])
            labels[list(first)] = unique_levels[0]
            labels[list(second)] = unique_levels[1]
            labels -= labels.mean()
            slopes.append(np.dot(values - values.mean(), labels) / np.square(labels).sum())
    slopes = np.asarray(slopes)
    p_value = np.mean(np.abs(slopes) >= abs(observed) - 1e-15)
    return observed, p_value, len(slopes)


def test_msc(signature):
    counts, metadata = read_gse255574()
    libraries = counts.sum(axis=0)
    logcpm = np.log2(counts.div(libraries, axis=1) * 1_000_000 + 1)
    score, n_positive, n_negative = score_samples(logcpm, signature)
    rows = []
    sample_rows = []
    for time in [24, 48]:
        subset = metadata.index[metadata.time_hours == time]
        stiffness = np.log2(metadata.loc[subset, "stiffness_pa"].to_numpy())
        values = score.loc[subset].to_numpy()
        slope, p_value, permutations = exact_three_group_slope(values, stiffness)
        rows.append(
            {
                "time_hours": time,
                "n_samples": len(subset),
                "positive_mapped": n_positive,
                "negative_mapped": n_negative,
                "score_slope_per_log2_pa": slope,
                "exact_permutation_p": p_value,
                "permutations": permutations,
            }
        )
        for sample in subset:
            sample_rows.append(
                {
                    "sample": sample,
                    "time_hours": time,
                    "stiffness_pa": metadata.loc[sample, "stiffness_pa"],
                    "score": score.loc[sample],
                }
            )
    results = pd.DataFrame(rows)
    results["BH_q"] = false_discovery_control(results.exact_permutation_p)
    results.to_csv(OUT / "MSC_transport_results.csv", index=False)
    pd.DataFrame(sample_rows).to_csv(OUT / "MSC_scores.csv", index=False)
    return results


def test_mouse(signature):
    pheno, logcpm, _ = prepare_data()
    homology = pd.read_csv("outputs/研究升级/one_to_one_homology.csv")
    human_to_mouse = homology.set_index("human").mouse
    mapped = {
        direction: human_to_mouse.reindex(genes).dropna().tolist()
        for direction, genes in signature.items()
    }
    mouse_signature = {
        direction: [gene for gene in genes if gene in logcpm.index]
        for direction, genes in mapped.items()
    }
    score, n_positive, n_negative = score_samples(logcpm, mouse_signature)
    pheno["signature_score"] = score.reindex(pheno.index)
    rows = []
    for outcome in ["maximum_load", "structural_stiffness"]:
        columns = ["signature_score", outcome, "sex_binary", "age", "weight", "generation"]
        data = pheno[columns].dropna()
        x = pd.DataFrame(index=data.index)
        for column in ["signature_score", "age", "weight"]:
            x[column] = (data[column] - data[column].mean()) / data[column].std(ddof=1)
        x["sex_binary"] = data.sex_binary
        x = pd.concat(
            [x, pd.get_dummies(data.generation, prefix="generation", drop_first=True, dtype=float)],
            axis=1,
        )
        x = sm.add_constant(x, has_constant="add")
        y = (data[outcome] - data[outcome].mean()) / data[outcome].std(ddof=1)
        fit = sm.OLS(y, x).fit(cov_type="HC3", use_t=True)
        ci = fit.conf_int().loc["signature_score"]
        rows.append(
            {
                "outcome": outcome,
                "n": len(data),
                "positive_mapped": n_positive,
                "negative_mapped": n_negative,
                "standardized_beta": fit.params.signature_score,
                "ci_low": ci.iloc[0],
                "ci_high": ci.iloc[1],
                "HC3_p": fit.pvalues.signature_score,
            }
        )
    results = pd.DataFrame(rows)
    results["BH_q"] = false_discovery_control(results.HC3_p)
    results.to_csv(OUT / "mouse_bone_results.csv", index=False)
    pheno[["generation", "signature_score", "maximum_load", "structural_stiffness"]].to_csv(
        OUT / "mouse_scores.csv"
    )
    return results


def main():
    discovery = read_gse22011()
    signature = discover_signature(discovery)
    msc = test_msc(signature)
    mouse = test_mouse(signature)
    manifest = {}
    paths = [
        DATA / "GSE22011_series_matrix.txt.gz",
        DATA / "GPL6244.annot.gz",
        DATA / "gencode.v38.annotation.gtf.gz",
    ] + sorted((DATA / "GSE255574_DMSO").glob("*.txt.gz"))
    for path in paths:
        manifest[str(path)] = {
            "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
    (OUT / "input_manifest.json").write_text(json.dumps(manifest, indent=2))
    summary = {
        "discovery_genes": discovery.shape[0],
        "signature_positive": len(signature["positive"]),
        "signature_negative": len(signature["negative"]),
        "MSC_transport_criterion_met": bool((msc.BH_q < 0.05).all() and (msc.score_slope_per_log2_pa > 0).all()),
        "mouse_any_adjusted_association": bool((mouse.BH_q < 0.05).any()),
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2))
    print("MSC transport")
    print(msc.to_string(index=False))
    print("Mouse bone")
    print(mouse.to_string(index=False))
    print(summary)


if __name__ == "__main__":
    main()
