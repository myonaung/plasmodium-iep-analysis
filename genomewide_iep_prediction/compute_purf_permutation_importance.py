#!/usr/bin/env python3
"""Permutation importance for SNP/site-level PURF features.

This is a practical, defensible analogue of Chou-style permutation importance
for the SNP-first model. It trains a supervised surrogate random forest using:
  positive class: known candidate IEP SNP/site rows
  background class: low-scoring unlabelled SNP/site rows

Feature importance is reported as the decrease in a held-out performance metric
after permuting each feature. With --scoring accuracy, this gives a practical
"mean decrease in accuracy" analogue. It is not exact Chou per-tree OOB mean
decrease in accuracy from the final tree-filtered PURF model.
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.metrics import accuracy_score, balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def parse_args():
    p = argparse.ArgumentParser(description="Compute permutation decrease in known-positive AUC.")
    p.add_argument("--scores", required=True, help="Site-level PURF <prefix>.all_gene_scores.tsv")
    p.add_argument("--feature-importance", required=True, help="PURF <prefix>.feature_importance.tsv")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--known-col", default="purf_reference_positive")
    p.add_argument("--score-col", default="purf_score")
    p.add_argument("--top-n-features", type=int, default=50)
    p.add_argument("--negative-score-max", type=float, default=0.20)
    p.add_argument("--negative-ratio", type=int, default=20)
    p.add_argument("--max-negatives", type=int, default=50000)
    p.add_argument("--n-estimators", type=int, default=500)
    p.add_argument("--n-repeats", type=int, default=10)
    p.add_argument("--test-size", type=float, default=0.30)
    p.add_argument(
        "--scoring",
        choices=["roc_auc", "accuracy", "balanced_accuracy"],
        default="roc_auc",
        help="Held-out metric to decrease by permutation.",
    )
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--random-state", type=int, default=7)
    return p.parse_args()


def boolish(series):
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).to_numpy(dtype=bool)
    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.notna().any():
        return numeric.fillna(0).to_numpy(dtype=float) != 0
    return series.astype(str).str.strip().str.lower().isin(["true", "t", "1", "yes", "y"]).to_numpy()


def feature_category(feature):
    f = str(feature).lower()
    if "tajima" in f or "snp" in f or "nonsyn" in f or "syn" in f or "nucleotide" in f:
        return "Genomic/population-genetic"
    if "epitope" in f or "b_cell" in f or "immunogenicity" in f or "tap" in f or "hla" in f or "mhc" in f:
        return "Immunological"
    if "rsa" in f or "solvent" in f or "surface" in f or "discotope" in f or "iupred" in f or "secondary" in f or "ss_is" in f:
        return "Structural"
    if "charge" in f or "polar" in f or "nonpolar" in f or "aromatic" in f or "hydropathy" in f or "amino" in f:
        return "Physicochemical"
    if "tm" in f or "transmembrane" in f or "signal" in f or "gpi" in f or "secret" in f:
        return "Proteomic/membrane-secretory"
    if str(feature).startswith("genelevel_"):
        return "Gene-level context"
    return "Other"


def read_top_features(feature_importance_path, top_n):
    fi = pd.read_csv(feature_importance_path, sep="\t")
    if not {"feature", "importance"}.issubset(fi.columns):
        raise ValueError("Feature-importance file must contain feature and importance columns.")
    fi["importance"] = pd.to_numeric(fi["importance"], errors="coerce").fillna(0)
    fi = fi.sort_values("importance", ascending=False)
    if top_n and top_n > 0:
        fi = fi.head(top_n)
    return fi["feature"].astype(str).tolist()


def load_score_subset(args, features):
    header = pd.read_csv(args.scores, sep="\t", nrows=0).columns.tolist()
    needed = [args.known_col, args.score_col] + [f for f in features if f in header]
    missing = [f for f in features if f not in header]
    if missing:
        print("WARNING: missing feature columns skipped: " + ", ".join(missing[:20]))
    if args.known_col not in header:
        raise ValueError(f"Known column not found: {args.known_col}")
    if args.score_col not in header:
        raise ValueError(f"Score column not found: {args.score_col}")
    if len(needed) <= 2:
        raise ValueError("No requested feature columns found in score table.")
    return pd.read_csv(args.scores, sep="\t", usecols=needed, low_memory=False), needed[2:]


def select_training_rows(df, args):
    known = boolish(df[args.known_col])
    score = pd.to_numeric(df[args.score_col], errors="coerce")

    pos_idx = np.flatnonzero(known)
    neg_idx = np.flatnonzero((~known) & (score <= args.negative_score_max))
    if len(pos_idx) < 2:
        raise ValueError("Need at least two known-positive rows.")
    if len(neg_idx) < 2:
        order = np.argsort(score.fillna(np.inf).to_numpy())
        neg_idx = np.array([i for i in order if not known[i]][: max(2, len(pos_idx) * args.negative_ratio)])

    rng = np.random.default_rng(args.random_state)
    n_neg = min(len(neg_idx), args.max_negatives, max(len(pos_idx) * args.negative_ratio, len(pos_idx)))
    neg_idx = rng.choice(neg_idx, size=n_neg, replace=False)

    idx = np.concatenate([pos_idx, neg_idx])
    y = np.concatenate([np.ones(len(pos_idx), dtype=int), np.zeros(len(neg_idx), dtype=int)])
    return idx, y


def metric_label(scoring):
    if scoring == "roc_auc":
        return "auc"
    if scoring == "balanced_accuracy":
        return "balanced_accuracy"
    return "accuracy"


def model_metric(model, x, y, scoring):
    if scoring == "roc_auc":
        return roc_auc_score(y, model.predict_proba(x)[:, 1])
    pred = model.predict(x)
    if scoring == "balanced_accuracy":
        return balanced_accuracy_score(y, pred)
    return accuracy_score(y, pred)


def group_permutation_importance(model, x_test, y_test, features, args):
    rng = np.random.default_rng(args.random_state + 1000)
    baseline = model_metric(model, x_test, y_test, args.scoring)
    label = metric_label(args.scoring)
    groups = {}
    for i, feature in enumerate(features):
        groups.setdefault(feature_category(feature), []).append(i)

    rows = []
    for group, cols in groups.items():
        drops = []
        for _ in range(args.n_repeats):
            xp = x_test.copy()
            for col in cols:
                xp[:, col] = rng.permutation(xp[:, col])
            score = model_metric(model, xp, y_test, args.scoring)
            drops.append(baseline - score)
        drops = np.asarray(drops, dtype=float)
        se = drops.std(ddof=1) / np.sqrt(len(drops)) if len(drops) > 1 else np.nan
        rows.append(
            {
                "category": group,
                "n_features": len(cols),
                f"mean_decrease_{label}": drops.mean(),
                f"se_decrease_{label}": se,
                f"scaled_decrease_{label}": drops.mean() / se if se and np.isfinite(se) and se > 0 else np.nan,
                f"baseline_{label}": baseline,
            }
        )
    return pd.DataFrame(rows).sort_values(f"mean_decrease_{label}", ascending=False)


def main():
    args = parse_args()
    features = read_top_features(args.feature_importance, args.top_n_features)
    df, features = load_score_subset(args, features)
    idx, y = select_training_rows(df, args)

    x_df = df.iloc[idx][features].apply(pd.to_numeric, errors="coerce")
    x = x_df.to_numpy(dtype=float)

    x_train, x_test, y_train, y_test = train_test_split(
        x,
        y,
        test_size=args.test_size,
        random_state=args.random_state,
        stratify=y,
    )

    model = make_pipeline(
        SimpleImputer(strategy="median"),
        StandardScaler(),
        RandomForestClassifier(
            n_estimators=args.n_estimators,
            max_features="sqrt",
            class_weight="balanced",
            n_jobs=args.threads,
            random_state=args.random_state,
        ),
    )
    model.fit(x_train, y_train)
    label = metric_label(args.scoring)
    baseline_score = model_metric(model, x_test, y_test, args.scoring)

    perm = permutation_importance(
        model,
        x_test,
        y_test,
        scoring=args.scoring,
        n_repeats=args.n_repeats,
        random_state=args.random_state + 1,
        n_jobs=args.threads,
    )

    se = perm.importances_std / np.sqrt(args.n_repeats)
    out = pd.DataFrame(
        {
            "feature": features,
            "category": [feature_category(f) for f in features],
            f"mean_decrease_{label}": perm.importances_mean,
            f"sd_decrease_{label}": perm.importances_std,
            f"se_decrease_{label}": se,
            f"scaled_decrease_{label}": np.where(se > 0, perm.importances_mean / se, np.nan),
            f"baseline_{label}": baseline_score,
            "n_known_positive": int(y.sum()),
            "n_pseudo_background": int((y == 0).sum()),
        }
    ).sort_values(f"mean_decrease_{label}", ascending=False)

    group = group_permutation_importance(model, x_test, y_test, features, args)
    key = pd.DataFrame(
        [
            {"metric": f"baseline_known_positive_{label}", "value": baseline_score},
            {"metric": "known_positive_rows", "value": int(y.sum())},
            {"metric": "pseudo_background_rows", "value": int((y == 0).sum())},
            {"metric": "features_tested", "value": len(features)},
            {"metric": "negative_score_max", "value": args.negative_score_max},
            {"metric": "n_estimators", "value": args.n_estimators},
            {"metric": "n_repeats", "value": args.n_repeats},
            {"metric": "scoring", "value": args.scoring},
        ]
    )

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    out_file = str(out_prefix) + f".permutation_decrease_known_positive_{label}.tsv"
    group_file = str(out_prefix) + f".group_permutation_decrease_known_positive_{label}.tsv"
    key_file = str(out_prefix) + f".permutation_{label}_key_numbers.tsv"
    out.to_csv(out_file, sep="\t", index=False)
    group.to_csv(group_file, sep="\t", index=False)
    key.to_csv(key_file, sep="\t", index=False)

    print(f"Baseline held-out known-positive {label}: {baseline_score:.4f}")
    print(f"Known-positive rows: {int(y.sum())}")
    print(f"Pseudo-background rows: {int((y == 0).sum())}")
    print(f"Wrote: {out_file}")
    print(f"Wrote: {group_file}")
    print(f"Wrote: {key_file}")


if __name__ == "__main__":
    main()
