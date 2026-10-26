#!/usr/bin/env python3

import argparse
import importlib
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import MinMaxScaler


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Run the RTChou/purf positive-unlabelled random forest backend for "
            "P. falciparum antigen prioritisation. This uses the package "
            "PURandomForestClassifier with its PU impurity criterion."
        )
    )
    p.add_argument("--input", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--gene-col", default="gene_id")
    p.add_argument("--label-col", default="gene_id")
    p.add_argument("--positive-unit", choices=["gene", "site"], default="gene")
    p.add_argument("--position-col", default="auto", help="AA/residue position column in --input when --positive-unit site")
    p.add_argument("--training-sites", required=True)
    p.add_argument("--training-gene-col", default="gene_id")
    p.add_argument("--training-positive-col", default="train_positive")
    p.add_argument("--training-site-col", default="position")
    p.add_argument("--training-strength-col", default="iep_evidence_score")
    p.add_argument("--min-positive-sites", type=int, default=1)
    p.add_argument("--reference-gene-files", default="")
    p.add_argument("--positive-gene-col", default="")
    p.add_argument("--features", default="")
    p.add_argument(
        "--exclude-regex",
        default=(
            r"(^gene_|_key$|positive|known_iep|n_residues|purf_|reference|"
            r"tree_filter|training|^antigen_label$|^sp_prob$|"
            r"^gpi_specificity$|signal|gpi|t[_\. -]?cell|tcell|"
            r"t_epitope|ctl_epitope|cd4|cd8|mhc|hla|netmhc|netctl|"
            r"net_charge|c?nonsyn_syn_snp_ratio)"
        ),
    )
    p.add_argument("--cor-cutoff", type=float, default=0.75)
    p.add_argument(
        "--prefer-feature-regex",
        default=(
            r"TajimaD_PNG__mean|TajimaD_PNG__max|TajimaD_PNG__prop_nonzero|"
            r"TajimaD_mean__mean|TajimaD_mean__max|TajimaD_mean__prop_nonzero|"
            r"rsa__mean|B_epitope__mean|IUPRED2_score__max|discotope__mean|"
            r"secondarystruct_Group2|aa_|amino|hydropathy|charge|"
            r"molecular_weight|nonpolar|polar|aromatic"
        ),
    )
    p.add_argument("--positive-levels", default="0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9")
    p.add_argument("--select-positive-level-metric", default="reference_auc")
    p.add_argument("--n-trees-tuning", type=int, default=5000)
    p.add_argument("--n-trees-final", type=int, default=100000)
    p.add_argument("--max-features", default="sqrt")
    p.add_argument("--min-samples-leaf", type=int, default=1)
    p.add_argument("--random-state", type=int, default=7)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--top-n", type=int, default=250)
    p.add_argument("--tree-filter", default="TRUE")
    p.add_argument("--min-trees-after-filter", type=int, default=100)
    p.add_argument(
        "--tree-filter-min-reference-recall",
        type=float,
        default=1.0,
        help=(
            "Minimum fraction of OOB reference positives a tree must predict as "
            "positive to be retained. Chou's strict rule is 1.0; use 0.8 or "
            "0.7 to relax it."
        ),
    )
    p.add_argument(
        "--tree-filter-fallback",
        choices=["all", "best"],
        default="all",
        help=(
            "What to do if fewer than --min-trees-after-filter pass the recall "
            "threshold. 'all' preserves the old fallback; 'best' keeps the best "
            "trees ranked by OOB reference recall and diagnostics."
        ),
    )
    p.add_argument("--variable-space-weighting", default="TRUE")
    p.add_argument("--max-positive-replicates", type=int, default=8)
    p.add_argument("--positive-replicate-power", type=float, default=1.0)
    p.add_argument("--write-model", default="")
    return p.parse_args()


def split_arg(x):
    if x is None or str(x).strip() == "":
        return []
    return [y.strip() for y in str(x).split(",") if y.strip()]


def boolish_value(x):
    if isinstance(x, bool):
        return x
    s = str(x).strip().lower()
    if s in {"true", "t", "1", "yes", "y"}:
        return True
    if s in {"false", "f", "0", "no", "n", "", "nan", "none"}:
        return False
    try:
        return float(s) != 0
    except ValueError:
        return False


def boolish_series(x):
    if pd.api.types.is_bool_dtype(x):
        return x.fillna(False).to_numpy(dtype=bool)
    numeric = pd.to_numeric(x, errors="coerce")
    if numeric.notna().any():
        return numeric.fillna(0).to_numpy(dtype=float) != 0
    return x.astype(str).str.strip().str.lower().isin(["true", "t", "1", "yes", "y"]).to_numpy()


def read_table(path):
    path = str(path)
    suffix = Path(path).suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    return pd.read_csv(path, sep="\t")


def normalise_gene_id(value):
    text = str(value).strip().upper()
    hit = re.search(r"PF3D7_?[0-9]{7}", text)
    if hit:
        out = hit.group(0)
        return re.sub(r"^PF3D7([0-9])", r"PF3D7_\1", out)
    text = re.sub(r":MRNA-P[0-9]+$", "", text)
    text = re.sub(r"\.?[0-9]+-P[0-9]+$", "", text)
    text = re.sub(r"\.[0-9]+$", "", text)
    return text


def normalise_gene_key(value):
    return re.sub(r"[^A-Z0-9]+", "", normalise_gene_id(value).upper())


def pick_column(df, requested, candidates, label):
    if requested and requested != "auto":
        if requested not in df.columns:
            raise ValueError(f"{label} column not found: {requested}")
        return requested
    for col in candidates:
        if col in df.columns:
            return col
    raise ValueError(f"Could not find {label} column. Tried: {', '.join(candidates)}")


def normalise_position(value):
    val = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(val):
        return ""
    if float(val).is_integer():
        return str(int(val))
    return str(val)


def site_key(gene, position):
    g = normalise_gene_id(gene)
    p = normalise_position(position)
    if not g or not p:
        return ""
    return f"{g}:{p}"


def extract_ids_from_text(text):
    ids = re.findall(r"PF3D7_?[0-9]{7}", text.upper())
    return sorted(set(normalise_gene_id(x) for x in ids))


def extract_ids_from_file(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    ext = path.suffix.lower()
    if ext in {".xlsx", ".docx"}:
        import zipfile
        with zipfile.ZipFile(path) as zf:
            text = "\n".join(
                zf.read(name).decode("utf-8", errors="ignore")
                for name in zf.namelist()
                if name.endswith(".xml")
            )
    else:
        text = path.read_text(errors="ignore")
    return extract_ids_from_text(text)


def add_training_positives(df, args):
    sites = read_table(args.training_sites)
    for col in [args.training_gene_col, args.training_positive_col]:
        if col not in sites.columns:
            raise ValueError(f"Training-sites table lacks required column: {col}")

    sites["_is_positive"] = boolish_series(sites[args.training_positive_col])
    sites = sites.loc[sites["_is_positive"]].copy()
    if sites.empty:
        raise ValueError("No positive rows in training-sites table.")

    site_col = args.training_site_col if args.training_site_col in sites.columns else None
    if args.positive_unit == "site" and site_col is None:
        raise ValueError("--positive-unit site requires --training-site-col in training-sites table.")

    sites["gene_id_normalised"] = sites[args.training_gene_col].map(normalise_gene_id)
    if site_col:
        sites["_site"] = sites[site_col].astype(str)
    else:
        sites["_site"] = np.arange(len(sites)).astype(str)

    if args.positive_unit == "site":
        input_site_col = pick_column(
            df,
            args.position_col,
            ["position_std", "position", "AA_position", "aa_position", "AA_pos", "residue_position", "site_index", "site"],
            "input site/position",
        )
        df["source_gene_id_normalised"] = df[args.gene_col].map(normalise_gene_id)
        df["purf_unit_id"] = [site_key(g, p) for g, p in zip(df[args.gene_col], df[input_site_col])]
        df["gene_id_normalised"] = df["purf_unit_id"]
        df["gene_key"] = df["purf_unit_id"].map(lambda x: re.sub(r"[^A-Z0-9]+", "", str(x).upper()))

        sites["purf_unit_id"] = [site_key(g, p) for g, p in zip(sites[args.training_gene_col], sites[site_col])]
        sites = sites.loc[sites["purf_unit_id"].ne("")].copy()
        counts = sites[["purf_unit_id"]].drop_duplicates().rename(columns={"purf_unit_id": "gene_id_normalised"})
        counts["known_iep_site_count"] = 1

        if args.training_strength_col in sites.columns:
            strength = pd.to_numeric(sites[args.training_strength_col], errors="coerce").fillna(0)
            sites["_strength"] = np.maximum(strength.to_numpy(), 0)
            strength_tbl = (
                sites.groupby("purf_unit_id")["_strength"]
                .max()
                .reset_index(name="purf_training_strength")
                .rename(columns={"purf_unit_id": "gene_id_normalised"})
            )
        else:
            strength_tbl = pd.DataFrame({"gene_id_normalised": counts["gene_id_normalised"], "purf_training_strength": 1.0})

        df = df.merge(counts, on="gene_id_normalised", how="left")
        df = df.merge(strength_tbl, on="gene_id_normalised", how="left")
        df["known_iep_site_count"] = df["known_iep_site_count"].fillna(0).astype(int)
        df["purf_training_strength"] = pd.to_numeric(df["purf_training_strength"], errors="coerce").fillna(0)
        df["purf_reference_positive"] = df["known_iep_site_count"] >= args.min_positive_sites
        return df

    counts = (
        sites.groupby("gene_id_normalised")["_site"]
        .nunique()
        .reset_index(name="known_iep_site_count")
    )

    if args.training_strength_col in sites.columns:
        strength = pd.to_numeric(sites[args.training_strength_col], errors="coerce").fillna(0)
        sites["_strength"] = np.maximum(strength.to_numpy(), 0)
        strength_tbl = (
            sites.groupby("gene_id_normalised")["_strength"]
            .max()
            .reset_index(name="purf_training_strength")
        )
    else:
        strength_tbl = pd.DataFrame({"gene_id_normalised": counts["gene_id_normalised"], "purf_training_strength": 1.0})

    df["gene_id_normalised"] = df[args.gene_col].map(normalise_gene_id)
    df["gene_key"] = df[args.gene_col].map(normalise_gene_key)
    df = df.merge(counts, on="gene_id_normalised", how="left")
    df = df.merge(strength_tbl, on="gene_id_normalised", how="left")
    df["known_iep_site_count"] = df["known_iep_site_count"].fillna(0).astype(int)
    df["purf_training_strength"] = pd.to_numeric(df["purf_training_strength"], errors="coerce").fillna(0)
    df["purf_reference_positive"] = df["known_iep_site_count"] >= args.min_positive_sites
    return df


def write_gene_summary_from_site_scores(out, args):
    if args.positive_unit != "site" or "source_gene_id_normalised" not in out.columns:
        return
    gene = (
        out.assign(
            high_scoring_site=lambda x: (~x["purf_reference_positive"]) & (x["purf_rank"] <= args.top_n),
            score_ge_080=lambda x: (~x["purf_reference_positive"]) & (x["purf_score"] >= 0.80),
        )
        .groupby("source_gene_id_normalised", as_index=False)
        .agg(
            n_site_rows=("purf_score", "size"),
            n_known_positive_iep_sites=("purf_reference_positive", "sum"),
            max_site_purf_score=("purf_score", "max"),
            mean_site_purf_score=("purf_score", "mean"),
            median_site_purf_score=("purf_score", "median"),
            n_top_ranked_nontraining_sites=("high_scoring_site", "sum"),
            n_nontraining_sites_score_ge_080=("score_ge_080", "sum"),
            best_site_rank=("purf_rank", "min"),
        )
        .sort_values(["max_site_purf_score", "n_nontraining_sites_score_ge_080"], ascending=[False, False])
    )
    gene.to_csv(f"{args.out_prefix}.gene_summary_from_site_scores.tsv", sep="\t", index=False)


def read_reference_ids(args, positive_ids):
    files = split_arg(args.reference_gene_files)
    if not files:
        return set(positive_ids)
    out = set()
    for path in files:
        out.update(extract_ids_from_file(path))
    return out


def choose_features(df, args):
    requested = split_arg(args.features)
    if requested:
        return [x for x in requested if x in df.columns]
    banned = re.compile(args.exclude_regex, flags=re.IGNORECASE)
    out = []
    for col in df.columns:
        if col in {args.gene_col, args.label_col}:
            continue
        if banned.search(col):
            continue
        if pd.api.types.is_numeric_dtype(df[col]) or pd.api.types.is_bool_dtype(df[col]):
            out.append(col)
    return out


def prune_correlated_features(x_df, args):
    if args.cor_cutoff <= 0 or args.cor_cutoff >= 1:
        return x_df, pd.DataFrame(columns=["removed_feature", "retained_feature", "abs_spearman"])
    if x_df.shape[1] < 2:
        return x_df, pd.DataFrame(columns=["removed_feature", "retained_feature", "abs_spearman"])

    prefer = re.compile(args.prefer_feature_regex, flags=re.IGNORECASE)
    cor = x_df.corr(method="spearman").abs()
    cols = list(cor.columns)
    keep = set(cols)
    removed = []

    for i, a in enumerate(cols):
        if a not in keep:
            continue
        for b in cols[(i + 1):]:
            if b not in keep:
                continue
            val = cor.loc[a, b]
            if not np.isfinite(val) or val < args.cor_cutoff:
                continue
            a_pref = bool(prefer.search(a))
            b_pref = bool(prefer.search(b))
            if a_pref and not b_pref:
                drop, retain = b, a
            elif b_pref and not a_pref:
                drop, retain = a, b
            else:
                drop, retain = (b, a) if x_df[a].isna().sum() <= x_df[b].isna().sum() else (a, b)
            if drop in keep:
                keep.remove(drop)
                removed.append({"removed_feature": drop, "retained_feature": retain, "abs_spearman": val})
            if a not in keep:
                break

    kept_cols = [c for c in cols if c in keep]
    return x_df[kept_cols], pd.DataFrame(removed)


def prepare_matrix(df, features):
    x_df = df[features].apply(pd.to_numeric, errors="coerce")
    non_constant = []
    removed = []
    for col in x_df.columns:
        vals = x_df[col].dropna()
        if len(vals) == 0 or vals.nunique() <= 1:
            removed.append({"removed_feature": col, "retained_feature": "", "abs_spearman": np.nan, "reason": "constant_or_all_missing"})
        else:
            non_constant.append(col)
    x_df = x_df[non_constant]
    return x_df, pd.DataFrame(removed)


def scaled_imputed_matrix(x_df):
    imputer = SimpleImputer(strategy="median")
    x = imputer.fit_transform(x_df)
    scaler = MinMaxScaler()
    return scaler.fit_transform(x).astype(np.float32), imputer, scaler


def positive_replicates(x, y_pos, args):
    reps = np.ones(x.shape[0], dtype=int)
    if not boolish_value(args.variable_space_weighting):
        return reps
    pos_idx = np.where(y_pos)[0]
    if len(pos_idx) < 3:
        return reps
    xp = x[pos_idx]
    dist = np.sqrt(((xp[:, None, :] - xp[None, :, :]) ** 2).sum(axis=2))
    np.fill_diagonal(dist, np.nan)
    nearest = np.nanmin(dist, axis=1)
    if not np.any(np.isfinite(nearest)) or np.nanmax(nearest) <= 0:
        return reps
    scaled = nearest / np.nanmax(nearest)
    scaled = np.power(np.maximum(scaled, 0), args.positive_replicate_power)
    pos_reps = 1 + np.rint(scaled * (args.max_positive_replicates - 1)).astype(int)
    reps[pos_idx] = np.maximum(1, pos_reps)
    return reps


def duplicate_by_replicates(x, y, reps, row_ids):
    idx = np.repeat(np.arange(len(y)), reps)
    return x[idx], y[idx], row_ids[idx], idx


def auc_rank(y_true, score):
    y_true = np.asarray(y_true, dtype=bool)
    score = np.asarray(score, dtype=float)
    ok = np.isfinite(score)
    if y_true[ok].sum() == 0 or (~y_true[ok]).sum() == 0:
        return np.nan
    return roc_auc_score(y_true[ok].astype(int), score[ok])


def score_percentile(score):
    s = pd.Series(score)
    return s.rank(method="average", pct=True).to_numpy()


def score_separation(unlabelled_score):
    s = np.asarray(unlabelled_score, dtype=float)
    s = s[np.isfinite(s)]
    if len(s) < 10:
        return np.nan, np.nan
    q90 = np.quantile(s, 0.90)
    q50 = np.quantile(s, 0.50)
    q10 = np.quantile(s, 0.10)
    sep = (q90 - q50) / (q50 - q10 + 1e-9)
    pseudo_auc = q90 - q10
    return float(sep), float(pseudo_auc)


def import_rtchou_purf():
    try:
        module = importlib.import_module("purf.pu_ensemble")
        return getattr(module, "PURandomForestClassifier")
    except Exception as exc:
        raise ImportError(
            "Could not import RTChou/purf. Install it in the active conda env: "
            "git clone https://github.com/RTChou/purf.git && cd purf && python3 setup.py install"
        ) from exc


def train_purf(PURandomForestClassifier, x, y, pos_level, n_trees, args, seed):
    model = PURandomForestClassifier(
        pos_level=float(pos_level),
        n_estimators=int(n_trees),
        criterion="pu_gini",
        oob_score=True,
        n_jobs=int(args.threads),
        random_state=int(seed),
        max_features=args.max_features,
        min_samples_leaf=int(args.min_samples_leaf),
    )
    model.fit(x, y.astype(int))
    return model


def model_score(model, x):
    proba = model.predict_proba(x)
    if isinstance(proba, list):
        proba = proba[0]
    return np.asarray(proba)[:, 1]


def oob_score(model, x):
    score = np.asarray(model.oob_decision_function_)[:, 1]
    missing = ~np.isfinite(score)
    if missing.any():
        score[missing] = model_score(model, x[missing])
    return score


def tune_positive_levels(PURandomForestClassifier, x_train, y_train, original_row_ids, y_pos_original, y_ref_original, args):
    rows = []
    for i, level in enumerate([float(x) for x in split_arg(args.positive_levels)]):
        model = train_purf(
            PURandomForestClassifier,
            x_train,
            y_train,
            level,
            args.n_trees_tuning,
            args,
            args.random_state + i,
        )
        score_train = oob_score(model, x_train)
        tmp = pd.DataFrame({"row_id": original_row_ids, "score": score_train})
        score = tmp.groupby("row_id")["score"].mean().reindex(np.arange(len(y_pos_original))).to_numpy()
        pct = score_percentile(score)
        sep, pseudo_auc = score_separation(score[~y_pos_original])
        rows.append(
            {
                "positive_level": level,
                "positive_auc": auc_rank(y_pos_original, score),
                "reference_auc": auc_rank(y_ref_original, score),
                "unlabelled_score_separation": sep,
                "pseudo_auc_from_unlabelled_mixture": pseudo_auc,
                "explicit_positive_recall": float(np.nanmean(score[y_pos_original] >= 0.5)),
                "mean_positive_score": float(np.nanmean(score[y_pos_original])),
                "mean_unlabelled_score": float(np.nanmean(score[~y_pos_original])),
                "median_positive_percentile": float(np.nanmedian(pct[y_pos_original])),
                "mean_positive_percentile": float(np.nanmean(pct[y_pos_original])),
            }
        )
    tuning = pd.DataFrame(rows)
    metric = args.select_positive_level_metric
    if metric not in tuning.columns:
        raise ValueError(f"--select-positive-level-metric not found in tuning table: {metric}")
    tuning = tuning.sort_values(
        [metric, "explicit_positive_recall", "positive_auc", "median_positive_percentile"],
        ascending=[False, False, False, False],
    )
    return tuning, float(tuning.iloc[0]["positive_level"])


def tree_filter_weights(model, x_train, x_original, y_ref_original, original_row_ids, args):
    try:
        from sklearn.ensemble._forest import _generate_unsampled_indices, _get_n_samples_bootstrap
    except Exception as exc:
        raise ImportError("Could not import sklearn forest OOB helper functions for Chou tree filtering.") from exc

    n_train = len(original_row_ids)
    n_samples_bootstrap = _get_n_samples_bootstrap(n_train, model.max_samples)
    weights = np.zeros(len(model.estimators_), dtype=float)
    stats = []
    for i, tree in enumerate(model.estimators_):
        unsampled = _generate_unsampled_indices(tree.random_state, n_train, n_samples_bootstrap)
        original_unsampled = original_row_ids[unsampled]
        ref_mask = y_ref_original[original_unsampled]
        n_ref_oob = int(ref_mask.sum())
        keep = False
        n_ref_pred_pos = 0
        reference_recall = np.nan
        if n_ref_oob > 0:
            pred = tree.predict(x_train[unsampled][ref_mask])
            n_ref_pred_pos = int(np.sum(np.asarray(pred).astype(int) == 1))
            reference_recall = n_ref_pred_pos / n_ref_oob
            keep = reference_recall >= args.tree_filter_min_reference_recall
        all_score = tree.predict_proba(x_original)[:, 1]
        ref_auc = auc_rank(y_ref_original, all_score)
        mean_ref = float(np.nanmean(all_score[y_ref_original])) if np.any(y_ref_original) else np.nan
        mean_other = float(np.nanmean(all_score[~y_ref_original])) if np.any(~y_ref_original) else np.nan
        weights[i] = 1.0 if keep else 0.0
        stats.append(
            {
                "tree": i + 1,
                "n_reference_oob": n_ref_oob,
                "n_reference_oob_predicted_positive": n_ref_pred_pos,
                "reference_recall": reference_recall,
                "tree_filter_min_reference_recall": args.tree_filter_min_reference_recall,
                "ref_auc": ref_auc,
                "mean_ref": mean_ref,
                "mean_other": mean_other,
                "mean_diff": mean_ref - mean_other,
                "keep": bool(keep),
            }
        )
    stats = pd.DataFrame(stats)
    if int(weights.sum()) < args.min_trees_after_filter:
        if args.tree_filter_fallback == "best":
            n_keep = min(args.min_trees_after_filter, len(stats))
            best_idx = (
                stats.assign(
                    reference_recall_rank_value=stats["reference_recall"].fillna(-1),
                    ref_auc_rank_value=stats["ref_auc"].fillna(-1),
                    mean_diff_rank_value=stats["mean_diff"].fillna(-np.inf),
                )
                .sort_values(
                    ["reference_recall_rank_value", "ref_auc_rank_value", "mean_diff_rank_value"],
                    ascending=[False, False, False],
                )
                .head(n_keep)
                .index
            )
            weights[:] = 0.0
            weights[best_idx] = 1.0
            stats["keep"] = False
            stats.loc[best_idx, "keep"] = True
            stats["tree_filter_fallback_used"] = True
            sys.stderr.write(
                f"WARNING: Relaxed tree filtering kept fewer than {args.min_trees_after_filter} trees; "
                f"using the best {n_keep} trees instead of all trees.\n"
            )
        else:
            sys.stderr.write(
                f"WARNING: Chou tree filtering kept only {int(weights.sum())} trees; "
                "using all trees instead. Lower --tree-filter-min-reference-recall or use "
                "--tree-filter-fallback best if filtering is desired.\n"
            )
            weights[:] = 1.0
            stats["keep"] = True
            stats["tree_filter_fallback_used"] = True
    else:
        stats["tree_filter_fallback_used"] = False
    return weights, stats


def weighted_model_score(model, x, weights):
    if hasattr(model, "predict_proba_with_weights"):
        proba = model.predict_proba_with_weights(x, weights=weights)
        if isinstance(proba, list):
            proba = proba[0]
        return np.asarray(proba)[:, 1]
    return model_score(model, x)


def main():
    args = parse_args()
    df = read_table(args.input)
    if args.gene_col not in df.columns:
        if "protein_id" in df.columns:
            args.gene_col = "protein_id"
        else:
            raise ValueError(f"Gene column not found: {args.gene_col}")
    if args.label_col not in df.columns:
        args.label_col = args.gene_col

    df = add_training_positives(df, args)
    positive_ids = set(df.loc[df["purf_reference_positive"], "gene_id_normalised"].dropna())
    reference_ids = read_reference_ids(args, positive_ids)
    df["purf_tree_filter_reference"] = df["gene_id_normalised"].isin(reference_ids)

    features = choose_features(df, args)
    x_df, constant_removed = prepare_matrix(df, features)
    x_df, corr_removed = prune_correlated_features(x_df, args)
    if x_df.shape[1] == 0:
        raise ValueError("No usable feature columns remain after filtering/pruning.")

    x, _, _ = scaled_imputed_matrix(x_df)
    y_pos = df["purf_reference_positive"].to_numpy(dtype=bool)
    y_ref = df["purf_tree_filter_reference"].to_numpy(dtype=bool)
    if y_pos.sum() < 2:
        unit_label = "positive sites" if args.positive_unit == "site" else "positive genes"
        raise ValueError(f"Fewer than two {unit_label} matched the feature table.")

    reps = positive_replicates(x, y_pos, args)
    x_train, y_train, original_row_ids, _ = duplicate_by_replicates(
        x,
        y_pos.astype(int),
        reps,
        np.arange(len(y_pos)),
    )

    PURandomForestClassifier = import_rtchou_purf()
    tuning, best_level = tune_positive_levels(
        PURandomForestClassifier,
        x_train,
        y_train,
        original_row_ids,
        y_pos,
        y_ref,
        args,
    )
    final_model = train_purf(
        PURandomForestClassifier,
        x_train,
        y_train,
        best_level,
        args.n_trees_final,
        args,
        args.random_state + 10000,
    )

    raw_score = model_score(final_model, x)
    tree_stats = pd.DataFrame()
    final_score = raw_score.copy()
    trees_used = args.n_trees_final
    if boolish_value(args.tree_filter):
        weights, tree_stats = tree_filter_weights(final_model, x_train, x, y_ref, original_row_ids, args)
        final_score = weighted_model_score(final_model, x, weights)
        trees_used = int(np.sum(weights > 0))

    out = df.copy()
    out["purf_score_raw"] = raw_score
    out["purf_score"] = final_score
    out["purf_rank"] = pd.Series(-out["purf_score"]).rank(method="dense").astype(int)
    out["purf_positive_level"] = best_level
    out["purf_training_mode"] = "rtchou_binary"
    out["purf_tree_filtered"] = boolish_value(args.tree_filter)
    out["purf_trees_used"] = trees_used
    out["purf_positive_replicates"] = reps
    out["novel_purf_candidate"] = (~out["purf_reference_positive"]) & (out["purf_rank"] <= args.top_n)
    out = out.sort_values(["purf_rank", args.gene_col])

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(f"{args.out_prefix}.all_gene_scores.tsv", sep="\t", index=False)
    out.loc[out["novel_purf_candidate"]].head(args.top_n).to_csv(
        f"{args.out_prefix}.top_hits.tsv", sep="\t", index=False
    )
    write_gene_summary_from_site_scores(out, args)
    tuning.sort_values("positive_level").to_csv(f"{args.out_prefix}.positive_level_tuning.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "feature": list(x_df.columns),
            "importance": getattr(final_model, "feature_importances_", np.repeat(np.nan, x_df.shape[1])),
        }
    ).sort_values("importance", ascending=False, na_position="last").to_csv(
        f"{args.out_prefix}.feature_importance.tsv", sep="\t", index=False
    )
    removed = pd.concat([constant_removed, corr_removed], ignore_index=True, sort=False)
    removed.to_csv(f"{args.out_prefix}.removed_correlated_features.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_id_normalised": out["gene_id_normalised"], "positive_replicates": reps}).to_csv(
        f"{args.out_prefix}.positive_variable_space_weights.tsv", sep="\t", index=False
    )
    if not tree_stats.empty:
        tree_stats.to_csv(f"{args.out_prefix}.tree_filtering.tsv", sep="\t", index=False)

    if args.write_model:
        import pickle
        with open(args.write_model, "wb") as handle:
            pickle.dump(final_model, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"RTChou/purf backend: PURandomForestClassifier")
    print(f"Selected positive level: {best_level}")
    unit_label = "sites" if args.positive_unit == "site" else "genes"
    print(f"Positive {unit_label}: {int(y_pos.sum())}")
    print(f"Reference {unit_label} for tree filtering: {int(y_ref.sum())}")
    print(f"Feature columns: {x_df.shape[1]}")
    print(f"Final forest trees: {args.n_trees_final}")
    print(f"Trees used after Chou OOB reference filtering: {trees_used}")
    print(f"Wrote: {args.out_prefix}.all_gene_scores.tsv")
    print(f"Wrote: {args.out_prefix}.top_hits.tsv")
    if args.positive_unit == "site":
        print(f"Wrote: {args.out_prefix}.gene_summary_from_site_scores.tsv")
    print(f"Wrote: {args.out_prefix}.positive_level_tuning.tsv")


if __name__ == "__main__":
    main()
