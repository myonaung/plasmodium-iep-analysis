#!/usr/bin/env python3
"""Attach VCF genomic coordinates to RTChou PURF site-level predictions.

Input is the site-first RTChou output, where each row is a residue/SNP site.
The output is a predicted IEP-like SNP table first, followed by the gene that
each high-scoring SNP belongs to.
"""

import argparse
import csv
import re
import subprocess
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Merge RTChou site scores with Pf7 VCF ANN SNP coordinates.")
    p.add_argument("--site-scores", required=True, help="RTChou <prefix>.all_gene_scores.tsv from --positive-unit site")
    p.add_argument("--vcf", required=True, help="Annotated Pf7 SNP VCF with INFO/ANN")
    p.add_argument("--out", required=True, help="Detailed output TSV")
    p.add_argument("--simple-out", default="", help="Optional simple GeneID/Position/AA_position table")
    p.add_argument("--score-threshold", type=float, default=0.80)
    p.add_argument("--top-n", type=int, default=250, help="Maximum rows to keep after score filtering; 0 keeps all")
    p.add_argument("--candidate-genes", default="", help="Optional candidate gene table/list to restrict output")
    p.add_argument("--candidate-gene-col", default="auto")
    p.add_argument("--include-training-positives", default="FALSE")
    p.add_argument("--gene-col", default="source_gene_id_normalised")
    p.add_argument("--position-col", default="position_std")
    p.add_argument("--score-col", default="purf_score")
    p.add_argument("--rank-col", default="purf_rank")
    p.add_argument("--known-col", default="purf_reference_positive")
    p.add_argument("--bcftools", default="bcftools")
    return p.parse_args()


def bool_arg(x):
    return str(x).strip().lower() in {"true", "t", "1", "yes", "y"}


def boolish_series(x):
    if pd.api.types.is_bool_dtype(x):
        return x.fillna(False)
    numeric = pd.to_numeric(x, errors="coerce")
    out = pd.Series(False, index=x.index)
    out.loc[numeric.notna()] = numeric.loc[numeric.notna()] != 0
    text = x.astype(str).str.strip().str.lower()
    out.loc[text.isin(["true", "t", "1", "yes", "y"])] = True
    return out


def normalise_gene_id(value):
    if pd.isna(value):
        return ""
    text = str(value).strip().upper()
    hit = re.search(r"PF3D7_?[0-9]{7}", text)
    if hit:
        return re.sub(r"^PF3D7([0-9])", r"PF3D7_\1", hit.group(0))
    text = re.sub(r":MRNA-P[0-9]+$", "", text)
    text = re.sub(r"\.?[0-9]+-P[0-9]+$", "", text)
    text = re.sub(r"\.[0-9]+$", "", text)
    return text


def pick_existing(df, requested, candidates, label):
    if requested and requested != "auto" and requested in df.columns:
        return requested
    for col in candidates:
        if col in df.columns:
            return col
    raise ValueError(f"Could not find {label} column. Tried: {', '.join([requested] + candidates)}")


def read_table_auto(path):
    path = Path(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, low_memory=False)
    return pd.read_csv(path, sep="\t", low_memory=False)


def load_candidate_genes(path, gene_col):
    table = read_table_auto(path)
    if gene_col == "auto":
        gene_col = pick_existing(
            table,
            "auto",
            [
                "Candidate gene PlasmoDB ID",
                "gene_id",
                "gene_ID",
                "GeneID",
                "Gene_ID",
                "gene_id_normalised",
                "Gene ID",
                "PlasmoDB ID",
                "PlasmoDB_ID",
                "plasmodb_id",
            ],
            "candidate gene",
        )
    elif gene_col not in table.columns:
        raise ValueError(f"Candidate gene column not found: {gene_col}")
    return {normalise_gene_id(x) for x in table[gene_col] if normalise_gene_id(x)}


def aa_position_from_ann(fields):
    if len(fields) > 13:
        hit = re.search(r"(\d+)", fields[13])
        if hit:
            return int(hit.group(1))
    if len(fields) > 10:
        hit = re.search(r"p\.[A-Za-z*]+(\d+)", fields[10])
        if hit:
            return int(hit.group(1))
    return None


def parse_ann(ann_text):
    for ann in str(ann_text).split(","):
        fields = ann.split("|")
        if len(fields) < 5:
            continue
        gene_id = normalise_gene_id(fields[4])
        aa_pos = aa_position_from_ann(fields)
        if not gene_id or aa_pos is None:
            continue
        yield {
            "GeneID": gene_id,
            "AA_position": aa_pos,
            "Effect": fields[1] if len(fields) > 1 else "",
            "Impact": fields[2] if len(fields) > 2 else "",
            "GeneName": fields[3] if len(fields) > 3 else "",
            "Transcript": fields[6] if len(fields) > 6 else "",
            "HGVS_c": fields[9] if len(fields) > 9 else "",
            "HGVS_p": fields[10] if len(fields) > 10 else "",
        }


def read_vcf_map(args, wanted_genes):
    view_cmd = [args.bcftools, "view", "-v", "snps", "-Ou", args.vcf]
    query_cmd = [args.bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n"]
    view = subprocess.Popen(view_cmd, stdout=subprocess.PIPE)
    query = subprocess.Popen(query_cmd, stdin=view.stdout, stdout=subprocess.PIPE, text=True)
    view.stdout.close()

    rows = []
    seen = set()
    for line in query.stdout:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5:
            continue
        chrom, pos, ref, alt, ann_text = parts[:5]
        if ann_text in {"", "."}:
            continue
        for rec in parse_ann(ann_text):
            if rec["GeneID"] not in wanted_genes:
                continue
            key = (rec["GeneID"], rec["AA_position"], chrom, pos, ref, alt, rec["Transcript"], rec["Effect"])
            if key in seen:
                continue
            seen.add(key)
            rec.update({"Chrom": chrom, "Position": int(pos), "Ref": ref, "Alt": alt})
            rows.append(rec)

    status = query.wait()
    view_status = view.wait()
    if status or view_status:
        raise RuntimeError("bcftools failed while reading VCF annotations.")
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    scores = pd.read_csv(args.site_scores, sep="\t", low_memory=False)
    gene_col = pick_existing(
        scores,
        args.gene_col,
        ["source_gene_id_normalised", "gene_id", "GeneID", "protein_id"],
        "site-score gene",
    )
    pos_col = pick_existing(
        scores,
        args.position_col,
        ["position_std", "position", "AA_position", "aa_position", "residue_position", "site"],
        "site-score amino-acid position",
    )
    for col in [args.score_col, args.rank_col]:
        if col not in scores.columns:
            raise ValueError(f"Required site-score column not found: {col}")

    scores = scores.copy()
    scores["GeneID"] = scores[gene_col].map(normalise_gene_id)
    scores["AA_position"] = pd.to_numeric(scores[pos_col], errors="coerce")
    scores[args.score_col] = pd.to_numeric(scores[args.score_col], errors="coerce")
    scores[args.rank_col] = pd.to_numeric(scores[args.rank_col], errors="coerce")
    scores = scores.loc[scores["GeneID"].ne("") & scores["AA_position"].notna()].copy()
    scores["AA_position"] = scores["AA_position"].astype(int)

    if args.known_col in scores.columns and not bool_arg(args.include_training_positives):
        scores = scores.loc[~boolish_series(scores[args.known_col])].copy()

    if args.candidate_genes:
        candidate_genes = load_candidate_genes(args.candidate_genes, args.candidate_gene_col)
        scores = scores.loc[scores["GeneID"].isin(candidate_genes)].copy()
        print(f"Candidate genes loaded: {len(candidate_genes)}")
        print(f"Site-score rows after candidate-gene filter: {len(scores)}")

    candidates = scores.loc[scores[args.score_col] >= args.score_threshold].copy()
    if candidates.empty:
        candidates = scores.sort_values(args.score_col, ascending=False).copy()
        if args.top_n and args.top_n > 0:
            candidates = candidates.head(args.top_n)
    else:
        candidates = candidates.sort_values([args.score_col, args.rank_col], ascending=[False, True])
        if args.top_n and args.top_n > 0:
            candidates = candidates.head(args.top_n)

    vcf_map = read_vcf_map(args, set(candidates["GeneID"]))
    merged = candidates.merge(vcf_map, on=["GeneID", "AA_position"], how="left")
    first_cols = [
        "GeneID", "Position", "AA_position", "Chrom", "Ref", "Alt",
        args.score_col, args.rank_col, "Effect", "Impact", "GeneName", "Transcript", "HGVS_c", "HGVS_p",
    ]
    cols = [c for c in first_cols if c in merged.columns] + [c for c in merged.columns if c not in first_cols]
    merged = merged[cols]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out, sep="\t", index=False)
    print(f"Wrote: {out}")
    print(f"Predicted IEP-like SNP rows: {len(merged)}")

    if args.simple_out:
        simple = merged[["GeneID", "Position", "AA_position"]].drop_duplicates()
        simple = simple.sort_values(["GeneID", "AA_position", "Position"], na_position="last")
        simple.to_csv(args.simple_out, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
        print(f"Wrote: {args.simple_out}")


if __name__ == "__main__":
    main()
