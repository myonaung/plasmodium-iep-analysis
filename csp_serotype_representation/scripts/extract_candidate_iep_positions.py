#!/usr/bin/env python3
import argparse
import csv
import math
import sys
from pathlib import Path


def norm_text(value: str) -> str:
    value = (value or "").strip().lower()
    if value.startswith("pf"):
        value = value[2:]
    return value


def parse_float(value):
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(out):
        return None
    return out


def parse_bool_like(value):
    if value is None:
        return None
    text = str(value).strip().lower()
    if text in {"ignore", "any", ""}:
        return None
    if text in {"true", "t", "1", "yes", "y"}:
        return True
    if text in {"false", "f", "0", "no", "n"}:
        return False
    return None


def parse_position(row):
    position = row.get("position")
    if position not in (None, ""):
      try:
          return int(float(position))
      except ValueError:
          pass

    site = row.get("site", "")
    if site.startswith("pos"):
        try:
            return int(site[3:])
        except ValueError:
            return None
    return None


def parse_site_filter(value):
    if value in (None, ""):
        return None
    out = set()
    for piece in str(value).split(","):
        piece = piece.strip()
        if not piece:
            continue
        try:
            out.add(int(piece))
        except ValueError as exc:
            raise SystemExit(f"Invalid --site value: {piece}") from exc
    return out if out else None


def row_matches(row, antigen, gene_id):
    antigen_ok = True
    gene_ok = True

    if antigen:
        antigen_ok = norm_text(row.get("antigen", "")) == norm_text(antigen)
    if gene_id:
        gene_ok = (row.get("gene_id", "") or "").strip() == gene_id.strip()

    return antigen_ok and gene_ok


def main():
    parser = argparse.ArgumentParser(description="Extract immune-escape positions from annotated site summary CSV.")
    parser.add_argument("--input", required=True, help="Annotated merged CSV")
    parser.add_argument("--antigen", default="", help="Antigen label to match from the antigen column")
    parser.add_argument("--gene-id", default="", help="Gene ID to match from the gene_id column")
    parser.add_argument(
        "--site",
        default="",
        help=(
            "Optional manual amino-acid site/position panel; can be comma-separated, "
            "e.g. 352,357. When supplied, these sites bypass the posterior/SNP/DAF "
            "filters and are used as the requested serotype-defining positions."
        ),
    )
    parser.add_argument("--prob-sym-gt-asym-threshold", "--prob_sym_gt_asym_threshold", dest="prob_sym_gt_asym_threshold", type=float, default=0.80)
    parser.add_argument("--prob-sym-gt-expected-threshold", "--prob_sym_gt_expected_threshold", dest="prob_sym_gt_expected_threshold", type=float, default=0.2)
    parser.add_argument("--require-high-confidence", "--require_high_confidence", dest="require_high_confidence", default="ignore", help="Require high_confidence_candidate == TRUE")
    parser.add_argument("--require-pf7-ns-snp", "--require_pf7_ns_snp", dest="require_pf7_ns_snp", default="ignore", help="Require pf7_ns_snp == 1")
    parser.add_argument("--require-pv4-snps", "--require_pv4_snps", dest="require_pv4_snps", default="ignore", help="Require Pv4_snps == 1")
    parser.add_argument("--require-snp-col", "--require_snp_col", dest="require_snp_col", default="", help="Require the named SNP column to equal 1, for example Pv4_snps")
    parser.add_argument("--max-dominant-allele-freq", "--max_dominant_allele_freq", dest="max_dominant_allele_freq", type=float, default=None, help="Require dominant_allele_freq <= this value")
    parser.add_argument("--require-sym-gt-expected-over-asym", "--require_prob_sym_gt_expected_over_asym", dest="require_sym_gt_expected_over_asym", default="true", help="Require prob_sym_gt_expected > prob_asym_gt_expected")
    parser.add_argument("--out", default="", help="Optional TSV output path for the matched rows")
    args = parser.parse_args()

    if not args.antigen and not args.gene_id:
        parser.error("Provide at least one of --antigen or --gene-id")

    input_path = Path(args.input)
    if not input_path.exists():
        raise SystemExit(f"Input file not found: {input_path}")

    require_high_conf = parse_bool_like(args.require_high_confidence)
    require_pf7_ns = parse_bool_like(args.require_pf7_ns_snp)
    require_pv4_snps = parse_bool_like(args.require_pv4_snps)
    require_sym_gt_expected_over_asym = parse_bool_like(args.require_sym_gt_expected_over_asym)
    site_filter = parse_site_filter(args.site)

    rows = []
    with input_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if not row_matches(row, args.antigen, args.gene_id):
                continue

            psa = parse_float(row.get("prob_sym_gt_asym"))
            pse = parse_float(row.get("prob_sym_gt_expected"))
            pae = parse_float(row.get("prob_asym_gt_expected"))
            position = parse_position(row)
            high_conf = parse_bool_like(row.get("high_confidence_candidate"))
            pf7_ns_snp = parse_float(row.get("pf7_ns_snp"))
            pv4_snps = parse_float(row.get("Pv4_snps", row.get("pv4_snps")))
            dominant_allele_freq = parse_float(row.get("dominant_allele_freq"))

            if position is None:
                continue

            if site_filter is not None and position not in site_filter:
                continue

            if site_filter is not None:
                row["position"] = position
                row["prob_sym_gt_asym"] = psa if psa is not None else ""
                row["prob_sym_gt_expected"] = pse if pse is not None else ""
                row["prob_asym_gt_expected"] = pae if pae is not None else ""
                row["_sort_prob_sym_gt_asym"] = psa if psa is not None else -math.inf
                row["_sort_prob_sym_gt_expected"] = pse if pse is not None else -math.inf
                rows.append(row)
                continue

            if psa is None or pse is None:
                continue

            if psa < args.prob_sym_gt_asym_threshold:
                continue
            if pse < args.prob_sym_gt_expected_threshold:
                continue

            if require_high_conf is True and high_conf is not True:
                continue
            if require_high_conf is False and high_conf is not False:
                continue

            if require_pf7_ns is True and pf7_ns_snp != 1:
                continue
            if require_pf7_ns is False and pf7_ns_snp == 1:
                continue

            if require_pv4_snps is True and pv4_snps != 1:
                continue
            if require_pv4_snps is False and pv4_snps == 1:
                continue

            if args.require_snp_col:
                required_snp_value = parse_float(row.get(args.require_snp_col))
                if required_snp_value != 1:
                    continue

            if args.max_dominant_allele_freq is not None:
                if dominant_allele_freq is None or dominant_allele_freq > args.max_dominant_allele_freq:
                    continue

            if require_sym_gt_expected_over_asym is True:
                if pae is None or not (pse > pae):
                    continue
            if require_sym_gt_expected_over_asym is False:
                if pae is not None and pse > pae:
                    continue

            row["position"] = position
            row["prob_sym_gt_asym"] = psa
            row["prob_sym_gt_expected"] = pse
            row["prob_asym_gt_expected"] = pae
            row["_sort_prob_sym_gt_asym"] = psa
            row["_sort_prob_sym_gt_expected"] = pse
            rows.append(row)

    if site_filter is not None:
        seen_positions = {row["position"] for row in rows}
        for position in sorted(site_filter - seen_positions):
            rows.append(
                {
                    "gene_id": args.gene_id,
                    "antigen": args.antigen,
                    "position": position,
                    "site": f"pos{position}",
                    "_sort_prob_sym_gt_asym": -math.inf,
                    "_sort_prob_sym_gt_expected": -math.inf,
                }
            )

    rows.sort(
        key=lambda x: (
            x["position"],
            -x["_sort_prob_sym_gt_asym"],
            -x["_sort_prob_sym_gt_expected"],
        )
    )

    unique_rows = []
    seen = set()
    for row in rows:
        pos = row["position"]
        if pos in seen:
            continue
        seen.add(pos)
        unique_rows.append(row)

    if args.out:
        out_path = Path(args.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = [
            "gene_id",
            "antigen",
            "position",
            "site",
            "prob_sym_gt_asym",
            "prob_sym_gt_expected",
            "prob_asym_gt_expected",
            "ranking_effect",
            "high_confidence_candidate",
            "pf7_ns_snp",
            "Pv4_snps",
            "dominant_allele_freq",
            "expected_mismatch",
            "sym_mean",
            "asym_mean",
        ]
        with out_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for row in unique_rows:
                writer.writerow({k: row.get(k, "") for k in fieldnames})

    positions = ",".join(str(row["position"]) for row in unique_rows)
    if not positions:
        print("", end="")
        return

    sys.stdout.write(positions)


if __name__ == "__main__":
    main()
