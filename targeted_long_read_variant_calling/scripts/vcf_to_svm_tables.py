#!/usr/bin/env python3
"""Convert Longshot and PEPPER VCFs to sm-SNIPER SVM input tables.

This is a lightweight replacement for the original GATK VariantsToTable +
longshot_tabular.sh/pepper_tabular.sh conversion path. It uses only the Python
standard library and writes the two text files consumed by unique_variants.R.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


MISSING = "NA"


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def parse_info(info_text: str) -> dict[str, str]:
    if info_text in ("", "."):
        return {}

    info: dict[str, str] = {}
    for item in info_text.rstrip(";").split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = "TRUE"
    return info


def parse_format(format_text: str, sample_text: str) -> dict[str, str]:
    if format_text in ("", ".") or sample_text in ("", "."):
        return {}

    keys = format_text.split(":")
    values = sample_text.split(":")
    return dict(zip(keys, values))


def first_sample(record: list[str], sample_name: str | None) -> str:
    if len(record) <= 9:
        return "."
    if sample_name is None:
        return record[9]
    return record[9]


def genotype_to_bases(gt: str, ref: str, alt: str) -> str:
    if gt in ("", ".", "./.", ".|."):
        return gt or "."

    separator = "|" if "|" in gt else "/"
    alleles = [ref] + alt.split(",")
    converted: list[str] = []
    for token in gt.replace("|", "/").split("/"):
        if token == ".":
            converted.append(".")
            continue
        try:
            converted.append(alleles[int(token)])
        except (ValueError, IndexError):
            converted.append(token)
    return separator.join(converted)


def qual_value(raw: str) -> str:
    if raw in ("", "."):
        return MISSING
    try:
        return f"{float(raw):.1f}"
    except ValueError:
        return raw


def value(mapping: dict[str, str], key: str) -> str:
    val = mapping.get(key, MISSING)
    return MISSING if val in ("", ".") else val


def iter_vcf_records(vcf: Path):
    with open_text(vcf) as handle:
        sample_names: list[str] = []
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.lstrip("#").split("\t")
                sample_names = header[9:]
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                continue
            yield fields, sample_names


def write_longshot_table(vcf: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "AC", "AQ", "SC", "PH", "GT", "DP", "GQ", "UQ"]

    with output.open("w") as out:
        out.write("\t".join(columns) + "\n")
        for record, sample_names in iter_vcf_records(vcf):
            chrom, pos, raw_id, ref, alt, qual, _flt, raw_info = record[:8]
            fmt = record[8] if len(record) > 8 else "."
            sample = first_sample(record, sample_names[0] if sample_names else None)
            info = parse_info(raw_info)
            fmt_values = parse_format(fmt, sample)
            gt = genotype_to_bases(value(fmt_values, "GT"), ref, alt)

            row = [
                chrom,
                pos,
                "longshot" if raw_id in ("", ".") else raw_id,
                ref,
                alt,
                qual_value(qual),
                value(info, "AC"),
                value(info, "AQ"),
                value(info, "SC"),
                value(info, "PH"),
                gt,
                value(info, "DP"),
                value(fmt_values, "GQ"),
                value(fmt_values, "UQ"),
            ]
            out.write("\t".join(row) + "\n")


def write_pepper_table(vcf: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "GT", "AP", "GQ", "DP", "AD", "VAF"]

    with output.open("w") as out:
        out.write("\t".join(columns) + "\n")
        for record, sample_names in iter_vcf_records(vcf):
            chrom, pos, raw_id, ref, alt, qual = record[:6]
            fmt = record[8] if len(record) > 8 else "."
            sample = first_sample(record, sample_names[0] if sample_names else None)
            fmt_values = parse_format(fmt, sample)
            gt = genotype_to_bases(value(fmt_values, "GT"), ref, alt)

            row = [
                chrom,
                pos,
                "PEPPER" if raw_id in ("", ".") else raw_id,
                ref,
                alt,
                qual_value(qual),
                gt,
                value(fmt_values, "AP"),
                value(fmt_values, "GQ"),
                value(fmt_values, "DP"),
                value(fmt_values, "AD"),
                value(fmt_values, "VAF"),
            ]
            out.write("\t".join(row) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create sm-SNIPER SVM input text tables from Longshot and PEPPER VCFs without GATK."
    )
    parser.add_argument("--sample", required=True, help="Sample name used in output filenames.")
    parser.add_argument("--longshot-vcf", required=True, type=Path, help="Longshot amplicon VCF.")
    parser.add_argument("--pepper-vcf", required=True, type=Path, help="PEPPER amplicon VCF.")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory for *_longshot.txt and *_PEPPER.txt.")
    args = parser.parse_args()

    write_longshot_table(args.longshot_vcf, args.outdir / f"{args.sample}_longshot.txt")
    write_pepper_table(args.pepper_vcf, args.outdir / f"{args.sample}_PEPPER.txt")


if __name__ == "__main__":
    main()
