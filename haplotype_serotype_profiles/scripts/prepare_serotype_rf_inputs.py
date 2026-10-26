#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Create RF-ready raw/alluvial flow inputs and an AA-position map "
            "restricted to antigens where the position map exactly matches the "
            "number of IEP tokens in immune_escape_haplotype."
        )
    )
    p.add_argument("--raw-flows", required=True, help="Raw flows TSV from the alluvial run")
    p.add_argument("--alluvial-flows", help="Displayed/alluvial flows TSV from the same run")
    p.add_argument("--position-map", required=True, help="Long TSV with species,gene_id,antigen,aa_position")
    p.add_argument("--out-prefix", required=True, help="Output prefix for cleaned files")
    return p.parse_args()


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        return reader.fieldnames, rows


def write_tsv(path, fieldnames, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def key(row):
    return (row["species"], row["gene_id"], row["antigen"])


def main():
    args = parse_args()

    map_fields, map_rows = read_tsv(args.position_map)
    raw_fields, raw_rows = read_tsv(args.raw_flows)

    required_map = {"species", "gene_id", "antigen", "aa_position"}
    required_flow = {"species", "gene_id", "antigen", "immune_escape_haplotype"}
    if not required_map.issubset(map_fields):
        missing = sorted(required_map - set(map_fields))
        raise SystemExit(f"Position map missing columns: {', '.join(missing)}")
    if not required_flow.issubset(raw_fields):
        missing = sorted(required_flow - set(raw_fields))
        raise SystemExit(f"Raw flows missing columns: {', '.join(missing)}")

    map_positions = {}
    for row in map_rows:
        map_positions.setdefault(key(row), set()).add(int(float(row["aa_position"])))

    flow_token_counts = {}
    for row in raw_rows:
        hap = row["immune_escape_haplotype"]
        if hap and hap not in {"__OTHER__", "__FILTERED_IEP__"}:
            flow_token_counts[key(row)] = max(
                flow_token_counts.get(key(row), 0),
                len(hap.split("|")),
            )

    all_keys = sorted(set(map_positions) | set(flow_token_counts))
    matched_keys = {
        k for k in all_keys
        if len(map_positions.get(k, set())) == flow_token_counts.get(k, 0)
        and flow_token_counts.get(k, 0) > 0
    }

    report_rows = []
    for k in all_keys:
        map_n = len(map_positions.get(k, set()))
        flow_n = flow_token_counts.get(k, 0)
        report_rows.append({
            "species": k[0],
            "gene_id": k[1],
            "antigen": k[2],
            "map_positions": map_n,
            "raw_flow_tokens": flow_n,
            "status": "OK" if k in matched_keys else "MISMATCH",
        })

    out_prefix = Path(args.out_prefix)
    out_map = f"{out_prefix}.iep_position_map.tsv"
    out_raw = f"{out_prefix}.raw_flows.tsv"
    out_report = f"{out_prefix}.match_report.tsv"

    write_tsv(out_map, map_fields, [r for r in map_rows if key(r) in matched_keys])
    write_tsv(out_raw, raw_fields, [r for r in raw_rows if key(r) in matched_keys])
    write_tsv(
        out_report,
        ["species", "gene_id", "antigen", "map_positions", "raw_flow_tokens", "status"],
        report_rows,
    )

    if args.alluvial_flows:
        alluvial_fields, alluvial_rows = read_tsv(args.alluvial_flows)
        required_alluvial = {"species", "gene_id", "antigen"}
        if not required_alluvial.issubset(alluvial_fields):
            missing = sorted(required_alluvial - set(alluvial_fields))
            raise SystemExit(f"Alluvial flows missing columns: {', '.join(missing)}")
        out_alluvial = f"{out_prefix}.alluvial_flows.tsv"
        write_tsv(
            out_alluvial,
            alluvial_fields,
            [r for r in alluvial_rows if key(r) in matched_keys],
        )
    else:
        out_alluvial = None

    print(f"Matched antigens: {len(matched_keys)}")
    print(f"Retained AA positions: {sum(1 for r in map_rows if key(r) in matched_keys)}")
    print(f"Wrote: {out_map}")
    print(f"Wrote: {out_raw}")
    if out_alluvial:
        print(f"Wrote: {out_alluvial}")
    print(f"Wrote: {out_report}")


if __name__ == "__main__":
    main()
