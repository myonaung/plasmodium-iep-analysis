#!/usr/bin/env python3

import argparse
import csv

import numpy as np


class UnionFind:
    def __init__(self, items):
        self.parent = {x: x for x in items}
        self.rank = {x: 0 for x in items}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def parse_args():
    p = argparse.ArgumentParser(description="Cluster similar samples from a genetic distance matrix.")
    p.add_argument("--matrix", required=True)
    p.add_argument("--sample-ids", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--out-clusters", required=True)
    p.add_argument("--out-pairs", required=True)
    p.add_argument("--distance-threshold", required=True, type=float)
    p.add_argument("--site-col", default="site")
    p.add_argument("--year-col", default="sample_year")
    p.add_argument("--year-threshold", type=int, default=0)
    return p.parse_args()


def read_ids(path):
    with open(path, encoding="utf-8") as fh:
        return [x.strip() for x in fh if x.strip()]


def read_metadata(path, site_col, year_col):
    with open(path, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    out = {}
    for row in rows:
        sample_id = row["sample_id"]
        year = row.get(year_col, "")
        out[sample_id] = {
            "site": row.get(site_col, ""),
            "sample_year": int(year) if year not in ("", None) else None,
        }
    return out


def main():
    args = parse_args()

    sample_ids = read_ids(args.sample_ids)
    metadata = read_metadata(args.metadata, args.site_col, args.year_col)
    matrix = np.load(args.matrix, mmap_mode="r")

    keep = [sample_id for sample_id in sample_ids if sample_id in metadata]
    idx = [sample_ids.index(sample_id) for sample_id in keep]
    matrix = matrix[np.ix_(idx, idx)]

    uf = UnionFind(keep)
    linked_pairs = []

    for i, sample_a in enumerate(keep):
        meta_a = metadata[sample_a]
        for j in range(i + 1, len(keep)):
            sample_b = keep[j]
            meta_b = metadata[sample_b]

            if not meta_a["site"] or not meta_b["site"] or meta_a["site"] != meta_b["site"]:
                continue

            year_a = meta_a["sample_year"]
            year_b = meta_b["sample_year"]
            if year_a is not None and year_b is not None and abs(year_a - year_b) > args.year_threshold:
                continue

            d = float(matrix[i, j])
            if np.isnan(d):
                continue
            if d <= args.distance_threshold:
                uf.union(sample_a, sample_b)
                linked_pairs.append((sample_a, sample_b, d, meta_a["site"], year_a, year_b))

    groups = {}
    for sample_id in keep:
        root = uf.find(sample_id)
        groups.setdefault(root, []).append(sample_id)

    rows = []
    cluster_num = 0
    for root, members in sorted(groups.items(), key=lambda kv: min(kv[1])):
        cluster_num += 1
        site = metadata[members[0]]["site"]
        cluster_id = f"{site}_cl{cluster_num}"
        for sample_id in sorted(members):
            rows.append(
                {
                    "sample_id": sample_id,
                    "cluster_id": cluster_id,
                    "cluster_size": len(members),
                    "site": metadata[sample_id]["site"],
                    "sample_year": metadata[sample_id]["sample_year"],
                }
            )

    with open(args.out_pairs, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["sample_a", "sample_b", "genetic_distance", "site", "year_a", "year_b"])
        writer.writerows(linked_pairs)

    with open(args.out_clusters, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["sample_id", "cluster_id", "cluster_size", "site", "sample_year"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
