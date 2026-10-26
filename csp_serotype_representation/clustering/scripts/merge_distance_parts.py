#!/usr/bin/env python3

import argparse
import glob
from pathlib import Path

import numpy as np


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Merge genetic-distance matrix parts. If part prefixes contain "
            "*.numerator.npy and *.counts.npy, these are summed before producing "
            "the final distance matrix. Otherwise, full *.npy distance matrices "
            "are merged by taking the mean of non-missing entries."
        )
    )
    p.add_argument(
        "--part-prefixes",
        default="",
        help="Comma-separated part prefixes, e.g. out/dist_part001,out/dist_part002",
    )
    p.add_argument(
        "--part-glob",
        default="",
        help="Glob for part prefixes or part matrices, e.g. 'out/dist_part*.npy'",
    )
    p.add_argument("--out-prefix", required=True, help="Output prefix for merged files")
    p.add_argument(
        "--sample-ids",
        default="",
        help="Optional sample ID file. If omitted, the first part .sample_ids.txt is used if present.",
    )
    p.add_argument("--keep-numerator", action="store_true", help="Keep merged numerator matrix")
    p.add_argument("--keep-counts", action="store_true", help="Keep merged callable-site count matrix")
    return p.parse_args()


def clean_prefix(path):
    path = str(path)
    for suffix in [".numerator.npy", ".counts.npy", ".npy"]:
        if path.endswith(suffix):
            return path[: -len(suffix)]
    return path


def get_part_prefixes(args):
    prefixes = []
    if args.part_prefixes.strip():
        prefixes.extend([x.strip() for x in args.part_prefixes.split(",") if x.strip()])
    if args.part_glob.strip():
        prefixes.extend([clean_prefix(x) for x in glob.glob(args.part_glob)])
    prefixes = sorted(dict.fromkeys(clean_prefix(x) for x in prefixes))
    if not prefixes:
        raise ValueError("Provide --part-prefixes or --part-glob.")
    return prefixes


def read_ids(path):
    with open(path, encoding="utf-8") as fh:
        return [x.strip() for x in fh if x.strip()]


def write_ids(path, ids):
    with open(path, "w", encoding="utf-8") as fh:
        for sample_id in ids:
            fh.write(f"{sample_id}\n")


def find_sample_ids(prefixes, requested):
    if requested:
        return read_ids(requested)
    for prefix in prefixes:
        candidate = Path(f"{prefix}.sample_ids.txt")
        if candidate.exists():
            return read_ids(candidate)
    return None


def all_have_numerator_counts(prefixes):
    return all(Path(f"{p}.numerator.npy").exists() and Path(f"{p}.counts.npy").exists() for p in prefixes)


def merge_numerator_counts(prefixes):
    numerator = None
    counts = None
    for prefix in prefixes:
        num = np.load(f"{prefix}.numerator.npy", mmap_mode="r")
        cnt = np.load(f"{prefix}.counts.npy", mmap_mode="r")
        if numerator is None:
            numerator = np.zeros(num.shape, dtype=np.float64)
            counts = np.zeros(cnt.shape, dtype=np.uint64)
        if num.shape != numerator.shape or cnt.shape != counts.shape:
            raise ValueError(f"Part has incompatible shape: {prefix}")
        numerator += np.asarray(num, dtype=np.float64)
        counts += np.asarray(cnt, dtype=np.uint64)
    out = np.full(numerator.shape, np.nan, dtype=np.float32)
    valid = counts > 0
    out[valid] = (numerator[valid] / counts[valid]).astype(np.float32)
    return out, numerator.astype(np.float32), counts


def merge_distance_matrices(prefixes):
    total = None
    counts = None
    for prefix in prefixes:
        matrix_path = f"{prefix}.npy"
        if not Path(matrix_path).exists():
            raise FileNotFoundError(matrix_path)
        mat = np.load(matrix_path, mmap_mode="r")
        if total is None:
            total = np.zeros(mat.shape, dtype=np.float64)
            counts = np.zeros(mat.shape, dtype=np.uint32)
        if mat.shape != total.shape:
            raise ValueError(f"Part has incompatible shape: {prefix}")
        arr = np.asarray(mat, dtype=np.float64)
        valid = np.isfinite(arr)
        total[valid] += arr[valid]
        counts[valid] += 1
    out = np.full(total.shape, np.nan, dtype=np.float32)
    valid = counts > 0
    out[valid] = (total[valid] / counts[valid]).astype(np.float32)
    return out


def main():
    args = parse_args()
    prefixes = get_part_prefixes(args)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    sample_ids = find_sample_ids(prefixes, args.sample_ids)

    if all_have_numerator_counts(prefixes):
        matrix, numerator, counts = merge_numerator_counts(prefixes)
        if args.keep_numerator:
            np.save(f"{out_prefix}.numerator.npy", numerator)
        if args.keep_counts:
            np.save(f"{out_prefix}.counts.npy", counts)
    else:
        matrix = merge_distance_matrices(prefixes)

    np.save(f"{out_prefix}.npy", matrix)
    if sample_ids is not None:
        if len(sample_ids) != matrix.shape[0]:
            raise ValueError("Sample ID count does not match merged matrix size.")
        write_ids(f"{out_prefix}.sample_ids.txt", sample_ids)

    with open(f"{out_prefix}.summary.txt", "w", encoding="utf-8") as fh:
        fh.write(f"n_parts\t{len(prefixes)}\n")
        fh.write(f"matrix_shape\t{matrix.shape}\n")
        fh.write(f"used_numerator_counts\t{all_have_numerator_counts(prefixes)}\n")
        fh.write(f"sample_ids_written\t{sample_ids is not None}\n")

    print(f"Wrote merged matrix: {out_prefix}.npy")
    if sample_ids is not None:
        print(f"Wrote sample IDs: {out_prefix}.sample_ids.txt")
    print(f"Wrote summary: {out_prefix}.summary.txt")


if __name__ == "__main__":
    main()
