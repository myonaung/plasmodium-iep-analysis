#!/usr/bin/env python3

import argparse
import os
import subprocess
import time
from pathlib import Path

import numpy as np


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Recalculate a PASS-only genetic distance matrix from a multi-sample VCF using "
            "read-depth-derived non-reference allele frequencies."
        )
    )
    p.add_argument("--vcf", required=True, help="PASS-filtered multi-sample VCF or VCF.GZ")
    p.add_argument("--out-prefix", required=True, help="Output prefix for matrix and metadata")
    p.add_argument("--bcftools", default="bcftools", help="bcftools executable")
    p.add_argument("--threads", type=int, default=4, help="bcftools threads")
    p.add_argument("--site-chunk", type=int, default=256, help="Number of variant sites to accumulate per chunk")
    p.add_argument("--sample-block", type=int, default=1024, help="Sample block size for matrix multiplication")
    p.add_argument("--min-total-depth", type=int, default=5, help="Minimum total depth to estimate f_i")
    p.add_argument("--min-allele-depth", type=int, default=2, help="Minimum allele depth before forcing 0 or 1")
    p.add_argument("--min-allele-frac", type=float, default=0.05, help="Minimum allele fraction before forcing 0 or 1")
    p.add_argument("--max-sites", type=int, default=None, help="Optional limit for testing on a small number of sites")
    p.add_argument("--keep-counts", action="store_true", help="Keep the pairwise valid-site count matrix")
    p.add_argument("--keep-numerator", action="store_true", help="Keep the pairwise distance numerator matrix")
    return p.parse_args()


def run_capture(cmd):
    return subprocess.check_output(cmd, text=True)


def get_sample_ids(vcf, bcftools):
    cmd = [bcftools, "query", "-l", vcf]
    return [x.strip() for x in run_capture(cmd).splitlines() if x.strip()]


def open_variant_stream(vcf, bcftools, threads):
    query_fmt = r"%REF\t%ALT[\t%AD]\n"
    cmd = [
        bcftools,
        "query",
        "--threads",
        str(threads),
        "-f",
        query_fmt,
        vcf,
    ]
    return subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1)


def parse_ad_site(line, n_samples, min_total_depth, min_allele_depth, min_allele_frac):
    parts = line.rstrip("\n").split("\t")
    if len(parts) < n_samples + 2:
        return None

    ref = parts[0]
    alt = parts[1]
    if len(ref) != 1 or len(alt) != 1 or "," in alt or alt == ".":
        return None

    f = np.zeros(n_samples, dtype=np.float32)
    m = np.zeros(n_samples, dtype=np.float32)

    for idx, ad in enumerate(parts[2:]):
        if not ad or ad == ".":
            continue
        vals = ad.split(",")
        if len(vals) < 2 or vals[0] == "." or vals[1] == ".":
            continue

        try:
            ref_depth = int(vals[0])
            alt_depth = int(vals[1])
        except ValueError:
            continue

        total_depth = ref_depth + alt_depth
        if total_depth < min_total_depth:
            continue

        alt_frac = alt_depth / total_depth
        ref_frac = ref_depth / total_depth

        if alt_depth < min_allele_depth or alt_frac < min_allele_frac:
            fi = 0.0
        elif ref_depth < min_allele_depth or ref_frac < min_allele_frac:
            fi = 1.0
        else:
            fi = alt_frac

        f[idx] = fi
        m[idx] = 1.0

    return f, m


def accumulate_chunk(numerator_mm, counts_mm, x_chunk, m_chunk, sample_block):
    n_samples = x_chunk.shape[0]
    for i0 in range(0, n_samples, sample_block):
        i1 = min(i0 + sample_block, n_samples)
        xi = x_chunk[i0:i1, :]
        mi = m_chunk[i0:i1, :]

        for j0 in range(i0, n_samples, sample_block):
            j1 = min(j0 + sample_block, n_samples)
            xj = x_chunk[j0:j1, :]
            mj = m_chunk[j0:j1, :]

            # Sum over chunk sites:
            # d_ij = f_i + f_j - 2 f_i f_j, only when both samples have callable data.
            add_num = (xi @ mj.T) + (mi @ xj.T) - (2.0 * (xi @ xj.T))
            add_cnt = np.rint(mi @ mj.T).astype(np.uint32)

            numerator_mm[i0:i1, j0:j1] += add_num.astype(np.float32, copy=False)
            counts_mm[i0:i1, j0:j1] += add_cnt

            if j0 != i0:
                numerator_mm[j0:j1, i0:i1] += add_num.T.astype(np.float32, copy=False)
                counts_mm[j0:j1, i0:i1] += add_cnt.T


def write_summary(summary_path, stats):
    with open(summary_path, "w", encoding="utf-8") as fh:
        for key, value in stats.items():
            fh.write(f"{key}\t{value}\n")


def main():
    args = parse_args()
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    sample_ids = get_sample_ids(args.vcf, args.bcftools)
    n_samples = len(sample_ids)
    if n_samples == 0:
        raise ValueError("No samples found in VCF.")

    sample_ids_path = Path(f"{out_prefix}.sample_ids.txt")
    with open(sample_ids_path, "w", encoding="utf-8") as fh:
        for sample_id in sample_ids:
            fh.write(f"{sample_id}\n")

    numerator_path = Path(f"{out_prefix}.numerator.npy")
    counts_path = Path(f"{out_prefix}.counts.npy")
    matrix_path = Path(f"{out_prefix}.npy")
    summary_path = Path(f"{out_prefix}.summary.txt")

    numerator_mm = np.lib.format.open_memmap(
        numerator_path,
        mode="w+",
        dtype=np.float32,
        shape=(n_samples, n_samples),
    )
    counts_mm = np.lib.format.open_memmap(
        counts_path,
        mode="w+",
        dtype=np.uint32,
        shape=(n_samples, n_samples),
    )
    numerator_mm[:] = 0
    counts_mm[:] = 0

    start = time.time()
    processed_sites = 0
    kept_sites = 0
    skipped_non_biallelic = 0
    callable_total = 0

    x_sites = []
    m_sites = []
    stopped_early = False

    proc = open_variant_stream(args.vcf, args.bcftools, args.threads)
    assert proc.stdout is not None

    try:
        for line in proc.stdout:
            processed_sites += 1
            parsed = parse_ad_site(
                line,
                n_samples=n_samples,
                min_total_depth=args.min_total_depth,
                min_allele_depth=args.min_allele_depth,
                min_allele_frac=args.min_allele_frac,
            )
            if parsed is None:
                skipped_non_biallelic += 1
                if args.max_sites is not None and processed_sites >= args.max_sites:
                    stopped_early = True
                    break
                continue

            f, m = parsed
            if np.any(m):
                kept_sites += 1
                callable_total += int(m.sum())
                x_sites.append(f)
                m_sites.append(m)

            if len(x_sites) >= args.site_chunk:
                x_chunk = np.stack(x_sites, axis=1)
                m_chunk = np.stack(m_sites, axis=1)
                accumulate_chunk(numerator_mm, counts_mm, x_chunk, m_chunk, args.sample_block)
                x_sites = []
                m_sites = []

            if args.max_sites is not None and processed_sites >= args.max_sites:
                stopped_early = True
                break

        if x_sites:
            x_chunk = np.stack(x_sites, axis=1)
            m_chunk = np.stack(m_sites, axis=1)
            accumulate_chunk(numerator_mm, counts_mm, x_chunk, m_chunk, args.sample_block)

    finally:
        if proc.stdout:
            proc.stdout.close()
        stderr = proc.stderr.read() if proc.stderr else ""
        if stopped_early:
            proc.terminate()
        return_code = proc.wait()
        if return_code != 0 and not stopped_early:
            raise RuntimeError(f"bcftools query failed with exit code {return_code}:\n{stderr}")

    matrix_mm = np.lib.format.open_memmap(
        matrix_path,
        mode="w+",
        dtype=np.float32,
        shape=(n_samples, n_samples),
    )
    matrix_mm[:] = np.nan

    for i0 in range(0, n_samples, args.sample_block):
        i1 = min(i0 + args.sample_block, n_samples)
        num_block = numerator_mm[i0:i1, :]
        cnt_block = counts_mm[i0:i1, :]
        out_block = np.full(num_block.shape, np.nan, dtype=np.float32)
        valid = cnt_block > 0
        out_block[valid] = num_block[valid] / cnt_block[valid]
        matrix_mm[i0:i1, :] = out_block

    matrix_mm.flush()
    numerator_mm.flush()
    counts_mm.flush()

    elapsed_sec = time.time() - start
    stats = {
        "vcf": args.vcf,
        "n_samples": n_samples,
        "processed_sites": processed_sites,
        "kept_sites": kept_sites,
        "skipped_non_biallelic_or_unparseable_sites": skipped_non_biallelic,
        "mean_callable_samples_per_kept_site": callable_total / kept_sites if kept_sites else 0,
        "site_chunk": args.site_chunk,
        "sample_block": args.sample_block,
        "min_total_depth": args.min_total_depth,
        "min_allele_depth": args.min_allele_depth,
        "min_allele_frac": args.min_allele_frac,
        "elapsed_seconds": round(elapsed_sec, 2),
    }
    write_summary(summary_path, stats)

    if not args.keep_numerator:
        try:
            os.remove(numerator_path)
        except OSError:
            pass
    if not args.keep_counts:
        try:
            os.remove(counts_path)
        except OSError:
            pass

    print(f"Wrote matrix: {matrix_path}")
    print(f"Wrote sample IDs: {sample_ids_path}")
    print(f"Wrote summary: {summary_path}")


if __name__ == "__main__":
    main()
