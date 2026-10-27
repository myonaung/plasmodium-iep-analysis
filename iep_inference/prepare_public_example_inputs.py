#!/usr/bin/env python3
"""Create a small, deidentified example dataset for candidate-IEP inference.

The script deliberately does not write an ID key. Participant and sample IDs are
replaced independently, sampling dates are converted to within-participant order,
the exposure covariate is replaced by a deterministic synthetic value, and only a
short segment of each aligned antigen sequence is retained.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--exposure", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--sample-column", default="samples")
    parser.add_argument("--participant-column", default="studyid")
    parser.add_argument("--time-column", default="study_day")
    parser.add_argument("--status-column", default="fever")
    parser.add_argument("--exposure-column", default="molFOI")
    parser.add_argument("--asymptomatic-value", default="A")
    parser.add_argument("--symptomatic-value", default="S")
    parser.add_argument("--participants", type=int, default=24)
    parser.add_argument("--max-samples-per-participant", type=int, default=3)
    parser.add_argument("--trim-start", type=int, default=150,
                        help="One-based first aligned residue to retain.")
    parser.add_argument("--trim-length", type=int, default=180)
    parser.add_argument("--salt", default="public-iep-example-v1")
    return parser.parse_args()


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header: str | None = None
    sequence: list[str] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header.split("_", 1)[0]] = "".join(sequence)
                header = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)
    if header is not None:
        records[header.split("_", 1)[0]] = "".join(sequence)
    return records


def stable_hash(value: str, salt: str) -> str:
    return hashlib.sha256(f"{salt}:{value}".encode()).hexdigest()


def numeric_time(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("inf")


def choose_rows(rows: list[dict[str, str]], status_col: str,
                max_samples: int) -> list[dict[str, str]]:
    rows = sorted(rows, key=lambda row: numeric_time(row["_time"]))
    chosen = [rows[0]]
    opposite = next(
        (row for row in rows[1:] if row[status_col] != rows[0][status_col]),
        None,
    )
    if opposite is not None:
        chosen.append(opposite)
    if rows[-1] not in chosen and len(chosen) < max_samples:
        chosen.append(rows[-1])
    for row in rows:
        if row not in chosen and len(chosen) < max_samples:
            chosen.append(row)
    return sorted(chosen, key=lambda row: numeric_time(row["_time"]))


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    sequences = read_fasta(args.fasta)
    with args.metadata.open(newline="") as handle:
        metadata = list(csv.DictReader(handle))
    with args.exposure.open(newline="") as handle:
        exposure_rows = list(csv.DictReader(handle))

    required_metadata = {
        args.sample_column,
        args.participant_column,
        args.time_column,
        args.status_column,
    }
    if not metadata or not required_metadata.issubset(metadata[0]):
        raise ValueError(f"Metadata must contain: {sorted(required_metadata)}")
    if not exposure_rows or not {
        args.participant_column,
        args.exposure_column,
    }.issubset(exposure_rows[0]):
        raise ValueError("Exposure file lacks the participant or exposure column")

    exposure_participants = {
        row[args.participant_column]
        for row in exposure_rows
        if row.get(args.participant_column)
    }
    allowed_status = {args.asymptomatic_value, args.symptomatic_value}
    by_participant: dict[str, list[dict[str, str]]] = defaultdict(list)
    seen_samples: set[str] = set()
    for row in metadata:
        sample = row.get(args.sample_column, "")
        participant = row.get(args.participant_column, "")
        status = row.get(args.status_column, "")
        if (
            sample in sequences
            and sample not in seen_samples
            and participant in exposure_participants
            and status in allowed_status
        ):
            clean_row = dict(row)
            clean_row["_time"] = row.get(args.time_column, "")
            by_participant[participant].append(clean_row)
            seen_samples.add(sample)

    eligible = [
        participant
        for participant, rows in by_participant.items()
        if len(rows) >= 2
        and allowed_status.issubset({row[args.status_column] for row in rows})
    ]
    eligible.sort(key=lambda value: stable_hash(value, args.salt))
    selected_participants = eligible[: args.participants]
    if len(selected_participants) < args.participants:
        raise ValueError(
            f"Only {len(selected_participants)} eligible participants were found"
        )

    start = args.trim_start - 1
    end = start + args.trim_length
    participant_ids = {
        original: f"P{index:03d}"
        for index, original in enumerate(selected_participants, start=1)
    }

    public_metadata: list[dict[str, object]] = []
    public_sequences: list[tuple[str, str]] = []
    sample_number = 1
    for original_participant in selected_participants:
        selected_rows = choose_rows(
            by_participant[original_participant],
            args.status_column,
            args.max_samples_per_participant,
        )
        for order, row in enumerate(selected_rows, start=1):
            original_sample = row[args.sample_column]
            public_sample = f"S{sample_number:03d}"
            sample_number += 1
            segment = sequences[original_sample][start:end]
            if len(segment) != args.trim_length:
                raise ValueError(
                    f"Sequence {original_sample} is too short for the requested segment"
                )
            public_metadata.append({
                "sample_id": public_sample,
                "participant_id": participant_ids[original_participant],
                "timepoint": order,
                "clinical_status": row[args.status_column],
            })
            public_sequences.append((public_sample, segment))

    metadata_path = args.outdir / "sample_metadata.csv"
    with metadata_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample_id", "participant_id", "timepoint", "clinical_status"],
        )
        writer.writeheader()
        writer.writerows(public_metadata)

    exposure_path = args.outdir / "exposure_covariate.csv"
    with exposure_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["participant_id", "molFOI"])
        writer.writeheader()
        for original_participant in selected_participants:
            digest = int(stable_hash(original_participant, args.salt + "-exposure")[:8], 16)
            synthetic_exposure = 1.0 + (digest % 401) / 100.0
            writer.writerow({
                "participant_id": participant_ids[original_participant],
                "molFOI": f"{synthetic_exposure:.2f}",
            })

    fasta_path = args.outdir / "ama1_example_alignment.fasta"
    with fasta_path.open("w") as handle:
        for sample, sequence in public_sequences:
            handle.write(f">{sample}_ama1_example\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset : offset + 80] + "\n")

    print(f"Wrote {len(public_sequences)} samples from "
          f"{len(selected_participants)} deidentified participants")
    print(f"Retained aligned residues {args.trim_start}-{end}")
    print(f"Output directory: {args.outdir}")


if __name__ == "__main__":
    main()
