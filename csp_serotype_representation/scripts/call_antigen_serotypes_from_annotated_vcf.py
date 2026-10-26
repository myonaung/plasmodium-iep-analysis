#!/usr/bin/env python3

import argparse
import csv
import gzip
import re
from collections import defaultdict


AA3_TO_1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Ter": "*",
    "Sec": "U",
    "Pyl": "O",
}

AA1 = {
    "A": "A", "R": "R", "N": "N", "D": "D", "C": "C",
    "Q": "Q", "E": "E", "G": "G", "H": "H", "I": "I",
    "L": "L", "K": "K", "M": "M", "F": "F", "P": "P",
    "S": "S", "T": "T", "W": "W", "Y": "Y", "V": "V",
    "*": "*", "X": "X",
}

CODON_TO_AA = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Call serotype-like sample classes from a SnpEff-annotated gene-region VCF, "
            "either from targeted amino-acid positions or from all variants in the gene."
        )
    )
    p.add_argument("--vcf", required=True, help="Annotated gene-region VCF(.gz)")
    p.add_argument("--gene", required=True, help="Gene identifier to match in ANN or SNPEFF_* annotations")
    p.add_argument("--mode", choices=["aa_positions", "all_variants"], required=True)
    p.add_argument("--aa-positions", help="Comma-separated amino-acid positions for aa_positions mode")
    p.add_argument(
        "--reference-aa",
        help="Optional reference residues like 352:K,359:E for targeted positions",
    )
    p.add_argument("--ann-field", default="ANN", help="INFO tag holding SnpEff ANN annotations when present")
    p.add_argument("--out", required=True, help="Output TSV with sample calls")
    p.add_argument("--variant-summary", help="Optional per-variant/per-position summary TSV")
    return p.parse_args()


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt", encoding="utf-8")


def parse_info(info_str):
    info = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        elif item:
            info[item] = True
    return info


def parse_reference_aa(arg):
    out = {}
    if not arg:
        return out
    for item in arg.split(","):
        item = item.strip()
        if not item:
            continue
        pos, aa = item.split(":", 1)
        out[int(pos)] = aa.strip()
    return out


def parse_target_positions(arg):
    if not arg:
        raise ValueError("--aa-positions is required in aa_positions mode")
    return sorted({int(x.strip()) for x in arg.split(",") if x.strip()})


def parse_gt(sample_field, gt_index):
    if sample_field in ("", ".", "./.", ".|."):
        return None
    pieces = sample_field.split(":")
    if gt_index >= len(pieces):
        return None
    gt = pieces[gt_index]
    if gt in ("", ".", "./.", ".|."):
        return None
    alleles = re.split(r"[\/|]", gt)
    if not alleles or any(a == "." for a in alleles):
        return None
    return alleles


def classify_genotype(alleles):
    if alleles is None:
        return ("MISSING", None)
    uniq = set(alleles)
    if uniq == {"0"}:
        return ("REF", None)
    if len(uniq) == 1 and "0" not in uniq:
        return ("ALT", next(iter(uniq)))
    return ("MIXED", None)


def ann_matches_gene(fields, gene):
    gene_tokens = {
        fields[3],
        fields[4],
        fields[6],
        fields[4].split(".", 1)[0] if fields[4] else fields[4],
        fields[6].split(".", 1)[0] if fields[6] else fields[6],
    }
    return gene in gene_tokens


def parse_hgvs_protein(hgvs_p):
    if not hgvs_p or hgvs_p in (".", ""):
        return None
    for pattern in (
        r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=|Ter|\*)$",
        r"p\.([A-Z\*])(\d+)([A-Z\*=])$",
        r"([A-Z\*])(\d+)([A-Z\*=])$",
        r"([A-Z])(\d+)$",
    ):
        m = re.match(pattern, hgvs_p)
        if not m:
            continue
        groups = m.groups()
        if len(groups) == 3:
            ref_token, pos, alt_token = groups
        else:
            ref_token, pos = groups
            alt_token = "="
        if len(ref_token) == 1:
            ref = AA1.get(ref_token)
            alt = ref if alt_token == "=" else AA1.get(alt_token)
        else:
            ref = AA3_TO_1.get(ref_token)
            alt = ref if alt_token == "=" else AA3_TO_1.get(alt_token, "*" if alt_token == "*" else None)
        if ref is None or alt is None:
            return None
        return int(pos), ref, alt
    return None


def parse_hgvs_coding(hgvs_c):
    if not hgvs_c or hgvs_c in (".", ""):
        return None
    m = re.match(r"c\.([1-9]\d*)([ACGT])>([ACGT])$", hgvs_c)
    if not m:
        return None
    cds_pos, ref_nt, alt_nt = m.groups()
    return int(cds_pos), ref_nt, alt_nt


def parse_snpeff_codon_change(codon_change):
    if not codon_change or codon_change in (".", "") or "/" not in codon_change:
        return None
    ref_codon, alt_codon = codon_change.split("/", 1)
    ref_codon = ref_codon.upper()
    alt_codon = alt_codon.upper()
    if len(ref_codon) != 3 or len(alt_codon) != 3:
        return None
    diffs = [(idx, ref_nt, alt_nt) for idx, (ref_nt, alt_nt) in enumerate(zip(ref_codon, alt_codon)) if ref_nt != alt_nt]
    if not diffs:
        return None
    return diffs


def parse_aa_pos_field(aa_field):
    if not aa_field or aa_field in (".", ""):
        return None
    value = aa_field.split("/", 1)[0]
    if value in ("", "."):
        return None
    try:
        return int(value)
    except ValueError:
        return None


def parse_effects_for_gene(info, ann_field, gene, alt_alleles, allow_ambiguous_multiallelic=False):
    records = []

    if ann_field in info:
        for ann in info[ann_field].split(","):
            fields = ann.split("|")
            if len(fields) < 14:
                continue
            if not ann_matches_gene(fields, gene):
                continue
            hgvs_c = fields[9]
            hgvs_p = fields[10]
            aa_pos = parse_aa_pos_field(fields[13])
            parsed = parse_hgvs_protein(hgvs_p)
            coding_change = parse_hgvs_coding(hgvs_c)
            if parsed is not None:
                aa_pos2, ref_aa, alt_aa = parsed
                aa_pos = aa_pos if aa_pos is not None else aa_pos2
            else:
                ref_aa = None
                alt_aa = None

            codon_idx = None
            ref_nt = None
            alt_nt = None
            if coding_change is not None:
                cds_pos, ref_nt, alt_nt = coding_change
                codon_idx = (cds_pos - 1) % 3

            records.append(
                {
                    "allele": fields[0],
                    "annotation": fields[1],
                    "hgvs_c": hgvs_c,
                    "hgvs_p": hgvs_p,
                    "aa_pos": aa_pos,
                    "ref_aa": ref_aa,
                    "alt_aa": alt_aa,
                    "ref_nt": ref_nt,
                    "alt_nt": alt_nt,
                    "codon_idx": codon_idx,
                }
            )
        return records

    gene_name = info.get("SNPEFF_GENE_NAME")
    transcript_id = info.get("SNPEFF_TRANSCRIPT_ID")
    gene_tokens = {
        gene_name,
        transcript_id,
        transcript_id.split(".", 1)[0] if transcript_id else transcript_id,
    }
    if gene not in gene_tokens:
        return []
    hgvs_p = info.get("SNPEFF_AMINO_ACID_CHANGE")
    parsed = parse_hgvs_protein(hgvs_p)
    if parsed is not None:
        aa_pos, ref_aa, alt_aa = parsed
    else:
        aa_pos = None
        ref_aa = None
        alt_aa = None

    codon_change = info.get("SNPEFF_CODON_CHANGE")
    diffs = parse_snpeff_codon_change(codon_change)
    if diffs and len(alt_alleles) == 1:
        codon_idx, ref_nt, alt_nt = diffs[0]
    else:
        codon_idx = None
        ref_nt = None
        alt_nt = None

    if len(alt_alleles) == 1:
        records.append(
            {
                "allele": alt_alleles[0],
                "annotation": info.get("SNPEFF_EFFECT", ""),
                "hgvs_c": "",
                "hgvs_p": hgvs_p or "",
                "aa_pos": aa_pos,
                "ref_aa": ref_aa,
                "alt_aa": alt_aa,
                "ref_nt": ref_nt,
                "alt_nt": alt_nt,
                "codon_idx": codon_idx,
            }
        )
        return records

    if not allow_ambiguous_multiallelic:
        return []

    records.append(
        {
            "allele": None,
            "annotation": info.get("SNPEFF_EFFECT", ""),
            "hgvs_c": "",
            "hgvs_p": hgvs_p or "",
            "aa_pos": aa_pos,
            "ref_aa": ref_aa,
            "alt_aa": alt_aa,
            "ref_nt": None,
            "alt_nt": None,
            "codon_idx": None,
        }
    )
    return records


def infer_alt_residue(ref_aa, ref_nts, alt_nts):
    if ref_aa is None or not alt_nts:
        return None

    candidates = []
    for codon, aa in CODON_TO_AA.items():
        if aa != ref_aa:
            continue
        if any(codon[idx] != nt for idx, nt in ref_nts.items()):
            continue
        candidates.append(codon)

    if not candidates:
        return None

    alt_aas = set()
    for codon in candidates:
        codon_list = list(codon)
        for idx, alt_nt in alt_nts.items():
            codon_list[idx] = alt_nt
        alt_codon = "".join(codon_list)
        alt_aa = CODON_TO_AA.get(alt_codon)
        if alt_aa is None:
            return None
        alt_aas.add(alt_aa)

    if len(alt_aas) == 1:
        return next(iter(alt_aas))
    return None


def aa_mode(vcf_path, gene, ann_field, target_positions, reference_aa):
    with open_text(vcf_path) as handle:
        sample_ids = []
        per_sample = None
        pos_ref = dict(reference_aa)
        variant_rows = []

        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_ids = header[9:]
                per_sample = {
                    sample_id: {
                        pos: {
                            "saw_ref": False,
                            "saw_missing": False,
                            "saw_mixed": False,
                            "residue": pos_ref.get(pos),
                            "ref_aa": pos_ref.get(pos),
                            "ref_nts": {},
                            "alt_nts": {},
                            "conflict": False,
                        }
                        for pos in target_positions
                    }
                    for sample_id in sample_ids
                }
                continue

            fields = line.rstrip("\n").split("\t")
            chrom, pos, _vid, ref, alt, _qual, _flt, info_str, fmt = fields[:9]
            info = parse_info(info_str)
            alt_alleles = alt.split(",")
            ann_records = parse_effects_for_gene(info, ann_field, gene, alt_alleles)
            if not ann_records:
                continue
            allele_effects = defaultdict(dict)
            for rec in ann_records:
                allele = rec["allele"]
                hgvs_c = rec["hgvs_c"]
                hgvs_p = rec["hgvs_p"]
                aa_pos = rec["aa_pos"]
                ref_aa = rec["ref_aa"]
                alt_aa = rec["alt_aa"]
                if aa_pos not in target_positions:
                    continue
                if aa_pos is None:
                    continue
                if aa_pos not in pos_ref and ref_aa is not None:
                    pos_ref[aa_pos] = ref_aa
                if allele not in alt_alleles:
                    continue
                alt_index = str(alt_alleles.index(allele) + 1)
                allele_effects[alt_index][aa_pos] = {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": allele,
                    "ref_aa": ref_aa,
                    "alt_aa": alt_aa,
                    "ref_nt": rec["ref_nt"],
                    "alt_nt": rec["alt_nt"],
                    "codon_idx": rec["codon_idx"],
                    "annotation": rec["annotation"],
                    "hgvs_c": hgvs_c,
                    "hgvs_p": hgvs_p,
                }

            if not allele_effects:
                continue

            format_fields = fmt.split(":")
            if "GT" not in format_fields:
                continue
            gt_index = format_fields.index("GT")
            affected_positions = sorted({p for pos_map in allele_effects.values() for p in pos_map.keys()})

            for sample_id, sample_field in zip(sample_ids, fields[9:]):
                state, alt_idx = classify_genotype(parse_gt(sample_field, gt_index))
                if state == "ALT" and alt_idx in allele_effects:
                    pos_iter = sorted(allele_effects[alt_idx].keys())
                else:
                    pos_iter = affected_positions
                for pos0 in pos_iter:
                    sample_pos = per_sample[sample_id][pos0]
                    if state == "REF":
                        sample_pos["saw_ref"] = True
                        if sample_pos["residue"] is None:
                            sample_pos["residue"] = pos_ref.get(pos0, sample_pos["residue"])
                        if sample_pos["ref_aa"] is None:
                            sample_pos["ref_aa"] = pos_ref.get(pos0)
                    elif state == "MISSING":
                        sample_pos["saw_missing"] = True
                    elif state == "MIXED":
                        sample_pos["saw_mixed"] = True
                        sample_pos["residue"] = None
                    elif state == "ALT" and alt_idx in allele_effects and pos0 in allele_effects[alt_idx]:
                        effect = allele_effects[alt_idx][pos0]
                        alt_aa = effect["alt_aa"]
                        if effect["ref_aa"] is not None:
                            sample_pos["ref_aa"] = effect["ref_aa"]
                            if pos0 not in pos_ref:
                                pos_ref[pos0] = effect["ref_aa"]
                        if effect["codon_idx"] is not None and effect["ref_nt"] is not None:
                            prev_ref_nt = sample_pos["ref_nts"].get(effect["codon_idx"])
                            if prev_ref_nt is not None and prev_ref_nt != effect["ref_nt"]:
                                sample_pos["conflict"] = True
                            sample_pos["ref_nts"][effect["codon_idx"]] = effect["ref_nt"]
                        if effect["codon_idx"] is not None and effect["alt_nt"] is not None:
                            prev_alt_nt = sample_pos["alt_nts"].get(effect["codon_idx"])
                            if prev_alt_nt is not None and prev_alt_nt != effect["alt_nt"]:
                                sample_pos["conflict"] = True
                            sample_pos["alt_nts"][effect["codon_idx"]] = effect["alt_nt"]
                        elif alt_aa is None:
                            sample_pos["saw_missing"] = True
                            sample_pos["residue"] = None

            for alt_idx, pos_map in allele_effects.items():
                for pos0, effect in pos_map.items():
                    variant_rows.append(
                        {
                            "analysis_mode": "aa_positions",
                            "chrom": effect["chrom"],
                            "pos": effect["pos"],
                            "ref": effect["ref"],
                            "alt": effect["alt"],
                            "aa_pos": pos0,
                            "ref_aa": effect["ref_aa"] or pos_ref.get(pos0),
                            "alt_aa": effect["alt_aa"],
                            "annotation": effect["annotation"],
                            "hgvs_c": effect["hgvs_c"],
                            "hgvs_p": effect["hgvs_p"],
                        }
                    )

        rows = []
        for sample_id in sample_ids:
            residues = []
            status = "CALLED"
            details = []
            for pos0 in target_positions:
                state_info = per_sample[sample_id][pos0]
                ref_aa = state_info["ref_aa"] or pos_ref.get(pos0)
                if state_info["conflict"] or state_info["saw_mixed"]:
                    state = "MIXED"
                    residue = None
                elif state_info["alt_nts"]:
                    residue = infer_alt_residue(ref_aa, state_info["ref_nts"], state_info["alt_nts"])
                    if residue is None:
                        state = "UNRESOLVED"
                    else:
                        state = "ALT"
                        state_info["residue"] = residue
                elif state_info["saw_ref"] or ref_aa is not None:
                    state = "REF"
                    residue = ref_aa
                    state_info["residue"] = residue
                elif state_info["saw_missing"]:
                    state = "MISSING"
                    residue = None
                else:
                    state = "UNRESOLVED"
                    residue = None

                if state == "MIXED":
                    status = "MIXED"
                    details.append(f"AA{pos0}=MIXED")
                elif state in ("MISSING", "UNRESOLVED"):
                    if status == "CALLED":
                        status = "UNRESOLVED"
                    if state == "MISSING":
                        details.append(f"AA{pos0}=MISSING")
                    else:
                        details.append(f"AA{pos0}=UNRESOLVED")
                else:
                    if residue is None:
                        if status == "CALLED":
                            status = "UNRESOLVED"
                        details.append(f"AA{pos0}=UNKNOWN")
                    else:
                        residues.append(f"{pos0}{residue}")
                        details.append(f"AA{pos0}={residue}")

            if status == "CALLED":
                residue_combo = "|".join(residues)
                reference_combo = []
                for pos0 in target_positions:
                    ref_aa = per_sample[sample_id][pos0]["ref_aa"] or pos_ref.get(pos0)
                    if ref_aa is None:
                        reference_combo = []
                        break
                    reference_combo.append(f"{pos0}{ref_aa}")
                assigned = "REF_HAPLOTYPE" if reference_combo and residue_combo == "|".join(reference_combo) else residue_combo
            else:
                assigned = status
            rows.append(
                {
                    "sample_id": sample_id,
                    "assigned_serotype": assigned,
                    "call_status": status,
                    "analysis_mode": "aa_positions",
                    "serotype_detail": ";".join(details),
                }
            )
        return rows, variant_rows


def all_variant_mode(vcf_path, gene, ann_field):
    with open_text(vcf_path) as handle:
        sample_ids = []
        sample_alts = None
        sample_missing = None
        sample_mixed = None
        variant_rows = []

        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_ids = header[9:]
                sample_alts = {sample_id: [] for sample_id in sample_ids}
                sample_missing = {sample_id: 0 for sample_id in sample_ids}
                sample_mixed = {sample_id: 0 for sample_id in sample_ids}
                continue

            fields = line.rstrip("\n").split("\t")
            chrom, pos, _vid, ref, alt, _qual, _flt, info_str, fmt = fields[:9]
            info = parse_info(info_str)
            alt_alleles = alt.split(",")
            ann_records = parse_effects_for_gene(
                info,
                ann_field,
                gene,
                alt_alleles,
                allow_ambiguous_multiallelic=True,
            )
            if not ann_records:
                continue
            if all(a == "*" for a in alt_alleles):
                continue

            format_fields = fmt.split(":")
            if "GT" not in format_fields:
                continue
            gt_index = format_fields.index("GT")
            site_id = f"{chrom}:{pos}:{ref}>{','.join(alt_alleles)}"

            variant_rows.append(
                {
                    "analysis_mode": "all_variants",
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": ",".join(alt_alleles),
                    "aa_pos": "|".join(
                        sorted(
                            {
                                str(rec["aa_pos"])
                                for rec in ann_records
                                if rec["aa_pos"] not in (None, "")
                            }
                        )
                    ),
                    "ref_aa": "|".join(
                        sorted(
                            {
                                rec["ref_aa"]
                                for rec in ann_records
                                if rec["ref_aa"] not in (None, "")
                            }
                        )
                    ),
                    "alt_aa": "|".join(
                        sorted(
                            {
                                rec["alt_aa"]
                                for rec in ann_records
                                if rec["alt_aa"] not in (None, "")
                            }
                        )
                    ),
                    "annotation": "|".join(sorted({rec["annotation"] for rec in ann_records if rec["annotation"]})),
                    "hgvs_c": "|".join(sorted({rec["hgvs_c"] for rec in ann_records if rec["hgvs_c"]})),
                    "hgvs_p": "|".join(sorted({rec["hgvs_p"] for rec in ann_records if rec["hgvs_p"]})),
                }
            )

            for sample_id, sample_field in zip(sample_ids, fields[9:]):
                state, alt_idx = classify_genotype(parse_gt(sample_field, gt_index))
                if state == "MISSING":
                    sample_missing[sample_id] += 1
                elif state == "MIXED":
                    sample_mixed[sample_id] += 1
                elif state == "ALT" and alt_idx is not None:
                    try:
                        alt_base = alt_alleles[int(alt_idx) - 1]
                    except (ValueError, IndexError):
                        alt_base = f"ALT{alt_idx}"
                    sample_alts[sample_id].append(f"{site_id}={alt_base}")

        rows = []
        for sample_id in sample_ids:
            if sample_missing[sample_id] > 0:
                status = "UNRESOLVED"
                assigned = "UNRESOLVED"
            elif sample_mixed[sample_id] > 0:
                status = "MIXED"
                assigned = "MIXED"
            else:
                status = "CALLED"
                assigned = ";".join(sample_alts[sample_id]) if sample_alts[sample_id] else "REF_HAPLOTYPE"

            rows.append(
                {
                    "sample_id": sample_id,
                    "assigned_serotype": assigned,
                    "call_status": status,
                    "analysis_mode": "all_variants",
                    "serotype_detail": (
                        f"n_alt_sites={len(sample_alts[sample_id])};"
                        f"n_missing_sites={sample_missing[sample_id]};"
                        f"n_mixed_sites={sample_mixed[sample_id]}"
                    ),
                }
            )
        return rows, variant_rows


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    if args.mode == "aa_positions":
        target_positions = parse_target_positions(args.aa_positions)
        reference_aa = parse_reference_aa(args.reference_aa)
        rows, variant_rows = aa_mode(args.vcf, args.gene, args.ann_field, target_positions, reference_aa)
    else:
        rows, variant_rows = all_variant_mode(args.vcf, args.gene, args.ann_field)

    write_tsv(
        args.out,
        rows,
        ["sample_id", "assigned_serotype", "call_status", "analysis_mode", "serotype_detail"],
    )
    if args.variant_summary:
        write_tsv(
            args.variant_summary,
            variant_rows,
            ["analysis_mode", "chrom", "pos", "ref", "alt", "aa_pos", "ref_aa", "alt_aa", "annotation", "hgvs_c", "hgvs_p"],
        )


if __name__ == "__main__":
    main()
