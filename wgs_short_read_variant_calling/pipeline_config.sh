#!/bin/bash
# =============================================================================
# pipeline_config.sh
#
# Central configuration for the WGS trim/align -> BQSR -> HaplotypeCaller ->
# VQSR variant-calling pipeline (01_trim_align.sh, 02_bqsr.sh,
# 03_haplotype_caller.sh, 04_vqsr.sh).
#
# Edit the defaults below or export variables with the same names before
# running. Environment values take precedence over these examples.
# =============================================================================

# ---- Alignment reference ----------------------------------------------------
# Used only by 01_trim_align.sh. This can be a combined/decoy reference
# (e.g. target organism + host genome) used to filter out contaminating
# reads before variant calling.
: "${ALIGN_REF_DIR:=/path/to/alignment_reference}"
: "${ALIGN_REFERENCE:=alignment_reference.fasta}"

# ---- Variant-calling reference ----------------------------------------------
# Used by 02_bqsr.sh, 03_haplotype_caller.sh, and 04_vqsr.sh.
: "${REF_DIR:=/path/to/reference}"
: "${REFERENCE:=reference_genome.fasta}"

# BED file defining the "core" / high-confidence callable regions of the
# reference genome (used to subset BAMs and restrict variant selection).
: "${CORE_REGIONS:=/path/to/core_regions.bed}"

# ---- Known-sites / training resources ---------------------------------------
# Used for BQSR (KNOWN_SITES1/2) and VQSR (TRAINING_SET, KNOWN_SITES1).
: "${KNOWN_SITES1:=/path/to/known_sites_1.vcf}"
: "${KNOWN_SITES2:=}"   # optional
: "${TRAINING_SET:=/path/to/vqsr_training_set.vcf.gz}"

# ---- Raw data ----------------------------------------------------------------
# Directory containing paired-end FASTQs named "<sample>_R1.fastq.gz" /
# "<sample>_R2.fastq.gz".
: "${RAW_FASTQ_DIR:=/path/to/raw_fastq}"

# ---- Trimmomatic --------------------------------------------------------------
: "${TRIMMOMATIC_JAR:=/path/to/trimmomatic-0.39.jar}"
: "${ADAPTER_FILE:=/path/to/adapters/NexteraPE-PE.fa}"

# ---- Working directory --------------------------------------------------------
# All intermediate and output files are written under here. Needs enough
# scratch space for BAMs/GVCFs for your sample set.
: "${WORK_DIR:=/path/to/work}"

# ---- Derived directories (created automatically; no need to edit) ------------
TRIMMED_DIR="${WORK_DIR}/trimmed_fastq"
BAM_DIR="${WORK_DIR}/bam"
RECAL_BAM_DIR="${WORK_DIR}/recal_bams"
INTERVAL_DIR="${WORK_DIR}/intervals"
GVCF_DIR="${WORK_DIR}/gvcf"
COMBINED_GVCF_DIR="${WORK_DIR}/gvcf_combined"
VCF_DIR="${WORK_DIR}/vcf"
LOG_DIR="${WORK_DIR}/logs"
LOCK_DIR="${WORK_DIR}/locks"
