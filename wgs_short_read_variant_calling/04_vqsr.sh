#!/bin/bash
#PBS -l walltime=48:00:00
#PBS -l mem=100gb
#PBS -l jobfs=5gb
#PBS -l ncpus=1
#PBS -q normal
#PBS -N VQSR
#PBS -o VQSR.out
#PBS -e VQSR.err

# =============================================================================
# 04_vqsr.sh
#
# Step 4 of 4: combines per-sample, per-chromosome GVCFs from
# 03_haplotype_caller.sh, joint-genotypes them, restricts to core-genome
# regions, and runs Variant Quality Score Recalibration (VQSR) separately
# for SNPs and INDELs. Produces the final recalibrated, core-region VCF.
#
# Configure paths in pipeline_config.sh or with environment variables.
# =============================================================================

echo "=== Starting Variant Recalibration (core regions only) ==="

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/pipeline_config.sh"

mkdir -p "$VCF_DIR" "$COMBINED_GVCF_DIR"

# Required tools must be available on PATH: Java, GATK, bcftools and tabix.

# ---------- Step 0: Sample Names ----------
echo "Extracting sample names..."
# The pattern below strips the "_<chrom>_v3.g.vcf.gz" suffix added by
# 03_haplotype_caller.sh. Adjust it if you changed the chromosome/interval
# naming scheme in that script.
if [[ ! -f "$GVCF_DIR/sample_names.txt" ]]; then
    find "$GVCF_DIR" -type f -name "*.g.vcf.gz" \
        | sed -E 's|.*/||' \
        | sed -E 's/_[^_]+_[0-9]+_v3\.g\.vcf\.gz$//' \
        | sort -u > "$GVCF_DIR/sample_names.txt"
fi

# ---------- Step 1: Concatenate GVCFs by Sample ----------
echo "Concatenating GVCFs by sample..."
while read -r sample; do
    merged_gvcf="${COMBINED_GVCF_DIR}/${sample}.g.vcf.gz"
    if [[ ! -f "$merged_gvcf" ]]; then
        files=$(find "$GVCF_DIR" -type f -name "${sample}_*_v3.g.vcf.gz" | sort)
        bcftools concat -a -Oz -o "$merged_gvcf" $files
        tabix -p vcf "$merged_gvcf"
    fi
done < "$GVCF_DIR/sample_names.txt"

# ---------- Step 2: Combine All GVCFs ----------
echo "Combining all sample GVCFs..."

if [[ ! -f "$VCF_DIR/combined.g.vcf.gz" ]]; then
    echo "Checking for missing .tbi indexes..."
    for gvcf in "$COMBINED_GVCF_DIR"/*.g.vcf.gz; do
        if [[ ! -f "${gvcf}.tbi" ]]; then
            echo "Indexing $gvcf"
            gatk --java-options "-Xmx4G" IndexFeatureFile -I "$gvcf"
        fi
    done

    find "$COMBINED_GVCF_DIR" -name "*.g.vcf.gz" > "$COMBINED_GVCF_DIR/gvcf_list.list"
    gatk --java-options "-Xmx50G" CombineGVCFs \
        -R "$REF_DIR/$REFERENCE" \
        --variant "$COMBINED_GVCF_DIR/gvcf_list.list" \
        -O "$VCF_DIR/combined.g.vcf.gz"
fi

# ---------- Step 3: Genotype ----------
echo "Genotyping combined GVCFs..."
if [[ ! -f "$VCF_DIR/raw_variants.vcf.gz" ]]; then
    gatk --java-options "-Xmx50G" GenotypeGVCFs \
        -R "$REF_DIR/$REFERENCE" \
        -V "$VCF_DIR/combined.g.vcf.gz" \
        -O "$VCF_DIR/raw_variants.vcf.gz"
fi

# ---------- Step 4: Select Variants from Core Regions ----------
echo "Selecting SNPs and INDELs from core..."
if [[ ! -f "$VCF_DIR/raw_snps_core.vcf.gz" ]]; then
    gatk --java-options "-Xmx16G" SelectVariants \
        -R "$REF_DIR/$REFERENCE" \
        -V "$VCF_DIR/raw_variants.vcf.gz" \
        -L "$CORE_REGIONS" \
        --select-type-to-include SNP \
        -O "$VCF_DIR/raw_snps_core.vcf.gz"
fi

if [[ ! -f "$VCF_DIR/raw_indels_core.vcf.gz" ]]; then
    gatk --java-options "-Xmx16G" SelectVariants \
        -R "$REF_DIR/$REFERENCE" \
        -V "$VCF_DIR/raw_variants.vcf.gz" \
        -L "$CORE_REGIONS" \
        --select-type-to-include INDEL \
        -O "$VCF_DIR/raw_indels_core.vcf.gz"
fi

# ---------- Step 5: VQSR - SNPs ----------
echo "Running SNP recalibration..."
if [[ ! -f "$VCF_DIR/snps.recal" ]]; then
    gatk --java-options "-Xmx32G" VariantRecalibrator \
        -R "$REF_DIR/$REFERENCE" \
        -V "$VCF_DIR/raw_snps_core.vcf.gz" \
        --resource:truth_training_known,training=true,truth=true,known=true,prior=15.0 "$TRAINING_SET" \
        --resource:known_sites,known=true,prior=2.0 "$KNOWN_SITES1" \
        -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
        --trust-all-polymorphic \
        --max-gaussians 8 \
        --mode SNP \
        -O "$VCF_DIR/snps.recal" \
        --tranches-file "$VCF_DIR/snps.tranches"
fi

if [[ ! -f "$VCF_DIR/recalibrated_snps_core.vcf.gz" ]]; then
    gatk --java-options "-Xmx16G" ApplyVQSR \
        -R "$REF_DIR/$REFERENCE" \
        -V "$VCF_DIR/raw_snps_core.vcf.gz" \
        --recal-file "$VCF_DIR/snps.recal" \
        --tranches-file "$VCF_DIR/snps.tranches" \
        --truth-sensitivity-filter-level 99.0 \
        -mode SNP \
        -O "$VCF_DIR/recalibrated_snps_core.vcf.gz"
fi

# ---------- Step 6: VQSR - INDELs ----------
echo "Running INDEL recalibration..."
if [[ ! -f "$VCF_DIR/indels.recal" ]]; then
    gatk --java-options "-Xmx32G" VariantRecalibrator \
        -R "$REF_DIR/$REFERENCE" \
        -V "$VCF_DIR/raw_indels_core.vcf.gz" \
        --resource:truth_training_known,training=true,truth=true,known=true,prior=15.0 "$TRAINING_SET" \
        -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
        --trust-all-polymorphic \
        --max-gaussians 8 \
        --mode INDEL \
        -O "$VCF_DIR/indels.recal" \
        --tranches-file "$VCF_DIR/indels.tranches"
fi

if [[ ! -f "$VCF_DIR/recalibrated_indels_core.vcf.gz" ]]; then
    gatk --java-options "-Xmx16G" ApplyVQSR \
        -R "$REF_DIR/$REFERENCE" \
        -V "$VCF_DIR/raw_indels_core.vcf.gz" \
        --recal-file "$VCF_DIR/indels.recal" \
        --tranches-file "$VCF_DIR/indels.tranches" \
        --truth-sensitivity-filter-level 99.0 \
        -mode INDEL \
        -O "$VCF_DIR/recalibrated_indels_core.vcf.gz"
fi

# ---------- Step 7: Merge Final VCF ----------
echo "Merging SNPs and INDELs into final core-only VCF..."
if [[ ! -f "$VCF_DIR/recalibrated_variants_core.vcf.gz" ]]; then
    gatk --java-options "-Xmx16G" MergeVcfs \
        -I "$VCF_DIR/recalibrated_snps_core.vcf.gz" \
        -I "$VCF_DIR/recalibrated_indels_core.vcf.gz" \
        -O "$VCF_DIR/recalibrated_variants_core.vcf.gz"
fi

echo "All done! Output: $VCF_DIR/recalibrated_variants_core.vcf.gz"
