# WGS variant-calling pipeline

Four PBS-ready scripts implementing a standard trim -> align -> BQSR ->
HaplotypeCaller -> VQSR short-read variant-calling workflow (BWA-MEM +
GATK 4).

## Inputs

- Paired Illumina FASTQ files named `SAMPLE_R1.fastq.gz` and
  `SAMPLE_R2.fastq.gz`.
- Alignment and variant-calling reference FASTA files with their required
  indices.
- A BED file defining callable core regions.
- Known-site and VQSR training VCF resources appropriate for the reference.
- Trimmomatic and its adapter FASTA.

Input locations are supplied through `pipeline_config.sh` or environment
variables; no institution-specific paths are embedded.

## Setup

1. Edit `pipeline_config.sh`, or export the corresponding variables, to set
   reference, resource, FASTQ and working-directory paths.
2. Make BWA, SAMtools, Java, GATK, R, BCFtools, HTSlib and Trimmomatic
   available on `PATH` as required by each stage.
3. Supply allocation and storage settings when submitting. For example on
   NCI Gadi: `qsub -P PROJECT -l storage=gdata/PROJECT 01_trim_align.sh`.
   Adjust the queue and resource directives for other schedulers.

## Run order

| Script | Purpose | Input | Output |
|---|---|---|---|
| `01_trim_align.sh` | Trim adapters, align, mark duplicates, subset to core regions | Raw paired FASTQs | `<sample>.sorted.dup.core.bam` |
| `02_bqsr.sh` | Base quality score recalibration | Core BAMs from step 1 | `<sample>.recal_reads.bam` |
| `03_haplotype_caller.sh` | Per-sample, per-chromosome variant calling (GVCF mode) | Recalibrated BAMs from step 2 | `<sample>_<chrom>.g.vcf.gz` |
| `04_vqsr.sh` | Joint genotyping + variant quality score recalibration | GVCFs from step 3 | `recalibrated_variants_core.vcf.gz` |

Each script is resumable: it checks for existing output files and lock
files before (re)processing a sample, so a failed/killed job can simply be
resubmitted.

Note: `03_haplotype_caller.sh` contains an example chromosome-interval list
for a specific reference build — replace it with the intervals matching
your own reference genome.
