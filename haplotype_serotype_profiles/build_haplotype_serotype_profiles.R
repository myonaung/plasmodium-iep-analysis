#!/usr/bin/env Rscript
# =============================================================================
# build_haplotype_serotype_profiles.R
#
# Builds an alluvial plot showing, per antigen, how "all-SNP" haplotypes
# (every polymorphic site) collapse down onto "immune-escape" haplotypes
# (only sites flagged as immune-escape candidates), for a chosen sample
# scope (a focal region by default, e.g. Papua New Guinea; or all of Asia;
# or global). Rare haplotypes on both sides are collapsed/merged so the
# alluvial stays readable.
#
# As a side effect, this script writes the two flow tables that
# downstream accumulation-curve analyses consume:
#   <output>_raw_flows.tsv   -> pass as --raw-flows
#   <output>_flows.tsv       -> pass as --alluvial-flows
# (The PDF/SVG alluvial plot itself is a separate, optional output — see
# --output/--output-svg below — and isn't needed just to produce those two
# TSVs.)
#
# Required inputs (via CLI flags; nothing is bundled)
# -----------------------------------------------------------------
# --pf-summary     Bayesian site-summary CSV for the first species.
# --pv-summary     Bayesian site-summary CSV for the second species.
# --pf-vcf-dir / --pv-vcf-dir     Directory of per-antigen VCFs (or use
#                  --pf-fasta-dir / --pv-fasta-dir for FASTA alignments
#                  instead — at least one source is required per species).
# --pf-metadata / --pv-metadata   Sample metadata table (required unless
#                  --pf-sample-scope/--pv-sample-scope is "global").
# --output         Output PDF path (required).
#
# Usage
# -----
# Rscript build_haplotype_serotype_profiles.R \
#   --pf-summary path/to/pf_site_summary.csv \
#   --pv-summary path/to/pv_site_summary.csv \
#   --pf-vcf-dir path/to/pf_vcf_dir \
#   --pv-vcf-dir path/to/pv_vcf_dir \
#   --pf-metadata path/to/pf_sample_metadata.tsv \
#   --pv-metadata path/to/pv_sample_metadata.tsv \
#   --output alluvial_plot.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(grid)
})

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
num <- function(x) suppressWarnings(as.numeric(x))

boolish <- function(x) {
  if (is.logical(x)) return(x)
  xn <- suppressWarnings(as.numeric(x))
  if (!all(is.na(xn))) return(xn == 1)
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

norm_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  gsub("[^a-z0-9]+", "", x)
}

split_arg <- function(x) {
  v <- trimws(unlist(strsplit(if (is.null(x) || identical(x, "")) "" else x, ",")))
  v[nzchar(v)]
}

pick_col <- function(df, cands) {
  h <- cands[cands %in% names(df)]
  if (length(h)) h[[1]] else NA_character_
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x)) || identical(x, "")) y else x
}

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parse_args <- function(args) {
  out <- list(
    pf_summary    = "",
    pv_summary    = "",
    pf_vcf_dir    = "",
    pv_vcf_dir    = "",
    pf_fasta_dir  = "",
    pv_fasta_dir  = "",
    fasta_filter_by_metadata = "FALSE",
    pf_metadata   = "",
    pv_metadata   = "",
    output              = NULL,
    output_svg          = NULL,
    summary_output      = NULL,
    flows_output        = NULL,
    raw_flows_output    = NULL,
    diagnostics_output  = NULL,
    pf_snp_col          = "pf7_ns_snp",
    pv_snp_col          = "Pv4_snps",
    metadata_sample_col = "",
    metadata_country_col= "",
    metadata_region_col = "",
    png_values          = "PNG,Papua New Guinea,Papua_New_Guinea",
    asia_values         = paste(
      c(
        "Asia", "South Asia", "Southeast Asia", "South-East Asia", "East Asia",
        "Afghanistan", "Bangladesh", "Bhutan", "Cambodia", "China", "India",
        "Indonesia", "Iran", "Laos",
        "Lao PDR", "Malaysia", "Myanmar", "Nepal", "Pakistan", "Philippines",
        "Sri Lanka", "Thailand", "Vietnam", "Viet Nam"
      ),
      collapse = ","
    ),
    sample_scope        = "png",
    pf_sample_scope     = "",
    pv_sample_scope     = "",
    pf_metadata_filter_col = "",
    pv_metadata_filter_col = "",
    pf_filter_values       = "",
    pv_filter_values       = "",
    prob_sym_gt_asym_threshold             = 0.80,
    prob_sym_gt_expected_threshold         = 0,
    require_prob_sym_gt_expected_over_asym = "TRUE",
    dominant_allele_freq_max               = 0.99,
    global_rare_filter                     = "TRUE",
    global_rare_filter_iep_fraction        = 0.50,
    global_rare_af_min                     = 0.10,
    global_rare_filter_class               = "significant",
    y_scale                       = "linear",
    min_png_samples               = 1,
    min_original_haplotype_percent = 1,
    min_reduced_haplotype_percent = 5,
    max_rare_original_haplotype_percent = 100,
    assign_rare_iep_to_nearest      = "TRUE",
    max_iep_hamming_distance        = 4,
    group_antigens_by_lifecycle     = "TRUE",
    top_n_per_species             = 5,
    force_pf_antigens             = "",
    force_pv_antigens             = "",
    exclude_pf_antigens           = "",
    exclude_pv_antigens           = "",
    facet_ncol                    = 5,
    # Both sides now use "H" prefix — "all-SNP haplotype" and
    # "immune-escape haplotype" are distinguished by axis label only
    left_prefix               = "H",
    right_prefix              = "H",
    ribbon_alpha              = 0.35,
    label_percent_threshold   = 0,
    label_min_height_pt       = 11,
    label_min_gap_pt          = 8,
    max_labels_per_side       = 8,
    max_left_labels_per_side  = 0,
    max_right_labels_per_side = 10,
    auto_max_labels_per_side  = 3,
    label_stagger_x           = 0,
    align_label_columns       = "TRUE",
    show_label_leaders        = "TRUE",
    show_filtered_iep_stub    = "TRUE",
    filtered_iep_stub_percent = 1,
    compact_narrow_labels     = "TRUE",
    use_ggrepel               = "TRUE",
    font_family               = "sans",
    # User-facing multiplier. Internally this is expanded below so
    # --font-scale 1 gives a manuscript-readable baseline.
    font_scale                = 1,
    pdf_width                 = 12,
    pdf_height                = 5
  )

  i <- 1L
  while (i <= length(args)) {
    if (!startsWith(args[[i]], "--") || i == length(args))
      stop("Bad arguments. Use --key value pairs.")
    k        <- gsub("-", "_", sub("^--", "", args[[i]]))
    out[[k]] <- args[[i + 1L]]
    i        <- i + 2L
  }

  if (is.null(out$output)) stop("Provide --output.")
  if (!nzchar(out$pf_summary)) stop("Provide --pf-summary (first species' Bayesian site-summary CSV).")
  if (!nzchar(out$pv_summary)) stop("Provide --pv-summary (second species' Bayesian site-summary CSV).")
  out$prob_sym_gt_asym_threshold             <- as.numeric(out$prob_sym_gt_asym_threshold)
  out$prob_sym_gt_expected_threshold         <- as.numeric(out$prob_sym_gt_expected_threshold)
  out$dominant_allele_freq_max               <- as.numeric(out$dominant_allele_freq_max)
  out$global_rare_filter                     <-
    tolower(out$global_rare_filter) %in% c("true", "t", "1", "yes", "y")
  out$global_rare_filter_iep_fraction        <- as.numeric(out$global_rare_filter_iep_fraction)
  if (!is.finite(out$global_rare_filter_iep_fraction)) out$global_rare_filter_iep_fraction <- 0.50
  out$global_rare_af_min                     <- as.numeric(out$global_rare_af_min)
  if (!is.finite(out$global_rare_af_min)) out$global_rare_af_min <- 0.10
  out$global_rare_filter_class               <- tolower(split_arg(out$global_rare_filter_class))
  out$min_png_samples                        <- as.integer(out$min_png_samples)
  out$min_original_haplotype_percent         <- as.numeric(out$min_original_haplotype_percent)
  out$min_reduced_haplotype_percent          <- as.numeric(out$min_reduced_haplotype_percent)
  out$max_rare_original_haplotype_percent    <- as.numeric(out$max_rare_original_haplotype_percent)
  if (!is.finite(out$max_rare_original_haplotype_percent) ||
      out$max_rare_original_haplotype_percent < 0 ||
      out$max_rare_original_haplotype_percent > 100) {
    stop("--max-rare-original-haplotype-percent must be between 0 and 100.")
  }
  out$assign_rare_iep_to_nearest <-
    tolower(out$assign_rare_iep_to_nearest) %in% c("true", "t", "1", "yes", "y")
  out$max_iep_hamming_distance <- as.integer(out$max_iep_hamming_distance)
  if (is.na(out$max_iep_hamming_distance) || out$max_iep_hamming_distance < 0) {
    stop("--max-iep-hamming-distance must be a non-negative integer.")
  }
  out$group_antigens_by_lifecycle <-
    tolower(out$group_antigens_by_lifecycle) %in% c("true", "t", "1", "yes", "y")
  out$top_n_per_species                      <- as.integer(out$top_n_per_species)
  out$force_pf_antigens                      <- split_arg(out$force_pf_antigens)
  out$force_pv_antigens                      <- split_arg(out$force_pv_antigens)
  out$exclude_pf_antigens                    <- split_arg(out$exclude_pf_antigens)
  out$exclude_pv_antigens                    <- split_arg(out$exclude_pv_antigens)
  out$facet_ncol                             <- as.integer(out$facet_ncol)
  out$ribbon_alpha                           <- as.numeric(out$ribbon_alpha)
  out$label_percent_threshold                <- as.numeric(out$label_percent_threshold)
  out$label_min_height_pt                    <- as.numeric(out$label_min_height_pt)
  out$label_min_gap_pt                       <- as.numeric(out$label_min_gap_pt)
  out$max_labels_per_side                    <- as.integer(out$max_labels_per_side)
  out$max_left_labels_per_side               <- as.integer(out$max_left_labels_per_side)
  out$max_right_labels_per_side              <- as.integer(out$max_right_labels_per_side)
  out$auto_max_labels_per_side               <- as.integer(out$auto_max_labels_per_side)
  if (is.na(out$auto_max_labels_per_side) || out$auto_max_labels_per_side < 1L) {
    out$auto_max_labels_per_side <- 3L
  }
  out$label_stagger_x                        <- as.numeric(out$label_stagger_x)
  out$align_label_columns                    <-
    tolower(out$align_label_columns) %in% c("true", "t", "1", "yes", "y")
  out$show_label_leaders                     <-
    tolower(out$show_label_leaders) %in% c("true", "t", "1", "yes", "y")
  out$show_filtered_iep_stub                 <-
    tolower(out$show_filtered_iep_stub) %in% c("true", "t", "1", "yes", "y")
  out$filtered_iep_stub_percent              <- as.numeric(out$filtered_iep_stub_percent)
  out$compact_narrow_labels                  <-
    tolower(out$compact_narrow_labels) %in% c("true", "t", "1", "yes", "y")
  out$use_ggrepel                            <-
    tolower(out$use_ggrepel) %in% c("true", "t", "1", "yes", "y")
  out$font_family <- as.character(out$font_family %||% "sans")
  if (!nzchar(out$font_family)) out$font_family <- "sans"
  out$font_scale <- as.numeric(out$font_scale)
  if (!is.finite(out$font_scale) || out$font_scale <= 0) out$font_scale <- 1
  out$font_scale <- out$font_scale * 2.5
  out$pdf_width                              <- as.numeric(out$pdf_width)
  out$pdf_height                             <- as.numeric(out$pdf_height)
  out$require_prob_sym_gt_expected_over_asym <-
    tolower(out$require_prob_sym_gt_expected_over_asym) %in% c("true", "t", "1", "yes", "y")
  out$fasta_filter_by_metadata <-
    tolower(out$fasta_filter_by_metadata) %in% c("true", "t", "1", "yes", "y")
  out$y_scale <- match.arg(out$y_scale, c("linear", "log10"))
  out$sample_scope <- match.arg(tolower(out$sample_scope), c("png", "asia", "global"))
  out$pf_sample_scope <- if (is.null(out$pf_sample_scope) || identical(out$pf_sample_scope, "")) {
    out$sample_scope
  } else {
    match.arg(tolower(out$pf_sample_scope), c("png", "asia", "global"))
  }
  out$pv_sample_scope <- if (is.null(out$pv_sample_scope) || identical(out$pv_sample_scope, "")) {
    out$sample_scope
  } else {
    match.arg(tolower(out$pv_sample_scope), c("png", "asia", "global"))
  }
  if (out$pf_sample_scope != "global" && !nzchar(out$pf_metadata)) {
    stop("Provide --pf-metadata, or run with --pf-sample-scope global.")
  }
  if (out$pv_sample_scope != "global" && !nzchar(out$pv_metadata)) {
    stop("Provide --pv-metadata, or run with --pv-sample-scope global.")
  }
  out
}

# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------
read_table_auto <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE))
    read_csv(path, show_col_types = FALSE, guess_max = 300000)
  else
    read_tsv(path, show_col_types = FALSE, guess_max = 300000)
}

metadata_samples_for_scope <- function(path, args, species_key, scope) {
  if (scope == "global") return(NULL)

  meta       <- read_table_auto(path)
  sample_col <- pick_col(meta, c(
    args$metadata_sample_col,
    "sample_id", "sample", "Sample", "Sample_ID", "sample_name",
    "SampleName", "id", "ID", "ox_code", "Sample_Name"
  ))
  if (is.na(sample_col)) stop("Could not find a sample ID column in metadata: ", path)

  filter_col_arg <- if (species_key == "pf") args$pf_metadata_filter_col else args$pv_metadata_filter_col
  filter_values_arg <- if (species_key == "pf") args$pf_filter_values else args$pv_filter_values

  if (scope == "png") {
    values <- if (nzchar(filter_values_arg)) filter_values_arg else args$png_values
    filter_col <- pick_col(meta, c(
      filter_col_arg,
      args$metadata_country_col,
      "country", "Country", "COUNTRY", "admin0", "Admin0",
      "country_name", "Country_name", "location", "Location"
    ))
    if (is.na(filter_col)) stop("Could not find a country/filter column in metadata: ", path)
  } else if (scope == "asia") {
    values <- if (nzchar(filter_values_arg)) filter_values_arg else args$asia_values
    filter_col <- pick_col(meta, c(
      filter_col_arg,
      args$metadata_region_col,
      "region", "Region", "REGION",
      "continent", "Continent", "CONTINENT",
      "subregion", "Subregion", "sub_region", "Sub_region",
      "geographic_region", "Geographic_region", "GeographicRegion",
      "population", "Population",
      args$metadata_country_col,
      "country", "Country", "COUNTRY", "admin0", "Admin0",
      "country_name", "Country_name", "location", "Location"
    ))
    if (is.na(filter_col)) stop("Could not find a region/country/filter column for Asia samples in metadata: ", path)
  } else {
    stop("Unknown sample scope: ", scope)
  }

  keep    <- norm_key(meta[[filter_col]]) %in% norm_key(split_arg(values))
  samples <- unique(as.character(meta[[sample_col]][keep]))
  samples <- samples[!is.na(samples) & nzchar(samples)]
  if (!length(samples)) {
    stop(
      "No ", toupper(scope), " samples found in metadata: ", path,
      ". Checked column '", filter_col, "' for values: ", values
    )
  }
  samples
}

# ---------------------------------------------------------------------------
# Escape classification
# ---------------------------------------------------------------------------
compute_escape <- function(df, snp_col, args) {
  if (!snp_col %in% names(df)) stop("SNP column not found: ", snp_col)
  psa     <- num(df$prob_sym_gt_asym)
  pse     <- num(df$prob_sym_gt_expected)
  pae     <- num(df$prob_asym_gt_expected)
  modeled <- if ("has_bayesian_summary" %in% names(df)) {
    boolish(df$has_bayesian_summary)
  } else {
    !is.na(psa) | !is.na(pse)
  }

  esc <- modeled &
    !is.na(psa) & psa >= args$prob_sym_gt_asym_threshold &
    !is.na(pse) & pse >= args$prob_sym_gt_expected_threshold &
    !is.na(num(df[[snp_col]])) & num(df[[snp_col]]) == 1

  if (isTRUE(args$require_prob_sym_gt_expected_over_asym))
    esc <- esc & !is.na(pae) & pse > pae

  if ("dominant_allele_freq" %in% names(df)) {
    daf <- num(df$dominant_allele_freq)
    esc <- esc & (is.na(daf) | daf <= args$dominant_allele_freq_max)
  }
  esc[is.na(esc)] <- FALSE
  esc
}

read_escape_positions <- function(path, snp_col, args) {
  df <- read_csv(path, show_col_types = FALSE, guess_max = 300000)
  if (!"position"         %in% names(df)) stop("Summary missing position column: ",         path)
  if (!"prob_sym_gt_asym" %in% names(df)) stop("Summary missing prob_sym_gt_asym column: ", path)

  antigen_col <- pick_col(df, c("antigen", "gene_name", "Gene_Name_annotated", "gene_id"))
  gene_col    <- pick_col(df, c("gene_id", "Gene_ID",   "Gene_Name_annotated", "gene_name"))
  if (is.na(antigen_col)) stop("Summary has no antigen/gene_name column: ", path)
  if (is.na(gene_col))    stop("Summary has no gene_id/gene_name column: ", path)

  df$immune_escape <- compute_escape(df, snp_col, args)
  df %>%
    transmute(
      gene_id      = as.character(.data[[gene_col]]),
      antigen      = as.character(.data[[antigen_col]]),
      gene_key     = norm_key(.data[[gene_col]]),
      antigen_key  = norm_key(.data[[antigen_col]]),
      aa_position  = as.integer(num(position)),
      immune_escape = immune_escape
    ) %>%
    filter(immune_escape, !is.na(aa_position)) %>%
    distinct(gene_key, antigen_key, aa_position, .keep_all = TRUE)
}

read_reference_aa_table <- function(path) {
  df <- read_csv(path, show_col_types = FALSE, guess_max = 300000)
  if (!"position" %in% names(df)) stop("Summary missing position column: ", path)

  antigen_col <- pick_col(df, c("antigen", "gene_name", "Gene_Name_annotated", "gene_id"))
  gene_col    <- pick_col(df, c("gene_id", "Gene_ID", "Gene_Name_annotated", "gene_name"))
  aa_col      <- pick_col(df, c(
    "ref_aa", "ref_aa_variant", "aminoacid", "AA",
    "reference_aa", "RefAA", "Reference_AA", "ref_amino_acid"
  ))
  if (is.na(antigen_col)) stop("Summary has no antigen/gene_name column: ", path)
  if (is.na(gene_col))    stop("Summary has no gene_id/gene_name column: ", path)
  if (is.na(aa_col)) {
    return(tibble(
      gene_id = character(),
      antigen = character(),
      gene_key = character(),
      antigen_key = character(),
      position = integer(),
      ref_aa = character()
    ))
  }

  df %>%
    transmute(
      gene_id     = as.character(.data[[gene_col]]),
      antigen     = as.character(.data[[antigen_col]]),
      gene_key    = norm_key(.data[[gene_col]]),
      antigen_key = norm_key(.data[[antigen_col]]),
      position    = as.integer(num(position)),
      ref_aa      = toupper(trimws(as.character(.data[[aa_col]])))
    ) %>%
    mutate(ref_aa = if_else(grepl("^[A-Z*]$", ref_aa), ref_aa, NA_character_)) %>%
    filter(!is.na(position), !is.na(ref_aa), nzchar(ref_aa)) %>%
    distinct(gene_key, antigen_key, position, .keep_all = TRUE)
}

read_escape_rank_table <- function(path, snp_col, args) {
  df <- read_csv(path, show_col_types = FALSE, guess_max = 300000)
  if (!"position" %in% names(df)) stop("Summary missing position column: ", path)
  if (!snp_col %in% names(df)) stop("SNP column not found: ", snp_col)

  antigen_col <- pick_col(df, c("antigen", "gene_name", "Gene_Name_annotated", "gene_id"))
  gene_col    <- pick_col(df, c("gene_id", "Gene_ID", "Gene_Name_annotated", "gene_name"))
  if (is.na(antigen_col)) stop("Summary has no antigen/gene_name column: ", path)
  if (is.na(gene_col))    stop("Summary has no gene_id/gene_name column: ", path)

  psa <- num(df$prob_sym_gt_asym)
  pse <- num(df$prob_sym_gt_expected)
  pae <- num(df$prob_asym_gt_expected)
  snp_ok <- !is.na(num(df[[snp_col]])) & num(df[[snp_col]]) == 1
  modeled <- if ("has_bayesian_summary" %in% names(df)) {
    boolish(df$has_bayesian_summary)
  } else {
    !is.na(psa) | !is.na(pse)
  }
  daf_ok <- rep(TRUE, nrow(df))
  if ("dominant_allele_freq" %in% names(df)) {
    daf <- num(df$dominant_allele_freq)
    daf_ok <- is.na(daf) | daf <= args$dominant_allele_freq_max
  }

  significant <- compute_escape(df, snp_col, args)
  potential <- modeled &
    snp_ok &
    daf_ok &
    !is.na(pse) & !is.na(pae) & pse > pae &
    !significant
  potential[is.na(potential)] <- FALSE

  df %>%
    mutate(
      gene_id = as.character(.data[[gene_col]]),
      antigen = as.character(.data[[antigen_col]]),
      gene_key = norm_key(.data[[gene_col]]),
      antigen_key = norm_key(.data[[antigen_col]]),
      aa_position = as.integer(num(position)),
      snp_ok = snp_ok,
      significant_iep = significant,
      potential_iep = potential
    ) %>%
    filter(!is.na(aa_position)) %>%
    distinct(gene_key, antigen_key, aa_position, .keep_all = TRUE) %>%
    group_by(gene_id, antigen, gene_key, antigen_key) %>%
    summarise(
      n_total_sites = n_distinct(aa_position),
      n_snp_sites = sum(snp_ok, na.rm = TRUE),
      n_significant_iep_sites = sum(significant_iep, na.rm = TRUE),
      n_potential_iep_sites = sum(potential_iep, na.rm = TRUE),
      weighted_iep_score = n_significant_iep_sites + 0.5 * n_potential_iep_sites,
      weighted_iep_density = if_else(n_snp_sites > 0, weighted_iep_score / n_snp_sites, 0),
      .groups = "drop"
    )
}

read_escape_site_table <- function(path, snp_col, args) {
  df <- read_csv(path, show_col_types = FALSE, guess_max = 300000)
  if (!"position" %in% names(df)) stop("Summary missing position column: ", path)
  if (!snp_col %in% names(df)) stop("SNP column not found: ", snp_col)

  antigen_col <- pick_col(df, c("antigen", "gene_name", "Gene_Name_annotated", "gene_id"))
  gene_col    <- pick_col(df, c("gene_id", "Gene_ID", "Gene_Name_annotated", "gene_name"))
  alt_col     <- pick_col(df, c("alt_aa", "alternative_aa", "variant_alt_aa", "AltAA", "mut_aa"))
  aa_change_col <- pick_col(df, c(
    "hgvs_p", "HGVS_P", "aa_change", "AA_change", "amino_acid_change",
    "variant_aa_change", "protein_change"
  ))
  if (is.na(antigen_col)) stop("Summary has no antigen/gene_name column: ", path)
  if (is.na(gene_col))    stop("Summary has no gene_id/gene_name column: ", path)

  psa <- num(df$prob_sym_gt_asym)
  pse <- num(df$prob_sym_gt_expected)
  pae <- num(df$prob_asym_gt_expected)
  snp_ok <- !is.na(num(df[[snp_col]])) & num(df[[snp_col]]) == 1
  modeled <- if ("has_bayesian_summary" %in% names(df)) {
    boolish(df$has_bayesian_summary)
  } else {
    !is.na(psa) | !is.na(pse)
  }
  daf_ok <- rep(TRUE, nrow(df))
  if ("dominant_allele_freq" %in% names(df)) {
    daf <- num(df$dominant_allele_freq)
    daf_ok <- is.na(daf) | daf <= args$dominant_allele_freq_max
  }

  significant <- compute_escape(df, snp_col, args)
  potential <- modeled &
    snp_ok &
    daf_ok &
    !is.na(pse) & !is.na(pae) & pse > pae &
    !significant
  potential[is.na(potential)] <- FALSE

  df %>%
    transmute(
      gene_id = as.character(.data[[gene_col]]),
      antigen = as.character(.data[[antigen_col]]),
      gene_key = norm_key(.data[[gene_col]]),
      antigen_key = norm_key(.data[[antigen_col]]),
      aa_position = as.integer(num(position)),
      alt_aa = if (!is.na(alt_col)) toupper(trimws(as.character(.data[[alt_col]]))) else NA_character_,
      aa_change_text = if (!is.na(aa_change_col)) as.character(.data[[aa_change_col]]) else NA_character_,
      snp_ok = snp_ok,
      significant_iep = significant,
      potential_iep = potential,
      global_vcf_af = NA_real_,
      global_rare_filtered = FALSE
    ) %>%
    mutate(
      alt_aa_from_change = vapply(aa_change_text, function(txt) {
        parsed <- parse_aa_change(txt)
        if (is.null(parsed) || !nrow(parsed)) NA_character_ else parsed$alt_aa[[1]]
      }, character(1)),
      alt_aa = if_else(grepl("^[A-Z*]$", alt_aa), alt_aa, alt_aa_from_change),
      alt_aa = if_else(grepl("^[A-Z*]$", alt_aa), alt_aa, NA_character_)
    ) %>%
    select(-aa_change_text, -alt_aa_from_change) %>%
    filter(!is.na(aa_position))
}

summarise_escape_rank <- function(site_table) {
  if (!nrow(site_table)) {
    return(tibble(
      gene_id = character(),
      antigen = character(),
      gene_key = character(),
      antigen_key = character(),
      n_total_sites = integer(),
      n_snp_sites = integer(),
      n_significant_iep_sites = integer(),
      n_potential_iep_sites = integer(),
      n_global_rare_filtered_iep_sites = integer(),
      weighted_iep_score = numeric(),
      weighted_iep_density = numeric()
    ))
  }

  site_table %>%
    group_by(gene_id, antigen, gene_key, antigen_key, aa_position) %>%
    summarise(
      snp_ok = any(snp_ok, na.rm = TRUE),
      significant_iep = any(significant_iep, na.rm = TRUE),
      potential_iep = any(potential_iep, na.rm = TRUE),
      global_rare_filtered = any(global_rare_filtered, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(gene_id, antigen, gene_key, antigen_key) %>%
    summarise(
      n_total_sites = n_distinct(aa_position),
      n_snp_sites = sum(snp_ok, na.rm = TRUE),
      n_significant_iep_sites = sum(significant_iep, na.rm = TRUE),
      n_potential_iep_sites = sum(potential_iep & !significant_iep, na.rm = TRUE),
      n_global_rare_filtered_iep_sites = sum(global_rare_filtered, na.rm = TRUE),
      weighted_iep_score = n_significant_iep_sites + 0.5 * n_potential_iep_sites,
      weighted_iep_density = if_else(n_snp_sites > 0, weighted_iep_score / n_snp_sites, 0),
      .groups = "drop"
    )
}

escape_positions_from_sites <- function(site_table) {
  site_table %>%
    filter(significant_iep, !is.na(aa_position)) %>%
    distinct(gene_key, antigen_key, aa_position, .keep_all = TRUE) %>%
    transmute(
      gene_id,
      antigen,
      gene_key,
      antigen_key,
      aa_position,
      immune_escape = TRUE
    )
}

# ---------------------------------------------------------------------------
# VCF parsing
# ---------------------------------------------------------------------------
parse_gene_antigen_from_path <- function(path) {
  base <- basename(path)
  base <- sub("\\.(vcf|fa|fasta|faa|fas|fna|aln|afa)(\\.gz)?$", "", base, ignore.case = TRUE)
  m    <- regexec("^([A-Za-z0-9]+_[0-9]+)_(.+)$", base)
  hit  <- regmatches(base, m)[[1]]
  if (length(hit) >= 3) list(gene_id = hit[[2]], antigen = hit[[3]])
  else                  list(gene_id = NA_character_, antigen = base)
}

parse_gene_antigen_from_vcf <- function(path) {
  parse_gene_antigen_from_path(path)
}

positions_for_vcf <- function(escape_positions, gene_id, antigen) {
  g <- norm_key(gene_id)
  a <- norm_key(antigen)
  escape_positions %>%
    filter(gene_key == g | antigen_key == a) %>%
    pull(aa_position) %>%
    unique()
}

infer_gene_antigen_from_path <- function(path, escape_positions) {
  parsed <- parse_gene_antigen_from_path(path)
  parsed_positions <- positions_for_vcf(escape_positions, parsed$gene_id, parsed$antigen)
  if (length(parsed_positions)) {
    g <- norm_key(parsed$gene_id)
    a <- norm_key(parsed$antigen)
    canonical <- escape_positions %>%
      distinct(gene_id, antigen, gene_key, antigen_key) %>%
      filter(gene_key == g | antigen_key == a) %>%
      arrange(desc(gene_key == g), desc(antigen_key == a), antigen)
    if (nrow(canonical)) {
      return(list(gene_id = canonical$gene_id[[1]], antigen = canonical$antigen[[1]]))
    }
    return(parsed)
  }

  path_key <- norm_key(paste(dirname(path), basename(path), sep = "/"))
  candidates <- escape_positions %>%
    distinct(gene_id, antigen, gene_key, antigen_key) %>%
    mutate(
      gene_match = vapply(gene_key, function(x) {
        !is.na(x) && nzchar(x) && grepl(x, path_key, fixed = TRUE)
      }, logical(1)),
      antigen_match = vapply(antigen_key, function(x) {
        !is.na(x) && nzchar(x) && grepl(x, path_key, fixed = TRUE)
      }, logical(1))
    ) %>%
    filter(gene_match | antigen_match) %>%
    mutate(match_len = pmax(nchar(gene_key), nchar(antigen_key), na.rm = TRUE)) %>%
    arrange(desc(match_len), antigen)

  if (nrow(candidates)) {
    list(gene_id = candidates$gene_id[[1]], antigen = candidates$antigen[[1]])
  } else {
    parsed
  }
}

extract_info_value <- function(info, key) {
  m <- regmatches(info, regexpr(paste0("(^|;)", key, "=[^;]+"), info))
  if (!length(m) || identical(m, character(0))) return(NA_character_)
  sub(paste0("^;?", key, "="), "", m)
}

aa_positions_from_text <- function(text) {
  if (is.na(text) || !nzchar(text) || text == ".") return(integer())
  text <- URLdecode(text)
  out  <- integer()
  patterns <- c(
    "p\\.[A-Z][a-z]{2}([0-9]+)([A-Z][a-z]{2}|=|Ter|\\*)",
    "p\\.([A-Z\\*])([0-9]+)([A-Z\\*=])",
    "(^|[^A-Za-z])([A-Z\\*])([0-9]+)([A-Z\\*=])([^A-Za-z]|$)",
    "(^|[^A-Za-z])([A-Z][a-z]{2})([0-9]+)([A-Z][a-z]{2}|=|Ter|\\*)([^A-Za-z]|$)"
  )
  for (pat in patterns) {
    hits <- gregexpr(pat, text, perl = TRUE)[[1]]
    if (hits[[1]] == -1) next
    vals <- regmatches(text, list(hits))[[1]]
    for (val in vals) {
      nums <- regmatches(val, gregexpr("[0-9]+", val, perl = TRUE))[[1]]
      if (length(nums)) out <- c(out, as.integer(nums[[1]]))
    }
  }
  unique(out[!is.na(out)])
}

ann_aa_positions <- function(info) {
  legacy <- c(
    extract_info_value(info, "SNPEFF_AMINO_ACID_CHANGE"),
    extract_info_value(info, "SNPEFF_HGVS_P"),
    extract_info_value(info, "HGVS_P")
  )
  legacy_pos <- unique(unlist(lapply(legacy, aa_positions_from_text)))
  m <- regmatches(info, regexpr("ANN=[^;]+", info))
  if (!length(m) || identical(m, character(0))) return(legacy_pos)
  ann     <- sub("^ANN=", "", m)
  records <- strsplit(ann, ",", fixed = TRUE)[[1]]
  pos     <- integer()
  for (rec in records) {
    fields <- strsplit(rec, "|", fixed = TRUE)[[1]]
    if (length(fields) >= 14 && nzchar(fields[[14]])) {
      p <- suppressWarnings(as.integer(sub("/.*$", "", fields[[14]])))
      if (!is.na(p)) pos <- c(pos, p)
    }
    if (length(fields) >= 11 && nzchar(fields[[11]]))
      pos <- c(pos, aa_positions_from_text(fields[[11]]))
  }
  unique(c(pos, legacy_pos))
}

aa3_to_1 <- function(x) {
  map <- c(
    Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C",
    Gln = "Q", Glu = "E", Gly = "G", His = "H", Ile = "I",
    Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P",
    Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V",
    Ter = "*", Sec = "U", Pyl = "O"
  )
  x <- as.character(x)
  ifelse(nchar(x) == 1, x, unname(map[x]))
}

parse_aa_change <- function(text) {
  if (is.na(text) || !nzchar(text) || text == ".") return(NULL)
  text <- URLdecode(text)

  m3 <- regexec("p\\.([A-Z][a-z]{2})([0-9]+)([A-Z][a-z]{2}|=|Ter|\\*)", text, perl = TRUE)
  h3 <- regmatches(text, m3)[[1]]
  if (length(h3) == 4) {
    ref <- aa3_to_1(h3[[2]])
    alt <- if (h3[[4]] == "=") ref else aa3_to_1(h3[[4]])
    if (!is.na(ref) && !is.na(alt)) return(tibble(pos = as.integer(h3[[3]]), ref_aa = ref, alt_aa = alt))
  }

  m1 <- regexec("p\\.([A-Z\\*])([0-9]+)([A-Z\\*=])", text, perl = TRUE)
  h1 <- regmatches(text, m1)[[1]]
  if (length(h1) == 4) {
    ref <- h1[[2]]
    alt <- if (h1[[4]] == "=") ref else h1[[4]]
    return(tibble(pos = as.integer(h1[[3]]), ref_aa = ref, alt_aa = alt))
  }

  m1_plain <- regexec("(^|[^A-Za-z])([A-Z\\*])([0-9]+)([A-Z\\*=])([^A-Za-z]|$)", text, perl = TRUE)
  h1_plain <- regmatches(text, m1_plain)[[1]]
  if (length(h1_plain) == 6) {
    ref <- h1_plain[[3]]
    alt <- if (h1_plain[[5]] == "=") ref else h1_plain[[5]]
    return(tibble(pos = as.integer(h1_plain[[4]]), ref_aa = ref, alt_aa = alt))
  }

  NULL
}

aa_changes_from_info <- function(info, alt_field) {
  alt_alleles <- strsplit(alt_field, ",", fixed = TRUE)[[1]]
  out <- list()

  m <- regmatches(info, regexpr("ANN=[^;]+", info))
  if (length(m) && !identical(m, character(0))) {
    ann <- sub("^ANN=", "", m)
    records <- strsplit(ann, ",", fixed = TRUE)[[1]]
    for (rec in records) {
      fields <- strsplit(rec, "|", fixed = TRUE)[[1]]
      allele <- if (length(fields) >= 1) fields[[1]] else NA_character_
      aa_txt <- if (length(fields) >= 11) fields[[11]] else NA_character_
      parsed <- parse_aa_change(aa_txt)
      if (is.null(parsed) || !nrow(parsed)) next
      allele_index <- match(allele, alt_alleles)
      if (is.na(allele_index)) next
      out[[length(out) + 1L]] <- parsed %>%
        mutate(allele_index = allele_index)
    }
  }

  if (!length(out)) {
    legacy <- c(
      extract_info_value(info, "SNPEFF_AMINO_ACID_CHANGE"),
      extract_info_value(info, "SNPEFF_HGVS_P"),
      extract_info_value(info, "HGVS_P")
    )
    for (txt in legacy) {
      parsed <- parse_aa_change(txt)
      if (!is.null(parsed) && nrow(parsed)) {
        out[[length(out) + 1L]] <- parsed %>% mutate(allele_index = 1L)
      }
    }
  }

  if (!length(out)) return(tibble())
  bind_rows(out) %>%
    filter(!is.na(pos), !is.na(ref_aa), !is.na(alt_aa), ref_aa != alt_aa) %>%
    distinct(pos, ref_aa, alt_aa, allele_index, .keep_all = TRUE)
}

normalize_gt <- function(x, gt_index) {
  vals <- vapply(strsplit(x, ":", fixed = TRUE), function(z) {
    if (length(z) >= gt_index) z[[gt_index]] else "."
  }, character(1))
  vals <- sub("\\|", "/", vals)
  vals[is.na(vals) | vals == "" | vals == "."] <- "./."
  vals
}

gt_allele_counts <- function(gt) {
  alleles <- unlist(strsplit(gt, "[/|]", perl = TRUE), use.names = FALSE)
  alleles <- alleles[nzchar(alleles) & alleles != "."]
  alleles <- suppressWarnings(as.integer(alleles))
  alleles <- alleles[!is.na(alleles)]
  list(total = length(alleles), counts = table(alleles))
}

gt_to_aa_state <- function(gt, ref_aa, alt_by_index) {
  vapply(gt, function(g) {
    if (is.na(g) || g %in% c(".", "./.", ".|.")) return(".")
    alleles <- unlist(strsplit(g, "[/|]", perl = TRUE))
    alleles <- alleles[nzchar(alleles) & alleles != "."]
    if (!length(alleles)) return(".")
    states <- vapply(alleles, function(a) {
      idx <- suppressWarnings(as.integer(a))
      if (is.na(idx)) return(".")
      if (idx == 0L) return(ref_aa)
      hit <- alt_by_index[[as.character(idx)]]
      if (is.null(hit) || is.na(hit) || !nzchar(hit)) "." else hit
    }, character(1))
    states <- unique(states[states != "."])
    if (!length(states)) "." else paste(sort(states), collapse = "/")
  }, character(1))
}

combine_position_states <- function(state_vectors, ref_aa = NA_character_) {
  mat <- do.call(cbind, state_vectors)
  apply(mat, 1, function(x) {
    vals <- unique(x[!is.na(x) & x != "."])
    if (!is.na(ref_aa) && nzchar(ref_aa) && any(vals != ref_aa)) {
      # Multiple nucleotide variants can map to the same AA coordinate.
      # A reference call at one record does not mean the sample is reference
      # when another record at the same AA position carries an alternate AA.
      vals <- vals[vals != ref_aa]
    }
    if (!length(vals)) "." else paste(sort(vals), collapse = "/")
  })
}

aa_state_variable_positions <- function(state_matrix) {
  keep <- vapply(seq_len(ncol(state_matrix)), function(j) {
    vals <- unique(state_matrix[, j])
    vals <- vals[!is.na(vals) & vals != "."]
    length(vals) > 1
  }, logical(1))
  which(keep)
}

state_matrix_to_long <- function(sample_ids, state_matrix) {
  if (is.null(state_matrix) || !length(state_matrix) ||
      !nrow(state_matrix) || !ncol(state_matrix)) {
    return(tibble(sample_id = character(), position = integer(), aa_state = character()))
  }
  rownames(state_matrix) <- seq_len(nrow(state_matrix))
  out <- as_tibble(as.table(state_matrix), .name_repair = "minimal")
  names(out) <- c("sample_i", "position", "aa_state")
  out %>%
    mutate(
      sample_i = as.integer(as.character(sample_i)),
      sample_id = sample_ids[sample_i],
      position = as.integer(as.character(position)),
      aa_state = as.character(aa_state)
    ) %>%
    filter(!is.na(position), !is.na(aa_state), aa_state != ".") %>%
    select(sample_id, position, aa_state)
}

vcf_haplotypes <- function(path, selected_samples, iep_positions) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  lines       <- readLines(con, warn = FALSE)
  header_line <- lines[startsWith(lines, "#CHROM")]
  if (!length(header_line)) stop("VCF has no #CHROM header: ", path)
  header      <- strsplit(header_line[[1]], "\t", fixed = TRUE)[[1]]
  if (length(header) < 10) stop("VCF has no sample columns: ", path)
  vcf_samples <- header[10:length(header)]
  keep_idx    <- if (is.null(selected_samples)) seq_along(vcf_samples) else which(vcf_samples %in% selected_samples)
  if (!length(keep_idx)) {
    return(list(sample_states = tibble(), sample_ids = character(),
                n_png_samples = 0,
                n_selected_samples = 0,
                n_all_snp_records = 0, n_iep_snp_records = 0))
  }
  kept_samples <- vcf_samples[keep_idx]

  record_lines  <- lines[!startsWith(lines, "#")]
  if (!length(record_lines)) {
    return(list(sample_states = tibble(), sample_ids = kept_samples,
                n_png_samples = length(keep_idx),
                n_selected_samples = length(keep_idx),
                n_all_snp_records = 0, n_iep_snp_records = 0))
  }

  states_by_pos <- list()
  ref_by_pos <- list()

  for (line in record_lines) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 10) next
    ref <- fields[[4]]; alt <- fields[[5]]
    if (nchar(ref) != 1 || !grepl("^[ACGT]", ref) || !grepl("^[ACGT,]+$", alt)) next
    aa_changes <- aa_changes_from_info(fields[[8]], alt)
    if (!nrow(aa_changes)) next
    fmt      <- strsplit(fields[[9]], ":", fixed = TRUE)[[1]]
    gt_index <- match("GT", fmt)
    if (is.na(gt_index)) next
    sample_fields <- fields[9 + keep_idx]
    gt            <- normalize_gt(sample_fields, gt_index)

    for (pos in unique(aa_changes$pos)) {
      pos_changes <- aa_changes %>% filter(pos == !!pos)
      ref_aa <- pos_changes$ref_aa[[1]]
      alt_by_index <- as.list(pos_changes$alt_aa)
      names(alt_by_index) <- as.character(pos_changes$allele_index)
      states <- gt_to_aa_state(gt, ref_aa, alt_by_index)
      key <- as.character(pos)
      if (is.null(ref_by_pos[[key]]) || is.na(ref_by_pos[[key]])) {
        ref_by_pos[[key]] <- ref_aa
      }
      states_by_pos[[key]] <- c(states_by_pos[[key]], list(states))
    }
  }

  if (!length(states_by_pos)) {
    return(list(sample_states = tibble(), sample_ids = kept_samples,
                n_png_samples = length(keep_idx),
                n_selected_samples = length(keep_idx),
                n_all_snp_records = 0, n_iep_snp_records = 0))
  }

  pos_names <- sort(as.integer(names(states_by_pos)))
  state_matrix <- do.call(cbind, lapply(as.character(pos_names), function(pos) {
    combine_position_states(states_by_pos[[pos]], ref_by_pos[[pos]])
  }))
  colnames(state_matrix) <- as.character(pos_names)

  variable_idx <- aa_state_variable_positions(state_matrix)
  if (!length(variable_idx)) {
    return(list(sample_states = tibble(), sample_ids = kept_samples,
                n_png_samples = length(keep_idx),
                n_selected_samples = length(keep_idx),
                n_all_snp_records = 0, n_iep_snp_records = 0))
  }

  state_matrix <- state_matrix[, variable_idx, drop = FALSE]
  all_positions <- as.integer(colnames(state_matrix))
  iep_cols <- which(all_positions %in% unique(as.integer(iep_positions)))

  list(
    sample_states = state_matrix_to_long(kept_samples, state_matrix),
    sample_ids          = kept_samples,
    n_png_samples        = length(kept_samples),
    n_selected_samples   = length(kept_samples),
    n_all_snp_records    = length(all_positions),
    n_iep_snp_records    = length(iep_cols)
  )
}

resolve_vcf_identity <- function(path, candidates) {
  parsed <- parse_gene_antigen_from_path(path)
  parsed_key <- list(gene_key = norm_key(parsed$gene_id), antigen_key = norm_key(parsed$antigen))
  hit <- candidates %>%
    distinct(gene_id, antigen, gene_key, antigen_key) %>%
    filter(gene_key == parsed_key$gene_key | antigen_key == parsed_key$antigen_key)
  if (nrow(hit)) {
    return(list(gene_id = hit$gene_id[[1]], antigen = hit$antigen[[1]]))
  }

  path_key <- norm_key(paste(dirname(path), basename(path), sep = "/"))
  hit <- candidates %>%
    distinct(gene_id, antigen, gene_key, antigen_key) %>%
    mutate(
      gene_match = vapply(gene_key, function(x) {
        !is.na(x) && nzchar(x) && grepl(x, path_key, fixed = TRUE)
      }, logical(1)),
      antigen_match = vapply(antigen_key, function(x) {
        !is.na(x) && nzchar(x) && grepl(x, path_key, fixed = TRUE)
      }, logical(1))
    ) %>%
    filter(gene_match | antigen_match) %>%
    mutate(match_len = pmax(nchar(gene_key), nchar(antigen_key), na.rm = TRUE)) %>%
    arrange(desc(match_len), antigen)
  if (nrow(hit)) {
    list(gene_id = hit$gene_id[[1]], antigen = hit$antigen[[1]])
  } else {
    parsed
  }
}

vcf_global_aa_frequencies <- function(path, species, gene_id, antigen) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)

  sample_count <- 0L
  out <- list()

  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break
    if (startsWith(line, "##")) next
    if (startsWith(line, "#CHROM")) {
      header <- strsplit(line, "\t", fixed = TRUE)[[1]]
      sample_count <- max(0L, length(header) - 9L)
      next
    }
    if (startsWith(line, "#")) next

    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 10) next
    ref <- fields[[4]]
    alt <- fields[[5]]
    if (nchar(ref) != 1 || !grepl("^[ACGT]", ref) || !grepl("^[ACGT,]+$", alt)) next

    aa_changes <- aa_changes_from_info(fields[[8]], alt)
    if (!nrow(aa_changes)) next

    fmt <- strsplit(fields[[9]], ":", fixed = TRUE)[[1]]
    gt_index <- match("GT", fmt)
    if (is.na(gt_index)) next
    gt <- normalize_gt(fields[10:length(fields)], gt_index)
    counted <- gt_allele_counts(gt)
    if (!counted$total) next

    for (i in seq_len(nrow(aa_changes))) {
      allele_index <- aa_changes$allele_index[[i]]
      allele_count <- as.integer(counted$counts[as.character(allele_index)] %||% 0L)
      out[[length(out) + 1L]] <- tibble(
        species = species,
        gene_id = gene_id,
        antigen = antigen,
        gene_key = norm_key(gene_id),
        antigen_key = norm_key(antigen),
        aa_position = aa_changes$pos[[i]],
        alt_aa = aa_changes$alt_aa[[i]],
        global_vcf_alt_count = allele_count,
        global_vcf_called_alleles = counted$total,
        global_vcf_af = allele_count / counted$total,
        global_vcf_samples = sample_count,
        global_vcf_file = basename(path)
      )
    }
  }

  if (!length(out)) return(tibble())
  bind_rows(out) %>%
    group_by(species, gene_id, antigen, gene_key, antigen_key, aa_position, alt_aa) %>%
    summarise(
      global_vcf_af = max(global_vcf_af, na.rm = TRUE),
      global_vcf_alt_count = sum(global_vcf_alt_count, na.rm = TRUE),
      global_vcf_called_alleles = max(global_vcf_called_alleles, na.rm = TRUE),
      global_vcf_samples = max(global_vcf_samples, na.rm = TRUE),
      global_vcf_file = paste(sort(unique(global_vcf_file)), collapse = ";"),
      .groups = "drop"
    )
}

load_global_vcf_frequencies <- function(vcf_files, species, site_table) {
  if (!length(vcf_files) || !nrow(site_table)) return(tibble())
  candidates <- site_table %>%
    distinct(gene_id, antigen, gene_key, antigen_key)

  bind_rows(lapply(vcf_files, function(path) {
    ids <- resolve_vcf_identity(path, candidates)
    if (is.na(ids$gene_id) && is.na(ids$antigen)) return(tibble())
    if (!any(candidates$gene_key == norm_key(ids$gene_id) |
             candidates$antigen_key == norm_key(ids$antigen), na.rm = TRUE)) {
      return(tibble())
    }
    vcf_global_aa_frequencies(path, species, ids$gene_id, ids$antigen)
  }))
}

apply_global_rare_filter_to_sites <- function(site_table, vcf_files, species, args) {
  site_table <- site_table %>%
    mutate(
      global_vcf_af = NA_real_,
      global_vcf_alt_count = NA_integer_,
      global_vcf_called_alleles = NA_integer_,
      global_vcf_samples = NA_integer_,
      global_vcf_file = NA_character_,
      global_rare_filtered = FALSE
    )

  if (!isTRUE(args$global_rare_filter) || !nrow(site_table) || !length(vcf_files)) {
    return(site_table)
  }

  enriched <- site_table %>%
    group_by(gene_id, antigen, gene_key, antigen_key, aa_position) %>%
    summarise(
      snp_ok = any(snp_ok, na.rm = TRUE),
      significant_iep = any(significant_iep, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(gene_id, antigen, gene_key, antigen_key) %>%
    summarise(
      n_snp_sites = sum(snp_ok, na.rm = TRUE),
      n_significant_iep_sites = sum(significant_iep, na.rm = TRUE),
      iep_fraction = if_else(n_snp_sites > 0, n_significant_iep_sites / n_snp_sites, 0),
      .groups = "drop"
    ) %>%
    filter(iep_fraction > args$global_rare_filter_iep_fraction)

  if (!nrow(enriched)) {
    message(
      species, " global rare filter: no antigens exceeded ",
      sprintf("%.0f%%", 100 * args$global_rare_filter_iep_fraction),
      " significant-IEP/SNP sites, so no global rare filtering was applied."
    )
    return(site_table)
  }
  message(
    species, " global rare filter: checking ", nrow(enriched),
    " enriched antigen(s): ",
    paste(head(paste(enriched$gene_id, enriched$antigen, sep = "/"), 8), collapse = ", "),
    if (nrow(enriched) > 8) ", ..." else ""
  )

  global_freqs <- load_global_vcf_frequencies(vcf_files, species, site_table)
  if (!nrow(global_freqs)) {
    message(
      species,
      " global rare filter: no amino-acid allele frequencies could be parsed from matching global VCF files."
    )
    return(site_table)
  }

  lookup <- global_freqs %>%
    transmute(
      gene_key,
      antigen_key,
      aa_position,
      alt_aa = toupper(trimws(alt_aa)),
      global_vcf_af,
      global_vcf_alt_count,
      global_vcf_called_alleles,
      global_vcf_samples,
      global_vcf_file
    )

  filter_significant <- "significant" %in% args$global_rare_filter_class
  filter_potential <- "potential" %in% args$global_rare_filter_class

  filtered <- site_table %>%
    select(-any_of(c("global_vcf_af", "global_vcf_alt_count",
                     "global_vcf_called_alleles", "global_vcf_samples",
                     "global_vcf_file", "global_rare_filtered"))) %>%
    left_join(
      enriched %>% select(gene_key, antigen_key, iep_fraction),
      by = c("gene_key", "antigen_key")
    ) %>%
    left_join(
      lookup,
      by = c("gene_key", "antigen_key", "aa_position", "alt_aa")
    ) %>%
    mutate(
      .class_selected = (filter_significant & significant_iep) |
        (filter_potential & potential_iep),
      global_rare_filtered = !is.na(iep_fraction) &
        .class_selected &
        !is.na(global_vcf_af) &
        global_vcf_af < args$global_rare_af_min,
      significant_iep = if_else(global_rare_filtered & filter_significant, FALSE, significant_iep),
      potential_iep = if_else(global_rare_filtered & filter_potential, FALSE, potential_iep)
    )

  n_candidates <- sum(filtered$.class_selected & !is.na(filtered$iep_fraction), na.rm = TRUE)
  n_matched <- sum(
    filtered$.class_selected & !is.na(filtered$iep_fraction) & !is.na(filtered$global_vcf_af),
    na.rm = TRUE
  )
  n_rare <- sum(filtered$global_rare_filtered, na.rm = TRUE)
  message(
    species, " global rare filter: candidate rows=", n_candidates,
    "; matched global AA AF=", n_matched,
    "; rare below ", sprintf("%.2f", args$global_rare_af_min), "=", n_rare, "."
  )

  filtered %>%
    select(-.class_selected, -iep_fraction)
}

read_filtered_escape_tables <- function(summary_path, snp_col, vcf_files, species, args) {
  sites <- read_escape_site_table(summary_path, snp_col, args)
  sites <- apply_global_rare_filter_to_sites(sites, vcf_files, species, args)
  list(
    positions = escape_positions_from_sites(sites),
    rank = summarise_escape_rank(sites),
    rare_filtered = sites %>% filter(global_rare_filtered)
  )
}

# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------
read_fasta_records <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  lines <- lines[nzchar(lines)]
  header_idx <- which(startsWith(lines, ">"))
  if (!length(header_idx)) stop("FASTA has no records: ", path)

  records <- vector("list", length(header_idx))
  for (i in seq_along(header_idx)) {
    start <- header_idx[[i]]
    end <- if (i < length(header_idx)) header_idx[[i + 1L]] - 1L else length(lines)
    header <- sub("^>", "", lines[[start]])
    seq_lines <- lines[(start + 1L):end]
    seq <- toupper(gsub("\\s+", "", paste(seq_lines, collapse = "")))
    records[[i]] <- tibble(header = header, sequence = seq)
  }
  bind_rows(records) %>%
    mutate(sample_id = vapply(header, fasta_primary_id, character(1)))
}

fasta_header_candidates <- function(header) {
  header <- sub("^>", "", as.character(header))
  token <- strsplit(header, "\\s+")[[1]][[1]]
  candidates <- c(
    header,
    token,
    unlist(strsplit(token, "[|]", perl = TRUE)),
    sub("^sample[:=]", "", grep("^sample[:=]", unlist(strsplit(header, "\\s+")), value = TRUE, ignore.case = TRUE))
  )
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  candidates
}

fasta_primary_id <- function(header) {
  fasta_header_candidates(header)[[1]]
}

fasta_keep_indices <- function(headers, selected_samples) {
  if (is.null(selected_samples)) return(seq_along(headers))
  selected_norm <- norm_key(selected_samples)
  keep <- vapply(headers, function(h) {
    any(norm_key(fasta_header_candidates(h)) %in% selected_norm)
  }, logical(1))
  which(keep)
}

sequence_matrix <- function(seqs) {
  chars <- strsplit(seqs, "", fixed = TRUE)
  max_len <- max(lengths(chars))
  mat <- matrix("-", nrow = length(chars), ncol = max_len)
  for (i in seq_along(chars)) {
    if (length(chars[[i]])) mat[i, seq_along(chars[[i]])] <- chars[[i]]
  }
  mat
}

variable_positions <- function(mat) {
  missing <- c("-", ".", "?", "X", "N", "*")
  keep <- vapply(seq_len(ncol(mat)), function(j) {
    vals <- unique(mat[, j])
    vals <- vals[!vals %in% missing]
    length(vals) > 1
  }, logical(1))
  which(keep)
}

hap_from_positions <- function(mat, positions) {
  if (!length(positions)) return(rep(NA_character_, nrow(mat)))
  positions <- positions[positions >= 1 & positions <= ncol(mat)]
  if (!length(positions)) return(rep(NA_character_, nrow(mat)))
  apply(mat[, positions, drop = FALSE], 1, paste, collapse = "|")
}

fasta_haplotypes <- function(path, selected_samples, iep_positions) {
  records <- read_fasta_records(path)
  keep_idx <- fasta_keep_indices(records$header, selected_samples)
  if (!length(keep_idx)) {
    return(list(sample_states = tibble(), sample_ids = character(),
                n_png_samples = 0,
                n_selected_samples = 0,
                n_all_snp_records = 0, n_iep_snp_records = 0))
  }

  kept <- records[keep_idx, , drop = FALSE]
  mat <- sequence_matrix(kept$sequence)
  all_var_pos <- variable_positions(mat)
  iep_var_pos <- intersect(all_var_pos, unique(as.integer(iep_positions)))
  state_pos <- sort(unique(c(all_var_pos, intersect(unique(as.integer(iep_positions)), seq_len(ncol(mat))))))

  if (!length(state_pos)) {
    return(list(
      sample_states = tibble(),
      sample_ids = kept$sample_id,
      n_png_samples = nrow(kept),
      n_selected_samples = nrow(kept),
      n_all_snp_records = length(all_var_pos),
      n_iep_snp_records = length(iep_var_pos)
    ))
  }

  state_matrix <- mat[, state_pos, drop = FALSE]
  colnames(state_matrix) <- as.character(state_pos)

  list(
    sample_states = state_matrix_to_long(kept$sample_id, state_matrix),
    sample_ids        = kept$sample_id,
    n_png_samples      = nrow(kept),
    n_selected_samples = nrow(kept),
    n_all_snp_records  = length(all_var_pos),
    n_iep_snp_records  = length(iep_var_pos)
  )
}

merged_states_to_haplotypes <- function(sample_table, state_table, ref_table, iep_positions) {
  if (!nrow(sample_table)) {
    return(list(
      sample_haps = tibble(),
      n_all_snp_records = 0L,
      n_iep_snp_records = 0L,
      all_snp_reference_haplotype = NA_character_,
      immune_escape_reference_haplotype = NA_character_
    ))
  }

  sample_table <- sample_table %>%
    mutate(sample_key = norm_key(sample_id)) %>%
    distinct(sample_key, .keep_all = TRUE)

  state_table <- state_table %>%
    mutate(
      sample_key = norm_key(sample_id),
      position = as.integer(position),
      aa_state = as.character(aa_state)
    ) %>%
    filter(sample_key %in% sample_table$sample_key,
           !is.na(position), !is.na(aa_state), aa_state != ".")

  iep_positions <- unique(as.integer(iep_positions))
  positions <- sort(unique(c(state_table$position, iep_positions)))
  positions <- positions[!is.na(positions)]
  if (!length(positions)) {
    return(list(
      sample_haps = tibble(),
      n_all_snp_records = 0L,
      n_iep_snp_records = 0L,
      all_snp_reference_haplotype = NA_character_,
      immune_escape_reference_haplotype = NA_character_
    ))
  }

  ref_lu <- ref_table %>%
    filter(position %in% positions) %>%
    distinct(position, .keep_all = TRUE) %>%
    select(position, ref_aa)

  state_collapsed <- state_table %>%
    group_by(sample_key, position) %>%
    summarise(
      aa_state = {
        vals <- unique(aa_state[!is.na(aa_state) & aa_state != "."])
        if (!length(vals)) "." else paste(sort(vals), collapse = "/")
      },
      .groups = "drop"
    )

  state_grid <- tidyr::expand_grid(
    sample_key = sample_table$sample_key,
    position = positions
  ) %>%
    left_join(state_collapsed, by = c("sample_key", "position")) %>%
    left_join(ref_lu, by = "position") %>%
    mutate(aa_state = coalesce(aa_state, ref_aa, ".")) %>%
    left_join(sample_table %>% select(sample_key, sample_id), by = "sample_key")

  state_matrix <- state_grid %>%
    select(sample_key, position, aa_state) %>%
    tidyr::pivot_wider(
      names_from = position,
      values_from = aa_state,
      values_fill = "."
    )

  matrix_sample_keys <- state_matrix$sample_key
  mat <- as.matrix(state_matrix[, setdiff(names(state_matrix), "sample_key"), drop = FALSE])
  storage.mode(mat) <- "character"
  rownames(mat) <- matrix_sample_keys

  variable_idx <- aa_state_variable_positions(mat)
  if (!length(variable_idx)) {
    return(list(
      sample_haps = tibble(),
      n_all_snp_records = 0L,
      n_iep_snp_records = 0L,
      all_snp_reference_haplotype = NA_character_,
      immune_escape_reference_haplotype = NA_character_
    ))
  }

  mat <- mat[, variable_idx, drop = FALSE]
  all_positions <- as.integer(colnames(mat))
  iep_cols <- which(all_positions %in% iep_positions)
  if (!length(iep_cols)) {
    return(list(
      sample_haps = tibble(),
      n_all_snp_records = length(all_positions),
      n_iep_snp_records = 0L,
      all_snp_reference_haplotype = NA_character_,
      immune_escape_reference_haplotype = NA_character_
    ))
  }

  reference_states <- tibble(position = all_positions) %>%
    left_join(ref_lu, by = "position") %>%
    mutate(ref_aa = coalesce(ref_aa, ".")) %>%
    pull(ref_aa)

  all_snp_reference_haplotype <- paste(reference_states, collapse = "|")
  immune_escape_reference_haplotype <- paste(reference_states[iep_cols], collapse = "|")

  sample_lookup <- sample_table %>%
    distinct(sample_key, sample_id)

  sample_haps <- tibble(
    sample_key = rownames(mat),
    all_snp_haplotype = apply(mat, 1, paste, collapse = "|"),
    immune_escape_haplotype = apply(mat[, iep_cols, drop = FALSE], 1, paste, collapse = "|")
  ) %>%
    left_join(sample_lookup, by = "sample_key") %>%
    select(sample_id, all_snp_haplotype, immune_escape_haplotype)

  list(
    sample_haps = sample_haps,
    n_all_snp_records = length(all_positions),
    n_iep_snp_records = length(iep_cols),
    all_snp_reference_haplotype = all_snp_reference_haplotype,
    immune_escape_reference_haplotype = immune_escape_reference_haplotype
  )
}

# ---------------------------------------------------------------------------
# Collapse/filter rare haplotypes:
#   1. Right (IEP): rare haps are removed from the displayed right side.
#   2. Left side keeps raw prevalence denominator from all selected samples.
#   3. Left (all-SNP): rare haps among all raw flows -> "Other".
# ---------------------------------------------------------------------------
OTHER_SENTINEL <- "__OTHER__"
FILTERED_IEP_SENTINEL <- "__FILTERED_IEP__"

iep_hamming_distance <- function(a, b) {
  aa <- strsplit(as.character(a), "|", fixed = TRUE)[[1]]
  bb <- strsplit(as.character(b), "|", fixed = TRUE)[[1]]
  max_len <- max(length(aa), length(bb))
  length(aa) <- max_len
  length(bb) <- max_len
  aa[is.na(aa)] <- "__MISSING__"
  bb[is.na(bb)] <- "__MISSING__"
  sum(aa != bb)
}

nearest_major_iep_map <- function(right_freq, max_distance) {
  split(right_freq, interaction(
    right_freq$species,
    right_freq$gene_id,
    right_freq$antigen,
    drop = TRUE
  )) %>%
    lapply(function(group_df) {
      major <- group_df %>%
        filter(!right_rare) %>%
        arrange(desc(right_n), immune_escape_haplotype)
      rare <- group_df %>% filter(right_rare)

      if (!nrow(rare) || !nrow(major)) {
        return(tibble(
          species = character(),
          gene_id = character(),
          antigen = character(),
          immune_escape_haplotype = character(),
          assigned_iep_haplotype = character(),
          nearest_iep_distance = integer()
        ))
      }

      bind_rows(lapply(seq_len(nrow(rare)), function(i) {
        distances <- vapply(
          major$immune_escape_haplotype,
          function(candidate) {
            iep_hamming_distance(
              rare$immune_escape_haplotype[[i]],
              candidate
            )
          },
          integer(1)
        )
        best_distance <- min(distances)
        best <- major[distances == best_distance, , drop = FALSE] %>%
          arrange(desc(right_n), immune_escape_haplotype) %>%
          slice(1)

        tibble(
          species = rare$species[[i]],
          gene_id = rare$gene_id[[i]],
          antigen = rare$antigen[[i]],
          immune_escape_haplotype = rare$immune_escape_haplotype[[i]],
          assigned_iep_haplotype = if (
            is.finite(best_distance) && best_distance <= max_distance
          ) best$immune_escape_haplotype[[1]] else FILTERED_IEP_SENTINEL,
          nearest_iep_distance = as.integer(best_distance)
        )
      }))
    }) %>%
    bind_rows()
}

collapse_rare_haps <- function(
    flows_df,
    left_threshold_pct,
    right_threshold_pct,
    assign_rare_to_nearest = TRUE,
    max_iep_hamming_distance = 4
) {
  right_freq <- flows_df %>%
    group_by(species, gene_id, antigen, immune_escape_haplotype) %>%
    summarise(right_n = sum(n_samples), .groups = "drop") %>%
    group_by(species, gene_id, antigen) %>%
    mutate(right_pct  = right_n / sum(right_n) * 100,
           right_rare = right_pct < right_threshold_pct) %>% ungroup()

  nearest_map <- if (isTRUE(assign_rare_to_nearest)) {
    nearest_major_iep_map(right_freq, max_iep_hamming_distance)
  } else {
    tibble(
      species = character(),
      gene_id = character(),
      antigen = character(),
      immune_escape_haplotype = character(),
      assigned_iep_haplotype = character(),
      nearest_iep_distance = integer()
    )
  }

  grouped <- flows_df %>%
    left_join(
      right_freq %>%
        select(species, gene_id, antigen, immune_escape_haplotype, right_rare),
              by = c("species", "gene_id", "antigen", "immune_escape_haplotype")) %>%
    left_join(
      nearest_map,
      by = c("species", "gene_id", "antigen", "immune_escape_haplotype")
    ) %>%
    mutate(
      immune_escape_haplotype = if_else(
        right_rare,
        coalesce(assigned_iep_haplotype, FILTERED_IEP_SENTINEL),
        immune_escape_haplotype
      )
    ) %>%
    select(-right_rare, -assigned_iep_haplotype, -nearest_iep_distance)

  if (!nrow(grouped)) return(grouped)

  left_freq <- grouped %>%
    group_by(species, gene_id, antigen, all_snp_haplotype) %>%
    summarise(
      left_n = sum(n_samples),
      raw_samples_denominator = first(raw_samples_denominator),
      .groups = "drop"
    ) %>%
    group_by(species, gene_id, antigen) %>%
    mutate(left_pct  = left_n / first(raw_samples_denominator) * 100,
           left_rare = left_pct < left_threshold_pct) %>% ungroup()

  grouped %>%
    left_join(left_freq %>% select(species, gene_id, antigen, all_snp_haplotype, left_rare),
              by = c("species", "gene_id", "antigen", "all_snp_haplotype")) %>%
    mutate(
      all_snp_haplotype = if_else(left_rare, OTHER_SENTINEL, all_snp_haplotype)
    ) %>%
    group_by(species, gene_id, antigen, all_snp_haplotype, immune_escape_haplotype) %>%
    summarise(
      n_samples = sum(n_samples),
      raw_samples_denominator = first(raw_samples_denominator),
      all_snp_reference_haplotype = first(all_snp_reference_haplotype),
      immune_escape_reference_haplotype = first(immune_escape_reference_haplotype),
      .groups = "drop"
    ) %>%
    group_by(species, gene_id, antigen) %>%
    mutate(
      total_samples_denominator = sum(n_samples[immune_escape_haplotype != FILTERED_IEP_SENTINEL]),
      raw_samples_denominator = first(raw_samples_denominator)
    ) %>%
    ungroup()
}

# ---------------------------------------------------------------------------
# Per-species processing
# ---------------------------------------------------------------------------
process_species <- function(species, vcf_dir, fasta_dir, metadata, summary_path, snp_col, sample_scope, args) {
  species_key <- if (species == "Pf7") "pf" else "pv"
  selected_samples <- metadata_samples_for_scope(metadata, args, species_key, sample_scope)
  vcf_files <- sort(c(
    list.files(vcf_dir, pattern = "\\.vcf\\.gz$", full.names = TRUE),
    list.files(vcf_dir, pattern = "\\.vcf$",      full.names = TRUE)
  ))
  escape_tables <- read_filtered_escape_tables(summary_path, snp_col, vcf_files, species, args)
  escape_positions <- escape_tables$positions
  escape_rank <- escape_tables$rank
  rare_filtered_escape <- escape_tables$rare_filtered
  if (isTRUE(args$global_rare_filter)) {
    message(
      species, " globally rare IEP rows removed: ",
      nrow(rare_filtered_escape),
      "; unique positions: ",
      n_distinct(rare_filtered_escape$gene_key, rare_filtered_escape$antigen_key,
                 rare_filtered_escape$aa_position)
    )
  }
  reference_aa <- read_reference_aa_table(summary_path)

  fasta_files <- character()
  if (!is.null(fasta_dir) && nzchar(fasta_dir)) {
    if (!dir.exists(fasta_dir)) stop("FASTA directory not found: ", fasta_dir)
    fasta_files <- sort(list.files(
      fasta_dir,
      pattern = "\\.(fa|fasta|faa|fas|fna|aln|afa)(\\.gz)?$",
      full.names = TRUE,
      ignore.case = TRUE,
      recursive = TRUE
    ))
  }

  if (!length(vcf_files) && !length(fasta_files))
    stop("No VCF or FASTA files found for ", species)

  summaries   <- list()
  flows       <- list()
  raw_flows   <- list()
  diagnostics <- list()
  source_states <- list()
  source_samples <- list()

  input_table <- bind_rows(
    tibble(input_file = vcf_files, input_type = "vcf"),
    tibble(input_file = fasta_files, input_type = "fasta")
  )

  for (i in seq_len(nrow(input_table))) {
    input_file <- input_table$input_file[[i]]
    input_type <- input_table$input_type[[i]]
    ids          <- infer_gene_antigen_from_path(input_file, escape_positions)
    iep_positions <- positions_for_vcf(escape_positions, ids$gene_id, ids$antigen)
    if (!length(iep_positions)) {
      parsed_ids <- parse_gene_antigen_from_path(input_file)
      diagnostics[[length(diagnostics) + 1L]] <- tibble(
        species                             = species,
        gene_id                             = parsed_ids$gene_id,
        antigen                             = parsed_ids$antigen,
        input_type                          = input_type,
        input_file                          = input_file,
        sample_scope                        = sample_scope,
        n_png_samples_in_vcf                = NA_integer_,
        n_selected_samples_in_vcf           = NA_integer_,
        n_immune_escape_aa_positions        = 0L,
        n_all_snp_records_in_vcf            = NA_integer_,
        n_immune_escape_snp_records_in_vcf  = NA_integer_,
        n_sample_haps_before_filter         = NA_integer_,
        n_reduced_haplotypes_before_filter  = NA_integer_,
        n_reduced_haplotypes_passing_filter = NA_integer_,
        n_sample_haps_after_filter          = NA_integer_,
        min_original_haplotype_percent      = args$min_original_haplotype_percent,
        min_reduced_haplotype_percent       = args$min_reduced_haplotype_percent,
        status                              = "skipped_no_escape_positions_matched_to_filename"
      )
      next
    }

    message("Reading ", input_type, " ", species, " ", ids$antigen, " (", ids$gene_id, ")")
    h <- if (input_type == "fasta") {
      fasta_selected <- if (isTRUE(args$fasta_filter_by_metadata)) selected_samples else NULL
      fasta_haplotypes(input_file, fasta_selected, iep_positions)
    } else {
      vcf_haplotypes(input_file, selected_samples, iep_positions)
    }

    diag_row <- tibble(
      species                             = species,
      gene_id                             = ids$gene_id,
      antigen                             = ids$antigen,
      input_type                          = input_type,
      input_file                          = input_file,
      sample_scope                        = sample_scope,
      n_png_samples_in_vcf                = h$n_png_samples,
      n_selected_samples_in_vcf           = h$n_selected_samples,
      n_immune_escape_aa_positions        = length(iep_positions),
      n_all_snp_records_in_vcf            = h$n_all_snp_records,
      n_immune_escape_snp_records_in_vcf  = h$n_iep_snp_records,
      n_sample_haps_before_filter         = NA_integer_,
      n_reduced_haplotypes_before_filter  = NA_integer_,
      n_reduced_haplotypes_passing_filter = NA_integer_,
      n_sample_haps_after_filter          = NA_integer_,
      min_original_haplotype_percent      = args$min_original_haplotype_percent,
      min_reduced_haplotype_percent       = args$min_reduced_haplotype_percent,
      status                              = "source_started"
    )
    if (h$n_selected_samples < args$min_png_samples) {
      diagnostics[[length(diagnostics) + 1L]] <- mutate(diag_row, status = "skipped_too_few_png_samples")
      next
    }
    if (!length(h$sample_ids)) {
      diagnostics[[length(diagnostics) + 1L]] <- mutate(diag_row, status = "skipped_no_selected_samples")
      next
    }

    diag_row$n_sample_haps_before_filter <- length(h$sample_ids)
    diagnostics[[length(diagnostics) + 1L]] <- mutate(
      diag_row,
      status = if_else(h$n_iep_snp_records > 0,
                       "source_loaded",
                       "source_loaded_no_iep_records")
    )

    source_samples[[length(source_samples) + 1L]] <- tibble(sample_id = h$sample_ids) %>%
      mutate(
        species = species,
        gene_id = ids$gene_id,
        antigen = ids$antigen,
        input_type = input_type,
        input_file = input_file,
        n_selected_samples_source = h$n_selected_samples,
        n_all_snp_records_source = h$n_all_snp_records,
        n_iep_snp_records_source = h$n_iep_snp_records,
        n_immune_escape_aa_positions = length(iep_positions)
      )

    if (nrow(h$sample_states)) {
      source_states[[length(source_states) + 1L]] <- h$sample_states %>%
        mutate(
          species = species,
          gene_id = ids$gene_id,
          antigen = ids$antigen,
          input_type = input_type,
          input_file = input_file
        )
    }
  }

  all_source_samples <- bind_rows(source_samples)
  all_source_states <- if (length(source_states)) {
    bind_rows(source_states)
  } else {
    tibble(
      sample_id = character(),
      position = integer(),
      aa_state = character(),
      species = character(),
      gene_id = character(),
      antigen = character(),
      input_type = character(),
      input_file = character()
    )
  }
  if (!nrow(all_source_samples)) {
    return(list(
      summary     = bind_rows(summaries),
      flows       = bind_rows(flows),
      raw_flows   = bind_rows(raw_flows),
      diagnostics = bind_rows(diagnostics)
    ))
  }

  for (key in split(seq_len(nrow(all_source_samples)), paste(all_source_samples$species, all_source_samples$gene_id, all_source_samples$antigen, sep = "\r"))) {
    antigen_samples_all <- all_source_samples[key, , drop = FALSE]
    antigen_samples <- antigen_samples_all %>%
      mutate(
        source_priority = if_else(input_type == "vcf", 1L, 2L),
        sample_key = norm_key(sample_id)
      ) %>%
      arrange(source_priority, input_file, sample_id) %>%
      distinct(sample_key, .keep_all = TRUE)

    ids <- list(
      gene_id = antigen_samples$gene_id[[1]],
      antigen = antigen_samples$antigen[[1]]
    )
    iep_positions <- positions_for_vcf(escape_positions, ids$gene_id, ids$antigen)
    rank_row <- escape_rank %>%
      filter(gene_key == norm_key(ids$gene_id) | antigen_key == norm_key(ids$antigen)) %>%
      arrange(desc(weighted_iep_density), desc(weighted_iep_score), antigen)
    if (!nrow(rank_row)) {
      rank_row <- tibble(
        n_total_sites = NA_integer_,
        n_snp_sites = NA_integer_,
        n_significant_iep_sites = NA_integer_,
        n_potential_iep_sites = NA_integer_,
        n_global_rare_filtered_iep_sites = NA_integer_,
        weighted_iep_score = NA_real_,
        weighted_iep_density = NA_real_
      )
    }
    input_types <- paste(sort(unique(antigen_samples_all$input_type)), collapse = "+")
    input_files <- paste(unique(antigen_samples_all$input_file), collapse = ";")
    source_stats <- antigen_samples_all %>%
      distinct(input_type, input_file, n_all_snp_records_source,
               n_iep_snp_records_source, n_immune_escape_aa_positions)
    message("Processing merged ", input_types, " ", species, " ", ids$antigen, " (", ids$gene_id, ")")

    antigen_states <- all_source_states %>%
      filter(species == !!species, gene_id == !!ids$gene_id, antigen == !!ids$antigen) %>%
      mutate(
        source_priority = if_else(input_type == "vcf", 1L, 2L),
        sample_key = norm_key(sample_id)
      ) %>%
      arrange(source_priority, input_file, sample_id, position) %>%
      distinct(sample_key, position, .keep_all = TRUE)

    antigen_ref <- reference_aa %>%
      filter(gene_key == norm_key(ids$gene_id) | antigen_key == norm_key(ids$antigen))

    merged_haps <- merged_states_to_haplotypes(
      antigen_samples %>% select(sample_id),
      antigen_states %>% select(sample_id, position, aa_state),
      antigen_ref,
      iep_positions
    )

    sample_haps <- merged_haps$sample_haps
    if (!nrow(sample_haps)) {
      diagnostics[[length(diagnostics) + 1L]] <- tibble(
        species                             = species,
        gene_id                             = ids$gene_id,
        antigen                             = ids$antigen,
        input_type                          = input_types,
        input_file                          = input_files,
        sample_scope                        = sample_scope,
        n_png_samples_in_vcf                = nrow(antigen_samples),
        n_selected_samples_in_vcf           = nrow(antigen_samples),
        n_immune_escape_aa_positions        = length(iep_positions),
        n_all_snp_records_in_vcf            = merged_haps$n_all_snp_records,
        n_immune_escape_snp_records_in_vcf  = merged_haps$n_iep_snp_records,
        n_sample_haps_before_filter         = nrow(antigen_samples),
        n_reduced_haplotypes_before_filter  = NA_integer_,
        n_reduced_haplotypes_passing_filter = NA_integer_,
        n_sample_haps_after_filter          = NA_integer_,
        min_original_haplotype_percent      = args$min_original_haplotype_percent,
        min_reduced_haplotype_percent       = args$min_reduced_haplotype_percent,
        status                              = "skipped_no_sample_haplotypes_after_source_merge"
      )
      next
    }

    # Count raw distinct haplotypes AFTER source states are merged onto one
    # shared amino-acid position matrix.
    n_raw_all_snp_haps <- n_distinct(sample_haps$all_snp_haplotype)
    n_raw_iep_haps     <- n_distinct(sample_haps$immune_escape_haplotype)

    # Count raw flows before collapsing
    flows_raw <- sample_haps %>%
      mutate(
        species = species,
        gene_id = ids$gene_id,
        antigen = ids$antigen,
        raw_samples_denominator = nrow(sample_haps),
        all_snp_reference_haplotype = merged_haps$all_snp_reference_haplotype,
        immune_escape_reference_haplotype = merged_haps$immune_escape_reference_haplotype
      ) %>%
      count(species, gene_id, antigen,
            all_snp_haplotype, immune_escape_haplotype,
            all_snp_reference_haplotype, immune_escape_reference_haplotype,
            raw_samples_denominator,
            name = "n_samples")
    raw_flows[[length(raw_flows) + 1L]] <- flows_raw

    diag_row <- tibble(
      species                             = species,
      gene_id                             = ids$gene_id,
      antigen                             = ids$antigen,
      input_type                          = input_types,
      input_file                          = input_files,
      sample_scope                        = sample_scope,
      n_png_samples_in_vcf                = nrow(sample_haps),
      n_selected_samples_in_vcf           = nrow(sample_haps),
      n_immune_escape_aa_positions        = max(source_stats$n_immune_escape_aa_positions, na.rm = TRUE),
      n_all_snp_records_in_vcf            = merged_haps$n_all_snp_records,
      n_immune_escape_snp_records_in_vcf  = merged_haps$n_iep_snp_records,
      n_sample_haps_before_filter         = nrow(sample_haps),
      n_reduced_haplotypes_before_filter  = n_raw_iep_haps,
      n_reduced_haplotypes_passing_filter = NA_integer_,
      n_sample_haps_after_filter          = NA_integer_,
      min_original_haplotype_percent      = args$min_original_haplotype_percent,
      min_reduced_haplotype_percent       = args$min_reduced_haplotype_percent,
      status                              = "started"
    )

    # Collapse rare original haplotypes into Other, but remove rare reduced
    # immune-escape haplotypes from the displayed alluvial.
    flows_c <- collapse_rare_haps(
      flows_raw,
      args$min_original_haplotype_percent,
      args$min_reduced_haplotype_percent,
      args$assign_rare_iep_to_nearest,
      args$max_iep_hamming_distance
    )

    diag_row$n_reduced_haplotypes_passing_filter <- n_distinct(
      flows_c$immune_escape_haplotype[
        !flows_c$immune_escape_haplotype %in% c(OTHER_SENTINEL, FILTERED_IEP_SENTINEL)
      ]
    )
    diag_row$n_sample_haps_after_filter <- sum(
      flows_c$n_samples[flows_c$immune_escape_haplotype != FILTERED_IEP_SENTINEL]
    )

    if (!nrow(flows_c)) {
      diagnostics[[length(diagnostics) + 1L]] <- mutate(diag_row, status = "skipped_no_haplotypes_after_collapse")
      next
    }
    diagnostics[[length(diagnostics) + 1L]] <- mutate(diag_row, status = "plotted")

    summaries[[length(summaries) + 1L]] <- tibble(
      species                             = species,
      gene_id                             = ids$gene_id,
      antigen                             = ids$antigen,
      input_type                          = input_types,
      input_file                          = input_files,
      sample_scope                        = sample_scope,
      n_png_samples                       = nrow(sample_haps),
      n_selected_samples                  = nrow(sample_haps),
      n_immune_escape_aa_positions        = diag_row$n_immune_escape_aa_positions,
      n_all_snp_records_in_vcf            = diag_row$n_all_snp_records_in_vcf,
      n_immune_escape_snp_records_in_vcf  = diag_row$n_immune_escape_snp_records_in_vcf,
      min_original_haplotype_percent      = args$min_original_haplotype_percent,
      min_reduced_haplotype_percent       = args$min_reduced_haplotype_percent,
      # Raw counts before any collapsing/dropping — used for antigen ranking
      n_all_snp_haplotypes_raw            = n_raw_all_snp_haps,
      n_iep_haplotypes_raw                = n_raw_iep_haps,
      n_total_sites_for_density           = rank_row$n_total_sites[[1]],
      n_snp_sites_for_density             = rank_row$n_snp_sites[[1]],
      n_significant_iep_sites             = rank_row$n_significant_iep_sites[[1]],
      n_potential_iep_sites               = rank_row$n_potential_iep_sites[[1]],
      n_global_rare_filtered_iep_sites    = rank_row$n_global_rare_filtered_iep_sites[[1]],
      weighted_iep_score                  = rank_row$weighted_iep_score[[1]],
      weighted_iep_density                = rank_row$weighted_iep_density[[1]],
      # Post-collapse displayed counts, excluding the merged Other block.
      n_all_snp_haplotypes                = n_distinct(flows_c$all_snp_haplotype[flows_c$all_snp_haplotype != OTHER_SENTINEL]),
      n_immune_escape_haplotypes          = n_distinct(
        flows_c$immune_escape_haplotype[
          !flows_c$immune_escape_haplotype %in% c(OTHER_SENTINEL, FILTERED_IEP_SENTINEL)
        ]
      ),
      n_png_samples_after_reduced_filter  = sum(
        flows_c$n_samples[flows_c$immune_escape_haplotype != FILTERED_IEP_SENTINEL]
      )
    )
    # Attach raw counts to flows so ranking can use them downstream
    flows_c <- flows_c %>%
      mutate(n_all_snp_haplotypes_raw = n_raw_all_snp_haps,
             n_iep_haplotypes_raw     = n_raw_iep_haps,
             weighted_iep_density     = rank_row$weighted_iep_density[[1]],
             weighted_iep_score       = rank_row$weighted_iep_score[[1]],
             n_significant_iep_sites  = rank_row$n_significant_iep_sites[[1]],
             n_potential_iep_sites    = rank_row$n_potential_iep_sites[[1]],
             n_global_rare_filtered_iep_sites = rank_row$n_global_rare_filtered_iep_sites[[1]],
             n_snp_sites_for_density  = rank_row$n_snp_sites[[1]])
    flows[[length(flows) + 1L]] <- flows_c
  }

  list(
    summary     = bind_rows(summaries),
    flows       = bind_rows(flows),
    raw_flows   = bind_rows(raw_flows),
    diagnostics = bind_rows(diagnostics)
  )
}

# ---------------------------------------------------------------------------
# Alluvial layout builders
# ---------------------------------------------------------------------------

# Build per-antigen Other counts lookup (for label annotation)
build_other_counts <- function(flows_df, hap_col) {
  flows_df %>%
    filter(.data[[hap_col]] == OTHER_SENTINEL) %>%
    group_by(species, gene_id, antigen) %>%
    summarise(n_other = sum(n_samples), .groups = "drop")
}

hap_distance_to_reference <- function(hap, ref_hap) {
  if (is.na(hap) || is.na(ref_hap) || !nzchar(hap) || !nzchar(ref_hap)) {
    return(NA_integer_)
  }
  aa <- strsplit(as.character(hap), "|", fixed = TRUE)[[1]]
  rr <- strsplit(as.character(ref_hap), "|", fixed = TRUE)[[1]]
  max_len <- max(length(aa), length(rr))
  length(aa) <- max_len
  length(rr) <- max_len
  missing <- c(NA_character_, ".", "-", "?", "X", "__MISSING__")
  keep <- !aa %in% missing & !rr %in% missing
  if (!any(keep)) return(NA_integer_)
  sum(aa[keep] != rr[keep])
}

rank_haplotype_labels <- function(flows_df, hap_col, prefix) {
  hap_col_sym <- rlang::ensym(hap_col)
  hap_col_name <- rlang::as_string(hap_col_sym)
  raw_total_col <- if (hap_col_name == "all_snp_haplotype") {
    "n_all_snp_haplotypes_raw"
  } else if (hap_col_name == "immune_escape_haplotype") {
    "n_iep_haplotypes_raw"
  } else {
    NA_character_
  }
  ref_col <- if (hap_col_name == "all_snp_haplotype") {
    "all_snp_reference_haplotype"
  } else if (hap_col_name == "immune_escape_haplotype") {
    "immune_escape_reference_haplotype"
  } else {
    NA_character_
  }

  flows_df %>%
    group_by(species, gene_id, antigen, hap = !!hap_col_sym) %>%
    summarise(
      n = sum(n_samples),
      raw_total = if (!is.na(raw_total_col) && raw_total_col %in% names(flows_df)) {
        first(.data[[raw_total_col]])
      } else {
        NA_integer_
      },
      reference_haplotype = if (!is.na(ref_col) && ref_col %in% names(flows_df)) {
        first(.data[[ref_col]])
      } else {
        NA_character_
      },
      .groups = "drop"
    ) %>%
    group_by(species, gene_id, antigen) %>%
    # Named haps sorted by frequency desc; sentinel groups always last.
    mutate(
      hap_sort = case_when(
        hap == OTHER_SENTINEL        ~ 2L,
        hap == FILTERED_IEP_SENTINEL ~ 1L,
        TRUE                         ~ 0L
      )
    ) %>%
    arrange(hap_sort, desc(n), hap, .by_group = TRUE) %>%
    mutate(
      rank       = cumsum(!(hap %in% c(OTHER_SENTINEL, FILTERED_IEP_SENTINEL))),
      named_n    = sum(!(hap %in% c(OTHER_SENTINEL, FILTERED_IEP_SENTINEL))),
      other_from = named_n + 1L,
      other_to   = if_else(!is.na(first(raw_total)) & first(raw_total) >= other_from,
                           as.integer(first(raw_total)), other_from),
      hap_label = case_when(
        hap == FILTERED_IEP_SENTINEL                   ~ "Rare IEP",
        hap == OTHER_SENTINEL & other_from == other_to ~ paste0(prefix, other_from),
        hap == OTHER_SENTINEL                          ~ paste0(prefix, other_from, "-", prefix, other_to),
        TRUE                                           ~ paste0(prefix, rank)
      ),
      reference_distance = vapply(
        seq_along(hap),
        function(i) {
          if (hap[[i]] %in% c(OTHER_SENTINEL, FILTERED_IEP_SENTINEL)) return(NA_integer_)
          hap_distance_to_reference(hap[[i]], reference_haplotype[[i]])
        },
        integer(1)
      ),
      reference_like = hap_col_name == "immune_escape_haplotype" &
        !(hap %in% c(OTHER_SENTINEL, FILTERED_IEP_SENTINEL)) &
        !is.na(reference_distance) &
        reference_distance == 0L,
      hap_label = if_else(
        reference_like,
        paste0(hap_label, "*"),
        hap_label
      )
    ) %>%
    select(-hap_sort, -named_n, -other_from, -other_to) %>%
    ungroup()
}

stack_strata <- function(strata_df) {
  strata_df %>%
    group_by(panel, side) %>%
    # Named haps sorted by frequency desc; Other and filtered remainder last.
    mutate(
      hap_sort = case_when(
        hap == FILTERED_IEP_SENTINEL ~ 2L,
        hap == OTHER_SENTINEL        ~ 1L,
        TRUE                         ~ 0L
      )
    ) %>%
    arrange(hap_sort, desc(n), hap_label, .by_group = TRUE) %>%
    mutate(
      prop = n / first(total_samples),
      ymax = 1 - cumsum(lag(prop, default = 0)),
      ymin = ymax - prop,
      ymid = (ymin + ymax) / 2,
      pct  = prop * 100
    ) %>%
    select(-hap_sort) %>%
    ungroup()
}

make_ribbon_polygon <- function(flow_row, n_points = 28) {
  t      <- seq(0, 1, length.out = n_points)
  smooth <- t * t * (3 - 2 * t)
  top_y  <- flow_row$left_ymax  + (flow_row$right_ymax  - flow_row$left_ymax)  * smooth
  bot_y  <- flow_row$left_ymin  + (flow_row$right_ymin  - flow_row$left_ymin)  * smooth
  tibble(
    flow_id    = flow_row$flow_id,
    panel      = flow_row$panel,
    species    = flow_row$species,
    antigen    = flow_row$antigen,
    left_label = flow_row$left_label,
    right_label= flow_row$right_label,
    n_samples  = flow_row$n_samples,
    x = c(t, rev(t)),
    y = c(top_y, rev(bot_y))
  )
}

build_alluvial_layers <- function(flows_df, args) {
  has_denominator <- "total_samples_denominator" %in% names(flows_df)
  has_raw_denominator <- "raw_samples_denominator" %in% names(flows_df)
  stub_prop <- if (isTRUE(args$show_filtered_iep_stub)) {
    pmin(pmax(args$filtered_iep_stub_percent / 100, 0), 0.20)
  } else {
    0
  }

  # Both sides use the same "H" prefix
  # Left = all-SNP haplotypes, Right = immune-escape haplotypes
  # Ranked independently so H1 on each side = most common haplotype on that side
  left_labels  <- rank_haplotype_labels(flows_df, all_snp_haplotype,       args$left_prefix) %>%
    rename(all_snp_haplotype = hap, left_label = hap_label, left_n = n) %>%
    select(species, gene_id, antigen, all_snp_haplotype, left_label, left_n)
  right_labels <- rank_haplotype_labels(flows_df, immune_escape_haplotype,  args$right_prefix) %>%
    rename(
      immune_escape_haplotype = hap,
      right_label = hap_label,
      right_n = n,
      right_reference_like = reference_like
    ) %>%
    select(
      species, gene_id, antigen, immune_escape_haplotype,
      right_label, right_n, right_reference_like
    )

  labelled <- flows_df %>%
    left_join(left_labels,  by = c("species", "gene_id", "antigen", "all_snp_haplotype")) %>%
    left_join(right_labels, by = c("species", "gene_id", "antigen", "immune_escape_haplotype")) %>%
    group_by(species, gene_id, antigen) %>%
    mutate(
      right_total_samples = if (has_denominator) first(total_samples_denominator) else sum(n_samples),
      left_total_samples = if (has_raw_denominator) first(raw_samples_denominator) else right_total_samples,
      has_filtered_iep = any(immune_escape_haplotype == FILTERED_IEP_SENTINEL),
      filtered_iep_total = sum(n_samples[immune_escape_haplotype == FILTERED_IEP_SENTINEL]),
      filtered_stub_prop = if_else(
        has_filtered_iep & stub_prop > 0 & right_total_samples > 0,
        stub_prop,
        0
      ),
      right_display_total_samples = if_else(
        filtered_stub_prop > 0,
        right_total_samples / (1 - filtered_stub_prop),
        right_total_samples
      )
    ) %>%
    ungroup() %>%
    mutate(
      panel   = paste(species, antigen, sep = ": "),
      flow_id = row_number()
    )

  left_strata  <- labelled %>%
    distinct(panel, species, antigen, side = "left",
             hap_label = left_label,  hap = all_snp_haplotype,      n = left_n,
             total_samples = left_total_samples) %>%
    mutate(reference_like = FALSE)
  labelled_visible <- labelled %>%
    filter(immune_escape_haplotype != FILTERED_IEP_SENTINEL)

  labelled_ribbons <- labelled %>%
    filter(immune_escape_haplotype != FILTERED_IEP_SENTINEL | filtered_stub_prop > 0)

  right_strata <- labelled_visible %>%
    distinct(panel, species, antigen, side = "right",
             hap_label = right_label, hap = immune_escape_haplotype, n = right_n,
             total_samples = right_display_total_samples,
             reference_like = right_reference_like)

  right_filtered_strata <- labelled %>%
    filter(immune_escape_haplotype == FILTERED_IEP_SENTINEL,
           filtered_stub_prop > 0) %>%
    distinct(panel, species, antigen, side = "right",
             hap_label = right_label, hap = immune_escape_haplotype,
             total_samples = right_display_total_samples,
             filtered_stub_prop) %>%
    mutate(
      n = total_samples * filtered_stub_prop,
      reference_like = FALSE
    )

  strata <- bind_rows(left_strata, right_strata, right_filtered_strata) %>%
    stack_strata() %>%
    mutate(
      x           = if_else(side == "left", 0, 1),
      xmin        = if_else(side == "left", -0.024,  0.976),
      xmax        = if_else(side == "left",  0.024,  1.024),
      label_x     = if_else(side == "left", -0.075,  1.075),
      label_hjust = if_else(side == "left",  1,       0),
      display_label = if_else(
        pct >= args$label_percent_threshold,
        paste0(hap_label, "\n", sprintf("%.1f%%", pct)),
        ""
      )
    )

  left_stack <- labelled_ribbons %>%
    arrange(panel, left_label, desc(n_samples), right_label) %>%
    group_by(panel, left_label) %>%
    mutate(
      left_total = sum(n_samples),
      left_prop  = n_samples / unique(left_total_samples),
      left_ymin  = first(strata$ymin[match(
        paste(panel, "left", left_label),
        paste(strata$panel, strata$side, strata$hap_label)
      )]) + cumsum(lag(left_prop, default = 0)),
      left_ymax  = left_ymin + left_prop
    ) %>%
    ungroup()

  right_stack <- left_stack %>%
    arrange(panel, right_label, desc(n_samples), left_label) %>%
    group_by(panel, right_label) %>%
    mutate(
      right_is_filtered = immune_escape_haplotype == FILTERED_IEP_SENTINEL,
      right_prop = if_else(
        right_is_filtered,
        n_samples / first(filtered_iep_total) * first(filtered_stub_prop),
        n_samples / unique(right_display_total_samples)
      ),
      right_ymin = first(strata$ymin[match(
        paste(panel, "right", right_label),
        paste(strata$panel, strata$side, strata$hap_label)
      )]) + cumsum(lag(right_prop, default = 0)),
      right_ymax = right_ymin + right_prop
    ) %>%
    select(-right_is_filtered) %>%
    ungroup()

  ribbons <- bind_rows(lapply(
    seq_len(nrow(right_stack)),
    function(i) make_ribbon_polygon(right_stack[i, ])
  ))

  list(flows = labelled, strata = strata, ribbons = ribbons)
}

# Colour original all-SNP haplotype blocks; rare merged blocks stay grey.
assign_node_colours <- function(strata) {
  palette <- c(
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
    "#56B4E9", "#F0E442", "#7F3C8D", "#11A579", "#3969AC",
    "#F2B701", "#E73F74", "#80BA5A", "#E68310", "#008695"
  )
  label_colours <- strata %>%
    filter(!hap %in% c(OTHER_SENTINEL, FILTERED_IEP_SENTINEL)) %>%
    distinct(side, hap_label) %>%
    group_by(side) %>%
    arrange(hap_label, .by_group = TRUE) %>%
    mutate(node_fill = rep(palette, length.out = n())) %>%
    ungroup()

  strata %>%
    left_join(label_colours, by = c("side", "hap_label")) %>%
    mutate(
      is_other = hap == OTHER_SENTINEL,
      is_filtered = hap == FILTERED_IEP_SENTINEL,
      node_fill = case_when(
        is_filtered ~ "#E6E6E6",
        is_other    ~ "#BDBDBD",
        TRUE        ~ node_fill
      ),
      node_fill = replace_na(node_fill, "white"),
      node_linetype = case_when(
        is_filtered ~ "dotted",
        is_other    ~ "dashed",
        TRUE        ~ "solid"
      )
    )
}

make_filtered_dots <- function(strata) {
  filtered <- strata %>% filter(is_filtered, ymax > ymin, xmax > xmin)
  if (!nrow(filtered)) {
    return(tibble(panel = character(), x = numeric(), y = numeric()))
  }

  bind_rows(lapply(seq_len(nrow(filtered)), function(i) {
    row <- filtered[i, ]
    y_step <- 0.026
    y_from <- row$ymin + y_step / 2
    y_to <- row$ymax - y_step / 2
    ys <- if (is.na(y_from) || is.na(y_to) || y_from > y_to) {
      row$ymid
    } else {
      seq(y_from, y_to, by = y_step)
    }
    xs <- seq(row$xmin + 0.006, row$xmax - 0.006, length.out = 3)
    expand.grid(x = xs, y = ys) %>%
      as_tibble() %>%
      mutate(panel = row$panel)
  }))
}

label_panel_geometry <- function(args, n_panels) {
  n_panels <- max(1L, as.integer(n_panels))
  ncol <- min(max(1L, args$facet_ncol), n_panels)
  nrow <- ceiling(n_panels / ncol)

  # Approximate drawable facet height. This keeps the label decision tied to
  # physical page size instead of a fixed percent cutoff.
  panel_height_in <- max(0.35, (args$pdf_height - 0.75) / nrow)
  panel_width_in  <- max(0.60, (args$pdf_width  - 0.65) / ncol)
  list(
    ncol = ncol,
    nrow = nrow,
    panel_height_pt = panel_height_in * 72,
    panel_width_pt  = panel_width_in  * 72
  )
}

adjust_label_group <- function(df, min_gap_y, y_min = 0, y_max = 1, x_stagger = 0) {
  if (!nrow(df)) return(df)

  df <- df %>% arrange(ymid, desc(pct), hap_label)
  y <- pmin(pmax(df$ymid, y_min + min_gap_y / 2), y_max - min_gap_y / 2)

  if (length(y) > 1L) {
    for (i in 2:length(y)) {
      y[[i]] <- max(y[[i]], y[[i - 1L]] + min_gap_y)
    }
    overflow <- y[[length(y)]] - (y_max - min_gap_y / 2)
    if (is.finite(overflow) && overflow > 0) y <- y - overflow
    y <- pmin(pmax(y, y_min + min_gap_y / 2), y_max - min_gap_y / 2)
  }

  side_sign <- ifelse(df$side == "left", -1, 1)
  stagger <- rep(c(0, 1, 2), length.out = nrow(df)) * x_stagger
  df$label_y <- y
  df$label_x_adj <- df$label_x + side_sign * stagger
  df
}

prepare_strata_labels <- function(strata, args, n_panels) {
  geom <- label_panel_geometry(args, n_panels)
  min_gap_y <- (args$label_min_gap_pt * args$font_scale * 1.35) / geom$panel_height_pt
  auto_label_cap <- max(
    1L,
    floor(geom$panel_height_pt / (args$label_min_gap_pt * args$font_scale))
  )
  auto_label_cap <- min(auto_label_cap, args$auto_max_labels_per_side)

  labelled <- strata %>%
    mutate(
      reference_like = replace_na(reference_like, FALSE),
      block_height_pt = pmax(0, ymax - ymin) * geom$panel_height_pt,
      label_has_room = block_height_pt >= args$label_min_height_pt * args$font_scale,
      base_show_label = is_other | is_filtered | pct >= args$label_percent_threshold
    ) %>%
    group_by(panel, side) %>%
    arrange(desc(is_other), desc(pct), hap_label, .by_group = TRUE) %>%
    mutate(
      label_rank = row_number(),
      side_label_cap = if_else(side == "left",
                               args$max_left_labels_per_side,
                               args$max_right_labels_per_side),
      side_label_cap = if_else(
        is.na(side_label_cap),
        auto_label_cap,
        pmin(side_label_cap, auto_label_cap)
      ),
      show_label = base_show_label &
        (is_other | is_filtered |
           reference_like |
           side_label_cap < 0L |
           (side_label_cap > 0L & label_rank <= side_label_cap))
    ) %>%
    ungroup() %>%
    mutate(
      display_label = if_else(
        show_label,
        hap_label,
        ""
      ),
      label_y = ymid,
      label_x_adj = label_x,
      label_size = case_when(
        is_other ~ 1.22 * args$font_scale,
        pct >= 10 ~ 1.28 * args$font_scale,
        !label_has_room ~ 0.95 * args$font_scale,
        TRUE ~ 1.08 * args$font_scale
      )
    )

  shown_input <- labelled %>% filter(show_label)
  shown <- if (nrow(shown_input)) {
    bind_rows(lapply(
      split(shown_input, list(shown_input$panel, shown_input$side), drop = TRUE),
      adjust_label_group,
      min_gap_y = min_gap_y,
      x_stagger = if (isTRUE(args$align_label_columns)) 0 else args$label_stagger_x
    )) %>%
      select(panel, side, hap_label, label_y, label_x_adj)
  } else {
    tibble(
      panel = strata$panel[0],
      side = character(),
      hap_label = character(),
      label_y = numeric(),
      label_x_adj = numeric()
    )
  }

  labelled %>%
    select(-label_y, -label_x_adj) %>%
    left_join(shown, by = c("panel", "side", "hap_label")) %>%
    mutate(
      label_y = coalesce(label_y, ymid),
      label_x_adj = coalesce(label_x_adj, label_x)
    )
}

make_label_layer <- function(strata, args) {
  label_df <- strata %>% filter(show_label)
  leader_df <- label_df %>%
    mutate(
      leader_x = if_else(side == "left", label_x_adj + 0.006, label_x_adj - 0.006),
      leader_xend = if_else(side == "left", xmin - 0.002, xmax + 0.002),
      leader_y = label_y,
      leader_yend = ymid
    )

  leader_layer <- if (isTRUE(args$show_label_leaders) && nrow(leader_df)) {
    geom_segment(
      data = leader_df,
      aes(x = leader_x, y = leader_y, xend = leader_xend, yend = leader_yend),
      inherit.aes = FALSE,
      colour = "#555555",
      linewidth = 0.12,
      lineend = "round",
      arrow = grid::arrow(length = unit(0.018, "inches"), type = "closed")
    )
  } else {
    NULL
  }

  text_layer <- if (isTRUE(args$use_ggrepel) && requireNamespace("ggrepel", quietly = TRUE)) {
    message("Using ggrepel for haplotype labels.")
    ggrepel::geom_text_repel(
      data = label_df,
      aes(
        x = label_x_adj, y = label_y, label = display_label,
        hjust = label_hjust, size = label_size
      ),
      direction = "y",
      seed = 11,
      force = 1.1,
      force_pull = 0.04,
      box.padding = 0.12,
      point.padding = 0.04,
      min.segment.length = 0,
      segment.size = 0,
      segment.colour = NA,
      max.overlaps = Inf,
      lineheight = 0.78,
      colour = "black",
      family = args$font_family
    )
  } else {
    if (isTRUE(args$use_ggrepel)) {
      warning("ggrepel is not installed; using staggered geom_text labels instead.")
    }
    geom_text(
      data = label_df,
      aes(
        x = label_x_adj, y = label_y, label = display_label,
        hjust = label_hjust, size = label_size
      ),
      lineheight = 0.78,
      colour = "black",
      family = args$font_family
    )
  }

  list(leader_layer, text_layer)
}

lifecycle_antigen_order <- tibble::tribble(
  ~species, ~lifecycle_stage, ~antigen_name, ~stage_order, ~antigen_order,
  "Pf7", "Pre-erythrocytic", "TRAP",     1L,  1L,
  "Pf7", "Pre-erythrocytic", "CSP",      1L,  2L,
  "Pf7", "Pre-erythrocytic", "LSA3",     1L,  3L,
  "Pf7", "Pre-erythrocytic", "CTRP",     1L,  4L,
  "Pf7", "Pre-erythrocytic", "STARP",    1L,  5L,
  "Pf7", "Pre-erythrocytic", "CelTOS",   1L,  6L,
  "Pf7", "Pre-erythrocytic", "TREP",     1L,  7L,
  "Pf7", "Gametocyte",       "Pfs25",    2L,  1L,
  "Pf7", "Gametocyte",       "P25",      2L,  1L,
  "Pf7", "Gametocyte",       "HAP2",     2L,  2L,
  "Pf7", "Gametocyte",       "Pfs230",   2L,  3L,
  "Pf7", "Gametocyte",       "P28",      2L,  4L,
  "Pf7", "Gametocyte",       "Pfs48/45", 2L,  5L,
  "Pf7", "Gametocyte",       "CHT1",     2L,  6L,
  "Pf7", "Blood-stage",      "AMA1",     3L,  1L,
  "Pf7", "Blood-stage",      "MSP6",     3L,  2L,
  "Pf7", "Blood-stage",      "MSP1",     3L,  3L,
  "Pf7", "Blood-stage",      "RALP1",    3L,  4L,
  "Pf7", "Blood-stage",      "RESA",     3L,  5L,
  "Pf7", "Blood-stage",      "RH4",      3L,  6L,
  "Pf7", "Blood-stage",      "RH5",      3L,  7L,
  "Pf7", "Blood-stage",      "RIPR",     3L,  8L,
  "Pf7", "Blood-stage",      "EBA175",   3L,  9L,
  "Pf7", "Blood-stage",      "CLAG8",    3L, 10L,
  "Pf7", "Blood-stage",      "DBLMSP",   3L, 11L,
  "Pf7", "Blood-stage",      "MSP9",     3L, 12L,
  "Pf7", "Blood-stage",      "SERA5",    3L, 13L,
  "Pf7", "Blood-stage",      "SERA8",    3L, 14L,
  "Pf7", "Blood-stage",      "CyRPA",    3L, 15L,
  "Pf7", "Blood-stage",      "MSP4",     3L, 16L,
  "Pf7", "Blood-stage",      "GLURP",    3L, 17L,
  "Pf7", "Blood-stage",      "MSP3",     3L, 18L,
  "Pf7", "Blood-stage",      "MSP2",     3L, 19L,
  "Pf7", "Blood-stage",      "EXP1",     3L, 20L,
  "Pf7", "Blood-stage",      "TRAMP",    3L, 21L,
  "Pv4", "Pre-erythrocytic", "TRAP",     1L,  1L,
  "Pv4", "Pre-erythrocytic", "CSP",      1L,  2L,
  "Pv4", "Gametocyte",       "P28",      2L,  1L,
  "Pv4", "Blood-stage",      "MSP1",     3L,  1L,
  "Pv4", "Blood-stage",      "CyRPA",    3L,  2L,
  "Pv4", "Blood-stage",      "MSP9",     3L,  3L,
  "Pv4", "Blood-stage",      "P41",      3L,  4L,
  "Pv4", "Blood-stage",      "AMA1",     3L,  5L,
  "Pv4", "Blood-stage",      "MSP8",     3L,  6L,
  "Pv4", "Blood-stage",      "RBP2b",    3L,  7L,
  "Pv4", "Blood-stage",      "RBP2a",    3L,  8L,
  "Pv4", "Blood-stage",      "RON2",     3L,  9L,
  "Pv4", "Blood-stage",      "DBP",      3L, 10L,
  "Pv4", "Blood-stage",      "TRAMP",    3L, 11L,
  "Pv4", "Blood-stage",      "RBP1a",    3L, 12L
) %>%
  mutate(antigen_key = norm_key(antigen_name))

add_lifecycle_group <- function(df, enabled = TRUE) {
  if (!isTRUE(enabled)) {
    return(df %>%
      mutate(
        lifecycle_stage = "All antigens",
        stage_order = 1L,
        antigen_order = species_rank
      ))
  }

  df %>%
    mutate(antigen_key = norm_key(antigen)) %>%
    left_join(
      lifecycle_antigen_order %>%
        select(species, antigen_key, lifecycle_stage, stage_order, antigen_order),
      by = c("species", "antigen_key")
    ) %>%
    mutate(
      lifecycle_stage = coalesce(lifecycle_stage, "Other"),
      stage_order = coalesce(stage_order, 4L),
      antigen_order = coalesce(antigen_order, 1000L + species_rank)
    )
}

panel_rank_table <- function(flows_df) {
  flows_df %>%
    mutate(panel = paste(species, antigen, sep = ": ")) %>%
    group_by(panel, species, antigen) %>%
    summarise(
      total = sum(n_samples),
      n_all_snp_haps_raw = first(n_all_snp_haplotypes_raw),
      n_iep_haps_raw = first(n_iep_haplotypes_raw),
      n_major_iep_serotypes = n_distinct(
        immune_escape_haplotype[
          !immune_escape_haplotype %in% c(
            OTHER_SENTINEL,
            FILTERED_IEP_SENTINEL
          )
        ]
      ),
      weighted_iep_density = first(weighted_iep_density),
      weighted_iep_score = first(weighted_iep_score),
      n_significant_iep_sites = first(n_significant_iep_sites),
      n_potential_iep_sites = first(n_potential_iep_sites),
      .groups = "drop"
    ) %>%
    group_by(species) %>%
    arrange(
      desc(replace_na(weighted_iep_density, -Inf)),
      desc(replace_na(weighted_iep_score, -Inf)),
      desc(replace_na(n_significant_iep_sites, -Inf)),
      desc(n_iep_haps_raw),
      antigen,
      .by_group = TRUE
    ) %>%
    mutate(species_rank = row_number()) %>%
    ungroup()
}

filter_top_antigens <- function(flows_df, args) {
  rare_original_stats <- flows_df %>%
    group_by(species, gene_id, antigen) %>%
    summarise(
      rare_original_n = sum(
        n_samples[all_snp_haplotype == OTHER_SENTINEL],
        na.rm = TRUE
      ),
      raw_samples_denominator = first(raw_samples_denominator),
      rare_original_haplotype_percent = if_else(
        is.finite(raw_samples_denominator) & raw_samples_denominator > 0,
        100 * rare_original_n / raw_samples_denominator,
        NA_real_
      ),
      .groups = "drop"
    )

  flows_df <- flows_df %>%
    left_join(
      rare_original_stats,
      by = c("species", "gene_id", "antigen"),
      suffix = c("", "_antigen_filter")
    ) %>%
    filter(
      !is.na(rare_original_haplotype_percent),
      rare_original_haplotype_percent <= args$max_rare_original_haplotype_percent
    )

  top_n <- args$top_n_per_species
  panel_rank <- panel_rank_table(flows_df) %>%
    mutate(antigen_key = norm_key(antigen))

  excluded <- bind_rows(
    tibble(species = "Pf7", antigen_key = norm_key(args$exclude_pf_antigens)),
    tibble(species = "Pv4", antigen_key = norm_key(args$exclude_pv_antigens))
  ) %>%
    filter(nzchar(antigen_key))

  forced <- bind_rows(
    tibble(species = "Pf7", antigen_key = norm_key(args$force_pf_antigens)),
    tibble(species = "Pv4", antigen_key = norm_key(args$force_pv_antigens))
  ) %>%
    filter(nzchar(antigen_key))

  if (nrow(excluded)) {
    panel_rank <- panel_rank %>%
      anti_join(excluded, by = c("species", "antigen_key"))
  }

  available_forced <- panel_rank %>%
    semi_join(forced, by = c("species", "antigen_key"))

  if (nrow(forced)) {
    missing_forced <- forced %>%
      anti_join(panel_rank, by = c("species", "antigen_key"))
    if (nrow(missing_forced)) {
      warning(
        "Forced antigen(s) not available after IEP/haplotype filtering: ",
        paste(paste(missing_forced$species, missing_forced$antigen_key, sep = ":"), collapse = ", ")
      )
    }
  }

  if (is.na(top_n) || top_n <= 0) {
    keep_panels <- panel_rank$panel
  } else {
    forced_counts <- available_forced %>%
      count(species, name = "n_forced")

    base_keep <- panel_rank %>%
      anti_join(available_forced %>% select(species, panel), by = c("species", "panel")) %>%
      left_join(forced_counts, by = "species") %>%
      mutate(n_forced = replace_na(n_forced, 0L),
             keep_n = pmax(top_n - n_forced, 0L)) %>%
      group_by(species) %>%
      arrange(species_rank, .by_group = TRUE) %>%
      filter(row_number() <= first(keep_n)) %>%
      ungroup()

    keep_panels <- unique(c(base_keep$panel, available_forced$panel))
  }

  flows_df %>%
    mutate(panel = paste(species, antigen, sep = ": ")) %>%
    filter(panel %in% keep_panels) %>%
    select(-panel)
}

format_panel_labels <- function(panel_rank, grouped = TRUE) {
  species_display <- c(
    "Pf7" = "Pf",
    "Pv4" = "Pv"
  )
  raw_labels <- paste0(species_display[panel_rank$species], panel_rank$antigen)
  setNames(raw_labels, panel_rank$panel)
}

make_combined_page <- function(flows_df, args) {
  flows_df <- filter_top_antigens(flows_df, args)
  if (!nrow(flows_df)) return(NULL)

  layers <- build_alluvial_layers(flows_df, args)
  panel_rank <- panel_rank_table(layers$flows) %>%
    add_lifecycle_group(args$group_antigens_by_lifecycle) %>%
    arrange(
      factor(species, levels = c("Pf7", "Pv4")),
      stage_order,
      desc(n_major_iep_serotypes),
      desc(n_all_snp_haps_raw),
      antigen
    )
  strip_labels <- format_panel_labels(
    panel_rank,
    args$group_antigens_by_lifecycle
  )

  ribbons <- layers$ribbons %>%
    mutate(panel = factor(strip_labels[panel], levels = strip_labels))
  strata <- layers$strata %>%
    mutate(panel = factor(strip_labels[panel], levels = strip_labels)) %>%
    assign_node_colours() %>%
    prepare_strata_labels(args, n_panels = length(strip_labels))

  ggplot() +
    geom_polygon(
      data = ribbons,
      aes(x = x, y = y, group = flow_id),
      fill = "#C8C8C8",
      alpha = args$ribbon_alpha,
      colour = "#8F8F8F",
      linewidth = 0.04
    ) +
    geom_rect(
      data = strata,
      aes(
        xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = node_fill, linetype = node_linetype
      ),
      colour = "black",
      linewidth = 0.18
    ) +
    make_label_layer(strata, args) +
    annotate("text", x = 0, y = -0.058, label = "All-SNP\nhaplotype",
             size = 1.65 * args$font_scale, family = args$font_family, hjust = 0.5) +
    annotate("text", x = 1, y = -0.058, label = "Major IEP\nhaplotype",
             size = 1.65 * args$font_scale, family = args$font_family, hjust = 0.5) +
    facet_wrap(~panel, ncol = args$facet_ncol) +
    scale_fill_identity() +
    scale_linetype_identity() +
    scale_size_identity() +
    scale_x_continuous(limits = c(-0.32, 1.32), breaks = NULL, expand = expansion(0)) +
    scale_y_continuous(
      limits = c(-0.065, 1.01),
      breaks = c(0, 0.5, 1),
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(0)
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 6 * args$font_scale, base_family = args$font_family) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.y = element_text(size = 4.5 * args$font_scale, colour = "black"),
      axis.title.y = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 5.8 * args$font_scale, face = "bold", colour = "black"),
      panel.spacing.x = unit(0.65, "lines"),
      panel.spacing.y = unit(0.35, "lines"),
      plot.margin = margin(2, 10, 2, 7)
    ) +
    coord_cartesian(clip = "off")
}

# ---------------------------------------------------------------------------
# Plot one species page — Nature-family aesthetic
# ---------------------------------------------------------------------------
make_species_page <- function(flows_df, species_id, args) {

  sp_flows <- flows_df %>% filter(species == species_id)
  if (!nrow(sp_flows)) return(NULL)

  # Build layers
  layers <- build_alluvial_layers(sp_flows, args)

  # Panel order: most distinct raw all-SNP haplotypes first
  # n_all_snp_haplotypes_raw is attached per-flow from process_species
  panel_rank <- layers$flows %>%
    group_by(panel, antigen) %>%
    summarise(
      total                = sum(n_samples),
      n_all_snp_haps_raw   = first(n_all_snp_haplotypes_raw),
      n_iep_haps_raw       = first(n_iep_haplotypes_raw),
      .groups = "drop"
    ) %>%
    arrange(desc(n_all_snp_haps_raw), desc(total), antigen)

  panel_levels <- panel_rank$panel

  # Clean panel label: strip species prefix.
  clean_label <- function(x) sub("^(Pf7|Pv4): ", "", x)
  strip_labels <- setNames(
    clean_label(panel_rank$panel),
    clean_label(panel_rank$panel)
  )

  ribbons <- layers$ribbons %>%
    mutate(panel = factor(
      strip_labels[clean_label(panel)],
      levels = strip_labels
    ))
  strata <- layers$strata %>%
    mutate(panel = factor(
      strip_labels[clean_label(panel)],
      levels = strip_labels
    )) %>%
    assign_node_colours() %>%
    prepare_strata_labels(args, n_panels = length(strip_labels))

  # Species display name
  sp_title <- if (species_id == "Pf7") "Plasmodium falciparum" else "Plasmodium vivax"
  scope_label <- if (args$sample_scope == "global") "Global" else "PNG"

  ggplot() +
    # Neutral ribbons: colour is reserved for original all-SNP haplotype blocks.
    geom_polygon(
      data    = ribbons,
      aes(x = x, y = y, group = flow_id),
      fill    = "#BDBDBD",
      alpha   = 0.42,
      colour  = "#8C8C8C",
      linewidth = 0.10
    ) +
    # Haplotype nodes: left side coloured by all-SNP haplotype; right side white.
    geom_rect(
      data = strata,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = strata$node_fill,
      colour = "black",
      linewidth = 0.25,
      linetype = ifelse(strata$is_other, "dashed", "solid")
    ) +
    # Labels — only shown above threshold; placed outside nodes
    make_label_layer(strata, args) +
    # Column axis headers per facet
    annotate("text", x = 0,   y = -0.09, label = "All-SNP\nhaplotype",
             size = 2.2 * args$font_scale, fontface = "plain", family = args$font_family, hjust = 0.5) +
    annotate("text", x = 1,   y = -0.09, label = "Immune-escape\nhaplotype",
             size = 2.2 * args$font_scale, fontface = "plain", family = args$font_family, hjust = 0.5) +
    facet_wrap(~panel, ncol = args$facet_ncol) +
    scale_size_identity() +
    scale_x_continuous(limits = c(-0.28, 1.28), breaks = NULL, expand = expansion(0)) +
    scale_y_continuous(
      limits = c(-0.14, 1.02),
      breaks = c(0, 0.25, 0.50, 0.75, 1),
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(0)
    ) +
    labs(
      title    = paste0(sp_title,
                        " -- ", scope_label, " haplotype collapse: all-SNP -> immune-escape"),
      subtitle = paste0(
        "Ranked by all-SNP haplotype count (high to low)  |  ",
        "Coloured left blocks = all-SNP haplotypes  |  ",
        "Grey left block = all-SNP haplotypes <", args$min_original_haplotype_percent, "% merged  |  ",
        "Grey right block = IEP haplotypes <", args$min_reduced_haplotype_percent, "% merged"
      ),
      x = NULL,
      y = paste0(scope_label, " sample proportion")
    ) +
    theme_classic(base_size = 7 * args$font_scale, base_family = args$font_family) +
    theme(
      plot.title      = element_text(size = 8 * args$font_scale, face = "italic", hjust = 0,
                                     margin = margin(b = 2)),
      plot.subtitle   = element_text(size = 6 * args$font_scale, face = "plain", colour = "grey40",
                                     hjust = 0, margin = margin(b = 6)),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      axis.line.x     = element_blank(),
      axis.text.y     = element_text(size = 5.5 * args$font_scale, colour = "black"),
      axis.title.y    = element_text(size = 6.5 * args$font_scale, colour = "black"),
      strip.background= element_blank(),
      strip.text      = element_text(size = 7 * args$font_scale, face = "bold", colour = "black"),
      panel.spacing.x = unit(1.2, "lines"),
      panel.spacing.y = unit(0.9, "lines"),
      plot.margin     = margin(6, 24, 6, 10)
    ) +
    coord_cartesian(clip = "off")
}

# ---------------------------------------------------------------------------
# Main — two pages: one per species
# ---------------------------------------------------------------------------
main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  pf <- process_species("Pf7", args$pf_vcf_dir, args$pf_fasta_dir, args$pf_metadata,
                        args$pf_summary, args$pf_snp_col, args$pf_sample_scope, args)
  pv <- process_species("Pv4", args$pv_vcf_dir, args$pv_fasta_dir, args$pv_metadata,
                        args$pv_summary, args$pv_snp_col, args$pv_sample_scope, args)

  summary_df     <- bind_rows(pf$summary,     pv$summary)
  flows_df       <- bind_rows(pf$flows,       pv$flows)
  raw_flows_df   <- bind_rows(pf$raw_flows,   pv$raw_flows)
  diagnostics_df <- bind_rows(pf$diagnostics, pv$diagnostics)

  if (!nrow(summary_df))
    stop("No antigens with immune escape SNP haplotypes were found.")

  plot_flows <- filter_top_antigens(flows_df, args)
  if (!nrow(plot_flows)) {
    stop("No flows remain after top-N antigen filtering.")
  }

  dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)

  p <- make_combined_page(flows_df, args)

  pdf(args$output, width = args$pdf_width, height = args$pdf_height,
      onefile = TRUE, useDingbats = FALSE)
  if (!is.null(p)) print(p)
  dev.off()

  if (!is.null(args$output_svg) && !identical(args$output_svg, "")) {
    dir.create(dirname(args$output_svg), recursive = TRUE, showWarnings = FALSE)
    svg(args$output_svg, width = args$pdf_width, height = args$pdf_height,
        family = args$font_family, onefile = TRUE)
    if (!is.null(p)) print(p)
    dev.off()
  }

  summary_output <- if (is.null(args$summary_output) || identical(args$summary_output, ""))
    paste0(tools::file_path_sans_ext(args$output), "_summary.tsv") else args$summary_output
  flows_output <- if (is.null(args$flows_output) || identical(args$flows_output, ""))
    paste0(tools::file_path_sans_ext(args$output), "_flows.tsv") else args$flows_output
  diagnostics_output <- if (is.null(args$diagnostics_output) || identical(args$diagnostics_output, ""))
    paste0(tools::file_path_sans_ext(args$output), "_diagnostics.tsv") else args$diagnostics_output
  raw_flows_output <- if (is.null(args$raw_flows_output) || identical(args$raw_flows_output, ""))
    paste0(tools::file_path_sans_ext(args$output), "_raw_flows.tsv") else args$raw_flows_output

  selected_panels <- plot_flows %>%
    distinct(species, gene_id, antigen)
  plot_summary <- summary_df %>%
    semi_join(selected_panels, by = c("species", "gene_id", "antigen")) %>%
    left_join(
      plot_flows %>%
        distinct(
          species, gene_id, antigen,
          rare_original_n,
          rare_original_haplotype_percent
        ),
      by = c("species", "gene_id", "antigen")
    )
  plot_raw_flows <- raw_flows_df %>%
    semi_join(selected_panels, by = c("species", "gene_id", "antigen"))

  for (f in c(summary_output, flows_output, raw_flows_output, diagnostics_output))
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)

  write_tsv(plot_summary,   summary_output)
  write_tsv(plot_flows,     flows_output)
  write_tsv(plot_raw_flows, raw_flows_output)
  write_tsv(diagnostics_df, diagnostics_output)

  message(sprintf("Pf7 antigens plotted : %d", n_distinct(plot_flows$antigen[plot_flows$species == "Pf7"])))
  message(sprintf("Pv4 antigens plotted : %d", n_distinct(plot_flows$antigen[plot_flows$species == "Pv4"])))
  message(sprintf("Wrote PDF            : %s", args$output))
  if (!is.null(args$output_svg) && !identical(args$output_svg, "")) {
    message(sprintf("Wrote SVG            : %s", args$output_svg))
  }
  message(sprintf("Wrote summary        : %s", summary_output))
  message(sprintf("Wrote flows          : %s", flows_output))
  message(sprintf("Wrote raw flows      : %s", raw_flows_output))
  message(sprintf("Wrote diagnostics    : %s", diagnostics_output))
}

main()
