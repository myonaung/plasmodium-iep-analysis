#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

parse_args <- function(args) {
  out <- list(
    pf_summary = "",
    pv_summary = "",
    pf_snp_col = "pf7_ns_snp",
    pv_snp_col = "Pv4_snps",
    prob_sym_gt_asym_threshold = 0.80,
    prob_sym_gt_expected_threshold = 0,
    require_prob_sym_gt_expected_over_asym = "TRUE",
    dominant_allele_freq_max = 0.99,
    out = "candidate_iep_position_map_psa080.tsv"
  )

  i <- 1L
  while (i <= length(args)) {
    if (!startsWith(args[[i]], "--") || i == length(args)) {
      stop("Bad argument: ", args[[i]])
    }
    key <- gsub("-", "_", sub("^--", "", args[[i]]))
    if (!key %in% names(out)) stop("Unknown option: --", key)
    out[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }

  out$prob_sym_gt_asym_threshold <- as.numeric(out$prob_sym_gt_asym_threshold)
  out$prob_sym_gt_expected_threshold <- as.numeric(out$prob_sym_gt_expected_threshold)
  out$dominant_allele_freq_max <- as.numeric(out$dominant_allele_freq_max)
  out$require_prob_sym_gt_expected_over_asym <-
    tolower(out$require_prob_sym_gt_expected_over_asym) %in%
    c("true", "t", "1", "yes", "y")
  if (!nzchar(out$pf_summary)) stop("--pf-summary is required.")
  if (!nzchar(out$pv_summary)) stop("--pv-summary is required.")
  out
}

num <- function(x) suppressWarnings(as.numeric(x))

pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) hit[[1]] else NA_character_
}

make_map <- function(path, species_code, species_label, snp_col, args) {
  df <- read_csv(path, show_col_types = FALSE, guess_max = 300000)
  if (!"position" %in% names(df)) stop("Summary missing position column: ", path)
  if (!snp_col %in% names(df)) stop("SNP column not found: ", snp_col)

  antigen_col <- pick_col(df, c("antigen", "gene_name", "Gene_Name_annotated", "gene_id"))
  gene_col <- pick_col(df, c("gene_id", "Gene_ID", "Gene_Name_annotated", "gene_name"))
  if (is.na(antigen_col)) stop("Summary has no antigen/gene_name column: ", path)
  if (is.na(gene_col)) stop("Summary has no gene_id/gene_name column: ", path)

  psa <- num(df$prob_sym_gt_asym)
  pse <- num(df$prob_sym_gt_expected)
  pae <- num(df$prob_asym_gt_expected)
  snp_ok <- !is.na(num(df[[snp_col]])) & num(df[[snp_col]]) == 1
  modeled <- if ("has_bayesian_summary" %in% names(df)) {
    tolower(as.character(df$has_bayesian_summary)) %in% c("true", "t", "1", "yes", "y")
  } else {
    !is.na(psa) | !is.na(pse)
  }

  iep <- modeled &
    snp_ok &
    !is.na(psa) & psa >= args$prob_sym_gt_asym_threshold &
    !is.na(pse) & pse >= args$prob_sym_gt_expected_threshold

  if (args$require_prob_sym_gt_expected_over_asym) {
    iep <- iep & !is.na(pae) & pse > pae
  }

  if ("dominant_allele_freq" %in% names(df)) {
    daf <- num(df$dominant_allele_freq)
    iep <- iep & (is.na(daf) | daf <= args$dominant_allele_freq_max)
  }
  iep[is.na(iep)] <- FALSE

  tibble(
    species = species_code,
    Species = species_label,
    antigen = as.character(df[[antigen_col]]),
    gene_id = as.character(df[[gene_col]]),
    aa_position = as.integer(num(df$position))
  ) %>%
    filter(iep, !is.na(aa_position)) %>%
    distinct(species, Species, antigen, gene_id, aa_position) %>%
    arrange(species, gene_id, antigen, aa_position)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  out <- bind_rows(
    make_map(args$pf_summary, "Pf7", "P. falciparum", args$pf_snp_col, args),
    make_map(args$pv_summary, "Pv4", "P. vivax", args$pv_snp_col, args)
  )
  write_tsv(out, args$out)

  summary <- out %>%
    count(species, Species, gene_id, antigen, name = "candidate_iep_positions_n") %>%
    arrange(species, gene_id, antigen)
  write_tsv(summary, sub("\\.tsv$", ".summary.tsv", args$out))

  message("Wrote: ", args$out)
  message("Wrote: ", sub("\\.tsv$", ".summary.tsv", args$out))
}

main()
