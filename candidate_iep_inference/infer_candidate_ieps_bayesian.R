#!/usr/bin/env Rscript
#
# infer_candidate_ieps_bayesian.R
#
# Purpose
# -------
# Detects antigen sites putatively under host-immune selection by comparing,
# within the same host over time, the rate of allelic mismatch between a
# later ("second") sample and an earlier sample, contrasting hosts whose
# later sample was symptomatic vs. asymptomatic. Each site's within-host
# mismatch is modelled with a Bayesian binomial (brms) model and compared
# against the background population-level expected mismatch (heterozygosity)
# at that site, so that sites flagged as candidates are unlikely to be
# explained by generic sequence diversity alone.
#
# Required inputs (all user-supplied, paths/columns are configurable)
# ---------------------------------------------------------------------
# 1. A FASTA alignment of the antigen locus, one sequence per sample, where
#    the sample ID is the part of the FASTA header before the first "_".
# 2. A metadata CSV with, at minimum: a sample ID column, a study/host ID
#    column (so repeat samples from the same host can be paired), a
#    study-day (or timepoint) column, and a binary status column
#    (e.g. symptomatic/asymptomatic, case/control, etc.).
# 3. A multiplicity-of-infection (MOI) or equivalent exposure-covariate CSV
#    keyed by the same study/host ID.
#
# Output files (written to `output_dir`, named "<Antigen>_<piece>.<ext>")
# ------------------------------------------------------------------------
#   <Antigen>_pairwise_comparisons.csv     - all within-host sample pairs
#   <Antigen>_model_data.csv               - aggregated per-site model input
#   <Antigen>_site_summary.csv             - posterior summary per site
#   <Antigen>_candidate_sites.csv          - high-confidence candidates only
#   <Antigen>_model_fit.rds                - saved brms model object
#   <Antigen>_summary_plot.jpeg            - effect-size summary plot
#
# Dependencies
# ------------
# install.packages(c("dplyr", "tidyr", "readr", "tibble", "seqinr",
#                     "stringr", "ggplot2", "ggrepel", "brms"))
# brms additionally requires a working Stan toolchain (rstan or cmdstanr).
#
# Usage
# -----
# source("infer_candidate_ieps_bayesian.R")
# result <- allele_specific_immunity_bayesian(
#   data_dir      = "path/to/data",
#   fasta_file    = "antigen_alignment.fasta",
#   metadata_file = "sample_metadata.csv",
#   molfoi_file   = "moi_estimates.csv",
#   Antigen       = "MyAntigen"
# )
#
# Notes
# -----
# This script does not change the R session's working directory; all file
# paths are resolved relative to `data_dir` / `output_dir`, which you supply.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(seqinr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(brms)
})

#' Elementwise mismatch indicator between two rows of allele calls
#'
#' @param a,b Data frames/vectors of equal length holding allele calls.
#' @return Integer vector of 0/1 (0 = match, 1 = mismatch).
compare_binary_mismatch <- function(a, b) {
  a_chr <- as.character(a)
  b_chr <- as.character(b)
  as.integer(a_chr != b_chr)
}

#' Build all pairwise within-host sample comparisons
#'
#' For every host ("studyid") with more than one sample, enumerates all
#' pairwise combinations of its samples and computes a per-site mismatch
#' indicator for each pair.
#'
#' @param seq_data Data frame of merged sequence + metadata, one row per
#'   sample, row names = sample IDs, must contain a `studyid` and `Status`
#'   column plus the columns named in `seq_cols`.
#' @param seq_cols Character vector of column names holding per-site allele
#'   calls.
#' @return A tibble with one row per sample pair and one column per site
#'   holding the mismatch indicator.
build_within_host_pairs <- function(seq_data, seq_cols) {
  sample_freq <- as.data.frame(table(seq_data$studyid), stringsAsFactors = FALSE)
  sample_freq <- sample_freq[sample_freq$Freq > 1, , drop = FALSE]
  if (!all(seq_cols %in% colnames(seq_data))) {
    missing_cols <- seq_cols[!seq_cols %in% colnames(seq_data)]
    stop("Sequence columns missing from paired data: ", paste(missing_cols, collapse = ", "))
  }

  pair_rows <- list()
  k <- 0L

  for (i in seq_len(nrow(sample_freq))) {
    study_id <- sample_freq[i, 1]
    sub <- subset(seq_data, seq_data$studyid == study_id)
    for (m in seq_len(nrow(sub) - 1L)) {
      samp1 <- sub[m, ]
      for (n in (m + 1L):nrow(sub)) {
        samp2 <- sub[n, ]
        k <- k + 1L
        mismatch <- compare_binary_mismatch(samp1[, seq_cols, drop = FALSE], samp2[, seq_cols, drop = FALSE])
        pair_rows[[k]] <- tibble(
          Comparison = paste0("Comp_", k),
          studyid = samp1$studyid,
          Samp_1_SampleID = rownames(samp1),
          Samp_2_SampleID = rownames(samp2),
          Samp_1_Status = as.character(samp1$Status),
          Samp_2_Status = as.character(samp2$Status),
          Transition = paste0(samp1$Status, "_", samp2$Status)
        ) %>%
          bind_cols(as_tibble(as.list(stats::setNames(mismatch, seq_cols))))
      }
    }
  }

  bind_rows(pair_rows)
}

#' Standardise a raw status column to "A" (control/asymptomatic-like) or
#' "S" (case/symptomatic-like)
#'
#' Recognises common spellings/abbreviations. Adjust the mapping here if
#' your status labels differ (e.g. "case"/"control", "positive"/"negative").
#'
#' @param x Character or factor vector of raw status labels.
#' @return Character vector of "A"/"S" (other values pass through unchanged
#'   and will be filtered out downstream).
standardise_status <- function(x) {
  out <- toupper(trimws(as.character(x)))
  out[out %in% c("ASYM", "ASYMPTOMATIC", "A")] <- "A"
  out[out %in% c("SYM", "SYMPTOMATIC", "S")] <- "S"
  out
}

#' Compute population-level expected mismatch (expected heterozygosity) at
#' each site
#'
#' @param seq_data Data frame containing one column per site.
#' @param polymorphic_sites Character vector of site column names.
#' @return A tibble with one row per site summarising allele diversity.
compute_site_expected_mismatch <- function(seq_data, polymorphic_sites) {
  bind_rows(lapply(polymorphic_sites, function(site) {
    vals <- as.character(seq_data[[site]])
    vals <- vals[!is.na(vals) & vals != "" & vals != "-"]

    if (length(vals) == 0) {
      return(tibble(
        site = site,
        expected_mismatch = NA_real_,
        callable_n = 0L,
        n_alleles = 0L,
        dominant_allele_freq = NA_real_
      ))
    }

    freqs <- prop.table(table(vals))
    tibble(
      site = site,
      expected_mismatch = 1 - sum(freqs^2),
      callable_n = length(vals),
      n_alleles = length(freqs),
      dominant_allele_freq = max(as.numeric(freqs))
    )
  }))
}

#' Load and merge the FASTA alignment, metadata, and MOI inputs
#'
#' @param data_dir Directory containing the input files.
#' @param source_path Optional path to an R script to `source()` before
#'   processing (e.g. for custom helper functions). Left NULL by default.
#' @param fasta_file FASTA alignment file name (relative to `data_dir`).
#' @param metadata_file Metadata CSV file name (relative to `data_dir`).
#' @param molfoi_file MOI CSV file name (relative to `data_dir`).
#' @param sample_col,studyid_col,timepoint_col,status_col,molfoi_col Column
#'   names in the respective input files; override if your files use
#'   different headers.
#' @return A list with the merged sequence data, metadata+MOI table,
#'   polymorphic site names, and per-site expected mismatch summary.
prepare_inputs <- function(
  data_dir,
  source_path = NULL,
  fasta_file,
  metadata_file = "sample_metadata.csv",
  molfoi_file = "moi_estimates.csv",
  sample_col = "samples",
  studyid_col = "studyid",
  timepoint_col = "timepoint",
  status_col = "status",
  molfoi_col = "molFOI"
) {
  set.seed(11111)

  if (!is.null(source_path) && nzchar(source_path) && file.exists(source_path)) {
    source(source_path)
  }

  fasta_path <- file.path(data_dir, fasta_file)
  metadata_path <- file.path(data_dir, metadata_file)
  molfoi_path <- file.path(data_dir, molfoi_file)

  fasta_seqs <- read.fasta(fasta_path, forceDNAtolower = FALSE)
  data <- as.data.frame(do.call(rbind, fasta_seqs))
  positions <- paste0("pos", seq_len(ncol(data)))
  colnames(data) <- positions

  data <- data %>% mutate(across(everything(), as.factor))
  data <- data[, sapply(data, nlevels) > 1, drop = FALSE]
  data <- tibble::rownames_to_column(data, "seq")
  data <- data.frame(lapply(data, function(x) gsub(".", "-", x, fixed = TRUE)))
  data <- data %>% mutate(samples = gsub("_.*$", "", seq)) %>% select(samples, everything())

  epi_data_raw <- read.csv(metadata_path, stringsAsFactors = FALSE)
  required_metadata <- c(sample_col, studyid_col, timepoint_col, status_col)
  if (!all(required_metadata %in% names(epi_data_raw))) {
    stop(
      "Metadata file must contain columns: ",
      paste(required_metadata, collapse = ", ")
    )
  }

  epi_data <- epi_data_raw %>%
    transmute(
      studyid = .data[[studyid_col]],
      samples = .data[[sample_col]],
      timepoint = .data[[timepoint_col]],
      Status = standardise_status(.data[[status_col]])
    ) %>%
    filter(!is.na(samples), !is.na(studyid), !is.na(timepoint), Status %in% c("A", "S"))

  allele_specific_seq <- merge(epi_data, data, by = "samples")
  if ("seq" %in% names(allele_specific_seq)) {
    seq_data <- allele_specific_seq %>% select(-seq)
  } else {
    seq_data <- allele_specific_seq
  }

  molFOI_raw <- read.csv(molfoi_path, stringsAsFactors = FALSE)
  if (!all(c(studyid_col, molfoi_col) %in% names(molFOI_raw))) {
    stop("MOI file must contain columns: ", studyid_col, ", ", molfoi_col)
  }
  molFOI <- molFOI_raw %>%
    transmute(
      studyid = .data[[studyid_col]],
      molFOI = suppressWarnings(as.numeric(.data[[molfoi_col]]))
    ) %>%
    filter(!is.na(studyid), is.finite(molFOI))

  epi_data_FOI <- merge(epi_data, molFOI, by = "studyid", all.x = TRUE)
  epi_data_FOI[is.na(epi_data_FOI)] <- 0

  seq_data <- distinct(seq_data, samples, .keep_all = TRUE)
  rownames(seq_data) <- seq_data$samples
  seq_data <- seq_data %>% select(-samples)
  seq_data <- seq_data[order(seq_data$timepoint), ]
  seq_data <- seq_data[order(seq_data$studyid), ]

  polymorphic_sites <- intersect(colnames(data[-c(1, 2)]), colnames(seq_data))
  if (length(polymorphic_sites) == 0) {
    stop("No polymorphic sequence columns remained after metadata merge")
  }
  expected_site_rates <- compute_site_expected_mismatch(seq_data, polymorphic_sites)

  list(
    seq_data = seq_data,
    epi_data_FOI = epi_data_FOI,
    polymorphic_sites = polymorphic_sites,
    expected_site_rates = expected_site_rates
  )
}

#' Reshape pairwise comparisons into a per-site binomial dataset for brms
#'
#' @param within_comp Output of `build_within_host_pairs()`.
#' @param epi_data_FOI Metadata + MOI table from `prepare_inputs()`.
#' @param polymorphic_sites Character vector of site column names.
#' @param expected_site_rates Per-site expected mismatch summary.
#' @return A list with the long-format pairwise data and the aggregated
#'   model-ready data frame.
build_bayesian_site_dataset <- function(within_comp, epi_data_FOI, polymorphic_sites, expected_site_rates) {
  long_pairs <- within_comp %>%
    rename(samples = Samp_1_SampleID) %>%
    pivot_longer(cols = all_of(polymorphic_sites), names_to = "site", values_to = "mismatch") %>%
    mutate(
      mismatch = as.integer(mismatch),
      second_status = ifelse(as.character(Samp_2_Status) == "S", 1L, 0L)
    ) %>%
    left_join(
      epi_data_FOI %>% select(studyid, molFOI) %>% distinct(),
      by = "studyid"
    ) %>%
    left_join(expected_site_rates, by = "site") %>%
    mutate(
      molFOI = ifelse(is.na(molFOI), 0, as.numeric(molFOI)),
      expected_mismatch = ifelse(is.na(expected_mismatch), 0, expected_mismatch),
      expected_mismatch_clipped = pmin(pmax(expected_mismatch, 1e-6), 1 - 1e-6),
      offset_logit_expected = qlogis(expected_mismatch_clipped)
    )

  model_data <- long_pairs %>%
    group_by(
      studyid,
      site,
      second_status,
      molFOI,
      expected_mismatch,
      expected_mismatch_clipped,
      offset_logit_expected,
      callable_n,
      n_alleles,
      dominant_allele_freq
    ) %>%
    summarise(
      event = sum(mismatch),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      site = factor(site),
      studyid = factor(studyid),
      molFOI_scaled = as.numeric(scale(molFOI))
    )

  list(long_pairs = long_pairs, model_data = model_data)
}

#' Fit the Bayesian binomial mismatch model
#'
#' @param model_data Output `model_data` from `build_bayesian_site_dataset()`.
#' @param draws,iter,warmup,chains,cores,adapt_delta,max_treedepth Standard
#'   brms/Stan sampling controls.
#' @param background_adjustment Either "posthoc" (default; background
#'   diversity is reported/compared after fitting) or "offset" (background
#'   diversity is included as a fixed offset in the linear predictor).
#' @return A fitted `brmsfit` object.
fit_allele_escape_model <- function(
  model_data,
  draws,
  iter,
  warmup,
  chains,
  cores,
  adapt_delta,
  max_treedepth,
  background_adjustment = c("posthoc", "offset")
) {
  background_adjustment <- match.arg(background_adjustment)
  priors <- c(
    set_prior("normal(0, 1.5)", class = "Intercept"),
    set_prior("normal(0, 1)", class = "b"),
    set_prior("student_t(3, 0, 1.5)", class = "sd")
  )

  model_formula <- if (background_adjustment == "offset") {
    bf(
      event | trials(n) ~ 1 + second_status + molFOI_scaled +
        offset(offset_logit_expected) +
        (1 | site) + (0 + second_status | site) + (1 | studyid)
    )
  } else {
    bf(
      event | trials(n) ~ 1 + second_status + molFOI_scaled +
        (1 | site) + (0 + second_status | site) + (1 | studyid)
    )
  }

  brm(
    formula = model_formula,
    data = model_data,
    family = stats::binomial(),
    prior = priors,
    iter = iter,
    warmup = warmup,
    chains = chains,
    cores = cores,
    seed = 1,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    refresh = 0
  )
}

#' Summarise posterior mismatch probabilities per site
#'
#' @param fit Fitted `brmsfit` from `fit_allele_escape_model()`.
#' @param model_data Model data used to fit `fit`.
#' @param draws Number of posterior draws to summarise over.
#' @param background_adjustment Must match the value used to fit `fit`.
#' @return A tibble with one row per site summarising posterior mismatch
#'   probabilities for the two status groups and their contrast.
summarise_site_posteriors <- function(fit, model_data, draws, background_adjustment = c("posthoc", "offset")) {
  background_adjustment <- match.arg(background_adjustment)
  sites <- levels(model_data$site)
  ref_studyid <- levels(model_data$studyid)[1]
  site_baseline <- model_data %>%
    select(site, expected_mismatch, expected_mismatch_clipped, offset_logit_expected, callable_n, n_alleles, dominant_allele_freq) %>%
    distinct()
  newdata <- tidyr::expand_grid(
    site = factor(sites, levels = sites),
    second_status = c(0L, 1L)
  ) %>%
    left_join(site_baseline, by = "site") %>%
    mutate(
      molFOI_scaled = 0,
      studyid = factor(ref_studyid, levels = levels(model_data$studyid)),
      n = 1L
    )

  posterior <- posterior_epred(
    fit,
    newdata = newdata,
    re_formula = ~(1 | site) + (0 + second_status | site),
    ndraws = draws
  )

  if (length(dim(posterior)) == 2) {
    post_mat <- posterior
  } else {
    post_mat <- posterior[, , 1, drop = FALSE]
    post_mat <- matrix(post_mat, nrow = dim(posterior)[1], ncol = dim(posterior)[2])
  }

  summary_rows <- lapply(seq_along(sites), function(i) {
    idx_asym <- (2 * i) - 1
    idx_sym <- 2 * i
    asym_vals <- post_mat[, idx_asym]
    sym_vals <- post_mat[, idx_sym]
    log_ratio_vals <- log2((sym_vals + 1e-6) / (asym_vals + 1e-6))
    expected_val <- site_baseline$expected_mismatch[match(sites[i], site_baseline$site)]
    asym_excess <- asym_vals - expected_val
    sym_excess <- sym_vals - expected_val

    tibble(
      site = sites[i],
      expected_mismatch = expected_val,
      callable_n = site_baseline$callable_n[match(sites[i], site_baseline$site)],
      n_alleles = site_baseline$n_alleles[match(sites[i], site_baseline$site)],
      dominant_allele_freq = site_baseline$dominant_allele_freq[match(sites[i], site_baseline$site)],
      asym_mean = mean(asym_vals),
      asym_lower_95 = quantile(asym_vals, 0.025),
      asym_upper_95 = quantile(asym_vals, 0.975),
      sym_mean = mean(sym_vals),
      sym_lower_95 = quantile(sym_vals, 0.025),
      sym_upper_95 = quantile(sym_vals, 0.975),
      diff_mean = mean(sym_vals - asym_vals),
      diff_lower_95 = quantile(sym_vals - asym_vals, 0.025),
      diff_upper_95 = quantile(sym_vals - asym_vals, 0.975),
      asym_excess_mean = mean(asym_excess),
      asym_excess_lower_95 = quantile(asym_excess, 0.025),
      asym_excess_upper_95 = quantile(asym_excess, 0.975),
      sym_excess_mean = mean(sym_excess),
      sym_excess_lower_95 = quantile(sym_excess, 0.025),
      sym_excess_upper_95 = quantile(sym_excess, 0.975),
      log2_ratio_mean = mean(log_ratio_vals),
      prob_sym_gt_asym = mean(sym_vals > asym_vals),
      prob_sym_gt_expected = mean(sym_vals > expected_val),
      prob_asym_gt_expected = mean(asym_vals > expected_val),
      background_adjustment = background_adjustment
    )
  })

  bind_rows(summary_rows)
}

#' Plot the site-level summary (effect size vs. background-adjusted excess)
#'
#' @param site_summary Output of the main pipeline's site summary table,
#'   must contain a `high_confidence_candidate` logical column.
#' @param antigen_name Label used in the plot title (name of the antigen /
#'   locus being analysed).
#' @param background_adjustment Must match the value used to build
#'   `site_summary`.
#' @return A ggplot object.
plot_site_summary <- function(site_summary, antigen_name, background_adjustment = c("posthoc", "offset")) {
  background_adjustment <- match.arg(background_adjustment)
  highlight_df <- site_summary %>%
    filter(high_confidence_candidate)

  if (background_adjustment == "offset") {
    ggplot(site_summary, aes(x = log2_ratio_mean, y = sym_excess_mean)) +
      geom_point(color = "deepskyblue3", alpha = 0.55) +
      geom_vline(xintercept = 0, color = "gray60") +
      geom_hline(yintercept = 0, color = "gray60") +
      geom_point(data = highlight_df, color = "red4", size = 3) +
      geom_text_repel(aes(label = site), size = 3.2) +
      theme_classic() +
      labs(
        title = paste("Bayesian allele-specific immunity model for", antigen_name),
        subtitle = "Hard background-diversity correction in the likelihood using site-specific expected mismatch",
        x = "Posterior mean log2(symptomatic mismatch probability / asymptomatic mismatch probability)",
        y = "Posterior mean symptomatic excess over expected mismatch"
      )
  } else {
    ggplot(site_summary, aes(x = log2_ratio_mean, y = diff_mean)) +
      geom_point(color = "deepskyblue3", alpha = 0.55) +
      geom_vline(xintercept = 0, color = "gray60") +
      geom_hline(yintercept = 0, color = "gray60") +
      geom_point(data = highlight_df, color = "red4", size = 3) +
      geom_text_repel(aes(label = site), size = 3.2) +
      theme_classic() +
      labs(
        title = paste("Bayesian allele-specific immunity model for", antigen_name),
        subtitle = "Primary ranking uses within-host symptomatic vs asymptomatic contrast; background mismatch is reported post hoc",
        x = "Posterior mean log2(symptomatic mismatch probability / asymptomatic mismatch probability)",
        y = "Posterior mean symptomatic minus asymptomatic mismatch probability"
      )
  }
}

#' Bayesian allele-specific immune-escape analysis pipeline
#'
#' End-to-end pipeline: loads a FASTA alignment and matching metadata/MOI
#' tables, builds within-host pairwise comparisons, fits a Bayesian
#' binomial mismatch model per site, ranks sites by evidence of elevated
#' mismatch in symptomatic (vs. asymptomatic) within-host transitions after
#' accounting for background population diversity, and writes all
#' intermediate/output tables, the fitted model, and a summary plot to
#' `output_dir`.
#'
#' @param data_dir Directory containing the input FASTA/metadata/MOI files.
#' @param output_dir Directory to write output files to. Defaults to
#'   `data_dir`. Created if it doesn't already exist.
#' @param source_path Optional helper R script to source before processing.
#' @param fasta_file FASTA alignment file name (relative to `data_dir`).
#' @param Antigen Character label for the antigen/locus being analysed;
#'   used in output file names and plot titles.
#' @param metadata_file,molfoi_file Input CSV file names (relative to
#'   `data_dir`).
#' @param sample_col,studyid_col,timepoint_col,status_col,molfoi_col Column
#'   names in the input files; override if your files use different
#'   headers.
#' @param draws Number of posterior draws used when summarising the model.
#' @param iter,warmup,chains,cores,adapt_delta,max_treedepth Stan/brms
#'   sampling controls, passed to `brms::brm()`.
#' @param min_callable_n,min_event_s,min_expected_mismatch,min_prob_sym_gt_asym
#'   Thresholds used to flag "high-confidence" candidate sites.
#' @param background_adjustment Either "posthoc" (default) or "offset"; see
#'   `fit_allele_escape_model()`.
#' @return (Invisibly) a list with the pairwise input data, model input
#'   data, site-level summary table, the fitted `brmsfit`, and the summary
#'   ggplot object. Also writes the following files to `output_dir`:
#'   `<Antigen>_<date>_pairwise_comparisons.csv`,
#'   `<Antigen>_<date>_model_data.csv`,
#'   `<Antigen>_<date>_site_summary.csv`,
#'   `<Antigen>_<date>_candidate_sites.csv`,
#'   `<Antigen>_<date>_model_fit.rds`, and
#'   `<Antigen>_<date>_summary_plot.jpeg`.
allele_specific_immunity_bayesian <- function(
  data_dir,
  output_dir = data_dir,
  source_path = NULL,
  fasta_file,
  Antigen,
  metadata_file = "sample_metadata.csv",
  molfoi_file = "moi_estimates.csv",
  sample_col = "samples",
  studyid_col = "studyid",
  timepoint_col = "timepoint",
  status_col = "status",
  molfoi_col = "molFOI",
  draws = 2000,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  adapt_delta = 0.99,
  max_treedepth = 15,
  min_callable_n = 20,
  min_event_s = 3,
  min_expected_mismatch = 0.01,
  min_prob_sym_gt_asym = 0.95,
  background_adjustment = c("posthoc", "offset")
) {
  background_adjustment <- match.arg(background_adjustment)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  # All output files share this prefix, e.g. "MyAntigen_2026-07-08_site_summary.csv"
  file_prefix <- paste0(Antigen, "_", Sys.Date())
  out_path <- function(suffix) file.path(output_dir, paste0(file_prefix, suffix))

  prepared <- prepare_inputs(
    data_dir = data_dir,
    source_path = source_path,
    fasta_file = fasta_file,
    metadata_file = metadata_file,
    molfoi_file = molfoi_file,
    sample_col = sample_col,
    studyid_col = studyid_col,
    timepoint_col = timepoint_col,
    status_col = status_col,
    molfoi_col = molfoi_col
  )

  within_comp <- build_within_host_pairs(prepared$seq_data, seq_cols = prepared$polymorphic_sites)
  write.csv(within_comp, out_path("_pairwise_comparisons.csv"), row.names = FALSE)

  model_inputs <- build_bayesian_site_dataset(
    within_comp,
    prepared$epi_data_FOI,
    prepared$polymorphic_sites,
    prepared$expected_site_rates
  )
  write.csv(model_inputs$model_data, out_path("_model_data.csv"), row.names = FALSE)

  fit <- fit_allele_escape_model(
    model_inputs$model_data,
    draws = draws,
    iter = iter,
    warmup = warmup,
    chains = chains,
    cores = cores,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    background_adjustment = background_adjustment
  )

  site_summary <- summarise_site_posteriors(
    fit,
    model_inputs$model_data,
    draws,
    background_adjustment = background_adjustment
  )

  raw_counts <- model_inputs$model_data %>%
    group_by(site, second_status) %>%
    summarise(event = sum(event), n = sum(n), .groups = "drop") %>%
    mutate(status_label = ifelse(second_status == 1L, "S", "A")) %>%
    select(-second_status) %>%
    pivot_wider(names_from = status_label, values_from = c(event, n), names_sep = "_")

  site_summary <- site_summary %>%
    left_join(raw_counts, by = "site") %>%
    mutate(
      event_A = ifelse(is.na(event_A), 0, event_A),
      event_S = ifelse(is.na(event_S), 0, event_S),
      n_A = ifelse(is.na(n_A), 0, n_A),
      n_S = ifelse(is.na(n_S), 0, n_S),
      site_index = as.numeric(gsub("pos", "", site)),
      ranking_effect = ifelse(background_adjustment == "offset", sym_excess_mean, diff_mean),
      immune_escape_flag = prob_sym_gt_asym >= min_prob_sym_gt_asym & ranking_effect > 0,
      high_confidence_candidate = prob_sym_gt_asym >= min_prob_sym_gt_asym &
        ranking_effect > 0 &
        callable_n >= min_callable_n &
        event_S >= min_event_s &
        expected_mismatch >= min_expected_mismatch
    ) %>%
    arrange(desc(prob_sym_gt_asym), desc(ranking_effect), site_index)

  write.csv(site_summary, out_path("_site_summary.csv"), row.names = FALSE)
  write.csv(
    site_summary %>%
      filter(high_confidence_candidate) %>%
      arrange(desc(prob_sym_gt_asym), desc(ranking_effect), site_index),
    out_path("_candidate_sites.csv"),
    row.names = FALSE
  )
  saveRDS(fit, out_path("_model_fit.rds"))

  plot_obj <- plot_site_summary(site_summary, Antigen, background_adjustment = background_adjustment)
  ggsave(
    filename = out_path("_summary_plot.jpeg"),
    plot = plot_obj,
    width = 15,
    height = 12,
    units = "in",
    dpi = 500
  )

  invisible(
    list(
      pairwise_input = within_comp,
      model_input = model_inputs$model_data,
      site_summary = site_summary,
      fit = fit,
      plot = plot_obj
    )
  )
}
