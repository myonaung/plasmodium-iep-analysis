#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(brms)
})

usage <- function() {
  cat(
    paste(
      "Usage:",
      "Rscript scripts/estimate_global_serotype_representation_bayesian.R",
      "--calls sample_calls.tsv",
      "--outdir results",
      "[--metadata sample_metadata.tsv]",
      "[--clusters csp_sample_clusters.tsv]",
      "[--top-n 20]",
      "[--draws 2000]",
      "[--iter 4000]",
      "[--warmup 2000]",
      "[--chains 4]",
      "[--cores 4]",
      "[--adapt-delta 0.99]",
      "[--max-treedepth 15]",
      sep = " "
    ),
    "\n"
  )
  quit(status = 1)
}

parse_args_simple <- function(args) {
  if (length(args) == 0) usage()

  out <- list(
    calls = NULL,
    metadata = NULL,
    clusters = NULL,
    outdir = NULL,
    top_n = 20L,
    draws = 2000L,
    iter = 4000L,
    warmup = 2000L,
    chains = 4L,
    cores = 4L,
    adapt_delta = 0.99,
    max_treedepth = 15L
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (!startsWith(arg, "--") || i == length(args)) stop("Bad arguments.")
    key <- gsub("-", "_", sub("^--", "", arg))
    value <- args[[i + 1L]]
    i <- i + 2L
    if (!key %in% names(out)) stop("Unknown option: --", gsub("_", "-", key))
    out[[key]] <- value
  }

  out$top_n <- as.integer(out$top_n)
  out$draws <- as.integer(out$draws)
  out$iter <- as.integer(out$iter)
  out$warmup <- as.integer(out$warmup)
  out$chains <- as.integer(out$chains)
  out$cores <- as.integer(out$cores)
  out$adapt_delta <- as.numeric(out$adapt_delta)
  out$max_treedepth <- as.integer(out$max_treedepth)

  if (is.null(out$calls) || is.null(out$outdir)) usage()
  out
}

prepare_metadata <- function(assignments, metadata_path = NULL, clusters_path = NULL) {
  metadata <- tibble(sample_id = assignments$sample_id)

  if (!is.null(metadata_path)) {
    metadata <- read_tsv(metadata_path, show_col_types = FALSE)
    if (!"sample_id" %in% names(metadata)) {
      stop("Metadata file must contain sample_id column.")
    }
  }

  assignments <- assignments %>%
    mutate(
      assigned_serotype = as.character(assigned_serotype),
      call_status = as.character(call_status)
    )

  merged <- assignments %>%
    left_join(metadata, by = "sample_id") %>%
    mutate(
      site = if ("site" %in% names(.)) as.character(site) else "GLOBAL_SITE",
      region = if ("region" %in% names(.)) as.character(region) else "GLOBAL"
    )

  merged$site[is.na(merged$site) | merged$site == ""] <- "GLOBAL_SITE"
  merged$region[is.na(merged$region) | merged$region == ""] <- "GLOBAL"

  if (!is.null(clusters_path)) {
    clusters <- read_tsv(clusters_path, show_col_types = FALSE)
    if (!all(c("sample_id", "cluster_id") %in% names(clusters))) {
      stop("Cluster file must contain sample_id and cluster_id columns.")
    }
    clusters <- clusters %>% select(sample_id, cluster_id) %>% distinct(sample_id, .keep_all = TRUE)
    merged <- merged %>% left_join(clusters, by = "sample_id")
  }

  merged %>% mutate(cluster_id = if ("cluster_id" %in% names(.)) as.character(cluster_id) else NA_character_)
}

build_raw_unit_table <- function(sample_data) {
  sample_data %>%
    transmute(
      unit_id = sample_id,
      sample_id = sample_id,
      site = site,
      region = region,
      assigned_serotype = as.character(assigned_serotype),
      call_status = as.character(call_status)
    )
}

build_cluster_unit_table <- function(sample_data) {
  if (all(is.na(sample_data$cluster_id))) return(NULL)

  sample_data %>%
    mutate(cluster_key = ifelse(is.na(cluster_id) | cluster_id == "", paste0("singleton:", sample_id), cluster_id)) %>%
    group_by(site, region, cluster_key) %>%
    summarise(
      sample_count = n(),
      called_serotypes = list(sort(unique(as.character(assigned_serotype[call_status == "CALLED"])))),
      any_mixed = any(call_status %in% c("MIXED", "AMBIGUOUS")),
      any_unresolved = any(call_status == "UNRESOLVED"),
      .groups = "drop"
    ) %>%
    mutate(
      assigned_serotype = case_when(
        lengths(called_serotypes) == 1L & !any_mixed ~ vapply(called_serotypes, `[`, character(1), 1L),
        lengths(called_serotypes) > 1L ~ "MIXED_CLUSTER",
        any_unresolved ~ "UNRESOLVED_CLUSTER",
        TRUE ~ "OTHER_CLUSTER"
      ),
      call_status = ifelse(assigned_serotype %in% c("MIXED_CLUSTER", "UNRESOLVED_CLUSTER", "OTHER_CLUSTER"), assigned_serotype, "CALLED"),
      unit_id = paste(site, cluster_key, sep = "::")
    ) %>%
    select(unit_id, site, region, assigned_serotype, call_status, sample_count)
}

dirichlet_draws <- function(category_counts, n_draws) {
  cats <- names(category_counts)
  alpha <- as.numeric(category_counts) + 1
  mat <- matrix(NA_real_, nrow = n_draws, ncol = length(alpha), dimnames = list(NULL, cats))
  for (i in seq_len(n_draws)) {
    g <- rgamma(length(alpha), shape = alpha, rate = 1)
    mat[i, ] <- g / sum(g)
  }
  mat
}

get_dirichlet_result <- function(unit_table, serotypes, n_draws, label = "raw_dirichlet") {
  called <- unit_table %>% filter(call_status == "CALLED", assigned_serotype %in% serotypes)
  total_units <- nrow(unit_table)
  counts <- setNames(rep(0L, length(serotypes)), serotypes)
  observed <- table(called$assigned_serotype)
  counts[names(observed)] <- as.integer(observed)
  other_count <- total_units - sum(counts)
  draw_mat <- dirichlet_draws(c(counts, `.OTHER` = other_count), n_draws)
  list(
    draws = draw_mat[, serotypes, drop = FALSE],
    model_method = label,
    counts = counts,
    other_count = other_count,
    total_units = total_units
  )
}

fit_bayesian_prevalence_model <- function(unit_table, serotypes, args, label) {
  called <- unit_table %>%
    filter(call_status == "CALLED", assigned_serotype %in% serotypes)

  total_units <- nrow(unit_table)
  counts <- setNames(rep(0L, length(serotypes)), serotypes)
  observed <- table(called$assigned_serotype)
  counts[names(observed)] <- as.integer(observed)
  other_count <- total_units - sum(counts)

  if (length(serotypes) == 1) {
    draw_mat <- dirichlet_draws(c(counts, `.OTHER` = other_count), args$draws)[, serotypes, drop = FALSE]
    return(list(
      draws = draw_mat,
      model_method = paste0("bayesian_fallback:single_serotype:", label),
      counts = counts,
      other_count = other_count,
      total_units = total_units,
      model = NULL
    ))
  }

  if (dplyr::n_distinct(unit_table$site) < 2) {
    fallback <- get_dirichlet_result(unit_table, serotypes, args$draws, "bayesian_fallback:single_site_dirichlet")
    fallback$model <- NULL
    return(fallback)
  }

  site_totals <- unit_table %>% count(site, name = "n")
  model_data <- tidyr::crossing(site = unique(unit_table$site), serotype = serotypes) %>%
    left_join(site_totals, by = "site") %>%
    left_join(
      unit_table %>%
        filter(call_status == "CALLED", assigned_serotype %in% serotypes) %>%
        count(site, assigned_serotype, name = "event"),
      by = c("site", "serotype" = "assigned_serotype")
    ) %>%
    mutate(event = ifelse(is.na(event), 0L, event))

  model_data$serotype <- factor(model_data$serotype, levels = serotypes)
  model_data$site <- factor(model_data$site)

  priors <- c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 2.5)", class = "sd")
  )

  fit <- brm(
    formula = bf(event | trials(n) ~ 0 + serotype + (1 | site)),
    data = model_data,
    family = binomial(),
    prior = priors,
    iter = args$iter,
    warmup = args$warmup,
    chains = args$chains,
    cores = args$cores,
    seed = 1,
    control = list(adapt_delta = args$adapt_delta, max_treedepth = args$max_treedepth),
    refresh = 0
  )

  newdata <- tibble(
    serotype = factor(serotypes, levels = serotypes),
    site = factor(rep(levels(model_data$site)[1], length(serotypes)), levels = levels(model_data$site)),
    n = 1L
  )

  posterior <- posterior_epred(
    fit,
    newdata = newdata,
    re_formula = NA,
    ndraws = args$draws
  )

  if (length(dim(posterior)) == 2) {
    draw_mat <- posterior
  } else {
    draw_mat <- posterior[, , 1, drop = FALSE]
    draw_mat <- matrix(draw_mat, nrow = dim(posterior)[1], ncol = dim(posterior)[2])
  }
  colnames(draw_mat) <- serotypes

  list(
    draws = draw_mat,
    model_method = paste0("brms_binomial_site:", label),
    counts = counts,
    other_count = other_count,
    total_units = total_units,
    model = fit
  )
}

summarise_draws <- function(draw_mat) {
  tibble(
    serotype = colnames(draw_mat),
    mean = colMeans(draw_mat),
    median = apply(draw_mat, 2, median),
    lower_95 = apply(draw_mat, 2, quantile, probs = 0.025),
    upper_95 = apply(draw_mat, 2, quantile, probs = 0.975)
  )
}

beta_posterior_draws <- function(success, total, n_draws) {
  rbeta(n_draws, success + 1, total - success + 1)
}

compute_cumulative_coverage <- function(ordered_serotypes, counts, total_units, top_n, n_draws) {
  ordered_serotypes <- ordered_serotypes[seq_len(min(length(ordered_serotypes), top_n))]
  bind_rows(lapply(seq_along(ordered_serotypes), function(i) {
    selected <- ordered_serotypes[seq_len(i)]
    success <- sum(as.integer(counts[selected]))
    vals <- beta_posterior_draws(success, total_units, n_draws)
    tibble(
      rank = i,
      antigen_count = i,
      included_serotypes = paste(selected, collapse = ","),
      coverage_mean = mean(vals),
      coverage_median = median(vals),
      coverage_lower_95 = quantile(vals, 0.025),
      coverage_upper_95 = quantile(vals, 0.975)
    )
  }))
}

main <- function() {
  args <- parse_args_simple(commandArgs(trailingOnly = TRUE))
  dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

  calls <- read_tsv(args$calls, show_col_types = FALSE)
  required <- c("sample_id", "assigned_serotype", "call_status")
  if (!all(required %in% names(calls))) {
    stop("Calls file must contain: sample_id, assigned_serotype, call_status")
  }

  sample_data <- prepare_metadata(calls, args$metadata, args$clusters)
  serotypes <- sample_data %>%
    filter(call_status == "CALLED") %>%
    pull(assigned_serotype) %>%
    unique() %>%
    sort()

  if (length(serotypes) == 0) {
    stop("No CALLED serotypes found in calls file.")
  }

  raw_units <- build_raw_unit_table(sample_data)
  cluster_units <- build_cluster_unit_table(sample_data)
  ranking_units <- if (!is.null(cluster_units)) cluster_units else raw_units
  ranking_label <- if (!is.null(cluster_units)) "cluster_adjusted" else "raw"

  ranking_result <- fit_bayesian_prevalence_model(ranking_units, serotypes, args, ranking_label)
  coverage_result <- get_dirichlet_result(raw_units, serotypes, args$draws, "bayesian_dirichlet_raw")

  if (!is.null(ranking_result$model)) {
    saveRDS(ranking_result$model, file.path(args$outdir, "ranking_brms_model.rds"))
  }

  ranking_summary <- summarise_draws(ranking_result$draws) %>%
    arrange(desc(mean), serotype) %>%
    mutate(rank_source = ranking_label, estimation_method = ranking_result$model_method)

  coverage_summary_by_serotype <- summarise_draws(coverage_result$draws) %>%
    rename(raw_mean = mean, raw_median = median, raw_lower_95 = lower_95, raw_upper_95 = upper_95)

  prevalence_summary <- ranking_summary %>%
    left_join(coverage_summary_by_serotype, by = "serotype") %>%
    left_join(
      tibble(
        serotype = serotypes,
        raw_count = as.integer(coverage_result$counts[serotypes]),
        raw_prevalence = as.integer(coverage_result$counts[serotypes]) / coverage_result$total_units,
        ranking_count = as.integer(ranking_result$counts[serotypes]),
        ranking_prevalence = as.integer(ranking_result$counts[serotypes]) / ranking_result$total_units
      ),
      by = "serotype"
    ) %>%
    arrange(desc(mean), serotype)

  cumulative <- compute_cumulative_coverage(
    prevalence_summary$serotype,
    coverage_result$counts,
    coverage_result$total_units,
    args$top_n,
    args$draws
  ) %>%
    mutate(
      ranking_source = ranking_label,
      coverage_source = "raw_bayesian_beta",
      ranking_estimation_method = ranking_result$model_method,
      coverage_estimation_method = "bayesian_beta_raw"
    )

  overview <- tibble(
    metric = c(
      "samples_total",
      "samples_called",
      "samples_mixed",
      "samples_unresolved",
      "serotypes_defined",
      "ranking_units",
      "ranking_unit_total",
      "coverage_unit_total",
      "ranking_estimation_method",
      "coverage_estimation_method",
      "iter",
      "warmup",
      "chains",
      "cores"
    ),
    value = c(
      nrow(sample_data),
      sum(sample_data$call_status == "CALLED"),
      sum(sample_data$call_status %in% c("MIXED", "AMBIGUOUS")),
      sum(sample_data$call_status == "UNRESOLVED"),
      length(serotypes),
      ranking_label,
      ranking_result$total_units,
      coverage_result$total_units,
      ranking_result$model_method,
      coverage_result$model_method,
      args$iter,
      args$warmup,
      args$chains,
      args$cores
    )
  )

  write_tsv(sample_data, file.path(args$outdir, "sample_serotype_calls.tsv"))
  write_tsv(prevalence_summary, file.path(args$outdir, "serotype_prevalence_summary.tsv"))
  write_tsv(cumulative, file.path(args$outdir, "cumulative_serotype_representation.tsv"))
  write_tsv(overview, file.path(args$outdir, "analysis_overview.tsv"))
  write_tsv(raw_units, file.path(args$outdir, "raw_units.tsv"))
  if (!is.null(cluster_units)) {
    write_tsv(cluster_units, file.path(args$outdir, "cluster_adjusted_units.tsv"))
  }
  message("Wrote outputs to: ", args$outdir)
}

main()
