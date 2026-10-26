#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(brms)
  library(ggplot2)
  library(scales)
})

usage <- function() {
  cat(
    paste(
      "Usage:",
      "Rscript scripts/estimate_country_serotype_representation_bayesian.R",
      "--calls sample_calls.tsv",
      "--metadata sample_metadata.tsv",
      "--out-prefix country_coverage_bayesian",
      "[--clusters csp_sample_clusters.tsv]",
      "[--country-col country]",
      "[--top-n 5]",
      "[--draws 2000]",
      "[--iter 4000]",
      "[--warmup 2000]",
      "[--chains 4]",
      "[--cores 4]",
      "[--adapt-delta 0.99]",
      "[--max-treedepth 15]",
      "[--estimation brms|dirichlet]",
      "[--save-models false]",
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
    out_prefix = NULL,
    country_col = "country",
    top_n = 5L,
    draws = 2000L,
    iter = 4000L,
    warmup = 2000L,
    chains = 4L,
    cores = 4L,
    adapt_delta = 0.99,
    max_treedepth = 15L,
    estimation = "brms",
    save_models = "false"
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
  out$estimation <- tolower(out$estimation)
  if (!out$estimation %in% c("brms", "dirichlet")) {
    stop("--estimation must be one of: brms, dirichlet")
  }
  out$save_models <- tolower(out$save_models) %in% c("true", "1", "yes")

  if (is.null(out$calls) || is.null(out$metadata) || is.null(out$out_prefix)) usage()
  out
}

prepare_sample_data <- function(calls_path, metadata_path, clusters_path = NULL, country_col = "country") {
  calls <- read_tsv(calls_path, show_col_types = FALSE)
  required_call_cols <- c("sample_id", "assigned_serotype", "call_status")
  if (!all(required_call_cols %in% names(calls))) {
    stop("Calls file must contain: sample_id, assigned_serotype, call_status")
  }
  calls <- calls %>%
    mutate(
      assigned_serotype = as.character(assigned_serotype),
      call_status = as.character(call_status)
    )

  metadata <- read_tsv(metadata_path, show_col_types = FALSE)
  if (!"sample_id" %in% names(metadata)) {
    stop("Metadata file must contain sample_id")
  }
  if (!country_col %in% names(metadata)) {
    stop("Metadata file must contain the country column: ", country_col)
  }

  sample_data <- calls %>%
    left_join(metadata, by = "sample_id") %>%
    mutate(
      country = .data[[country_col]],
      site = if ("site" %in% names(.)) as.character(site) else "GLOBAL_SITE",
      region = if ("region" %in% names(.)) as.character(region) else "GLOBAL"
    )

  if (!is.null(clusters_path)) {
    clusters <- read_tsv(clusters_path, show_col_types = FALSE)
    if (!all(c("sample_id", "cluster_id") %in% names(clusters))) {
      stop("Cluster file must contain sample_id and cluster_id")
    }
    clusters <- clusters %>% select(sample_id, cluster_id) %>% distinct()
    sample_data <- sample_data %>% left_join(clusters, by = "sample_id")
  } else {
    sample_data$cluster_id <- NA_character_
  }

  sample_data %>%
    mutate(
      country = ifelse(is.na(country) | country == "", "Unknown", as.character(country)),
      site = ifelse(is.na(site) | site == "", "GLOBAL_SITE", as.character(site)),
      region = ifelse(is.na(region) | region == "", "GLOBAL", as.character(region)),
      cluster_id = as.character(cluster_id)
    )
}

build_raw_unit_table <- function(sample_data) {
  sample_data %>%
    transmute(
      unit_id = sample_id,
      sample_id = sample_id,
      country = country,
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
    group_by(country, site, region, cluster_key) %>%
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
      call_status = case_when(
        assigned_serotype == "MIXED_CLUSTER" ~ "MIXED_CLUSTER",
        assigned_serotype == "UNRESOLVED_CLUSTER" ~ "UNRESOLVED_CLUSTER",
        assigned_serotype == "OTHER_CLUSTER" ~ "OTHER_CLUSTER",
        TRUE ~ "CALLED"
      ),
      unit_id = paste(site, cluster_key, sep = "::")
    ) %>%
    select(unit_id, country, site, region, assigned_serotype, call_status, sample_count)
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

get_dirichlet_result <- function(unit_table, serotypes, n_draws, label) {
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

  if (identical(args$estimation, "dirichlet")) {
    fallback <- get_dirichlet_result(unit_table, serotypes, args$draws, paste0("bayesian_dirichlet:", label))
    fallback$model <- NULL
    return(fallback)
  }

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
    fallback <- get_dirichlet_result(unit_table, serotypes, args$draws, paste0("bayesian_fallback:single_site:", label))
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

beta_interval <- function(success, total) {
  if (total <= 0) return(c(NA_real_, NA_real_))
  c(qbeta(0.025, success + 1, total - success + 1), qbeta(0.975, success + 1, total - success + 1))
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

summarise_country <- function(country_name, sample_data, cluster_units, args) {
  raw_country <- build_raw_unit_table(sample_data %>% filter(country == country_name))
  if (nrow(raw_country) == 0) return(NULL)

  ranking_country <- if (!is.null(cluster_units)) {
    cluster_units %>% filter(country == country_name)
  } else {
    raw_country
  }

  if (nrow(ranking_country) == 0) return(NULL)

  serotypes <- ranking_country %>%
    filter(call_status == "CALLED") %>%
    pull(assigned_serotype) %>%
    unique() %>%
    sort()

  if (length(serotypes) == 0) return(NULL)

  ranking_label <- if (!is.null(cluster_units)) paste0("cluster_adjusted:", country_name) else paste0("raw:", country_name)
  ranking_result <- fit_bayesian_prevalence_model(ranking_country, serotypes, args, ranking_label)
  coverage_result <- get_dirichlet_result(raw_country, serotypes, args$draws, paste0("country_dirichlet_raw:", country_name))

  ranking_summary <- summarise_draws(ranking_result$draws) %>%
    arrange(desc(mean), serotype) %>%
    mutate(
      country = country_name,
      rank_source = if (!is.null(cluster_units)) "cluster_adjusted" else "raw",
      estimation_method = ranking_result$model_method
    )

  coverage_summary_by_serotype <- summarise_draws(coverage_result$draws) %>%
    rename(raw_mean = mean, raw_median = median, raw_lower_95 = lower_95, raw_upper_95 = upper_95)

  serotype_summary <- ranking_summary %>%
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
    serotype_summary$serotype,
    coverage_result$counts,
    coverage_result$total_units,
    args$top_n,
    args$draws
  ) %>%
    mutate(
      country = country_name,
      ranking_source = if (!is.null(cluster_units)) "cluster_adjusted" else "raw",
      coverage_source = "raw_bayesian_beta",
      ranking_estimation_method = ranking_result$model_method,
      coverage_estimation_method = "bayesian_beta_raw"
    )

  overview <- tibble(
    country = country_name,
    metric = c(
      "samples_total",
      "samples_called",
      "samples_mixed",
      "samples_unresolved",
      "serotypes_defined",
      "ranking_unit_total",
      "coverage_unit_total",
      "ranking_estimation_method",
      "coverage_estimation_method"
    ),
    value = c(
      nrow(sample_data %>% filter(country == country_name)),
      sum((sample_data %>% filter(country == country_name))$call_status == "CALLED"),
      sum((sample_data %>% filter(country == country_name))$call_status %in% c("MIXED", "AMBIGUOUS")),
      sum((sample_data %>% filter(country == country_name))$call_status == "UNRESOLVED"),
      length(serotypes),
      ranking_result$total_units,
      coverage_result$total_units,
      ranking_result$model_method,
      coverage_result$model_method
    )
  )

  list(
    serotypes = serotype_summary,
    cumulative = cumulative,
    overview = overview,
    model = ranking_result$model
  )
}

main <- function() {
  args <- parse_args_simple(commandArgs(trailingOnly = TRUE))
  dir.create(dirname(args$out_prefix), recursive = TRUE, showWarnings = FALSE)

  sample_data <- prepare_sample_data(args$calls, args$metadata, args$clusters, args$country_col)
  cluster_units <- build_cluster_unit_table(sample_data)

  countries <- sort(unique(sample_data$country))
  pieces <- lapply(countries, summarise_country, sample_data = sample_data, cluster_units = cluster_units, args = args)
  keep <- !vapply(pieces, is.null, logical(1))
  pieces <- pieces[keep]
  countries <- countries[keep]

  if (length(pieces) == 0) {
    stop("No countries had callable serotype data.")
  }

  serotype_summary <- bind_rows(lapply(pieces, `[[`, "serotypes"))
  cumulative <- bind_rows(lapply(pieces, `[[`, "cumulative"))
  overview <- bind_rows(lapply(pieces, `[[`, "overview"))

  if (args$save_models) {
    model_dir <- paste0(args$out_prefix, ".brms_models")
    dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
    for (i in seq_along(pieces)) {
      model <- pieces[[i]]$model
      if (!is.null(model)) {
        safe_country <- gsub("[^A-Za-z0-9._-]+", "_", countries[[i]])
        saveRDS(model, file.path(model_dir, paste0(safe_country, ".rds")))
      }
    }
  }

  coverage_heatmap <- cumulative %>%
    mutate(
      country = reorder(country, -coverage_mean, FUN = max),
      antigen_count = factor(antigen_count, levels = sort(unique(antigen_count)))
    ) %>%
    ggplot(aes(x = antigen_count, y = country, fill = coverage_mean)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = percent(coverage_mean, accuracy = 1)), size = 2.8) +
    scale_fill_gradientn(
      colours = c("lightblue", "skyblue2", "gold", "orange2", "red3", "red4"),
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = percent_format(accuracy = 1),
      oob = scales::squish,
      na.value = "#d0d0d0",
      name = "Representation"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(face = "bold")
    ) +
    labs(
      title = "Country-level serotype representation (Bayesian ranking)",
      subtitle = "Serotypes ranked by Bayesian cluster-adjusted prevalence; cells show raw sample representation",
      caption = "Columns are the number of serotypes included in the hypothetical combination."
    )

  write_tsv(serotype_summary, paste0(args$out_prefix, ".country_serotype_summary.tsv"))
  write_tsv(cumulative, paste0(args$out_prefix, ".country_cumulative_representation.tsv"))
  write_tsv(overview, paste0(args$out_prefix, ".country_analysis_overview.tsv"))
  ggsave(paste0(args$out_prefix, ".representation_heatmap.png"), coverage_heatmap, width = 11, height = 7, dpi = 300)
  ggsave(paste0(args$out_prefix, ".representation_heatmap.pdf"), coverage_heatmap, width = 11, height = 7)

  message("Wrote: ", paste0(args$out_prefix, ".country_serotype_summary.tsv"))
  message("Wrote: ", paste0(args$out_prefix, ".country_cumulative_representation.tsv"))
  message("Wrote: ", paste0(args$out_prefix, ".country_analysis_overview.tsv"))
  message("Wrote: ", paste0(args$out_prefix, ".representation_heatmap.png"))
  message("Wrote: ", paste0(args$out_prefix, ".representation_heatmap.pdf"))
  if (args$save_models) {
    message("Wrote per-country models under: ", paste0(args$out_prefix, ".brms_models"))
  }
}

main()
