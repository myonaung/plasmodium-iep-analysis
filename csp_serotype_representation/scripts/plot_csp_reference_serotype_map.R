#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
  library(scales)
  library(tibble)
})

usage <- function() {
  cat(
    paste(
      "Usage:",
      "Rscript scripts/plot_csp_reference_serotype_map.R",
      "--calls sample_calls.tsv",
      "--metadata sample_metadata.tsv",
      "--out-prefix csp_reference_serotype_map",
      "[--clusters csp_sample_clusters.tsv]",
      "[--variant-summary sample_variant_summary.tsv]",
      "[--nonsynonymous-only true|false]",
      "[--annotation-keep NON_SYNONYMOUS_CODING,missense_variant]",
      "[--bayesian true|false]",
      "[--draws 2000]",
      "[--iter 4000]",
      "[--warmup 2000]",
      "[--chains 4]",
      "[--cores 4]",
      "[--adapt-delta 0.99]",
      "[--max-treedepth 15]",
      "[--reference-serotype PROFILE]",
      "[--country-col country]",
      "[--title '3D7-like CSP putative-serotype representation']",
      sep = " "
    ),
    "\n"
  )
  quit(status = 1)
}

parse_args_simple <- function(args) {
  if (length(args) == 0) {
    usage()
  }

  out <- list(
    calls = NULL,
    metadata = NULL,
    clusters = NULL,
    variant_summary = NULL,
    out_prefix = NULL,
    nonsynonymous_only = "false",
    annotation_keep = "NON_SYNONYMOUS_CODING",
    bayesian = "false",
    draws = 2000L,
    iter = 4000L,
    warmup = 2000L,
    chains = 4L,
    cores = 4L,
    adapt_delta = 0.99,
    max_treedepth = 15L,
    reference_serotype = NULL,
    country_col = "country",
    title = "3D7-like CSP putative-serotype representation by country"
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (!startsWith(arg, "--") || i == length(args)) {
      stop("Bad arguments.")
    }
    key <- gsub("-", "_", sub("^--", "", arg))
    value <- args[[i + 1L]]
    i <- i + 2L
    if (!key %in% names(out)) {
      stop("Unknown option: --", gsub("_", "-", key))
    }
    out[[key]] <- value
  }

  if (is.null(out$calls) || is.null(out$metadata) || is.null(out$out_prefix)) {
    usage()
  }
  out$nonsynonymous_only <- tolower(as.character(out$nonsynonymous_only))
  if (!out$nonsynonymous_only %in% c("true", "false")) {
    stop("--nonsynonymous-only must be true or false")
  }
  out$bayesian <- tolower(as.character(out$bayesian))
  if (!out$bayesian %in% c("true", "false")) {
    stop("--bayesian must be true or false")
  }
  out$annotation_keep <- as.character(out$annotation_keep)
  out$draws <- as.integer(out$draws)
  out$iter <- as.integer(out$iter)
  out$warmup <- as.integer(out$warmup)
  out$chains <- as.integer(out$chains)
  out$cores <- as.integer(out$cores)
  out$adapt_delta <- as.numeric(out$adapt_delta)
  out$max_treedepth <- as.integer(out$max_treedepth)
  out
}

parse_annotation_keep <- function(annotation_keep) {
  values <- trimws(unlist(strsplit(annotation_keep, ",", fixed = TRUE), use.names = FALSE))
  values <- values[nzchar(values)]
  unique(values)
}

build_annotation_key_set <- function(variant_summary_path, annotation_keep) {
  variant_summary <- read_tsv(variant_summary_path, show_col_types = FALSE)
  required_cols <- c("analysis_mode", "chrom", "pos", "annotation")
  if (!all(required_cols %in% names(variant_summary))) {
    stop("Variant summary must contain: analysis_mode, chrom, pos, annotation")
  }

  keep_values <- parse_annotation_keep(annotation_keep)
  if (length(keep_values) == 0) {
    stop("--annotation-keep must contain at least one annotation label")
  }

  variant_summary %>%
    filter(
      analysis_mode == "all_variants",
      Reduce(
        `|`,
        lapply(keep_values, function(label) grepl(label, annotation, fixed = TRUE))
      )
    ) %>%
    transmute(site_key = paste(chrom, pos, sep = ":")) %>%
    distinct() %>%
    pull(site_key)
}

filter_all_variant_haplotype <- function(haplotype, nonsynonymous_keys) {
  if (is.na(haplotype) || haplotype %in% c("", "REF_HAPLOTYPE", "MIXED", "UNRESOLVED")) {
    return(haplotype)
  }

  site_entries <- unlist(strsplit(haplotype, ";", fixed = TRUE), use.names = FALSE)
  if (length(site_entries) == 0) {
    return("REF_HAPLOTYPE")
  }

  keep <- vapply(
    site_entries,
    function(entry) {
      if (!grepl("=", entry, fixed = TRUE)) {
        return(FALSE)
      }
      lhs <- strsplit(entry, "=", fixed = TRUE)[[1]][1]
      parts <- strsplit(lhs, ":", fixed = TRUE)[[1]]
      if (length(parts) < 2) {
        return(FALSE)
      }
      site_key <- paste(parts[1], parts[2], sep = ":")
      site_key %in% nonsynonymous_keys
    },
    logical(1)
  )

  kept_entries <- site_entries[keep]
  if (length(kept_entries) == 0) {
    "REF_HAPLOTYPE"
  } else {
    paste(kept_entries, collapse = ";")
  }
}

apply_annotation_filter <- function(calls, variant_summary_path, annotation_keep, reference_serotype = NULL) {
  retained_keys <- build_annotation_key_set(variant_summary_path, annotation_keep)

  if (!"analysis_mode" %in% names(calls)) {
    calls$analysis_mode <- NA_character_
  }

  calls <- calls %>%
    mutate(
      assigned_serotype = ifelse(
        analysis_mode == "all_variants" & call_status == "CALLED",
        vapply(assigned_serotype, filter_all_variant_haplotype, character(1), nonsynonymous_keys = retained_keys),
        assigned_serotype
      )
    )

  if (is.null(reference_serotype) || !nzchar(reference_serotype)) {
    reference_filtered <- reference_serotype
  } else if (reference_serotype %in% c("REF_HAPLOTYPE", "MIXED", "UNRESOLVED")) {
    reference_filtered <- reference_serotype
  } else if (grepl("=", reference_serotype, fixed = TRUE)) {
    reference_filtered <- filter_all_variant_haplotype(reference_serotype, retained_keys)
  } else {
    reference_filtered <- reference_serotype
  }

  list(calls = calls, reference_serotype = reference_filtered)
}

prepare_sample_data <- function(
  calls_path,
  metadata_path,
  clusters_path = NULL,
  country_col = "country",
  variant_summary_path = NULL,
  nonsynonymous_only = FALSE,
  annotation_keep = "NON_SYNONYMOUS_CODING",
  reference_serotype = NULL
) {
  calls <- read_tsv(calls_path, show_col_types = FALSE)
  required_call_cols <- c("sample_id", "assigned_serotype", "call_status")
  if (!all(required_call_cols %in% names(calls))) {
    stop("Calls file must contain: sample_id, assigned_serotype, call_status")
  }

  filtered_reference <- reference_serotype
  if (nonsynonymous_only) {
    if (is.null(variant_summary_path)) {
      stop("--variant-summary is required when --nonsynonymous-only true")
    }
    filtered <- apply_annotation_filter(calls, variant_summary_path, annotation_keep, reference_serotype)
    calls <- filtered$calls
    filtered_reference <- filtered$reference_serotype
  }

  metadata <- read_tsv(metadata_path, show_col_types = FALSE)
  if (!"sample_id" %in% names(metadata)) {
    stop("Metadata file must contain sample_id")
  }
  if (!country_col %in% names(metadata)) {
    stop("Metadata file must contain country column: ", country_col)
  }

  sample_data <- calls %>%
    left_join(metadata, by = "sample_id") %>%
    mutate(
      country = .data[[country_col]],
      site = if ("site" %in% names(.)) as.character(site) else "GLOBAL_SITE"
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

  sample_data <- sample_data %>%
    mutate(
      country = ifelse(is.na(country) | country == "", "Unknown", as.character(country)),
      site = ifelse(is.na(site) | site == "", "GLOBAL_SITE", as.character(site)),
      cluster_id = as.character(cluster_id)
    )

  list(sample_data = sample_data, reference_serotype = filtered_reference)
}

build_cluster_adjusted_units <- function(sample_data) {
  sample_data %>%
    mutate(cluster_key = ifelse(is.na(cluster_id) | cluster_id == "", paste0("singleton:", sample_id), cluster_id)) %>%
    group_by(country, site, cluster_key) %>%
    summarise(
      sample_count = n(),
      called_serotypes = list(sort(unique(assigned_serotype[call_status == "CALLED"]))),
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
    )
}

beta_interval <- function(success, total) {
  if (total <= 0) {
    return(c(NA_real_, NA_real_))
  }
  c(qbeta(0.025, success + 1, total - success + 1), qbeta(0.975, success + 1, total - success + 1))
}

beta_posterior_draws <- function(success, total, n_draws) {
  rbeta(n_draws, success + 1, total - success + 1)
}

summarise_reference_country_bayesian <- function(cluster_units, reference_serotype, args) {
  countries <- sort(unique(cluster_units$country))

  bind_rows(lapply(countries, function(country_name) {
    country_units <- cluster_units %>% filter(country == country_name)
    if (nrow(country_units) == 0) {
      return(NULL)
    }

    success <- sum(country_units$call_status == "CALLED" & country_units$assigned_serotype == reference_serotype)
    total <- nrow(country_units)

    if (dplyr::n_distinct(country_units$site) < 2) {
      vals <- beta_posterior_draws(success, total, args$draws)
      return(tibble(
        country = country_name,
        cluster_adjusted_total = total,
        cluster_adjusted_called = sum(country_units$call_status == "CALLED"),
        cluster_adjusted_reference = success,
        cluster_adjusted_fraction_all = mean(vals),
        cluster_adjusted_fraction_called = ifelse(
          sum(country_units$call_status == "CALLED") > 0,
          success / sum(country_units$call_status == "CALLED"),
          NA_real_
        ),
        lower_95_all = quantile(vals, 0.025),
        upper_95_all = quantile(vals, 0.975),
        estimation_method = "bayesian_beta_single_site"
      ))
    }

    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("brms is required for --bayesian true")
    }

    site_data <- country_units %>%
      group_by(site) %>%
      summarise(
        event = sum(call_status == "CALLED" & assigned_serotype == reference_serotype),
        n = n(),
        .groups = "drop"
      )
    site_data$site <- factor(site_data$site)

    fit <- brms::brm(
      formula = brms::bf(event | trials(n) ~ 1 + (1 | site)),
      data = site_data,
      family = stats::binomial(),
      prior = c(
        brms::set_prior("normal(0, 2.5)", class = "Intercept"),
        brms::set_prior("student_t(3, 0, 2.5)", class = "sd")
      ),
      iter = args$iter,
      warmup = args$warmup,
      chains = args$chains,
      cores = args$cores,
      seed = 1,
      control = list(adapt_delta = args$adapt_delta, max_treedepth = args$max_treedepth),
      refresh = 0
    )

    newdata <- tibble(site = factor(levels(site_data$site)[1], levels = levels(site_data$site)), n = 1L)
    posterior <- brms::posterior_epred(
      fit,
      newdata = newdata,
      re_formula = NA,
      ndraws = args$draws
    )

    if (length(dim(posterior)) == 2) {
      vals <- posterior[, 1]
    } else {
      vals <- posterior[, 1, 1]
    }

    tibble(
      country = country_name,
      cluster_adjusted_total = total,
      cluster_adjusted_called = sum(country_units$call_status == "CALLED"),
      cluster_adjusted_reference = success,
      cluster_adjusted_fraction_all = mean(vals),
      cluster_adjusted_fraction_called = ifelse(
        sum(country_units$call_status == "CALLED") > 0,
        success / sum(country_units$call_status == "CALLED"),
        NA_real_
      ),
      lower_95_all = quantile(vals, 0.025),
      upper_95_all = quantile(vals, 0.975),
      estimation_method = "brms_binomial_site_reference"
    )
  }))
}

normalise_country_names <- function(x) {
  recode(
    x,
    "Democratic Republic of the Congo" = "Dem. Rep. Congo",
    "Republic of the Congo" = "Congo",
    "Cote d'Ivoire" = "Ivory Coast",
    "Côte d’Ivoire" = "Ivory Coast",
    "The Gambia" = "Gambia",
    "United Republic of Tanzania" = "Tanzania",
    "Myanmar (Burma)" = "Myanmar",
    "Lao PDR" = "Laos",
    "Papua New Guinea" = "Papua New Guinea",
    .default = x
  )
}

infer_reference_serotype <- function(sample_data, user_value = NULL) {
  if (!is.null(user_value) && nzchar(user_value)) {
    return(user_value)
  }
  called_serotypes <- sort(unique(sample_data$assigned_serotype[sample_data$call_status == "CALLED"]))
  if ("REF_HAPLOTYPE" %in% called_serotypes) {
    return("REF_HAPLOTYPE")
  }
  stop(
    "Could not infer the reference serotype automatically. ",
    "Please pass --reference-serotype using a profile present in the calls table."
  )
}

main <- function() {
  args <- parse_args_simple(commandArgs(trailingOnly = TRUE))
  out_prefix <- args$out_prefix
  nonsynonymous_only <- identical(args$nonsynonymous_only, "true")
  bayesian <- identical(args$bayesian, "true")

  prepared <- prepare_sample_data(
    args$calls,
    args$metadata,
    args$clusters,
    args$country_col,
    args$variant_summary,
    nonsynonymous_only,
    args$annotation_keep,
    args$reference_serotype
  )
  sample_data <- prepared$sample_data
  reference_serotype <- infer_reference_serotype(sample_data, prepared$reference_serotype)

  cluster_units <- build_cluster_adjusted_units(sample_data)

  if (bayesian) {
    country_summary <- summarise_reference_country_bayesian(cluster_units, reference_serotype, args)
  } else {
    country_summary <- cluster_units %>%
      group_by(country) %>%
      summarise(
        cluster_adjusted_total = n(),
        cluster_adjusted_called = sum(call_status == "CALLED"),
        cluster_adjusted_reference = sum(call_status == "CALLED" & assigned_serotype == reference_serotype),
        cluster_adjusted_fraction_all = cluster_adjusted_reference / cluster_adjusted_total,
        cluster_adjusted_fraction_called = ifelse(cluster_adjusted_called > 0, cluster_adjusted_reference / cluster_adjusted_called, NA_real_),
        .groups = "drop"
      )

    intervals <- t(mapply(beta_interval, country_summary$cluster_adjusted_reference, country_summary$cluster_adjusted_total))
    country_summary$lower_95_all <- intervals[, 1]
    country_summary$upper_95_all <- intervals[, 2]
    country_summary$estimation_method <- "cluster_adjusted_beta_interval"
  }
  country_summary$reference_serotype <- reference_serotype

  raw_summary <- sample_data %>%
    group_by(country) %>%
    summarise(
      raw_total = n(),
      raw_called = sum(call_status == "CALLED"),
      raw_reference = sum(call_status == "CALLED" & assigned_serotype == reference_serotype),
      raw_fraction_all = raw_reference / raw_total,
      raw_fraction_called = ifelse(raw_called > 0, raw_reference / raw_called, NA_real_),
      .groups = "drop"
    )

  country_summary <- country_summary %>%
    left_join(raw_summary, by = "country") %>%
    arrange(desc(cluster_adjusted_fraction_all), country)

  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
    mutate(country_join = normalise_country_names(admin))

  map_data <- world %>%
    left_join(
      country_summary %>% mutate(country_join = normalise_country_names(country)),
      by = "country_join"
    )

  legend_breaks <- seq(0, 1, by = 0.2)

  p <- ggplot(map_data) +
    geom_sf(aes(fill = cluster_adjusted_fraction_all), color = "#6F6F6F", linewidth = 0.38) +
    geom_sf(fill = NA, color = "#3A3A3A", linewidth = 0.24) +
    scale_fill_gradientn(
      name = "Reference serotype\nfraction",
      colours = c("#B2182B", "#EF8A62", "#F7F7F7", "#67A9CF", "#2166AC"),
      breaks = legend_breaks,
      limits = c(0, 1),
      na.value = "#d0d0d0",
      labels = percent_format(accuracy = 1),
      oob = scales::squish
    ) +
    coord_sf(expand = FALSE) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.background = element_rect(fill = "#F8F8F8", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.key.height = grid::unit(0.7, "cm"),
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    ) +
    labs(
      title = args$title,
      subtitle = paste("Reference putative serotype:", reference_serotype),
      caption = if (bayesian && nonsynonymous_only) {
        paste0(
          "Bayesian site-adjusted prevalence on cluster-adjusted units; retained annotations: ",
          paste(parse_annotation_keep(args$annotation_keep), collapse = ", "),
          "; legend fixed from 0% to 100%."
        )
      } else if (bayesian) {
        "Bayesian site-adjusted prevalence on cluster-adjusted units; legend fixed from 0% to 100%."
      } else if (nonsynonymous_only) {
        paste0(
          "Fraction uses cluster-adjusted units; retained annotations: ",
          paste(parse_annotation_keep(args$annotation_keep), collapse = ", "),
          "; legend fixed from 0% to 100%."
        )
      } else {
        "Fraction uses cluster-adjusted units; legend fixed from 0% to 100%."
      }
    )

  write_tsv(country_summary, paste0(out_prefix, ".country_summary.tsv"))
  write_tsv(cluster_units, paste0(out_prefix, ".cluster_adjusted_units.tsv"))
  ggsave(paste0(out_prefix, ".png"), p, width = 12, height = 6, dpi = 300)
  ggsave(paste0(out_prefix, ".pdf"), p, width = 12, height = 6)

  message("Reference serotype: ", reference_serotype)
  message("Wrote: ", paste0(out_prefix, ".country_summary.tsv"))
  message("Wrote: ", paste0(out_prefix, ".cluster_adjusted_units.tsv"))
  message("Wrote: ", paste0(out_prefix, ".png"))
  message("Wrote: ", paste0(out_prefix, ".pdf"))
}

main()
