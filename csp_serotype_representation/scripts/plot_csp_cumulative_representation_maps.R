#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
  library(scales)
})

usage <- function() {
  cat(
    paste(
      "Usage:",
      "Rscript scripts/plot_csp_cumulative_representation_maps.R",
      "[--coverage country_cumulative_representation.tsv]",
      "[--serotypes country_serotype_summary.tsv]",
      "--out-prefix csp_cumulative_representation_maps",
      "[--max-rank 4]",
      "[--coverage-measure raw_prevalence|ranking_prevalence]",
      "[--title 'CSP putative-serotype representation by rank']",
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
    coverage = NULL,
    serotypes = NULL,
    out_prefix = NULL,
    max_rank = 4L,
    coverage_measure = "raw_prevalence",
    title = "CSP putative-serotype representation by rank"
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

  out$max_rank <- as.integer(out$max_rank)
  if (is.null(out$out_prefix)) {
    usage()
  }
  if (is.null(out$coverage) && is.null(out$serotypes)) {
    stop("Provide either --coverage or --serotypes.")
  }
  if (!is.null(out$coverage) && !is.null(out$serotypes)) {
    stop("Use either --coverage or --serotypes, not both.")
  }
  if (!out$coverage_measure %in% c("raw_prevalence", "ranking_prevalence")) {
    stop("--coverage-measure must be raw_prevalence or ranking_prevalence")
  }
  out
}

normalise_country_names <- function(x) {
  dplyr::recode(
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

wrap_panel_label <- function(x, width = 34) {
  vapply(
    x,
    function(one) paste(strwrap(one, width = width), collapse = "\n"),
    character(1)
  )
}

build_from_coverage <- function(coverage_path, max_rank) {
  coverage <- read_tsv(coverage_path, show_col_types = FALSE)
  required <- c("country", "rank", "coverage_mean")
  if (!all(required %in% names(coverage))) {
    stop("Coverage file must contain country, rank, and coverage_mean")
  }

  available_ranks <- sort(unique(coverage$rank))
  wanted_ranks <- available_ranks[available_ranks <= max_rank]
  if (length(wanted_ranks) == 0) {
    stop("No coverage rows available for the requested ranks.")
  }

  coverage_subset <- coverage %>%
    filter(rank %in% wanted_ranks) %>%
    mutate(
      country_join = normalise_country_names(country),
      rank_label = paste0("Rank ", rank)
    )

  summary_table <- coverage_subset %>%
    select(country, rank, included_serotypes, coverage_mean, coverage_lower_95, coverage_upper_95) %>%
    arrange(rank, desc(coverage_mean), country)

  list(
    coverage_subset = coverage_subset,
    summary_table = summary_table,
    subtitle = paste0("Maps show cumulative coverage using the top 1 to top ", max(wanted_ranks), " serotypes"),
    caption = "Representation is calculated from raw sample fractions; serotype ranking may be cluster-adjusted upstream.",
    extra_tables = NULL
  )
}

build_from_serotypes <- function(serotypes_path, max_rank, coverage_measure) {
  serotypes <- read_tsv(serotypes_path, show_col_types = FALSE)
  required <- c("country", "serotype", "ranking_count", "ranking_prevalence", "raw_prevalence")
  if (!all(required %in% names(serotypes))) {
    stop("Serotype summary must contain: country, serotype, ranking_count, ranking_prevalence, raw_prevalence")
  }

  global_ranked <- serotypes %>%
    group_by(serotype) %>%
    summarise(
      global_count = sum(ranking_count, na.rm = TRUE),
      global_measure = sum(.data[[coverage_measure]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(global_count), desc(global_measure), serotype)

  all_serotypes <- global_ranked$serotype
  selected_serotypes <- character(0)
  selection_rows <- list()

  if (length(all_serotypes) > 0) {
    first_serotype <- global_ranked$serotype[[1]]
    selected_serotypes <- first_serotype
    selection_rows[[1]] <- tibble(
      rank = 1L,
      selected_serotype = first_serotype,
      objective_score = global_ranked$global_measure[global_ranked$serotype == first_serotype],
      selection_method = "global_most_common"
    )
  }

  score_combo <- function(included) {
    serotypes %>%
      filter(serotype %in% included) %>%
      group_by(country) %>%
      summarise(country_coverage = sum(.data[[coverage_measure]], na.rm = TRUE), .groups = "drop") %>%
      summarise(total_score = sum(country_coverage, na.rm = TRUE)) %>%
      pull(total_score)
  }

  while (length(selected_serotypes) < min(max_rank, length(all_serotypes))) {
    remaining <- setdiff(all_serotypes, selected_serotypes)
    candidate_scores <- lapply(remaining, function(candidate) {
      combo <- c(selected_serotypes, candidate)
      tibble(
        serotype = candidate,
        objective_score = score_combo(combo)
      )
    }) %>%
      bind_rows() %>%
      left_join(global_ranked, by = "serotype") %>%
      arrange(desc(objective_score), desc(global_count), desc(global_measure), serotype)

    best <- candidate_scores[1, ]
    selected_serotypes <- c(selected_serotypes, best$serotype)
    selection_rows[[length(selection_rows) + 1L]] <- tibble(
      rank = length(selected_serotypes),
      selected_serotype = best$serotype,
      objective_score = best$objective_score,
      selection_method = "greedy_max_added_coverage"
    )
  }

  if (length(selected_serotypes) == 0) {
    stop("No globally ranked serotypes available to plot.")
  }

  countries <- sort(unique(serotypes$country))
  coverage_subset <- lapply(countries, function(country_name) {
    country_df <- serotypes %>%
      filter(country == country_name)

    lapply(seq_along(selected_serotypes), function(i) {
      included <- selected_serotypes[seq_len(i)]
      coverage_value <- country_df %>%
        filter(serotype %in% included) %>%
        summarise(val = sum(.data[[coverage_measure]], na.rm = TRUE)) %>%
        pull(val)

      tibble(
        country = country_name,
        rank = i,
        included_serotypes = paste(included, collapse = ","),
        coverage_mean = coverage_value
      )
    }) %>%
      bind_rows()
  }) %>%
    bind_rows() %>%
    mutate(
      country_join = normalise_country_names(country),
      rank_label = paste0("Rank ", rank)
    )

  summary_table <- coverage_subset %>%
    arrange(rank, desc(coverage_mean), country)

  selected_summary <- bind_rows(selection_rows) %>%
    mutate(
      included_serotypes = vapply(
        seq_len(n()),
        function(i) paste(selected_serotypes[seq_len(i)], collapse = ","),
        character(1)
      )
    )

  measure_label <- if (coverage_measure == "raw_prevalence") {
    "raw"
  } else {
    "cluster-adjusted"
  }

  list(
    coverage_subset = coverage_subset,
    summary_table = summary_table,
    subtitle = paste0(
      "Maps show cumulative coverage using a global greedy serotype ordering; coverage uses ",
      measure_label,
      " fractions"
    ),
    caption = "Rank 1 uses the globally most common serotype; later ranks add the serotype that maximises extra cross-country coverage.",
    extra_tables = list(
      global_serotype_ranking = global_ranked,
      selected_serotype_order = selected_summary
    )
  )
}

main <- function() {
  args <- parse_args_simple(commandArgs(trailingOnly = TRUE))

  built <- if (!is.null(args$coverage)) {
    build_from_coverage(args$coverage, args$max_rank)
  } else {
    build_from_serotypes(args$serotypes, args$max_rank, args$coverage_measure)
  }

  coverage_subset <- built$coverage_subset
  wanted_ranks <- sort(unique(coverage_subset$rank))

  panel_labels <- coverage_subset %>%
    group_by(rank) %>%
    summarise(
      n_combo = dplyr::n_distinct(included_serotypes),
      panel_label = if (first(n_combo) == 1L) {
        first(included_serotypes)
      } else {
        paste0("Rank ", first(rank), "\n(country-specific combinations)")
      },
      .groups = "drop"
    ) %>%
    mutate(
      rank_label = paste0("Rank ", rank),
      panel_label = wrap_panel_label(panel_label)
    )

  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
    mutate(country_join = normalise_country_names(admin))

  rank_maps <- lapply(wanted_ranks, function(this_rank) {
    world %>%
      left_join(
        coverage_subset %>%
          filter(rank == this_rank) %>%
          select(country_join, coverage_mean, included_serotypes, country, rank_label),
        by = "country_join"
      ) %>%
      mutate(rank = this_rank, rank_label = paste0("Rank ", this_rank))
  })

  map_data <- do.call(rbind, rank_maps)
  map_data <- map_data %>%
    left_join(
      panel_labels %>% select(rank_label, panel_label),
      by = "rank_label"
    )
  map_data$panel_label <- factor(
    map_data$panel_label,
    levels = panel_labels$panel_label[match(paste0("Rank ", wanted_ranks), panel_labels$rank_label)]
  )

  p <- ggplot(map_data) +
    geom_sf(aes(fill = coverage_mean), color = "#6F6F6F", linewidth = 0.34) +
    geom_sf(fill = NA, color = "#3A3A3A", linewidth = 0.22) +
    scale_fill_gradientn(
      name = "Coverage",
      colours = c("#B2182B", "#EF8A62", "#F7F7F7", "#67A9CF", "#2166AC"),
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = percent_format(accuracy = 1),
      oob = scales::squish,
      na.value = "#d0d0d0"
    ) +
    coord_sf(expand = FALSE) +
    facet_wrap(~ panel_label, ncol = 2) +
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
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      title = args$title,
      subtitle = built$subtitle,
      caption = built$caption
    )

  write_tsv(built$summary_table, paste0(args$out_prefix, ".rank_summary.tsv"))
  if (!is.null(built$extra_tables$global_serotype_ranking)) {
    write_tsv(built$extra_tables$global_serotype_ranking, paste0(args$out_prefix, ".global_serotype_ranking.tsv"))
  }
  if (!is.null(built$extra_tables$selected_serotype_order)) {
    write_tsv(built$extra_tables$selected_serotype_order, paste0(args$out_prefix, ".selected_serotype_order.tsv"))
  }
  ggsave(paste0(args$out_prefix, ".png"), p, width = 14, height = 9, dpi = 300)
  ggsave(paste0(args$out_prefix, ".pdf"), p, width = 14, height = 9)

  message("Wrote: ", paste0(args$out_prefix, ".rank_summary.tsv"))
  if (!is.null(built$extra_tables$global_serotype_ranking)) {
    message("Wrote: ", paste0(args$out_prefix, ".global_serotype_ranking.tsv"))
  }
  if (!is.null(built$extra_tables$selected_serotype_order)) {
    message("Wrote: ", paste0(args$out_prefix, ".selected_serotype_order.tsv"))
  }
  message("Wrote: ", paste0(args$out_prefix, ".png"))
  message("Wrote: ", paste0(args$out_prefix, ".pdf"))
}

main()
