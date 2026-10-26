#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

parse_args <- function(args) {
  out <- list(
    raw_flows = "",
    alluvial_flows = "",
    output = "pf_pv_haplotype_serotype_accumulation.pdf",
    output_svg = "",
    summary_output = "",
    permutations = 2000,
    seed = 1,
    facet_ncol = 5,
    min_original_haplotype_percent = 1,
    min_reduced_haplotype_percent = 2,
    min_serotypes_exclusive = 1,
    max_serotypes_exclusive = 11,
    max_rare_original_haplotype_percent = 65,
    max_iep_hamming_distance = 4,
    pf_antigens = "",
    pv_antigens = "",
    point_size = 0,
    ribbon_type = "se",
    show_iqr = "TRUE",
    show_confidence_interval = "",
    font_family = "Arial",
    panel_aspect_ratio = 1,
    pdf_width = 10,
    pdf_height = 6.5
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

  if (!nzchar(out$raw_flows)) stop("--raw-flows is required.")
  if (!nzchar(out$alluvial_flows)) stop("--alluvial-flows is required.")
  if (!nzchar(out$output)) stop("--output is required.")
  out$permutations <- as.integer(out$permutations)
  out$seed <- as.integer(out$seed)
  out$facet_ncol <- as.integer(out$facet_ncol)
  out$min_original_haplotype_percent <- as.numeric(out$min_original_haplotype_percent)
  out$min_reduced_haplotype_percent <- as.numeric(out$min_reduced_haplotype_percent)
  out$min_serotypes_exclusive <- as.numeric(out$min_serotypes_exclusive)
  out$max_serotypes_exclusive <- as.numeric(out$max_serotypes_exclusive)
  out$max_rare_original_haplotype_percent <-
    as.numeric(out$max_rare_original_haplotype_percent)
  if (
    !is.finite(out$max_rare_original_haplotype_percent) ||
      out$max_rare_original_haplotype_percent < 0 ||
      out$max_rare_original_haplotype_percent > 100
  ) {
    stop("--max-rare-original-haplotype-percent must be between 0 and 100.")
  }
  out$max_iep_hamming_distance <- as.integer(out$max_iep_hamming_distance)
  if (is.na(out$max_iep_hamming_distance) || out$max_iep_hamming_distance < 0) {
    stop("--max-iep-hamming-distance must be a non-negative integer.")
  }
  out$pf_antigens <- trimws(strsplit(out$pf_antigens, ",", fixed = TRUE)[[1]])
  out$pf_antigens <- out$pf_antigens[nzchar(out$pf_antigens)]
  out$pv_antigens <- trimws(strsplit(out$pv_antigens, ",", fixed = TRUE)[[1]])
  out$pv_antigens <- out$pv_antigens[nzchar(out$pv_antigens)]
  out$point_size <- as.numeric(out$point_size)
  out$ribbon_type <- match.arg(
    tolower(out$ribbon_type),
    c("none", "se", "sd", "iqr", "ci95")
  )
  out$show_iqr <-
    tolower(out$show_iqr) %in% c("true", "t", "1", "yes", "y")
  if (nzchar(out$show_confidence_interval)) {
    warning(
      "--show-confidence-interval is retained for compatibility; ",
      "use --show-iqr instead."
    )
    out$ribbon_type <- if (
      tolower(out$show_confidence_interval) %in% c("true", "t", "1", "yes", "y")
    ) "ci95" else "none"
  } else if (!out$show_iqr && identical(out$ribbon_type, "iqr")) {
    out$ribbon_type <- "none"
  }
  out$pdf_width <- as.numeric(out$pdf_width)
  out$pdf_height <- as.numeric(out$pdf_height)
  out$panel_aspect_ratio <- as.numeric(out$panel_aspect_ratio)
  if (!is.finite(out$panel_aspect_ratio) || out$panel_aspect_ratio <= 0) {
    stop("--panel-aspect-ratio must be greater than zero.")
  }
  if (out$permutations < 100) stop("--permutations must be at least 100.")
  out
}

norm_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  gsub("[^a-z0-9]+", "", x)
}

lifecycle_stage <- function(species, antigen) {
  key <- norm_key(antigen)

  pf_pre <- norm_key(c("TRAP", "CSP", "LSA3", "CTRP", "STARP", "CelTOS", "TREP"))
  pf_gametocyte <- norm_key(c("Pfs25", "P25", "HAP2", "Pfs230", "P28", "Pfs48/45", "CHT1"))
  pv_pre <- norm_key(c("TRAP", "CSP"))
  pv_gametocyte <- norm_key(c("P28"))

  case_when(
    species == "Pf7" & key %in% pf_pre ~ "Pre-erythrocytic",
    species == "Pf7" & key %in% pf_gametocyte ~ "Gametocyte",
    species == "Pv4" & key %in% pv_pre ~ "Pre-erythrocytic",
    species == "Pv4" & key %in% pv_gametocyte ~ "Gametocyte",
    TRUE ~ "Blood-stage"
  )
}

stage_number <- function(stage) {
  match(stage, c("Pre-erythrocytic", "Gametocyte", "Blood-stage"))
}

accumulate_mapping <- function(mapping, permutations, seed, central = "mean") {
  n_haplotypes <- nrow(mapping)
  set.seed(seed)
  accumulation <- matrix(0L, nrow = permutations, ncol = n_haplotypes)

  for (i in seq_len(permutations)) {
    ordered_serotypes <- mapping$major_iep_serotype[
      sample.int(n_haplotypes, n_haplotypes, replace = FALSE)
    ]
    accumulation[i, ] <- cumsum(
      !is.na(ordered_serotypes) & !duplicated(ordered_serotypes)
    )
  }

  total_serotypes <- dplyr::n_distinct(mapping$major_iep_serotype, na.rm = TRUE)
  curve <- tibble(
    haplotypes_sampled = seq_len(n_haplotypes),
    mean_serotypes = colMeans(accumulation),
    sd_serotypes = apply(accumulation, 2, sd),
    se_serotypes = sd_serotypes / sqrt(permutations),
    lower_se = pmax(0, mean_serotypes - se_serotypes),
    upper_se = pmin(total_serotypes, mean_serotypes + se_serotypes),
    lower_sd = pmax(0, mean_serotypes - sd_serotypes),
    upper_sd = pmin(total_serotypes, mean_serotypes + sd_serotypes),
    median_serotypes = apply(accumulation, 2, median),
    lower_iqr = apply(accumulation, 2, quantile, probs = 0.25),
    upper_iqr = apply(accumulation, 2, quantile, probs = 0.75),
    lower_95 = apply(accumulation, 2, quantile, probs = 0.025),
    upper_95 = apply(accumulation, 2, quantile, probs = 0.975)
  )

  central_column <- if (identical(central, "median")) {
    "median_serotypes"
  } else {
    "mean_serotypes"
  }
  n_for_90 <- curve %>%
    filter(.data[[central_column]] >= 0.90 * total_serotypes) %>%
    summarise(value = min(haplotypes_sampled)) %>%
    pull(value)
  if (!length(n_for_90) || !is.finite(n_for_90)) n_for_90 <- NA_integer_

  list(
    curve = curve,
    n_haplotypes = n_haplotypes,
    n_serotypes = total_serotypes,
    n_for_90 = n_for_90
  )
}

integer_breaks <- function(limits) {
  lower <- ceiling(min(limits, na.rm = TRUE))
  upper <- floor(max(limits, na.rm = TRUE))
  if (!is.finite(lower) || !is.finite(upper) || lower > upper) return(numeric())

  candidates <- unique(as.integer(round(pretty(c(lower, upper), n = 5))))
  candidates[candidates >= lower & candidates <= upper]
}

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

build_nearest_serotype_map <- function(flows, major_serotypes, max_distance) {
  raw_profiles <- flows %>%
    group_by(species, gene_id, antigen, immune_escape_haplotype) %>%
    summarise(profile_n = sum(n_samples), .groups = "drop")

  major_profiles <- raw_profiles %>%
    inner_join(
      major_serotypes,
      by = c("species", "gene_id", "antigen", "immune_escape_haplotype")
    )

  split(
    raw_profiles,
    interaction(
      raw_profiles$species,
      raw_profiles$gene_id,
      raw_profiles$antigen,
      drop = TRUE
    )
  ) %>%
    lapply(function(group_df) {
      candidates <- major_profiles %>%
        filter(
          species == group_df$species[[1]],
          gene_id == group_df$gene_id[[1]],
          antigen == group_df$antigen[[1]]
        ) %>%
        arrange(desc(profile_n), immune_escape_haplotype)

      if (!nrow(candidates)) return(NULL)

      bind_rows(lapply(seq_len(nrow(group_df)), function(i) {
        profile <- group_df$immune_escape_haplotype[[i]]
        distances <- vapply(
          candidates$immune_escape_haplotype,
          function(candidate) iep_hamming_distance(profile, candidate),
          integer(1)
        )
        best_distance <- min(distances)
        best <- candidates[distances == best_distance, , drop = FALSE] %>%
          arrange(desc(profile_n), immune_escape_haplotype) %>%
          slice(1)

        tibble(
          species = group_df$species[[i]],
          gene_id = group_df$gene_id[[i]],
          antigen = group_df$antigen[[i]],
          immune_escape_haplotype = profile,
          major_iep_serotype = if (
            is.finite(best_distance) && best_distance <= max_distance
          ) best$immune_escape_haplotype[[1]] else NA_character_,
          nearest_iep_distance = as.integer(best_distance)
        )
      }))
    }) %>%
    bind_rows()
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  flows <- read_tsv(args$raw_flows, show_col_types = FALSE)
  displayed_flows <- read_tsv(args$alluvial_flows, show_col_types = FALSE)
  required <- c(
    "species", "gene_id", "antigen",
    "all_snp_haplotype", "immune_escape_haplotype", "n_samples"
  )
  missing <- setdiff(required, names(flows))
  if (length(missing)) stop("Raw flow table is missing: ", paste(missing, collapse = ", "))

  displayed_required <- c(
    "species", "gene_id", "antigen",
    "all_snp_haplotype", "immune_escape_haplotype",
    "n_all_snp_haplotypes_raw"
  )
  displayed_missing <- setdiff(displayed_required, names(displayed_flows))
  if (length(displayed_missing)) {
    stop("Alluvial flow table is missing: ", paste(displayed_missing, collapse = ", "))
  }

  selected_panels <- displayed_flows %>%
    distinct(species, gene_id, antigen)

  requested_panels <- bind_rows(
    if (length(args$pf_antigens)) {
      tibble(
        species = rep("Pf7", length(args$pf_antigens)),
        antigen = args$pf_antigens,
        requested_order = seq_along(args$pf_antigens)
      )
    },
    if (length(args$pv_antigens)) {
      tibble(
        species = rep("Pv4", length(args$pv_antigens)),
        antigen = args$pv_antigens,
        requested_order = seq_along(args$pv_antigens)
      )
    }
  )
  if (nrow(requested_panels)) {
    missing_requested <- requested_panels %>%
      anti_join(selected_panels, by = c("species", "antigen"))
    if (nrow(missing_requested)) {
      stop(
        "Requested antigen(s) absent from the alluvial flow table: ",
        paste(paste(missing_requested$species, missing_requested$antigen, sep = ":"), collapse = ", ")
      )
    }
    selected_panels <- selected_panels %>%
      semi_join(requested_panels, by = c("species", "antigen"))
  }

  alluvial_counts <- displayed_flows %>%
    group_by(species, gene_id, antigen) %>%
    summarise(
      alluvial_haplotype_count = first(n_all_snp_haplotypes_raw),
      alluvial_serotype_count = n_distinct(
        immune_escape_haplotype[
          !immune_escape_haplotype %in% c("__FILTERED_IEP__", "__OTHER__")
        ]
      ),
      rare_original_n = sum(
        n_samples[all_snp_haplotype == "__OTHER__"],
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

  eligible_panels <- alluvial_counts %>%
    filter(
      alluvial_serotype_count > args$min_serotypes_exclusive,
      alluvial_serotype_count < args$max_serotypes_exclusive,
      !is.na(rare_original_haplotype_percent),
      rare_original_haplotype_percent <=
        args$max_rare_original_haplotype_percent
    ) %>%
    select(species, gene_id, antigen)

  if (!nrow(eligible_panels)) {
    stop(
      "No antigens have major-serotype counts > ",
      args$min_serotypes_exclusive,
      " and < ",
      args$max_serotypes_exclusive,
      "."
    )
  }

  major_serotypes <- displayed_flows %>%
    semi_join(eligible_panels, by = c("species", "gene_id", "antigen")) %>%
    filter(!immune_escape_haplotype %in% c("__FILTERED_IEP__", "__OTHER__")) %>%
    distinct(species, gene_id, antigen, immune_escape_haplotype) %>%
    mutate(major_iep_serotype = immune_escape_haplotype)

  nearest_serotype_map <- build_nearest_serotype_map(
    flows %>%
      semi_join(selected_panels, by = c("species", "gene_id", "antigen")) %>%
      semi_join(eligible_panels, by = c("species", "gene_id", "antigen")),
    major_serotypes,
    args$max_iep_hamming_distance
  )

  mapping <- flows %>%
    semi_join(selected_panels, by = c("species", "gene_id", "antigen")) %>%
    semi_join(eligible_panels, by = c("species", "gene_id", "antigen")) %>%
    filter(
      !is.na(all_snp_haplotype),
      nzchar(all_snp_haplotype),
      !is.na(immune_escape_haplotype),
      nzchar(immune_escape_haplotype)
    ) %>%
    left_join(
      nearest_serotype_map,
      by = c("species", "gene_id", "antigen", "immune_escape_haplotype")
    ) %>%
    group_by(species, gene_id, antigen, all_snp_haplotype) %>%
    arrange(desc(n_samples), immune_escape_haplotype, .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(
      species, gene_id, antigen, all_snp_haplotype,
      immune_escape_haplotype, major_iep_serotype, nearest_iep_distance
    ) %>%
    left_join(
      alluvial_counts,
      by = c("species", "gene_id", "antigen")
    )

  groups <- split(
    mapping,
    interaction(mapping$species, mapping$gene_id, mapping$antigen, drop = TRUE)
  )
  if (!length(groups)) stop("No antigen mappings were available.")

  results <- lapply(seq_along(groups), function(i) {
    d <- groups[[i]]
    expected_haplotypes <- unique(d$alluvial_haplotype_count)
    expected_serotypes <- unique(d$alluvial_serotype_count)
    if (
      length(expected_haplotypes) != 1 ||
        n_distinct(d$all_snp_haplotype) != expected_haplotypes
    ) {
      stop(
        "Raw/alluvial haplotype count mismatch for ",
        unique(d$species), " ", unique(d$antigen), "."
      )
    }
    result <- accumulate_mapping(
      d,
      args$permutations,
      args$seed + i,
      central = if (identical(args$ribbon_type, "iqr")) "median" else "mean"
    )
    if (length(expected_serotypes) != 1 || result$n_serotypes != expected_serotypes) {
      stop(
        "Raw/alluvial serotype count mismatch for ",
        unique(d$species), " ", unique(d$antigen), "."
      )
    }
    panel <- paste(unique(d$species), unique(d$antigen), sep = " ")

    curve <- result$curve %>%
      mutate(
        species = unique(d$species),
        gene_id = unique(d$gene_id),
        antigen = unique(d$antigen),
        panel = panel,
        total_haplotypes = result$n_haplotypes,
        total_serotypes = result$n_serotypes,
        haplotypes_for_90_percent_serotypes = result$n_for_90
      )

    summary <- tibble(
      species = unique(d$species),
      gene_id = unique(d$gene_id),
      antigen = unique(d$antigen),
      total_haplotypes = result$n_haplotypes,
      total_serotypes = result$n_serotypes,
      haplotype_to_serotype_reduction_percent =
        100 * (1 - result$n_serotypes / result$n_haplotypes),
      haplotypes_for_90_percent_serotypes = result$n_for_90,
      min_original_haplotype_percent = args$min_original_haplotype_percent,
      min_reduced_haplotype_percent = args$min_reduced_haplotype_percent,
      max_iep_hamming_distance = args$max_iep_hamming_distance,
      rare_original_haplotype_percent =
        unique(d$rare_original_haplotype_percent),
      permutations = args$permutations
    )
    list(curve = curve, summary = summary)
  })

  curve_df <- bind_rows(lapply(results, `[[`, "curve"))
  summary_df <- bind_rows(lapply(results, `[[`, "summary"))

  panel_order <- summary_df %>%
    mutate(
      lifecycle_stage = lifecycle_stage(species, antigen),
      stage_order = stage_number(lifecycle_stage),
      species_order = match(species, c("Pf7", "Pv4")),
      panel = paste(species, antigen),
      panel_label = paste0(
        if_else(species == "Pf7", "Pf", "Pv"),
        antigen
      )
    ) %>%
    arrange(
      species_order,
      stage_order,
      desc(total_serotypes),
      desc(total_haplotypes),
      antigen
    )

  panel_labels <- setNames(panel_order$panel_label, panel_order$panel)
  curve_df$panel <- factor(
    panel_labels[curve_df$panel],
    levels = panel_order$panel_label
  )

  summary_df <- summary_df %>%
    left_join(
      panel_order %>%
        select(species, gene_id, antigen, lifecycle_stage, stage_order),
      by = c("species", "gene_id", "antigen")
    ) %>%
    arrange(
      match(species, c("Pf7", "Pv4")),
      stage_order,
      desc(total_serotypes),
      desc(total_haplotypes),
      antigen
    )

  species_colours <- c("Pf7" = "#B2182B", "Pv4" = "#2166AC")
  common_x_limits <- c(0, max(curve_df$haplotypes_sampled, na.rm = TRUE))
  common_y_limits <- c(0, max(curve_df$total_serotypes, na.rm = TRUE))
  common_x_breaks <- integer_breaks(common_x_limits)
  common_y_breaks <- integer_breaks(common_y_limits)

  make_antigen_plot <- function(plot_data) {
    species_value <- unique(plot_data$species)
    panel_title <- unique(as.character(plot_data$panel))
    central_value <- if (identical(args$ribbon_type, "iqr")) {
      "median_serotypes"
    } else {
      "mean_serotypes"
    }

    plot_out <- ggplot(
      plot_data,
      aes(x = haplotypes_sampled, y = .data[[central_value]])
    )

    if (!identical(args$ribbon_type, "none")) {
      ribbon_limits <- switch(
        args$ribbon_type,
        se = c("lower_se", "upper_se"),
        sd = c("lower_sd", "upper_sd"),
        iqr = c("lower_iqr", "upper_iqr"),
        ci95 = c("lower_95", "upper_95")
      )
      plot_out <- plot_out + geom_ribbon(
        aes(
          ymin = .data[[ribbon_limits[[1]]]],
          ymax = .data[[ribbon_limits[[2]]]],
          fill = species
        ),
        alpha = if (identical(args$ribbon_type, "se")) 0.16 else 0.22,
        colour = NA
      )
    }

    plot_out <- plot_out +
      geom_line(aes(colour = species), linewidth = 0.75) +
      geom_hline(
        aes(yintercept = total_serotypes),
        linewidth = 0.28,
        linetype = "dashed",
        colour = "grey45"
      )

    if (args$point_size > 0) {
      plot_out <- plot_out + geom_point(
        aes(colour = species),
        size = args$point_size,
        alpha = 0.7
      )
    }

    plot_out +
      scale_colour_manual(values = species_colours, guide = "none") +
      scale_fill_manual(values = species_colours, guide = "none") +
      scale_x_continuous(
        limits = common_x_limits,
        breaks = common_x_breaks,
        labels = function(x) format(as.integer(x), scientific = FALSE)
      ) +
      scale_y_continuous(
        limits = common_y_limits,
        breaks = common_y_breaks,
        labels = function(x) format(as.integer(x), scientific = FALSE)
      ) +
      labs(
        x = "Number of haplotypes",
        y = "Number of putative serotypes",
        title = panel_title
      ) +
      theme_classic(base_family = args$font_family, base_size = 6) +
      theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.35),
        axis.text = element_text(size = 6, colour = "black"),
        axis.title = element_text(size = 7, colour = "black"),
        plot.title = element_text(
          face = "bold",
          size = 7,
          colour = "black",
          hjust = 0.5,
          margin = margin(b = 3)
        ),
        aspect.ratio = args$panel_aspect_ratio,
        plot.margin = margin(2, 2, 2, 2)
      )
  }

  species_values <- c("Pf7", "Pv4")
  species_values <- species_values[
    vapply(
      species_values,
      function(species_value) any(curve_df$species == species_value),
      logical(1)
    )
  ]
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "The patchwork R package is required to place Pv on a new row. ",
      "Install it with install.packages('patchwork')."
    )
  }

  panel_plots <- list()
  for (species_value in species_values) {
    species_panels <- panel_order %>%
      filter(species == species_value) %>%
      pull(panel_label)

    for (panel_value in species_panels) {
      panel_data <- curve_df %>%
        filter(as.character(panel) == panel_value)
      panel_plots[[length(panel_plots) + 1L]] <- make_antigen_plot(panel_data)
    }

    n_padding <- (
      args$facet_ncol - length(species_panels) %% args$facet_ncol
    ) %% args$facet_ncol
    if (n_padding > 0) {
      panel_plots <- c(
        panel_plots,
        rep(list(patchwork::plot_spacer()), n_padding)
      )
    }
  }

  p <- patchwork::wrap_plots(
    panel_plots,
    ncol = args$facet_ncol,
    widths = rep(1, args$facet_ncol),
    axes = "keep",
    axis_titles = "collect"
  )

  ggsave(
    args$output,
    p,
    width = args$pdf_width,
    height = args$pdf_height,
    device = cairo_pdf
  )
  message("Wrote: ", args$output)

  if (nzchar(args$output_svg)) {
    ggsave(args$output_svg, p, width = args$pdf_width, height = args$pdf_height)
    message("Wrote: ", args$output_svg)
  }

  summary_output <- if (nzchar(args$summary_output)) {
    args$summary_output
  } else {
    paste0(tools::file_path_sans_ext(args$output), "_summary.tsv")
  }
  write_tsv(summary_df, summary_output)
  message("Wrote: ", summary_output)

  curve_output <- paste0(tools::file_path_sans_ext(summary_output), "_curves.tsv")
  write_tsv(curve_df, curve_output)
  message("Wrote: ", curve_output)
}

main()
