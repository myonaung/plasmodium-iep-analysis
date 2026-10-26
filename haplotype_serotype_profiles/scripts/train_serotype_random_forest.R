#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ranger)
  library(grid)
})

parse_args <- function(args) {
  out <- list(
    raw_flows = "",
    alluvial_flows = "",
    iep_position_map = "",
    out_prefix = "iep_serotype_rf",
    max_iep_hamming_distance = 4,
    min_serotypes_exclusive = 1,
    max_serotypes_exclusive = 999,
    max_rare_original_haplotype_percent = 100,
    probability_threshold = 0.70,
    trees = 1000,
    seed = 1,
    top_n_plot = 5,
    plot_width = NA,
    plot_height = NA,
    plot_facet_width = 2.8,
    plot_facet_height = 2.8,
    plot_pv_height_ratio = 1,
    font_family = "Arial",
    facet_ncol = 4,
    plot_panels = paste(
      c(
        "PfTRAP", "PfCSP", "PfRH4", "PfAMA1",
        "PfMSP6", "PfRESA", "PfRIPR",
        "PvTRAP", "PvAMA1", "PvCyRPA", "PvMSP9",
        "PvMSP8", "PvP41"
      ),
      collapse = ","
    ),
    plot_x_max = NA
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
  if (!nzchar(out$out_prefix)) stop("--out-prefix is required.")

  out$max_iep_hamming_distance <- as.integer(out$max_iep_hamming_distance)
  out$min_serotypes_exclusive <- as.numeric(out$min_serotypes_exclusive)
  out$max_serotypes_exclusive <- as.numeric(out$max_serotypes_exclusive)
  out$max_rare_original_haplotype_percent <-
    as.numeric(out$max_rare_original_haplotype_percent)
  out$probability_threshold <- as.numeric(out$probability_threshold)
  out$trees <- as.integer(out$trees)
  out$seed <- as.integer(out$seed)
  out$top_n_plot <- as.integer(out$top_n_plot)
  out$plot_width <- as.numeric(out$plot_width)
  out$plot_height <- as.numeric(out$plot_height)
  out$plot_facet_width <- as.numeric(out$plot_facet_width)
  out$plot_facet_height <- as.numeric(out$plot_facet_height)
  out$plot_pv_height_ratio <- as.numeric(out$plot_pv_height_ratio)
  out$font_family <- as.character(out$font_family)
  out$facet_ncol <- as.integer(out$facet_ncol)
  out$plot_panels <- trimws(strsplit(out$plot_panels, ",", fixed = TRUE)[[1]])
  out$plot_panels <- out$plot_panels[nzchar(out$plot_panels)]
  out$plot_x_max <- suppressWarnings(as.numeric(out$plot_x_max))
  out
}

norm_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  gsub("[^a-z0-9]+", "", x)
}

read_delim_auto <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    read_csv(path, show_col_types = FALSE, guess_max = 300000)
  } else {
    read_tsv(path, show_col_types = FALSE, guess_max = 300000)
  }
}

pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) hit[[1]] else NA_character_
}

read_iep_position_map <- function(path) {
  if (!nzchar(path)) return(tibble())

  df <- read_delim_auto(path)
  position_col <- pick_col(df, c("aa_position", "position", "site", "pos"))
  gene_col <- pick_col(df, c("gene_id", "Gene_ID", "gene", "plasmodb_id"))
  antigen_col <- pick_col(df, c("antigen", "gene_name", "Gene_Name_annotated", "name"))
  species_col <- pick_col(df, c("species", "Species"))
  order_col <- pick_col(df, c("position_order", "iep_index", "site_order"))
  feature_col <- pick_col(df, c("feature", "iep_feature"))

  if (is.na(position_col)) {
    stop("--iep-position-map must contain an aa_position/position/site/pos column.")
  }
  if (is.na(gene_col) && is.na(antigen_col)) {
    stop("--iep-position-map must contain gene_id and/or antigen.")
  }

  out <- df %>%
    transmute(
      species = if (!is.na(species_col)) as.character(.data[[species_col]]) else NA_character_,
      gene_id = if (!is.na(gene_col)) as.character(.data[[gene_col]]) else NA_character_,
      antigen = if (!is.na(antigen_col)) as.character(.data[[antigen_col]]) else NA_character_,
      gene_key = if (!is.na(gene_col)) norm_key(.data[[gene_col]]) else NA_character_,
      antigen_key = if (!is.na(antigen_col)) norm_key(.data[[antigen_col]]) else NA_character_,
      aa_position = as.integer(as.numeric(.data[[position_col]])),
      iep_index = if (!is.na(order_col)) as.integer(as.numeric(.data[[order_col]])) else NA_integer_,
      feature = if (!is.na(feature_col)) as.character(.data[[feature_col]]) else NA_character_
    ) %>%
    filter(!is.na(aa_position)) %>%
    distinct(species, gene_id, antigen, gene_key, antigen_key, aa_position, iep_index, feature) %>%
    arrange(species, gene_key, antigen_key, coalesce(iep_index, aa_position), aa_position) %>%
    group_by(species, gene_key, antigen_key) %>%
    mutate(
      iep_index = coalesce(iep_index, row_number()),
      feature = if_else(is.na(feature) | !nzchar(feature), sprintf("iep_%03d", iep_index), feature),
      feature_label = paste0("aa_", aa_position)
    ) %>%
    ungroup()

  out
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
    summarise(profile_n = sum(n_samples, na.rm = TRUE), .groups = "drop")

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

split_iep_profile <- function(x) {
  pieces <- strsplit(as.character(x), "|", fixed = TRUE)
  max_len <- max(lengths(pieces))
  mat <- do.call(rbind, lapply(pieces, function(v) {
    length(v) <- max_len
    v[is.na(v) | !nzchar(v)] <- "__MISSING__"
    v
  }))
  out <- as_tibble(mat, .name_repair = "minimal")
  names(out) <- sprintf("iep_%03d", seq_len(ncol(out)))
  mutate(out, across(everything(), factor))
}

lookup_positions_for_group <- function(feature_lookup, group_df) {
  if (!nrow(feature_lookup)) return(tibble())

  group_species <- first(group_df$species)
  group_gene_key <- norm_key(first(group_df$gene_id))
  group_antigen_key <- norm_key(first(group_df$antigen))

  hit <- feature_lookup %>%
    filter(
      (is.na(species) | species == group_species),
      (is.na(gene_key) | gene_key == group_gene_key),
      (is.na(antigen_key) | antigen_key == group_antigen_key)
    ) %>%
    distinct(feature, aa_position, feature_label)

  if (!nrow(hit)) {
    hit <- feature_lookup %>%
      filter(
        (is.na(species) | species == group_species),
        gene_key == group_gene_key | antigen_key == group_antigen_key
      ) %>%
      distinct(feature, aa_position, feature_label)
  }

  hit
}

fit_one_antigen <- function(group_df, trees, seed, probability_threshold, feature_lookup) {
  group_df <- group_df %>%
    filter(
      !is.na(major_iep_serotype),
      nzchar(major_iep_serotype),
      !major_iep_serotype %in% c("__FILTERED_IEP__", "__OTHER__")
    )

  n_serotypes <- n_distinct(group_df$major_iep_serotype)
  if (nrow(group_df) < 3 || n_serotypes < 2) {
    empty_assignments <- group_df %>%
      mutate(
        rf_predicted_serotype = NA_character_,
        rf_max_probability = NA_real_,
        rf_assignment_status = "RF_NOT_FIT"
      )
    empty_perf <- tibble(
      species = first(group_df$species),
      gene_id = first(group_df$gene_id),
      antigen = first(group_df$antigen),
      n_haplotypes = nrow(group_df),
      n_serotypes = n_serotypes,
      rf_oob_accuracy = NA_real_,
      rf_status = "not_fit_fewer_than_two_serotypes"
    )
    return(list(assignments = empty_assignments, performance = empty_perf, importance = tibble()))
  }

  features <- split_iep_profile(group_df$immune_escape_haplotype)
  model_df <- bind_cols(
    tibble(major_iep_serotype = factor(group_df$major_iep_serotype)),
    features
  )

  set.seed(seed)
  fit <- ranger(
    major_iep_serotype ~ .,
    data = model_df,
    num.trees = trees,
    probability = TRUE,
    importance = "impurity",
    case.weights = pmax(1, group_df$n_samples),
    classification = TRUE,
    seed = seed
  )

  pred <- predict(fit, data = model_df)$predictions
  predicted_serotype <- colnames(pred)[max.col(pred, ties.method = "first")]
  max_probability <- apply(pred, 1, max)

  assignments <- group_df %>%
    mutate(
      rf_predicted_serotype = predicted_serotype,
      rf_max_probability = max_probability,
      rf_assignment_status = case_when(
        rf_predicted_serotype == major_iep_serotype &
          rf_max_probability >= probability_threshold ~ "RF_MATCH_HIGH_CONFIDENCE",
        rf_predicted_serotype == major_iep_serotype ~ "RF_MATCH_LOW_CONFIDENCE",
        rf_max_probability >= probability_threshold ~ "RF_DIFFERENT_HIGH_CONFIDENCE",
        TRUE ~ "RF_DIFFERENT_LOW_CONFIDENCE"
      )
    )

  performance <- assignments %>%
    summarise(
      species = first(species),
      gene_id = first(gene_id),
      antigen = first(antigen),
      n_haplotypes = n(),
      n_serotypes = n_distinct(major_iep_serotype),
      rf_oob_accuracy = 1 - fit$prediction.error,
      rf_status = "fit",
      .groups = "drop"
    )

  position_lookup <- lookup_positions_for_group(feature_lookup, group_df)

  importance <- tibble(
    species = first(group_df$species),
    gene_id = first(group_df$gene_id),
    antigen = first(group_df$antigen),
    feature = names(fit$variable.importance),
    importance = as.numeric(fit$variable.importance)
  ) %>%
    left_join(position_lookup, by = "feature") %>%
    mutate(
      feature_label = coalesce(feature_label, feature)
    ) %>%
    arrange(desc(importance), feature)

  list(assignments = assignments, performance = performance, importance = importance)
}

preferred_panel_levels <- function(labels) {
  preferred <- c(
    "PfTRAP", "PfCSP", "PfRH4", "PfAMA1",
    "PfMSP6", "PfRESA", "PfRIPR",
    "PvTRAP", "PvAMA1", "PvCyRPA", "PvMSP9",
    "PvMSP8", "PvP41"
  )
  c(preferred[preferred %in% labels], sort(setdiff(labels, preferred)))
}

make_importance_block <- function(
  plot_df, panels, x_max, facet_ncol, font_family,
  title = NULL, x_label = NULL, y_label = NULL
) {
  block_df <- plot_df %>%
    filter(as.character(panel) %in% panels) %>%
    mutate(panel = factor(as.character(panel), levels = panels))
  if (!nrow(block_df)) return(NULL)

  panel_labels <- setNames(panels, panels)
  panel_labels[startsWith(names(panel_labels), "__blank_panel__")] <- ""
  fixed_x <- is.finite(x_max) && x_max > 0
  facet_scales <- if (fixed_x) "free_y" else "free"
  coord_layer <- if (fixed_x) {
    coord_cartesian(xlim = c(0, x_max), expand = FALSE)
  } else {
    coord_cartesian(expand = FALSE)
  }
  x_scale <- if (fixed_x) {
    scale_x_continuous(
      breaks = pretty(c(0, x_max), n = 4),
      expand = expansion(mult = c(0, 0.02))
    )
  } else {
    scale_x_continuous(expand = expansion(mult = c(0, 0.04)))
  }

  ggplot(block_df, aes(x = importance, y = y_key, fill = species_label)) +
    geom_blank(data = block_df %>% mutate(importance = 0)) +
    geom_col(
      data = block_df %>% filter(!is_placeholder),
      width = 0.48,
      colour = "black",
      linewidth = 0.18
    ) +
    facet_wrap(
      ~ panel,
      scales = facet_scales,
      ncol = facet_ncol,
      drop = FALSE,
      labeller = as_labeller(panel_labels)
    ) +
    scale_fill_manual(
      values = c("Pf7" = "#C51B3A", "Pv4" = "#2166AC"),
      guide = "none"
    ) +
    scale_y_discrete(
      labels = function(x) {
        label_lookup <- setNames(as.character(block_df$position_label), as.character(block_df$y_key))
        unname(label_lookup[x])
      }
    ) +
    coord_layer +
    x_scale +
    labs(x = x_label, y = y_label, title = title) +
    theme_classic(base_size = 12, base_family = font_family) +
    theme(
      plot.title = element_text(face = "bold", size = 24, hjust = 0),
      axis.title.x = element_text(size = 14, colour = "black", margin = margin(t = 4)),
      axis.title.y = element_text(size = 14, colour = "black", margin = margin(r = 5)),
      axis.text.x = element_text(size = 12, colour = "black"),
      axis.text.y = element_text(size = 12, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "grey35"),
      axis.line = element_line(linewidth = 0.25, colour = "grey35"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 14, colour = "black", margin = margin(0, 0, 3, 0)),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.35),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.42, "lines"),
      panel.spacing.y = unit(0.38, "lines"),
      plot.margin = margin(2, 2, 10, 2)
    )
}

pad_panels <- function(panels, facet_ncol, prefix) {
  n_padding <- (facet_ncol - length(panels) %% facet_ncol) %% facet_ncol
  if (n_padding <= 0) return(panels)
  c(panels, paste0("__blank_panel__", prefix, "_", seq_len(n_padding)))
}

make_importance_panel <- function(panel_df, panel_id, x_max, font_family) {
  fixed_x <- is.finite(x_max) && x_max > 0
  coord_layer <- if (fixed_x) {
    coord_cartesian(xlim = c(0, x_max), expand = FALSE)
  } else {
    coord_cartesian(expand = FALSE)
  }
  x_scale <- if (fixed_x) {
    scale_x_continuous(
      breaks = pretty(c(0, x_max), n = 4),
      expand = expansion(mult = c(0, 0.02))
    )
  } else {
    scale_x_continuous(expand = expansion(mult = c(0, 0.04)))
  }

  panel_df <- panel_df %>%
    mutate(
      y_slot = factor(slot, levels = rev(sort(unique(slot))))
    )
  label_lookup <- setNames(as.character(panel_df$position_label), as.character(panel_df$y_slot))

  ggplot(panel_df, aes(x = importance, y = y_slot, fill = species_label)) +
    geom_blank(data = panel_df %>% mutate(importance = 0)) +
    geom_col(
      data = panel_df %>% filter(!is_placeholder),
      width = 0.48,
      colour = "black",
      linewidth = 0.18
    ) +
    scale_fill_manual(
      values = c("Pf7" = "#C51B3A", "Pv4" = "#2166AC"),
      guide = "none"
    ) +
    scale_y_discrete(labels = function(x) unname(label_lookup[x])) +
    coord_layer +
    x_scale +
    labs(
      x = "Random forest importance",
      y = "IEP amino-acid position",
      title = panel_id
    ) +
    theme_classic(base_size = 12, base_family = font_family) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.35),
      axis.text = element_text(size = 12, colour = "black"),
      axis.title = element_text(size = 14, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "grey35"),
      axis.line = element_line(linewidth = 0.25, colour = "grey35"),
      plot.title = element_text(
        face = "bold",
        size = 14,
        colour = "black",
        hjust = 0.5,
        margin = margin(b = 3)
      ),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      aspect.ratio = 1,
      plot.margin = margin(2, 2, 2, 2)
    )
}

save_importance_patchwork_plot <- function(
  plot_df, out_path, top_n, width, height, facet_ncol, x_max,
  facet_width, facet_height, font_family
) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "The patchwork R package is required for aligned square RF panels. ",
      "Install it with install.packages('patchwork')."
    )
  }

  fixed_x <- is.finite(x_max) && x_max > 0
  if (!fixed_x) x_max <- NA_real_

  all_panels <- levels(plot_df$panel)
  pf_panels <- all_panels[startsWith(all_panels, "Pf")]
  pv_panels <- all_panels[startsWith(all_panels, "Pv")]

  build_panel_plots <- function(panels) {
    out <- lapply(panels, function(panel_id) {
      make_importance_panel(
        plot_df %>% filter(as.character(panel) == panel_id),
        panel_id,
        x_max,
        font_family
      )
    })
    n_padding <- (facet_ncol - length(panels) %% facet_ncol) %% facet_ncol
    if (n_padding > 0) {
      out <- c(out, rep(list(patchwork::plot_spacer()), n_padding))
    }
    out
  }

  panel_plots <- c(build_panel_plots(pf_panels), build_panel_plots(pv_panels))
  n_rows <- ceiling(length(pf_panels) / facet_ncol) +
    ceiling(length(pv_panels) / facet_ncol)

  if (!is.finite(width) || width <= 0) {
    width <- 0.6 + facet_ncol * facet_width
  }
  if (!is.finite(height) || height <= 0) {
    height <- 0.9 + n_rows * facet_height
  }

  p <- patchwork::wrap_plots(
    panel_plots,
    ncol = facet_ncol,
    widths = rep(1, facet_ncol),
    axes = "keep",
    axis_titles = "collect"
  ) +
    patchwork::plot_annotation(
      title = "IEP positions discriminating putative serotypes",
      theme = theme(
        plot.title = element_text(
          family = font_family,
          face = "bold",
          size = 24,
          hjust = 0
        )
      )
    )

  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(
    paste0(out_path, ".pdf"),
    p,
    width = width,
    height = height,
    device = pdf_device,
    family = if (capabilities("cairo")) font_family else "sans"
  )

  png_device <- function(filename, width, height, ...) {
    grDevices::png(
      filename,
      width = width,
      height = height,
      units = "in",
      res = 300,
      type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"),
      ...
    )
  }
  ggplot2::ggsave(
    paste0(out_path, ".png"),
    p,
    width = width,
    height = height,
    device = png_device
  )

  ggplot2::ggsave(
    paste0(out_path, ".svg"),
    p,
    width = width,
    height = height,
    device = grDevices::svg
  )
}

draw_merged_importance_plot <- function(pf_plot, pv_plot, pf_rows, pv_rows, pv_height_ratio) {
  pv_width_ratio <- 1
  if (!is.null(attr(pv_plot, "row_width_ratio"))) {
    pv_width_ratio <- attr(pv_plot, "row_width_ratio")
  }
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 2,
    ncol = 1,
    heights = unit(c(pf_rows, pv_rows * pv_height_ratio), "null")
  )))
  if (!is.null(pf_plot)) {
    print(pf_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  }
  if (!is.null(pv_plot)) {
    print(
      pv_plot,
      vp = viewport(
        layout.pos.row = 2,
        layout.pos.col = 1,
        x = 0,
        width = pv_width_ratio,
        just = "left"
      )
    )
  }
  popViewport()
}

save_importance_plot <- function(
  plot_df, out_path, top_n, width, height, facet_ncol, x_max,
  facet_width, facet_height, pv_height_ratio, font_family
) {
  if (!nrow(plot_df)) return(invisible(NULL))
  fixed_x <- is.finite(x_max) && x_max > 0
  if (!fixed_x) x_max <- NA_real_

  pf_panels <- levels(plot_df$panel)[startsWith(levels(plot_df$panel), "Pf")]
  pv_panels <- levels(plot_df$panel)[startsWith(levels(plot_df$panel), "Pv")]
  pf_panels_plot <- pf_panels
  pv_panels_plot <- pv_panels
  pf_rows <- max(1, ceiling(length(pf_panels_plot) / facet_ncol))
  pv_rows <- max(1, ceiling(length(pv_panels_plot) / facet_ncol))
  used_cols <- max(1, min(facet_ncol, max(length(pf_panels), length(pv_panels))))
  used_rows <- pf_rows + pv_rows
  pv_width_ratio <- min(1, length(pv_panels_plot) / max(1, min(facet_ncol, length(pf_panels_plot))))

  if (!is.finite(width) || width <= 0) {
    width <- 0.6 + used_cols * facet_width
  }
  if (!is.finite(height) || height <= 0) {
    height <- 0.9 + (pf_rows * facet_height) + (pv_rows * facet_height * pv_height_ratio)
  }

  pf_plot <- make_importance_block(
    plot_df, pf_panels_plot, x_max, facet_ncol, font_family,
    title = "IEP positions discriminating putative serotypes",
    x_label = NULL,
    y_label = "IEP amino-acid position"
  )
  pv_plot <- make_importance_block(
    plot_df, pv_panels_plot, x_max, facet_ncol, font_family,
    title = NULL,
    x_label = "Random forest importance",
    y_label = "IEP amino-acid position"
  )
  if (!is.null(pv_plot)) {
    attr(pv_plot, "row_width_ratio") <- pv_width_ratio
  }

  if (capabilities("cairo")) {
    grDevices::cairo_pdf(
      paste0(out_path, ".pdf"),
      width = width,
      height = height,
      family = font_family,
      onefile = TRUE
    )
  } else {
    pdf(
      paste0(out_path, ".pdf"),
      width = width,
      height = height,
      onefile = TRUE,
      useDingbats = FALSE,
      family = "sans"
    )
  }
  draw_merged_importance_plot(pf_plot, pv_plot, pf_rows, pv_rows, pv_height_ratio)
  dev.off()

  png(
    paste0(out_path, ".png"),
    width = width,
    height = height,
    units = "in",
    res = 300,
    type = if (capabilities("cairo")) "cairo" else getOption("bitmapType")
  )
  draw_merged_importance_plot(pf_plot, pv_plot, pf_rows, pv_rows, pv_height_ratio)
  dev.off()

  svg(paste0(out_path, ".svg"), width = width, height = height)
  draw_merged_importance_plot(pf_plot, pv_plot, pf_rows, pv_rows, pv_height_ratio)
  dev.off()
}

plot_importance <- function(
  importance, performance, out_prefix, top_n, width, height, facet_ncol,
  plot_panels, x_max, facet_width, facet_height, pv_height_ratio, font_family
) {
  if (!nrow(importance)) return(invisible(NULL))

  real_df <- importance %>%
    left_join(
      performance %>% select(species, gene_id, antigen, n_serotypes, rf_oob_accuracy),
      by = c("species", "gene_id", "antigen")
    ) %>%
    mutate(
      species_label = if_else(species == "Pf7", "Pf7", "Pv4"),
      species_prefix = if_else(species == "Pf7", "Pf", "Pv"),
      panel = paste0(species_prefix, antigen),
      position_label = if_else(
        !is.na(aa_position),
        paste0("aa ", aa_position),
        feature
      )
    ) %>%
    filter(panel %in% plot_panels) %>%
    mutate(panel = factor(panel, levels = plot_panels)) %>%
    group_by(panel) %>%
    arrange(desc(importance), position_label, .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    mutate(slot = row_number()) %>%
    ungroup()

  if (!nrow(real_df)) return(invisible(NULL))

  slot_grid <- tidyr::expand_grid(
    panel = factor(plot_panels, levels = plot_panels),
    slot = seq_len(top_n)
  ) %>%
    mutate(
      species_label = case_when(
        startsWith(as.character(panel), "Pf") ~ "Pf7",
        startsWith(as.character(panel), "Pv") ~ "Pv4",
        TRUE ~ "Pf7"
      )
    )

  y_levels <- unlist(
    lapply(plot_panels, function(panel_id) {
      paste(panel_id, rev(seq_len(top_n)), sep = "___")
    }),
    use.names = FALSE
  )

  plot_df <- slot_grid %>%
    left_join(
      real_df %>%
        select(panel, slot, importance, position_label),
      by = c("panel", "slot")
    ) %>%
    mutate(
      is_placeholder = is.na(importance),
      position_label = if_else(is_placeholder, "", position_label),
      y_key = factor(
        paste(as.character(panel), slot, sep = "___"),
        levels = y_levels
      )
    )

  base <- paste0(out_prefix, ".rf_iep_position_importance_top", top_n)
  save_importance_patchwork_plot(
    plot_df, base, top_n, width, height, facet_ncol, x_max,
    facet_width, facet_height, font_family
  )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(dirname(args$out_prefix), recursive = TRUE, showWarnings = FALSE)

  flows <- read_tsv(args$raw_flows, show_col_types = FALSE)
  displayed_flows <- read_tsv(args$alluvial_flows, show_col_types = FALSE)
  feature_lookup <- read_iep_position_map(args$iep_position_map)

  required <- c(
    "species", "gene_id", "antigen",
    "all_snp_haplotype", "immune_escape_haplotype", "n_samples"
  )
  missing <- setdiff(required, names(flows))
  if (length(missing)) {
    stop("Raw flow table is missing: ", paste(missing, collapse = ", "))
  }

  displayed_required <- c(
    "species", "gene_id", "antigen",
    "all_snp_haplotype", "immune_escape_haplotype",
    "n_all_snp_haplotypes_raw", "raw_samples_denominator"
  )
  displayed_missing <- setdiff(displayed_required, names(displayed_flows))
  if (length(displayed_missing)) {
    stop("Alluvial flow table is missing: ", paste(displayed_missing, collapse = ", "))
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

  raw_counts <- flows %>%
    group_by(species, gene_id, antigen) %>%
    summarise(
      raw_flow_rows = n(),
      raw_haplotype_count = n_distinct(all_snp_haplotype),
      raw_serotype_count = n_distinct(
        immune_escape_haplotype[
          !immune_escape_haplotype %in% c("__FILTERED_IEP__", "__OTHER__")
        ]
      ),
      raw_n_samples = sum(n_samples, na.rm = TRUE),
      .groups = "drop"
    )

  map_counts <- feature_lookup %>%
    group_by(species, gene_id, antigen) %>%
    summarise(n_iep_positions_in_map = n_distinct(aa_position), .groups = "drop")

  requested_panels <- tibble(
    requested_panel = args$plot_panels,
    species = if_else(startsWith(args$plot_panels, "Pf"), "Pf7", "Pv4"),
    antigen = sub("^(Pf|Pv)", "", args$plot_panels)
  ) %>%
    filter(grepl("^(Pf|Pv)", requested_panel))

  panel_audit_base <- bind_rows(
    map_counts %>% mutate(source = "position_map"),
    raw_counts %>%
      select(species, gene_id, antigen) %>%
      mutate(n_iep_positions_in_map = NA_integer_, source = "raw_flows"),
    alluvial_counts %>%
      select(species, gene_id, antigen) %>%
      mutate(n_iep_positions_in_map = NA_integer_, source = "alluvial_flows")
  ) %>%
    group_by(species, gene_id, antigen) %>%
    summarise(
      n_iep_positions_in_map = max(n_iep_positions_in_map, na.rm = TRUE),
      sources = paste(sort(unique(source)), collapse = ";"),
      .groups = "drop"
    ) %>%
    mutate(n_iep_positions_in_map = if_else(
      is.infinite(n_iep_positions_in_map),
      NA_integer_,
      as.integer(n_iep_positions_in_map)
    )) %>%
    full_join(
      requested_panels %>% select(species, antigen, requested_panel),
      by = c("species", "antigen")
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
    stop("No eligible antigens after serotype-count and rare-haplotype filters.")
  }

  major_serotypes <- displayed_flows %>%
    semi_join(eligible_panels, by = c("species", "gene_id", "antigen")) %>%
    filter(!immune_escape_haplotype %in% c("__FILTERED_IEP__", "__OTHER__")) %>%
    distinct(species, gene_id, antigen, immune_escape_haplotype) %>%
    mutate(major_iep_serotype = immune_escape_haplotype)

  nearest_serotype_map <- build_nearest_serotype_map(
    flows %>% semi_join(eligible_panels, by = c("species", "gene_id", "antigen")),
    major_serotypes,
    args$max_iep_hamming_distance
  )

  mapping <- flows %>%
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
    left_join(alluvial_counts, by = c("species", "gene_id", "antigen")) %>%
    select(
      species, gene_id, antigen, all_snp_haplotype,
      immune_escape_haplotype, major_iep_serotype, nearest_iep_distance,
      n_samples, alluvial_haplotype_count, alluvial_serotype_count,
      rare_original_haplotype_percent
    )

  write_tsv(mapping, paste0(args$out_prefix, ".iep_hamming_serotype_mapping.tsv"))

  fitted <- mapping %>%
    group_split(species, gene_id, antigen) %>%
    map(fit_one_antigen, trees = args$trees, seed = args$seed,
        probability_threshold = args$probability_threshold,
        feature_lookup = feature_lookup)

  assignments <- bind_rows(map(fitted, "assignments"))
  performance <- bind_rows(map(fitted, "performance"))
  importance <- bind_rows(map(fitted, "importance"))

  write_tsv(assignments, paste0(args$out_prefix, ".rf_serotype_assignments.tsv"))
  write_tsv(performance, paste0(args$out_prefix, ".rf_serotype_performance.tsv"))
  write_tsv(importance, paste0(args$out_prefix, ".rf_iep_position_importance.tsv"))

  panel_audit <- panel_audit_base %>%
    left_join(raw_counts, by = c("species", "gene_id", "antigen")) %>%
    left_join(alluvial_counts, by = c("species", "gene_id", "antigen")) %>%
    left_join(
      performance %>%
        select(
          species, gene_id, antigen,
          rf_n_haplotypes = n_haplotypes,
          rf_n_serotypes = n_serotypes,
          rf_oob_accuracy,
          rf_status
        ),
      by = c("species", "gene_id", "antigen")
    ) %>%
    mutate(
      eligibility_status = case_when(
        is.na(raw_flow_rows) & is.na(alluvial_serotype_count) ~ "missing_from_raw_and_alluvial_flows",
        is.na(raw_flow_rows) ~ "missing_from_raw_flows",
        is.na(alluvial_serotype_count) ~ "missing_from_alluvial_flows",
        alluvial_serotype_count <= args$min_serotypes_exclusive ~ "excluded_fewer_than_two_serotypes",
        alluvial_serotype_count >= args$max_serotypes_exclusive ~ "excluded_too_many_serotypes",
        is.na(rare_original_haplotype_percent) ~ "excluded_missing_rare_original_haplotype_percent",
        rare_original_haplotype_percent > args$max_rare_original_haplotype_percent ~ "excluded_rare_original_haplotype_percent",
        is.na(rf_status) ~ "eligible_but_not_fit",
        TRUE ~ rf_status
      ),
      panel = paste0(if_else(species == "Pf7", "Pf", "Pv"), antigen),
      requested_for_plot = !is.na(requested_panel)
    ) %>%
    arrange(factor(panel, levels = args$plot_panels), species, antigen)

  write_tsv(panel_audit, paste0(args$out_prefix, ".rf_panel_eligibility.tsv"))

  requested_panels <- tibble(
    panel = args$plot_panels,
    species = if_else(startsWith(args$plot_panels, "Pf"), "Pf7", "Pv4"),
    antigen = sub("^(Pf|Pv)", "", args$plot_panels)
  ) %>%
    filter(grepl("^(Pf|Pv)", panel))
  missing_plot_panels <- requested_panels %>%
    anti_join(
      importance %>%
        mutate(panel = paste0(if_else(species == "Pf7", "Pf", "Pv"), antigen)) %>%
        distinct(panel),
      by = "panel"
    )
  if (nrow(missing_plot_panels)) {
    write_tsv(
      missing_plot_panels,
      paste0(args$out_prefix, ".rf_missing_requested_plot_panels.tsv")
    )
    warning(
      "Some requested plot panels have no RF importance rows: ",
      paste(missing_plot_panels$panel, collapse = ", ")
    )
  }

  plot_importance(
    importance,
    performance,
    args$out_prefix,
    args$top_n_plot,
    args$plot_width,
    args$plot_height,
    args$facet_ncol,
    args$plot_panels,
    args$plot_x_max,
    args$plot_facet_width,
    args$plot_facet_height,
    args$plot_pv_height_ratio,
    args$font_family
  )

  message("Wrote: ", args$out_prefix, ".iep_hamming_serotype_mapping.tsv")
  message("Wrote: ", args$out_prefix, ".rf_serotype_assignments.tsv")
  message("Wrote: ", args$out_prefix, ".rf_serotype_performance.tsv")
  message("Wrote: ", args$out_prefix, ".rf_iep_position_importance.tsv")
  message(
    "Wrote: ", args$out_prefix, ".rf_iep_position_importance_top",
    args$top_n_plot, ".pdf/.png/.svg"
  )
}

main()
