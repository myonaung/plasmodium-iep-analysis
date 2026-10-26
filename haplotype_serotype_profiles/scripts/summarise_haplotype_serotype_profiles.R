#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

split_arg <- function(x) {
  x <- trimws(unlist(strsplit(if (is.null(x) || identical(x, "")) "" else x, ",")))
  x[nzchar(x)]
}

parse_args <- function(args) {
  out <- list(
    summary = "",
    flows = "",
    raw_flows = "",
    pf_site_summary = "",
    pv_site_summary = "",
    pf_snp_col = "pf7_ns_snp",
    pv_snp_col = "Pv4_snps",
    prob_sym_gt_asym_threshold = 0.80,
    prob_sym_gt_expected_threshold = 0,
    require_prob_sym_gt_expected_over_asym = "TRUE",
    dominant_allele_freq_max = 0.99,
    output = "haplotype_serotype_summary.tsv",
    left_prefix = "H",
    right_prefix = "H"
  )

  i <- 1L
  while (i <= length(args)) {
    if (!startsWith(args[[i]], "--") || i == length(args)) {
      stop("Bad arguments. Use --key value pairs.")
    }
    key <- gsub("-", "_", sub("^--", "", args[[i]]))
    out[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }

  if (!nzchar(out$summary)) stop("Provide --summary.")
  if (!nzchar(out$flows)) stop("Provide --flows.")
  if (!nzchar(out$output)) stop("Provide --output.")
  out$prob_sym_gt_asym_threshold <- as.numeric(out$prob_sym_gt_asym_threshold)
  out$prob_sym_gt_expected_threshold <- as.numeric(out$prob_sym_gt_expected_threshold)
  out$dominant_allele_freq_max <- as.numeric(out$dominant_allele_freq_max)
  out$require_prob_sym_gt_expected_over_asym <-
    tolower(out$require_prob_sym_gt_expected_over_asym) %in% c("true", "t", "1", "yes", "y")
  out
}

num <- function(x) suppressWarnings(as.numeric(x))

boolish <- function(x) {
  if (is.logical(x)) return(x)
  xn <- suppressWarnings(as.numeric(x))
  if (!all(is.na(xn))) return(xn == 1)
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

pick_col <- function(df, cands) {
  hit <- cands[cands %in% names(df)]
  if (length(hit)) hit[[1]] else NA_character_
}

first_existing <- function(df, cols) {
  cols <- cols[cols %in% names(df)]
  if (!length(cols)) return(rep(NA, nrow(df)))
  df[[cols[[1]]]]
}

first_non_missing <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x)) x[[1]] else NA_character_
}

as_out <- function(x) {
  x <- as.character(x)
  ifelse(is.na(x) | !nzchar(x), "XXX", x)
}

norm_antigen <- function(species, antigen) {
  antigen <- as.character(antigen)
  x <- toupper(gsub("[^A-Za-z0-9]+", "", antigen))
  ifelse(
    species == "Pf7" & x == "PFS25", "P25",
    ifelse(
      species == "Pf7" & x == "PFS230", "P230",
      ifelse(species == "Pf7" & x == "PFS4845", "P4845", x)
    )
  )
}

ordered_antigens <- function() {
  tribble(
    ~species, ~Species, ~`Life-cycle stage`, ~antigen, ~gene_id,
    "Pf7", "P. falciparum", "Pre-erythrocytic", "CELTOS",  "PF3D7_1216600",
    "Pf7", "P. falciparum", "Pre-erythrocytic", "CSP",     "PF3D7_0304600",
    "Pf7", "P. falciparum", "Pre-erythrocytic", "CTRP",    "PF3D7_0315200",
    "Pf7", "P. falciparum", "Pre-erythrocytic", "EXP1",    "PF3D7_1121600",
    "Pf7", "P. falciparum", "Pre-erythrocytic", "LSA3",    "PF3D7_0220000",
    "Pf7", "P. falciparum", "Pre-erythrocytic", "STARP",   "PF3D7_0702300",
    "Pf7", "P. falciparum", "Pre-erythrocytic", "TRAP",    "PF3D7_1335900",
    "Pf7", "P. falciparum", "Pre-erythrocytic", "TREP",    "PF3D7_1442600",
    "Pf7", "P. falciparum", "Blood-stage",      "AMA1",    "PF3D7_1133400",
    "Pf7", "P. falciparum", "Blood-stage",      "CLAG8",   "PF3D7_0831600",
    "Pf7", "P. falciparum", "Blood-stage",      "CYRPA",   "PF3D7_0423800",
    "Pf7", "P. falciparum", "Blood-stage",      "DBLMSP",  "PF3D7_1035700",
    "Pf7", "P. falciparum", "Blood-stage",      "EBA175",  "PF3D7_0731500",
    "Pf7", "P. falciparum", "Blood-stage",      "GLURP",   "PF3D7_1035300",
    "Pf7", "P. falciparum", "Blood-stage",      "MSP1",    "PF3D7_0930300",
    "Pf7", "P. falciparum", "Blood-stage",      "MSP2",    "PF3D7_0206800",
    "Pf7", "P. falciparum", "Blood-stage",      "MSP3",    "PF3D7_1035400",
    "Pf7", "P. falciparum", "Blood-stage",      "MSP4",    "PF3D7_0207000",
    "Pf7", "P. falciparum", "Blood-stage",      "MSP6",    "PF3D7_1035500",
    "Pf7", "P. falciparum", "Blood-stage",      "MSP9",    "PF3D7_1228600",
    "Pf7", "P. falciparum", "Blood-stage",      "RALP1",   "PF3D7_0722200",
    "Pf7", "P. falciparum", "Blood-stage",      "RESA",    "PF3D7_0102200",
    "Pf7", "P. falciparum", "Blood-stage",      "RH4",     "PF3D7_0424200",
    "Pf7", "P. falciparum", "Blood-stage",      "RH5",     "PF3D7_0424100",
    "Pf7", "P. falciparum", "Blood-stage",      "RIPR",    "PF3D7_0323400",
    "Pf7", "P. falciparum", "Blood-stage",      "SERA5",   "PF3D7_0207600",
    "Pf7", "P. falciparum", "Blood-stage",      "SERA8",   "PF3D7_0207300",
    "Pf7", "P. falciparum", "Blood-stage",      "TRAMP",   "PF3D7_1218000",
    "Pf7", "P. falciparum", "Gametocyte",       "CHT1",    "PF3D7_1252200",
    "Pf7", "P. falciparum", "Gametocyte",       "HAP2",    "PF3D7_1014200",
    "Pf7", "P. falciparum", "Gametocyte",       "P28",     "PF3D7_1030900",
    "Pf7", "P. falciparum", "Gametocyte",       "Pfs25",   "PF3D7_1031000",
    "Pf7", "P. falciparum", "Gametocyte",       "Pfs230",  "PF3D7_0209000",
    "Pf7", "P. falciparum", "Gametocyte",       "Pfs48/45","PF3D7_1346700",
    "Pv4", "P. vivax",      "Pre-erythrocytic", "CSP",     "PVP01_0835600",
    "Pv4", "P. vivax",      "Pre-erythrocytic", "TRAP",    "PVP01_1218700",
    "Pv4", "P. vivax",      "Blood-stage",      "AMA1",    "PVP01_0934200",
    "Pv4", "P. vivax",      "Blood-stage",      "CyRPA",   "PVP01_0532400",
    "Pv4", "P. vivax",      "Blood-stage",      "P41",     "PVP01_0304300",
    "Pv4", "P. vivax",      "Blood-stage",      "MSP1",    "PVP01_0728900",
    "Pv4", "P. vivax",      "Blood-stage",      "MSP8",    "PVP01_1032600",
    "Pv4", "P. vivax",      "Blood-stage",      "MSP9",    "PVP01_1446800",
    "Pv4", "P. vivax",      "Blood-stage",      "RBP1a",   "PVP01_0701200",
    "Pv4", "P. vivax",      "Blood-stage",      "RBP2a",   "PVP01_1402400",
    "Pv4", "P. vivax",      "Blood-stage",      "RBP2b",   "PVP01_0800700",
    "Pv4", "P. vivax",      "Blood-stage",      "RON2",    "PVP01_1255000",
    "Pv4", "P. vivax",      "Blood-stage",      "DBP",     "PVP01_0623800",
    "Pv4", "P. vivax",      "Blood-stage",      "TRAMP",   "PVP01_1436800",
    "Pv4", "P. vivax",      "Gametocyte",       "P28",     "PVP01_0616000"
  ) %>%
    mutate(
      antigen_key = norm_antigen(species, antigen),
      row_order = row_number()
    )
}

not_displayed_haplotype <- function(x) {
  x <- as.character(x)
  is.na(x) |
    x %in% c("__OTHER__", "__FILTERED_IEP__") |
    str_detect(x, regex("other|filtered", ignore_case = TRUE))
}

haplotype_diversity <- function(n) {
  n <- as.numeric(n)
  n <- n[is.finite(n) & n > 0]
  total <- sum(n)
  if (!is.finite(total) || total <= 1) return(NA_real_)
  p <- n / total
  total / (total - 1) * (1 - sum(p^2))
}

haplotype_stats <- function(flows) {
  required <- c("species", "antigen", "all_snp_haplotype", "n_samples")
  missing <- setdiff(required, names(flows))
  if (length(missing)) {
    stop("Missing required columns for haplotype diversity: ", paste(missing, collapse = ", "))
  }

  flows %>%
    filter(!not_displayed_haplotype(all_snp_haplotype)) %>%
    mutate(antigen_key = norm_antigen(species, antigen)) %>%
    group_by(species, antigen_key, all_snp_haplotype) %>%
    summarise(n_samples = sum(n_samples, na.rm = TRUE), .groups = "drop") %>%
    group_by(species, antigen_key) %>%
    summarise(
      `Haplotype number` = n_distinct(all_snp_haplotype),
      `Haplotype diversity (Hd)` = haplotype_diversity(n_samples),
      .groups = "drop"
    )
}

deduplicated_site_counts <- function(path, species, snp_col, args) {
  if (!nzchar(path) || !file.exists(path)) {
    return(tibble(
      species = character(),
      antigen_key = character(),
      `Non-synonymous polymorphic sites (n)` = integer(),
      `Candidate IEPs (n)` = integer()
    ))
  }

  df <- readr::read_csv(path, show_col_types = FALSE)
  antigen_col <- pick_col(df, c("antigen", "gene_name", "Gene_Name_annotated", "variant_gene_name"))
  gene_col <- pick_col(df, c("gene_id", "Gene_ID", "PlasmoDB_ID"))
  pos_col <- pick_col(df, c("position", "site", "variant_aa_pos", "aa_position", "pos"))

  if (is.na(antigen_col) || is.na(pos_col)) {
    warning("Could not deduplicate site counts for ", species, ": missing antigen or position column in ", path)
    return(tibble(
      species = character(),
      antigen_key = character(),
      `Non-synonymous polymorphic sites (n)` = integer(),
      `Candidate IEPs (n)` = integer()
    ))
  }

  if (!snp_col %in% names(df)) {
    warning("Could not find SNP column ", snp_col, " in ", path, "; SNP-site deduplication skipped for ", species)
    snp_flag <- rep(TRUE, nrow(df))
  } else {
    snp_flag <- boolish(df[[snp_col]])
  }

  if (!all(c("prob_sym_gt_asym", "prob_sym_gt_expected") %in% names(df))) {
    warning(
      "Could not calculate candidate IEP counts for ", species,
      ": missing posterior probability columns in ", path,
      ". Candidate IEP counts will be zero."
    )
    candidate_flag <- rep(FALSE, nrow(df))
  } else {
    prob_sa <- num(df$prob_sym_gt_asym)
    prob_se <- num(df$prob_sym_gt_expected)
    candidate_flag <- is.finite(prob_sa) &
      prob_sa >= args$prob_sym_gt_asym_threshold &
      is.finite(prob_se) &
      prob_se > args$prob_sym_gt_expected_threshold
    if (isTRUE(args$require_prob_sym_gt_expected_over_asym) &&
        "prob_asym_gt_expected" %in% names(df)) {
      candidate_flag <- candidate_flag & prob_se > num(df$prob_asym_gt_expected)
    }
  }

  if ("dominant_allele_freq" %in% names(df)) {
    candidate_flag <- candidate_flag &
      (is.na(num(df$dominant_allele_freq)) | num(df$dominant_allele_freq) <= args$dominant_allele_freq_max)
  }
  if ("global_rare_filtered" %in% names(df)) {
    candidate_flag <- candidate_flag & !boolish(df$global_rare_filtered)
  }

  site_df <- df %>%
    transmute(
      species = species,
      gene_id = if (!is.na(gene_col)) as.character(.data[[gene_col]]) else NA_character_,
      antigen = as.character(.data[[antigen_col]]),
      antigen_key = norm_antigen(species, antigen),
      position = num(.data[[pos_col]]),
      snp_flag = snp_flag,
      candidate_flag = candidate_flag
    ) %>%
    filter(!is.na(position), nzchar(antigen_key))

  snp_counts <- site_df %>%
    filter(snp_flag) %>%
    distinct(species, antigen_key, position) %>%
    count(species, antigen_key, name = "Non-synonymous polymorphic sites (n)")

  iep_counts <- site_df %>%
    filter(candidate_flag) %>%
    distinct(species, antigen_key, position) %>%
    count(species, antigen_key, name = "Candidate IEPs (n)")

  full_join(snp_counts, iep_counts, by = c("species", "antigen_key")) %>%
    mutate(
      `Non-synonymous polymorphic sites (n)` = coalesce(`Non-synonymous polymorphic sites (n)`, 0L),
      `Candidate IEPs (n)` = coalesce(`Candidate IEPs (n)`, 0L)
    )
}

rank_haplotypes <- function(flows, hap_col, ref_col, prefix) {
  required <- c("species", "gene_id", "antigen", hap_col, "n_samples")
  missing <- setdiff(required, names(flows))
  if (length(missing)) {
    stop("Missing required columns in flows file: ", paste(missing, collapse = ", "))
  }

  if (!ref_col %in% names(flows)) {
    flows[[ref_col]] <- NA_character_
  }

  flows %>%
    filter(!not_displayed_haplotype(.data[[hap_col]])) %>%
    mutate(antigen_key = norm_antigen(species, antigen)) %>%
    group_by(species, antigen_key, gene_id, antigen, hap = .data[[hap_col]]) %>%
    summarise(
      n_samples = sum(n_samples, na.rm = TRUE),
      reference_haplotype = first_non_missing(.data[[ref_col]]),
      .groups = "drop"
    ) %>%
    group_by(species, antigen_key) %>%
    arrange(desc(n_samples), hap, .by_group = TRUE) %>%
    mutate(
      haplotype_number = row_number(),
      haplotype_label = paste0(prefix, haplotype_number),
      reference_like = !is.na(reference_haplotype) &
        nzchar(reference_haplotype) &
        hap == reference_haplotype
    ) %>%
    ungroup()
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  summary_df <- read_tsv(args$summary, show_col_types = FALSE)
  flows_df <- read_tsv(args$flows, show_col_types = FALSE)
  raw_flows_df <- if (nzchar(args$raw_flows)) {
    read_tsv(args$raw_flows, show_col_types = FALSE)
  } else {
    flows_df
  }

  antigen_order <- ordered_antigens()
  hap_stats_df <- haplotype_stats(raw_flows_df)
  site_counts_df <- bind_rows(
    deduplicated_site_counts(args$pf_site_summary, "Pf7", args$pf_snp_col, args),
    deduplicated_site_counts(args$pv_site_summary, "Pv4", args$pv_snp_col, args)
  ) %>%
    rename(
      `Non-synonymous polymorphic sites (n)_dedup` = `Non-synonymous polymorphic sites (n)`,
      `Candidate IEPs (n)_dedup` = `Candidate IEPs (n)`
  )

  right_rank_df <- rank_haplotypes(
    flows_df,
    hap_col = "immune_escape_haplotype",
    ref_col = "immune_escape_reference_haplotype",
    prefix = args$right_prefix
  )

  left_rank_df <- rank_haplotypes(
    flows_df,
    hap_col = "all_snp_haplotype",
    ref_col = "all_snp_reference_haplotype",
    prefix = args$left_prefix
  )

  serotype_summary <- right_rank_df %>%
    group_by(species, antigen_key) %>%
    summarise(
      putative_serotypes_n = n_distinct(hap),
      total_serotype_samples = sum(n_samples, na.rm = TRUE),
      reference_like_putative_serotype_found = any(reference_like),
      reference_like_putative_serotype_samples = sum(n_samples[reference_like], na.rm = TRUE),
      reference_like_putative_serotype_percent = if_else(
        reference_like_putative_serotype_found & total_serotype_samples > 0,
        100 * reference_like_putative_serotype_samples / total_serotype_samples,
        NA_real_
      ),
      reference_like_putative_serotype = {
        x <- haplotype_label[reference_like]
        if (length(x)) paste0(x[[1]], "*") else "XXX"
      },
      reference_like_putative_serotype_number = {
        x <- haplotype_number[reference_like]
        if (length(x)) as.character(x[[1]]) else "XXX"
      },
      reference_like_iep_profile = {
        x <- hap[reference_like]
        if (length(x)) x[[1]] else "XXX"
      },
      .groups = "drop"
    )

  all_snp_reference_summary <- left_rank_df %>%
    group_by(species, antigen_key) %>%
    summarise(
      reference_all_snp_haplotype = {
        x <- haplotype_label[reference_like]
        if (length(x)) paste0(x[[1]], "*") else "XXX"
      },
      reference_all_snp_haplotype_number = {
        x <- haplotype_number[reference_like]
        if (length(x)) as.character(x[[1]]) else "XXX"
      },
      reference_all_snp_profile = {
        x <- hap[reference_like]
        if (length(x)) x[[1]] else "XXX"
      },
      .groups = "drop"
    )

  summary_cols <- c(
    "species", "gene_id", "antigen",
    "n_all_snp_haplotypes_raw",
    "n_all_snp_haplotypes",
    "n_immune_escape_haplotypes",
    "n_snp_sites_for_density",
    "n_significant_iep_sites",
    "n_potential_iep_sites",
    "n_global_rare_filtered_iep_sites"
  )
  keep_cols <- intersect(summary_cols, names(summary_df))

  summary_prepped <- summary_df %>%
    select(all_of(keep_cols)) %>%
    mutate(antigen_key = norm_antigen(species, antigen)) %>%
    distinct(species, antigen_key, .keep_all = TRUE)

  out <- antigen_order %>%
    left_join(
      summary_prepped,
      by = c("species", "antigen_key"),
      suffix = c("", "_summary")
    ) %>%
    left_join(hap_stats_df, by = c("species", "antigen_key")) %>%
    left_join(site_counts_df, by = c("species", "antigen_key")) %>%
    left_join(serotype_summary, by = c("species", "antigen_key")) %>%
    left_join(all_snp_reference_summary, by = c("species", "antigen_key")) %>%
    mutate(
      gene_id = coalesce(gene_id_summary, gene_id),
      haplotype_number_all_snp_raw = first_existing(
        pick(everything()),
        c("n_all_snp_haplotypes_raw", "n_all_snp_haplotypes")
      ),
      `Non-synonymous polymorphic sites (n)` = first_existing(
        pick(everything()),
        c("n_snp_sites_for_density")
      ),
      `Non-synonymous polymorphic sites (n)` = coalesce(
        `Non-synonymous polymorphic sites (n)_dedup`,
        `Non-synonymous polymorphic sites (n)`
      ),
      `Candidate IEPs (n)` = first_existing(
        pick(everything()),
        c("n_significant_iep_sites")
      ),
      `Candidate IEPs (n)` = coalesce(
        `Candidate IEPs (n)_dedup`,
        `Candidate IEPs (n)`
      ),
      `Candidate IEPs (n)` = coalesce(`Candidate IEPs (n)`, 0),
      `Haplotype number` = coalesce(`Haplotype number`, haplotype_number_all_snp_raw),
      `Putative serotypes (n)` = coalesce(
        putative_serotypes_n,
        first_existing(pick(everything()), c("n_immune_escape_haplotypes"))
      ),
      `Putative serotypes (n)` = if_else(
        suppressWarnings(as.numeric(`Candidate IEPs (n)`)) == 0 &
          (is.na(`Putative serotypes (n)`) | !is.finite(suppressWarnings(as.numeric(`Putative serotypes (n)`)))),
        1,
        `Putative serotypes (n)`
      ),
      `Reference-like putative serotype number` = coalesce(reference_like_putative_serotype_number, "XXX"),
      `Reference-like putative serotype label` = coalesce(reference_like_putative_serotype, "XXX"),
      `Reference-like putative serotype number` = if_else(
        suppressWarnings(as.numeric(`Candidate IEPs (n)`)) == 0,
        "1",
        as.character(`Reference-like putative serotype number`)
      ),
      `Reference-like putative serotype label` = if_else(
        suppressWarnings(as.numeric(`Candidate IEPs (n)`)) == 0,
        paste0(args$right_prefix, "1*"),
        as.character(`Reference-like putative serotype label`)
      ),
      `Reference-like putative serotype (%)` = round(reference_like_putative_serotype_percent, 2),
      `Reference-like putative serotype (%)` = if_else(
        suppressWarnings(as.numeric(`Candidate IEPs (n)`)) == 0 &
          (is.na(`Reference-like putative serotype (%)`) |
             !is.finite(`Reference-like putative serotype (%)`)),
        100,
        `Reference-like putative serotype (%)`
      ),
      `Reference all-SNP haplotype number` = coalesce(reference_all_snp_haplotype_number, "XXX"),
      `Reference all-SNP haplotype label` = coalesce(reference_all_snp_haplotype, "XXX")
    ) %>%
    select(
      row_order,
      Species,
      `Life-cycle stage`,
      gene_id,
      antigen,
      `Non-synonymous polymorphic sites (n)`,
      `Candidate IEPs (n)`,
      `Haplotype number`,
      `Haplotype diversity (Hd)`,
      `Reference all-SNP haplotype number`,
      `Reference all-SNP haplotype label`,
      `Putative serotypes (n)`,
      `Reference-like putative serotype number`,
      `Reference-like putative serotype label`,
      `Reference-like putative serotype (%)`,
      reference_all_snp_profile,
      reference_like_iep_profile,
      any_of(c(
        "n_potential_iep_sites",
        "n_global_rare_filtered_iep_sites",
        "n_all_snp_haplotypes",
        "n_immune_escape_haplotypes"
      ))
    ) %>%
    arrange(row_order) %>%
    mutate(across(-row_order, as_out)) %>%
    select(-any_of("row_order"))

  dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
  write_tsv(out, args$output)
  message("Wrote: ", args$output)
  message("Rows: ", nrow(out))
  out %>%
    mutate(
      candidate_iep_n = num(`Candidate IEPs (n)`),
      nonsyn_site_n = num(`Non-synonymous polymorphic sites (n)`)
    ) %>%
    group_by(Species) %>%
    summarise(
      antigen_rows = n(),
      candidate_iep_sites_sum = sum(candidate_iep_n, na.rm = TRUE),
      nonsyn_sites_sum = sum(nonsyn_site_n, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(message_text = paste0(
      Species, ": antigen rows=", antigen_rows,
      "; summed candidate IEP sites=", candidate_iep_sites_sum,
      "; summed non-synonymous polymorphic sites=", nonsyn_sites_sum
    )) %>%
    pull(message_text) %>%
    lapply(message)
}

main()
