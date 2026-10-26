#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

usage <- function() {
  cat(
    paste(
      "Usage:",
      "Rscript scripts/prepare_pf7_metadata.R",
      "--pf7-metadata Pf7_samples.txt",
      "--out pf7_sample_metadata.tsv",
      "[--sample-col Sample]",
      "[--study-col Study]",
      "[--country-col Country]",
      "[--admin1-col 'Admin level 1']",
      "[--region-col Population]",
      "[--year-col Year]",
      "[--qc-col 'QC pass']",
      "[--qc-keep TRUE]",
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
    pf7_metadata = NULL,
    out = NULL,
    sample_col = "Sample",
    study_col = "Study",
    country_col = "Country",
    admin1_col = "Admin level 1",
    region_col = "Population",
    year_col = "Year",
    qc_col = "QC pass",
    qc_keep = "TRUE"
  )

  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (!startsWith(args[[i]], "--") || i == length(args)) {
      stop("Bad arguments.")
    }
    value <- args[[i + 1L]]
    i <- i + 2L
    key2 <- gsub("-", "_", key)
    if (!key2 %in% names(out)) {
      stop("Unknown option: --", key)
    }
    out[[key2]] <- value
  }

  if (is.null(out$pf7_metadata) || is.null(out$out)) {
    usage()
  }

  out
}

main <- function() {
  args <- parse_args_simple(commandArgs(trailingOnly = TRUE))

  x <- read_tsv(args$pf7_metadata, show_col_types = FALSE)
  required <- c(args$sample_col, args$study_col, args$country_col, args$admin1_col, args$region_col, args$year_col)
  missing_cols <- setdiff(required, names(x))
  if (length(missing_cols) > 0) {
    stop("PF7 metadata file missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (args$qc_col %in% names(x)) {
    keep_value <- toupper(args$qc_keep)
    x <- x %>%
      mutate(.qc_value = toupper(as.character(.data[[args$qc_col]]))) %>%
      filter(.qc_value == keep_value)
  }

  out <- x %>%
    transmute(
      sample_id = as.character(.data[[args$sample_col]]),
      study = as.character(.data[[args$study_col]]),
      country = as.character(.data[[args$country_col]]),
      admin_level_1 = as.character(.data[[args$admin1_col]]),
      region = as.character(.data[[args$region_col]]),
      sample_year = suppressWarnings(as.integer(.data[[args$year_col]])),
      site = str_replace_all(
        paste(
          coalesce(as.character(.data[[args$study_col]]), "NA"),
          coalesce(as.character(.data[[args$country_col]]), "NA"),
          coalesce(as.character(.data[[args$admin1_col]]), "NA"),
          sep = "|"
        ),
        "\\s+",
        "_"
      )
    ) %>%
    distinct(sample_id, .keep_all = TRUE) %>%
    arrange(sample_id)

  write_tsv(out, args$out)
  message("Wrote PF7-ready metadata to: ", args$out)
}

main()
