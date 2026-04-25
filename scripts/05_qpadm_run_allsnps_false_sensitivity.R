#!/usr/bin/env Rscript
# 05: Sensitivity check for the core 1-source screen under qpadm's default
# SNP-intersection mode, equivalent to allsnps=FALSE.

suppressMessages({
  library(admixtools)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
})

local({
  ns <- asNamespace("admixtools")
  shim_env <- new.env(parent = ns)
  shim_env$read_table2 <- readr::read_table
  for (nm in ls(ns)) {
    obj <- get(nm, envir = ns, inherits = FALSE)
    if (is.function(obj) && any(grepl("read_table2", deparse(body(obj)), fixed = TRUE))) {
      environment(obj) <- shim_env
      assignInNamespace(nm, obj, "admixtools")
    }
  }
})

local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- args[grep("^--file=", args)]
  d <- if (length(fa)) normalizePath(dirname(sub("^--file=", "", fa[1])))
       else normalizePath(getwd())
  sys.source(file.path(d, "_paths.R"), envir = globalenv())
})

RES <- RESULTS_DIR

RIGHT <- c(
  "OldAfrica", "Steppe_BA", "EHG", "Iron_Gates_HG", "Anatolia_N",
  "Iran_N", "Iberia_IA", "Greece_Minoan", "CroatiaMLBA_SloveniaIA",
  "Netherlands_MBA_IA", "Steppe_IA", "SoutheastTurkey_Byzantine",
  "Baltic_BA"
)

TARGETS <- c(
  "Croatian",
  "Serbian_Serb",
  "Romanian",
  "Bulgarian",
  "Albanian",
  "Greek",
  "Cypriot"
)

PUBLISHED_SOURCES <- c(
  "CroatiaSerbia_RomanLocal",
  "Aegean_BA_IA",
  "Albania_BA_IA",
  "Croatia_IA",
  "Bulgaria_IA",
  "Serbia_BA"
)

EXTENDED_SOURCES <- c(PUBLISHED_SOURCES, "NorthMacedonia_IA")

cat("[*] Running strict AADR-label target panel with allsnps=FALSE.\n")

run_oneway_grid <- function(source_order) {
  out_rows <- list()

  for (tgt in TARGETS) {
    for (src in source_order) {
      cat(sprintf("%-14s ~ %s\n", tgt, src))
      fit <- tryCatch(
        qpadm(
          data = PREF_MOD,
          left = src,
          right = RIGHT,
          target = tgt,
          verbose = FALSE
        ),
        error = function(e) {
          cat("   !! ", conditionMessage(e), "\n")
          NULL
        }
      )

      out_rows[[paste(tgt, src, sep = "|")]] <- tibble(
        target = tgt,
        source = src,
        pvalue = if (is.null(fit)) NA_real_ else fit$rankdrop$p[1],
        pass = if (is.null(fit)) NA else fit$rankdrop$p[1] > 0.05,
        error = is.null(fit)
      )
    }
  }

  bind_rows(out_rows) %>%
    mutate(
      target = factor(target, levels = TARGETS),
      source = factor(source, levels = source_order)
    ) %>%
    arrange(target, source)
}

write_outputs <- function(df) {
  long_path <- file.path(RES, "table8_oneway_with_northmac_allsnps_false_long.tsv")
  wide_path <- file.path(RES, "table8_oneway_with_northmac_allsnps_false_wide.tsv")

  write_tsv(
    df %>% mutate(target = as.character(target), source = as.character(source)),
    long_path
  )

  wide <- df %>%
    mutate(
      pvalue_fmt = case_when(
        is.na(pvalue) ~ "NA",
        pvalue >= 0.001 ~ sprintf("%.3f", pvalue),
        TRUE ~ sprintf("%.2e", pvalue)
      ),
      cell = case_when(
        is.na(pass) ~ "ERROR",
        pass ~ paste0("PASS ", pvalue_fmt),
        TRUE ~ paste0("fail ", pvalue_fmt)
      )
    ) %>%
    select(target, source, cell) %>%
    pivot_wider(names_from = source, values_from = cell)

  write_tsv(wide, wide_path)
  cat("[*] Wrote", basename(long_path), "and", basename(wide_path), "\n")
}

extended_grid <- run_oneway_grid(EXTENDED_SOURCES)
write_outputs(extended_grid)

summary_df <- extended_grid %>%
  mutate(source = as.character(source)) %>%
  group_by(target) %>%
  summarize(
    best_published_source = source[which.max(if_else(source %in% PUBLISHED_SOURCES, pvalue, -Inf))],
    best_published_p = max(if_else(source %in% PUBLISHED_SOURCES, pvalue, NA_real_), na.rm = TRUE),
    northmac_p = pvalue[source == "NorthMacedonia_IA"][1],
    northmac_pass = northmac_p > 0.05,
    .groups = "drop"
  ) %>%
  arrange(factor(target, levels = TARGETS))

write_tsv(summary_df, file.path(RES, "table8_oneway_summary_allsnps_false.tsv"))

if (file.exists(file.path(RES, "table8_oneway_summary.tsv"))) {
  true_summary <- read_tsv(file.path(RES, "table8_oneway_summary.tsv"), show_col_types = FALSE) %>%
    transmute(
      target,
      allsnps_true_best_published_p = best_published_p,
      allsnps_true_northmac_p = northmac_p,
      allsnps_true_northmac_pass = northmac_pass
    )

  comparison <- summary_df %>%
    transmute(
      target,
      allsnps_false_best_published_source = best_published_source,
      allsnps_false_best_published_p = best_published_p,
      allsnps_false_northmac_p = northmac_p,
      allsnps_false_northmac_pass = northmac_pass
    ) %>%
    left_join(true_summary, by = "target") %>%
    select(target, starts_with("allsnps_true"), starts_with("allsnps_false"))

  write_tsv(comparison, file.path(RES, "table8_oneway_allsnps_true_false_comparison.tsv"))
}

cat("[*] Wrote table8_oneway_summary_allsnps_false.tsv\n")
