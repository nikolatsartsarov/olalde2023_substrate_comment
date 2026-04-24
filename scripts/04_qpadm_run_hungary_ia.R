#!/usr/bin/env Rscript
# 04: Test Hungary Iron Age candidates on the same strict 1-source qpAdm
# screen used for the NorthMacedonia_IA extension analysis.

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
T9_PATH <- file.path(RES, "T9_allsnps_sensitivity.tsv")

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

INLINE_HUNGARY_IA_SOURCES <- character()

requested_hungary_sources <- unique(c(
  INLINE_HUNGARY_IA_SOURCES,
  read_population_list(
    path = HUNGARY_IA_CANDIDATES,
    env_name = "OLALDE_HUNGARY_IA_SOURCES"
  )
))

if (!length(requested_hungary_sources)) {
  stop(
    "No Hungary IA candidates configured.\n",
    "Provide direct AADR population labels via either:\n",
    "  - data/hungary_ia_candidates.tsv with a 'source' column, or\n",
    "  - OLALDE_HUNGARY_IA_SOURCES as a comma/newline-separated list.\n",
    "You may also hard-code labels in INLINE_HUNGARY_IA_SOURCES if preferred."
  )
}

if (!file.exists(IND_FILE)) {
  stop(
    "Modified .ind file not found at:\n  ", IND_FILE, "\n",
    "Run scripts/02_build_ind_and_f2cache.R first."
  )
}

ind <- read.table(
  IND_FILE,
  col.names = c("id", "sex", "pop"),
  stringsAsFactors = FALSE
)
pop_sizes <- ind %>% count(pop, name = "sample_n")
available_pops <- pop_sizes$pop

present_hungary_sources <- requested_hungary_sources[requested_hungary_sources %in% available_pops]
missing_hungary_sources <- setdiff(requested_hungary_sources, present_hungary_sources)

source_panel_df <- tibble(source = requested_hungary_sources) %>%
  left_join(pop_sizes, by = c("source" = "pop")) %>%
  mutate(
    sample_n = coalesce(sample_n, 0L),
    present = sample_n > 0L
  )

write_tsv(source_panel_df, file.path(RES, "hungary_ia_candidate_panel.tsv"))

cat("[*] Requested Hungary IA candidates:\n")
print(requested_hungary_sources)
cat("[*] Candidate presence written to ", file.path(RES, "hungary_ia_candidate_panel.tsv"), "\n", sep = "")

if (length(missing_hungary_sources) > 0) {
  cat("\n!!! Missing Hungary IA candidate labels in ", IND_FILE, ":\n", sep = "")
  print(missing_hungary_sources)
}

if (!length(present_hungary_sources)) {
  stop(
    "None of the requested Hungary IA candidates are present in ", IND_FILE, ".\n",
    "Either supply direct AADR labels that already exist in the modified .ind,\n",
    "or extend scripts/02_build_ind_and_f2cache.R / data/cluster_definitions.tsv\n",
    "to relabel or pool the desired Hungary IA individuals."
  )
}

ALL_SOURCES <- c(PUBLISHED_SOURCES, present_hungary_sources)

cat("[*] Using strict AADR-label target panel with allsnps=TRUE.\n")

run_oneway_grid_live <- function(source_order, targets = TARGETS, prefix = "") {
  out_rows <- list()

  for (tgt in targets) {
    for (src in source_order) {
      cat(sprintf("%s%-14s ~ %s\n", prefix, tgt, src))
      fit <- tryCatch(
        qpadm(
          data = PREF_MOD,
          left = src,
          right = RIGHT,
          target = tgt,
          allsnps = TRUE,
          verbose = FALSE
        ),
        error = function(e) {
          cat("   !! ", conditionMessage(e), "\n")
          NULL
        }
      )
      if (is.null(fit)) next

      out_rows[[paste(tgt, src, sep = "|")]] <- tibble(
        target = tgt,
        source = src,
        pvalue = fit$rankdrop$p[1],
        pass = fit$rankdrop$p[1] > 0.05
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

persist_live_fills <- function(filled) {
  if (nrow(filled) == 0) return(invisible(NULL))

  if (file.exists(T9_PATH)) {
    t9_full <- read_tsv(T9_PATH, show_col_types = FALSE)
    template_cols <- names(t9_full)
  } else {
    t9_full <- tibble()
    template_cols <- c("model", "target", "n_sources", "pvalue")
  }

  cache_rows <- filled %>%
    transmute(
      model = paste0("1src:", source),
      target,
      n_sources = 1,
      pvalue = as.numeric(pvalue)
    )

  extra_cols <- setdiff(template_cols, names(cache_rows))
  for (col in extra_cols) {
    cache_rows[[col]] <- NA
  }
  cache_rows <- cache_rows[, unique(c(template_cols, names(cache_rows))), drop = FALSE]

  if (nrow(t9_full) > 0) {
    t9_full <- rows_upsert(t9_full, cache_rows, by = c("model", "target"))
  } else {
    t9_full <- cache_rows
  }

  write_tsv(t9_full, T9_PATH, na = "NA")
  cat("[*] Updated cache at ", T9_PATH, "\n", sep = "")
}

run_oneway_grid <- function(source_order) {
  if (file.exists(T9_PATH)) {
    cat("[*] Reusing existing allsnps results from ", T9_PATH, "\n", sep = "")
    t9 <- read_tsv(T9_PATH, show_col_types = FALSE) %>%
      filter(n_sources == 1, target %in% TARGETS) %>%
      transmute(
        target,
        source = sub("^1src:", "", model),
        pvalue = as.numeric(pvalue),
        pass = pvalue > 0.05
      )

    grid <- expand_grid(target = TARGETS, source = source_order) %>%
      left_join(t9, by = c("target", "source")) %>%
      mutate(
        target = factor(target, levels = TARGETS),
        source = factor(source, levels = source_order)
      ) %>%
      arrange(target, source)

    missing_pairs <- grid %>%
      filter(is.na(pvalue) | is.na(pass)) %>%
      distinct(target, source) %>%
      mutate(target = as.character(target), source = as.character(source))

    if (nrow(missing_pairs) > 0) {
      cat("[*] Filling ", nrow(missing_pairs), " missing cache cells live.\n", sep = "")
      filled <- run_oneway_grid_live(
        source_order = unique(missing_pairs$source),
        targets = unique(missing_pairs$target),
        prefix = "  [fill] "
      ) %>%
        mutate(target = as.character(target), source = as.character(source))

      persist_live_fills(filled)

      grid <- grid %>%
        mutate(target = as.character(target), source = as.character(source)) %>%
        rows_update(filled, by = c("target", "source")) %>%
        mutate(
          target = factor(target, levels = TARGETS),
          source = factor(source, levels = source_order)
        ) %>%
        arrange(target, source)
    }

    return(grid)
  }

  cat("[*] No T9 cache found; running live qpadm(allsnps=TRUE) from ", PREF_MOD, " ...\n", sep = "")
  run_oneway_grid_live(source_order)
}

write_outputs <- function(df, stem) {
  long_path <- file.path(RES, paste0(stem, "_long.tsv"))
  wide_path <- file.path(RES, paste0(stem, "_wide.tsv"))

  write_tsv(
    df %>% mutate(target = as.character(target), source = as.character(source)),
    long_path
  )

  wide <- df %>%
    mutate(
      pvalue_fmt = if_else(
        pvalue >= 0.001,
        sprintf("%.3f", pvalue),
        sprintf("%.2e", pvalue)
      ),
      cell = if_else(pass, paste0("PASS ", pvalue_fmt), paste0("fail ", pvalue_fmt))
    ) %>%
    select(target, source, cell) %>%
    pivot_wider(names_from = source, values_from = cell)

  write_tsv(wide, wide_path)
  cat("[*] Wrote", basename(long_path), "and", basename(wide_path), "\n")
}

cat("\n========== Published screen + Hungary IA candidates ==========\n")
combined_grid <- run_oneway_grid(ALL_SOURCES)
write_outputs(combined_grid, "hungary_ia_oneway")

summary_df <- combined_grid %>%
  mutate(source = as.character(source)) %>%
  group_by(target) %>%
  summarize(
    best_published_source = source[which.max(if_else(source %in% PUBLISHED_SOURCES, pvalue, -Inf))],
    best_published_p = max(if_else(source %in% PUBLISHED_SOURCES, pvalue, NA_real_), na.rm = TRUE),
    best_hungary_source = source[which.max(if_else(source %in% present_hungary_sources, pvalue, -Inf))],
    best_hungary_p = max(if_else(source %in% present_hungary_sources, pvalue, NA_real_), na.rm = TRUE),
    any_hungary_pass = any(if_else(source %in% present_hungary_sources, pass, FALSE), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(factor(target, levels = TARGETS))

write_tsv(summary_df, file.path(RES, "hungary_ia_oneway_summary.tsv"))

cat("\n=== Summary ===\n")
print(as.data.frame(summary_df), row.names = FALSE, right = FALSE)

cat("\n[*] DONE. Results in", RES, "\n")
