#!/usr/bin/env Rscript
# 11: Generate manuscript tables for the strict 1-source AADR-label analysis.

suppressMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

`%+%` <- function(a, b) paste0(a, b)
ROW_END <- "\\\\"

local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- args[grep("^--file=", args)]
  d <- if (length(fa)) normalizePath(dirname(sub("^--file=", "", fa[1])))
       else normalizePath(getwd())
  sys.source(file.path(d, "_paths.R"), envir = globalenv())
})

RES <- RESULTS_DIR
TBL <- TABLES_DIR

fmt_p <- function(p) {
  ifelse(is.na(p), "--",
         ifelse(p >= 0.001, sprintf("%.3f", p), sprintf("%.1e", p)))
}

fmt_p_digits <- function(p, digits = 3) {
  ifelse(is.na(p), "--", sprintf(paste0("%.", digits, "f"), p))
}

escape_tex <- function(x) gsub("_", "\\_", x, fixed = TRUE)

oneway <- read_tsv(
  file.path(RES, "table8_oneway_with_northmac_long.tsv"),
  show_col_types = FALSE
)
summary_df <- read_tsv(
  file.path(RES, "table8_oneway_summary.tsv"),
  show_col_types = FALSE
)
panel_df <- read_tsv(
  file.path(RES, "table8_target_panel.tsv"),
  show_col_types = FALSE
)

source_order <- c(
  "CroatiaSerbia_RomanLocal",
  "Aegean_BA_IA",
  "Albania_BA_IA",
  "Croatia_IA",
  "Bulgaria_IA",
  "Serbia_BA",
  "NorthMacedonia_IA"
)
target_order <- c(
  "Croatian",
  "Serbian_Serb",
  "Romanian",
  "Bulgarian",
  "Albanian",
  "Greek",
  "Cypriot"
)

summary_df <- summary_df %>%
  mutate(
    target = factor(target, levels = target_order),
    verdict = if_else(northmac_pass, "passes", "fails")
  ) %>%
  arrange(target)

write_tsv(summary_df, file.path(TBL, "table1_main.tsv"))

best_min_row <- summary_df %>%
  slice_min(best_published_p, n = 1, with_ties = FALSE)
best_max_row <- summary_df %>%
  slice_max(best_published_p, n = 1, with_ties = FALSE)
northmac_min_row <- summary_df %>%
  slice_min(northmac_p, n = 1, with_ties = FALSE)
northmac_max_row <- summary_df %>%
  slice_max(northmac_p, n = 1, with_ties = FALSE)

scalars_tex <- c(
  sprintf("\\newcommand{\\BestPublishedMinP}{%s}", fmt_p_digits(best_min_row$best_published_p, 4)),
  sprintf("\\newcommand{\\BestPublishedMinTarget}{%s}", escape_tex(as.character(best_min_row$target))),
  sprintf("\\newcommand{\\BestPublishedMaxP}{%s}", fmt_p_digits(best_max_row$best_published_p, 4)),
  sprintf("\\newcommand{\\BestPublishedMaxTarget}{%s}", escape_tex(as.character(best_max_row$target))),
  sprintf("\\newcommand{\\NorthmacMinP}{%s}", fmt_p_digits(northmac_min_row$northmac_p, 3)),
  sprintf("\\newcommand{\\NorthmacMinTarget}{%s}", escape_tex(as.character(northmac_min_row$target))),
  sprintf("\\newcommand{\\NorthmacMaxP}{%s}", fmt_p_digits(northmac_max_row$northmac_p, 3)),
  sprintf("\\newcommand{\\NorthmacMaxTarget}{%s}", escape_tex(as.character(northmac_max_row$target)))
)
writeLines(scalars_tex, file.path(TBL, "scalars.tex"))
cat("[*] Wrote scalars.tex\n")

tex_main <- c(
  "\\begin{table}[H]",
  "\\centering",
  "\\caption{Strict AADR-label 1-way screen versus the added \\texttt{NorthMacedonia\\_IA} column.}",
  "\\label{tab:main}",
  "\\small",
  "\\begin{tabular}{l l r r l}",
  "\\toprule",
  paste0("Target & Best published source & Best published $p$ & \\texttt{NorthMacedonia\\_IA} $p$ & Verdict ", ROW_END),
  "\\midrule"
)

for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  tex_main <- c(
    tex_main,
    sprintf("%s & %s & %s & %s & %s %s",
            escape_tex(as.character(r$target)),
            escape_tex(r$best_published_source),
            fmt_p(r$best_published_p),
            fmt_p(r$northmac_p),
            r$verdict,
            ROW_END)
  )
}

tex_main <- c(tex_main, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex_main, file.path(TBL, "table1_main.tex"))
cat("[*] Wrote table1_main.tex / .tsv\n")

oneway_full <- oneway %>%
  mutate(
    target = factor(target, levels = target_order),
    source = factor(source, levels = source_order),
    cell = if_else(pass, paste0("\\textbf{", fmt_p(pvalue), "}"), fmt_p(pvalue))
  ) %>%
  arrange(target, source)

oneway_wide <- oneway_full %>%
  select(target, source, cell) %>%
  pivot_wider(names_from = source, values_from = cell)

write_tsv(oneway_wide, file.path(TBL, "tableS1_oneway_full.tsv"))

tex_s1 <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{\\textbf{Full 1-way qpAdm matrix for the strict AADR-label target panel plus \\texttt{NorthMacedonia\\_IA}.} Rows are the seven exact AADR-HO target labels retained in this comment. Columns are Olalde's six published 1-way substrates plus the added \\texttt{NorthMacedonia\\_IA} column. Entries are qpAdm $p$-values; bold indicates $p > 0.05$.}",
  "\\label{tab:oneway-full}",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{3pt}",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{l " %+% paste(rep("l", length(source_order)), collapse = " ") %+% "}",
  "\\toprule",
  paste0("Target & ", paste(escape_tex(source_order), collapse = " & "), " ", ROW_END),
  "\\midrule"
)

for (tg in target_order) {
  row <- oneway_wide %>% filter(target == tg)
  if (nrow(row) == 0) next
  cells <- vapply(source_order, function(s) {
    val <- row[[s]]
    if (is.null(val) || is.na(val)) "--" else val
  }, character(1))
  tex_s1 <- c(tex_s1, paste0(escape_tex(tg), " & ", paste(cells, collapse = " & "), " ", ROW_END))
}

tex_s1 <- c(tex_s1, "\\bottomrule", "\\end{tabular}", "}", "\\end{table}")
writeLines(tex_s1, file.path(TBL, "tableS1_oneway_full.tex"))
cat("[*] Wrote tableS1_oneway_full.tex / .tsv\n")

write_tsv(panel_df, file.path(TBL, "tableS2_target_panel.tsv"))

tex_s2 <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Exact public AADR Human Origins target labels retained in the strict public-label panel used in this comment.}",
  "\\label{tab:target-panel}",
  "\\small",
  "\\begin{tabular}{l}",
  "\\toprule",
  paste0("Target label ", ROW_END),
  "\\midrule"
)

for (i in seq_len(nrow(panel_df))) {
  r <- panel_df[i, ]
  tex_s2 <- c(tex_s2, sprintf("%s %s", escape_tex(r$target), ROW_END))
}

tex_s2 <- c(tex_s2, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex_s2, file.path(TBL, "tableS2_target_panel.tex"))
cat("[*] Wrote tableS2_target_panel.tex / .tsv\n")

cat("\n[*] DONE. All tables written to ", TBL, "\n")
