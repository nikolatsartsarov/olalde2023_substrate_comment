#!/usr/bin/env Rscript
# 10: Generate publication figures for the strict 1-source AADR-label analysis.

suppressMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- args[grep("^--file=", args)]
  d <- if (length(fa)) normalizePath(dirname(sub("^--file=", "", fa[1])))
       else normalizePath(getwd())
  sys.source(file.path(d, "_paths.R"), envir = globalenv())
})

RES <- RESULTS_DIR
FIG <- FIGS_DIR

oneway <- read_tsv(
  file.path(RES, "table8_oneway_with_northmac_long.tsv"),
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
source_labels <- c(
  "RomanLocal",
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

oneway_plot <- oneway %>%
  mutate(
    source = factor(source, levels = source_order, labels = source_labels),
    target = factor(target, levels = rev(target_order)),
    verdict = if_else(pass, "pass", "fail"),
    p_label = if_else(pvalue >= 0.001, sprintf("%.3f", pvalue), sprintf("%.1e", pvalue))
  )

p1 <- ggplot(oneway_plot, aes(x = source, y = target, fill = verdict)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = p_label), size = 2.4) +
  scale_fill_manual(values = c("fail" = "#d73027", "pass" = "#1a9850"), name = NULL) +
  labs(
    x = NULL,
    y = NULL,
    title = "Strict AADR-label 1-way qpAdm screen plus NorthMacedonia_IA",
    subtitle = "Seven exact AADR-HO target labels from the public panel; the added NorthMacedonia_IA column passes all seven."
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

ggsave(file.path(FIG, "fig1_table8_oneway.pdf"), p1,
       width = 8.6, height = 4.4, units = "in")
ggsave(file.path(FIG, "fig1_table8_oneway.png"), p1,
       width = 8.6, height = 4.4, units = "in", dpi = 300)
cat("[*] Wrote fig1_table8_oneway.{pdf,png}\n")
cat("\n[*] DONE.\n")
