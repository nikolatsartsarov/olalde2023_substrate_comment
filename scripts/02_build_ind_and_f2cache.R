#!/usr/bin/env Rscript
# 02: Build the modified genotype prefix used by the strict 7-label pipeline.

suppressMessages({
  library(admixtools); library(dplyr); library(readr)
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

# Central path config (env vars override defaults; see scripts/_paths.R).
local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- args[grep("^--file=", args)]
  d <- if (length(fa)) normalizePath(dirname(sub("^--file=", "", fa[1])))
       else normalizePath(getwd())
  sys.source(file.path(d, "_paths.R"), envir = globalenv())
})
require_aadr()

dir.create(dirname(PREF_MOD), recursive = TRUE, showWarnings = FALSE)

cat("[1] Reading .ind file ...\n")
ind <- read.table(paste0(PREF_ORIG, ".ind"),
                  col.names = c("id", "sex", "pop"),
                  stringsAsFactors = FALSE)
cat("  Total rows:", nrow(ind), "\n")

clu_path <- if (file.exists(CLUSTER_DEF)) CLUSTER_DEF else CLU_TSV
cat("  Using cluster definition file:", clu_path, "\n")
clu <- read_tsv(clu_path, show_col_types = FALSE)
clu <- as.data.frame(clu)
cat("  Olalde-cluster IDs:", nrow(clu), "\n")

# Strip suffixes from AADR ids: ".AG", ".SG", ".DG", ".HO", "_d", "_published",
# trailing dots, etc. The Olalde IDs in Table 3/6 are bare (e.g. "I5723", "AV1").
strip_id <- function(x) {
  x <- sub("\\.[A-Za-z0-9]+$", "", x)         # drop ".AG", ".SG", ...
  x <- sub("_d$", "", x)                       # drop "_d"
  x <- sub("_published$", "", x)               # drop "_published"
  x
}
ind$base_id <- strip_id(ind$id)

cat("[2] Applying Olalde cluster relabeling by ID match ...\n")
match_idx <- match(ind$base_id, clu$id)
to_relabel <- !is.na(match_idx)
ind$pop[to_relabel] <- clu$olalde_cluster[match_idx[to_relabel]]
cat("  Rows relabeled:", sum(to_relabel), "\n")
print(table(ind$pop[to_relabel]))

# CRITICAL FIX 1: drop _d duplicates if base ID is also present.
cat("[3a] Dropping _d duplicates whose base ID is present ...\n")
is_d_id <- grepl("_d\\.[A-Za-z0-9]+$", ind$id)
base_of_d <- sub("_d(\\.[A-Za-z0-9]+)$", "\\1", ind$id[is_d_id])
all_ids <- ind$id
drop_d <- is_d_id
drop_d[is_d_id] <- base_of_d %in% all_ids
cat(sprintf("  Marking %d _d duplicates as Ignore_*\n", sum(drop_d)))
ind$pop[drop_d] <- paste0("Ignore_", ind$pop[drop_d])

# CRITICAL FIX 2: 3 trailing IDs trigger malloc bug in extract_f2.
cat("[3b] Excluding 3 trailing IDs known to crash extract_f2 ...\n")
bad_ids <- c("I13519.AG", "I13536.AG", "I23495.TW")
bad_idx <- which(ind$id %in% bad_ids)
ind$pop[bad_idx] <- paste0("Ignore_", ind$pop[bad_idx])

# Pools we still need for the EXTENSION analyses (not Olalde clusters):
extension_pools <- list(
  Albania_BA_IA_pool = c("Albania_Cinamak_IA", "Albania_Cinamak_BA_IA")
)
cat("[4] Applying extension-source pooling ...\n")
for (newpop in names(extension_pools)) {
  members <- extension_pools[[newpop]]
  hits <- ind$pop %in% members
  ind$pop[hits] <- newpop
  cat(sprintf("  %-30s n=%d\n", newpop, sum(hits)))
}

# ----- Population sets -----
modern_targets <- c(
  "Bulgarian", "Romanian", "Albanian", "Croatian", "Greek",
  "Serbian_Serb", "Bosnian_Serb", "Cypriot", "Greek_WGA",
  "Gagauz", "Italian_North"
)

# Olalde's right-set (outgroups) for present-day Balkan analysis:
olalde_right_set <- c(
  "OldAfrica", "Steppe_BA", "EHG", "Iron_Gates_HG", "Anatolia_N",
  "Iran_N", "Iberia_IA", "Greece_Minoan", "CroatiaMLBA_SloveniaIA",
  "Netherlands_MBA_IA", "Steppe_IA", "SoutheastTurkey_Byzantine",
  "Baltic_BA"
)

# Olalde's source candidates for Table 8:
olalde_sources <- c(
  "CroatiaSerbia_RomanLocal", "CroatiaSerbia_RomanAnatolian",
  "CEE_EarlyMedieval", "WestAnatolia_Ottoman",
  "Aegean_BA_IA", "Albania_BA_IA", "Croatia_IA",
  "Bulgaria_IA", "Serbia_BA"
)

# Optional Hungary IA direct labels supplied via file/env:
hungary_ia_sources <- read_population_list(
  path = HUNGARY_IA_CANDIDATES,
  env_name = "OLALDE_HUNGARY_IA_SOURCES"
)
if (length(hungary_ia_sources) > 0) {
  cat("[4b] Hungary IA candidate labels requested:\n")
  print(hungary_ia_sources)
}

# Extension sources we want to test:
extension_sources <- c(
  "NorthMacedonia_IA",            # AADR direct label
  "Albania_BA_IA_pool"            # pool of Cinamak_IA + Cinamak_BA_IA
)
extension_sources <- unique(c(extension_sources, hungary_ia_sources))

needed_pops <- unique(c(
  modern_targets, olalde_right_set, olalde_sources, extension_sources
))
needed_pops <- intersect(needed_pops, unique(ind$pop))

cat("\n=== Population sample sizes (post-relabeling) ===\n")
size_tbl <- ind %>% filter(pop %in% needed_pops) %>% count(pop, sort = TRUE)
print(as.data.frame(size_tbl), row.names = FALSE)

missing_pops <- setdiff(
  unique(c(modern_targets, olalde_right_set, olalde_sources, extension_sources)),
  ind$pop)
if (length(missing_pops) > 0) {
  cat("\n!!! MISSING populations (will be skipped):\n")
  print(missing_pops)
}

cat("\n[5] Writing modified .ind ...\n")
write.table(ind[, c("id","sex","pop")], paste0(PREF_MOD, ".ind"),
            quote = FALSE, row.names = FALSE, col.names = FALSE,
            sep = "\t")
for (ext in c(".geno", ".snp")) {
  src <- paste0(PREF_ORIG, ext)
  dst <- paste0(PREF_MOD, ext)
  if (!file.exists(dst)) file.symlink(src, dst)
}

cat("\n[6] Modified genotype prefix ready for allsnps=TRUE runs.\n")
cat("  Prefix:", PREF_MOD, "\n")
cat("DONE.\n")
