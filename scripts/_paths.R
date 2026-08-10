## scripts/_paths.R
##
## Central path / configuration file for the Olalde 2023 substrate-comment
## pipeline. Every numbered script sources this file; users override
## defaults via environment variables. No script in this repo should
## contain hard-coded absolute paths besides this one.
##
## Required environment variables (set before running any script):
##   OLALDE_AADR_HO_PREFIX  Path prefix to AADR v66.0 Human Origins genotype
##                          files, WITHOUT the .geno/.snp/.ind extension.
##                          Example: /data/aadr/v66.HO.aadr.PUB
##                          (the three files .geno, .snp, .ind must all
##                          exist with this prefix).
##
## Optional environment variables:
##   OLALDE_WORK_DIR        Working directory for intermediate files (the
##                          modified .ind, the f2 cache, and Olalde
##                          ID->cluster TSVs). Large; can live anywhere.
##                          Default: <REPO_DIR>/work/
##   OLALDE_HUNGARY_IA_CANDIDATES
##                          Optional path to a TSV or plain-text file listing
##                          Hungary Iron Age candidate population labels.
##                          Default:
##                          <REPO_DIR>/data/hungary_ia_candidates.tsv
##   OLALDE_HUNGARY_IA_SOURCES
##                          Optional comma/newline/semicolon-separated list of
##                          Hungary Iron Age population labels. These labels
##                          are appended to any file-based candidate list.
##   OLALDE_S2_XLSX         Path to Olalde et al. 2023 Supplementary
##                          Data S2 (NIHMS1944038-supplement-Data_S2.xlsx).
##                          Required only by scripts 01 and 14.
##   OLALDE_REPO_DIR        Override repository root. Default: parent of
##                          the directory containing this _paths.R file.
##
## All other paths in the pipeline are derived from these.

# ----------------------------------------------------------------------
# Helper: read env var, fall back to default if unset/empty
# ----------------------------------------------------------------------
.getenv_or <- function(name, default = NULL) {
  v <- Sys.getenv(name, unset = "")
  if (nzchar(v)) v else default
}

.split_list_env <- function(value) {
  if (!nzchar(value)) return(character())
  parts <- unlist(strsplit(value, "[,\n;]", perl = TRUE), use.names = FALSE)
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

read_population_list <- function(path = NULL, env_name = NULL, default = character()) {
  vals <- default

  if (!is.null(env_name) && nzchar(env_name)) {
    vals <- c(vals, .split_list_env(Sys.getenv(env_name, unset = "")))
  }

  if (!is.null(path) && nzchar(path) && file.exists(path)) {
    file_vals <- tryCatch({
      tab <- utils::read.delim(
        path,
        sep = "\t",
        header = TRUE,
        stringsAsFactors = FALSE,
        comment.char = ""
      )
      if ("source" %in% names(tab)) {
        tab$source
      } else if (ncol(tab) >= 1L) {
        tab[[1L]]
      } else {
        character()
      }
    }, error = function(e) {
      lines <- readLines(path, warn = FALSE)
      lines <- trimws(lines)
      lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
      if (length(lines) && identical(tolower(lines[[1L]]), "source")) {
        lines[-1L]
      } else {
        lines
      }
    })
    vals <- c(vals, file_vals)
  }

  vals <- trimws(as.character(vals))
  vals <- vals[nzchar(vals)]
  unique(vals)
}

# ----------------------------------------------------------------------
# Locate this script's own directory (works under Rscript and interactive)
# ----------------------------------------------------------------------
.this_file <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  fa <- cmd_args[grep("^--file=", cmd_args)]
  if (length(fa) > 0L) {
    return(normalizePath(sub("^--file=", "", fa[1L])))
  }
  ofile <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(normalizePath(ofile))
  NA_character_
}

.scripts_dir <- tryCatch({
  tf <- .this_file()
  if (!is.na(tf)) normalizePath(dirname(tf)) else normalizePath(getwd())
}, error = function(e) normalizePath(getwd()))

# ----------------------------------------------------------------------
# REPO_DIR  = root of the olalde2023_substrate_comment repository
# ----------------------------------------------------------------------
REPO_DIR <- .getenv_or(
  "OLALDE_REPO_DIR",
  normalizePath(file.path(.scripts_dir, ".."), mustWork = FALSE)
)

# ----------------------------------------------------------------------
# WORK_DIR  = intermediate files (modified .ind, f2 cache, ID clusters)
# ----------------------------------------------------------------------
WORK_DIR <- .getenv_or(
  "OLALDE_WORK_DIR",
  file.path(REPO_DIR, "work")
)
dir.create(file.path(WORK_DIR, "intermediate"),
           recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------
# AADR genotype prefix (required only by script 02; downstream scripts
# read the modified .ind / f2 cache produced by script 02 and do not
# touch the original AADR files).
# ----------------------------------------------------------------------
PREF_ORIG <- .getenv_or("OLALDE_AADR_HO_PREFIX", "")
require_aadr <- function() {
  if (!nzchar(PREF_ORIG)) {
    stop(
      "OLALDE_AADR_HO_PREFIX is not set.\n",
      "Provide the path PREFIX to AADR v66.0 Human Origins genotype files\n",
      "(without .geno/.snp/.ind extension), e.g.\n",
      "  export OLALDE_AADR_HO_PREFIX=/data/aadr/v66.HO.aadr.PUB\n",
      "before running scripts that read the raw AADR files (script 02).\n",
      "See scripts/_paths.R for details."
    )
  }
  invisible(TRUE)
}

# ----------------------------------------------------------------------
# Olalde Supplementary Data S2 (.xlsx) -- required for scripts 01 and 14
# ----------------------------------------------------------------------
OLALDE_S2_XLSX <- .getenv_or("OLALDE_S2_XLSX", "")

# ----------------------------------------------------------------------
# Derived paths used throughout the pipeline
# ----------------------------------------------------------------------
HERE        <- WORK_DIR                                          # back-compat
PREF_MOD    <- file.path(WORK_DIR, "intermediate", "v66.HO.olalde_repro")
F2_DIR      <- file.path(WORK_DIR, "intermediate", "f2_blocks")
CLU_TSV     <- file.path(WORK_DIR, "intermediate", "olalde_id_clusters.tsv")
IND_FILE    <- paste0(PREF_MOD, ".ind")

# Repo-relative outputs (results, figures, tables, data)
RESULTS_DIR <- file.path(REPO_DIR, "results")
FIGS_DIR    <- file.path(REPO_DIR, "manuscript", "figs")
TABLES_DIR  <- file.path(REPO_DIR, "manuscript", "tables")
DATA_DIR    <- file.path(REPO_DIR, "data")
PAPER_DIR   <- REPO_DIR
CLUSTER_DEF <- file.path(DATA_DIR, "cluster_definitions.tsv")
HUNGARY_IA_CANDIDATES <- .getenv_or(
  "OLALDE_HUNGARY_IA_CANDIDATES",
  file.path(DATA_DIR, "hungary_ia_candidates.tsv")
)
OUT_SUMMARY <- file.path(DATA_DIR, "id_reconciliation.tsv")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGS_DIR,    recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR,    recursive = TRUE, showWarnings = FALSE)

invisible(NULL)
