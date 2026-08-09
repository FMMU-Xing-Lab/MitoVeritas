#!/usr/bin/env Rscript
# Prepare the list of candidate variant sites for ML feature extraction.
#
# Usage:
#   Rscript get_output_true_file.R <sample_dir> <cutoff> <vaf>
#
#   sample_dir : directory that contains sample_name.txt and one sub-directory per sample
#   cutoff     : maximum mutation frequency (%) used to keep low-frequency variants (0-100).
#                By default the filter is NOT applied (to reproduce the original behaviour);
#                set the environment variable APPLY_CUTOFF=1 to enable it.
#   vaf        : frequency-threshold tag used in the hetro/homo file names (e.g. 0.1)
#
# For each sample this script merges the hetro/homo mutation files produced by mtDNApipe
# and writes an "output_true_<vaf>" file with columns:
#   chrM, position, major allele, minor allele, sample name, mutation frequency

options(stringsAsFactors = FALSE)

args <- commandArgs(TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript get_output_true_file.R <sample_dir> <cutoff> <vaf>")
}

path   <- normalizePath(args[1], mustWork = TRUE)
cutoff <- as.numeric(args[2])
vaf    <- args[3]
setwd(path)

samples <- as.matrix(read.table(file = "sample_name.txt", sep = "\t", header = FALSE))
if (nrow(samples) == 0) {
  stop("sample_name.txt is empty")
}
samples <- samples[samples[, 1] != "", , drop = FALSE]

for (i in seq_len(nrow(samples))) {
  sample_name <- samples[i, 1]
  sample_path <- file.path(path, sample_name)
  if (!dir.exists(sample_path)) {
    warning(sprintf("Sample directory not found: %s", sample_path))
    next
  }
  setwd(sample_path)

  hetro.filename <- grep(paste0(".hetro_", vaf, ".txt"), dir(), value = TRUE)
  homo.filename  <- grep(paste0(".homo_", vaf, ".txt"), dir(), value = TRUE)
  if (length(hetro.filename) > 1) {
    warning(sprintf("Multiple hetro files found for %s, using the first one", sample_name))
    hetro.filename <- hetro.filename[1]
  }
  if (length(homo.filename) > 1) {
    warning(sprintf("Multiple homo files found for %s, using the first one", sample_name))
    homo.filename <- homo.filename[1]
  }

  if (length(hetro.filename) == 0 & length(homo.filename) != 0) {
    print(paste0("The hetro/homo files do not exist!!!(", sample_name, ")"))
  }

  hetro <- NULL
  homo  <- NULL
  if (length(hetro.filename) != 0 && file.exists(hetro.filename) && file.info(hetro.filename)$size > 0) {
    hetro <- read.table(file = hetro.filename, sep = "\t", header = FALSE, colClasses = "character")
  }
  if (length(homo.filename) != 0 && file.exists(homo.filename) && file.info(homo.filename)$size > 0) {
    homo <- read.table(file = homo.filename, sep = "\t", header = FALSE, colClasses = "character")
  }

  # drop a trailing all-empty column if present (files exported by mtDNApipe may end with a tab)
  drop_trailing_empty <- function(dt) {
    if (is.null(dt) || ncol(dt) == 0) return(dt)
    if (all(dt[, ncol(dt)] == "")) return(dt[, -ncol(dt), drop = FALSE])
    dt
  }
  hetro <- drop_trailing_empty(hetro)
  homo  <- drop_trailing_empty(homo)

  table <- NULL
  if (!is.null(hetro) && !is.null(homo)) {
    table <- rbind(hetro, homo)
  } else if (!is.null(hetro)) {
    table <- hetro
  } else if (!is.null(homo)) {
    table <- homo
  }

  if (is.null(table) || nrow(table) == 0) {
    print(paste0("No variants for sample ", sample_name,
                 "; writing empty output_true_", vaf))
    write.table(data.frame(), file = paste0("output_true_", vaf), sep = "\t",
                col.names = FALSE, row.names = FALSE, quote = FALSE)
    setwd(path)
    next
  }

  output <- cbind("chrM", table[, 1:3], sample_name, table[, ncol(table)])

  # Optional filter: keep only variants with frequency <= cutoff (%).
  # Disabled by default to reproduce the original behaviour.
  if (identical(Sys.getenv("APPLY_CUTOFF"), "1")) {
    output <- output[which(as.numeric(output[, 6]) <= as.numeric(cutoff) * 0.01), , drop = FALSE]
  }

  write.table(output, file = paste0("output_true_", vaf), sep = "\t",
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  setwd(path)
}

print(paste0("get_output_true_file.R finished: ", nrow(samples), " sample(s) processed"))
