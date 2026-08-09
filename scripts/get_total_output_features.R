#!/usr/bin/env Rscript
# Aggregate the per-sample output_feature_<vaf> files into one table.
#
# Usage:
#   Rscript get_total_output_features.R <sample_dir> <cutoff>
#
# Output:
#   <sample_dir>/Total_output_feature_<vaf>.txt

options(stringsAsFactors = FALSE)

args <- commandArgs(TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript get_total_output_features.R <sample_dir> <cutoff>")
}

path   <- normalizePath(args[1], mustWork = TRUE)
cutoff <- as.numeric(args[2])
setwd(path)

samples <- as.matrix(read.table(file = "sample_name.txt", sep = "\t", header = FALSE))
if (nrow(samples) == 0) {
  stop("sample_name.txt is empty")
}
samples <- samples[samples[, 1] != "", , drop = FALSE]
feature <- list()

for (i in seq_len(nrow(samples))) {
  sample_name <- samples[i, 1]
  sample_path <- file.path(path, sample_name)
  setwd(sample_path)

  # the python feature extractor writes output_feature_<vaf>.tmp first, then the
  # cleaned output_feature_<vaf>; aggregate the cleaned file
  files <- list.files(pattern = "^output_feature_.*")
  files <- files[!grepl("\\.tmp$", files)]

  if (length(files) == 0) {
    print(paste0("The output_feature file does not exist!!!(", sample_name, ")"))
    setwd(path)
    next
  }
  filename <- files[1]

  if (!file.exists(filename) || file.info(filename)$size == 0) {
    print(paste0("The file ", filename, " is empty or does not exist!!!(", sample_name, ")"))
    setwd(path)
    next
  }

  feature[[i]] <- tryCatch(
    read.table(file = filename, sep = "\t", header = TRUE, check.names = FALSE),
    error = function(e) {
      print(paste0("Failed to read ", filename, " (", sample_name, "): ", conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(feature[[i]])) {
    # sanity check: number of rows should match output_true_<vaf> (when the cutoff
    # filter is applied this only reports a warning)
    true_file <- list.files(pattern = "^output_true_")
    var.table <- NULL
    if (length(true_file) > 0) {
      var.table <- tryCatch(
        read.table(file = true_file[1], sep = "\t", header = FALSE, check.names = FALSE),
        error = function(e) NULL
      )
    }
    if (!is.null(var.table) && nrow(var.table) > 0 && ncol(var.table) >= 6) {
      var.table <- var.table[which(as.numeric(var.table[, 6]) <= as.numeric(cutoff) * 0.01), , drop = FALSE]
      if (nrow(var.table) != nrow(feature[[i]])) {
        print(paste0("Some variants may be missing in the output_feature file!!!(", sample_name, ")"))
      }
    }
  }
  setwd(path)
}

if (length(feature) == 0 || all(vapply(feature, is.null, logical(1)))) {
  stop("No per-sample output_feature files could be aggregated; check the warning messages above")
}
result <- do.call(rbind, feature)
setwd(path)
write.table(result, file = paste0("Total_", filename, ".txt"), sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)

print(paste0("get_total_output_features.R finished: ", nrow(result),
             " variant rows written to Total_", filename, ".txt"))
