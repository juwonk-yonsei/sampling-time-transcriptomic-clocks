#!/usr/bin/env Rscript
# Part F is run only for development datasets whose Part D decision is pass.

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
root <- normalizePath(file.path(dirname(file_arg), ".."))
dec_path <- file.path(root, "results", "stage1", "positive_control_decision.tsv")
dec <- read.delim(dec_path, stringsAsFactors = FALSE)
passed <- dec$dataset[dec$decision == "pass"]
status <- data.frame(
  dataset = dec$dataset,
  decision = dec$decision,
  part_f = ifelse(dec$decision == "pass", "run", "not_run"),
  stringsAsFactors = FALSE
)
out <- file.path(root, "results", "stage1", "partF_status.tsv")
write.table(status, out, sep = "\t", quote = FALSE, row.names = FALSE)
if (length(passed) == 0) {
  message("Part F not run: no development dataset passed the positive control.")
  quit(status = 0)
}
stop("Part F implementation is required for: ", paste(passed, collapse = ", "))
