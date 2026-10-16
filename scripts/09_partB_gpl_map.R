#!/usr/bin/env Rscript
# Part B. Map GPL24324 probeset_id to SYMBOL. Do not read GSE122725 expression.

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(clariomshumanhttranscriptcluster.db)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
root <- normalizePath(file.path(dirname(file_arg), ".."))
out_dir <- file.path(root, "results", "stage2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pkg <- "clariomshumanhttranscriptcluster.db"
desc <- packageDescription(pkg)
desc_path <- system.file("DESCRIPTION", package = pkg)

platform_path <- file.path(root, "data", "raw", "gpl", "GPL24324_platform.txt")
con <- file(platform_path, open = "rt")
on.exit(close(con), add = TRUE)
repeat {
  line <- readLines(con, n = 1, warn = FALSE)
  if (length(line) == 0) stop("platform table not found")
  if (startsWith(line, "!platform_table_begin")) break
}
header <- strsplit(readLines(con, n = 1, warn = FALSE), "\t", fixed = TRUE)[[1]]
idx <- match("probeset_id", header)
if (is.na(idx)) stop("probeset_id column missing")
id_list <- vector("list", 50000)
n_rows <- 0L
repeat {
  line <- readLines(con, n = 1, warn = FALSE)
  if (length(line) == 0 || startsWith(line, "!platform_table_end")) break
  n_rows <- n_rows + 1L
  if (n_rows > length(id_list)) {
    id_list <- c(id_list, vector("list", length(id_list)))
  }
  parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
  id_list[[n_rows]] <- if (length(parts) >= idx) parts[[idx]] else NA_character_
}
ids <- unlist(id_list[seq_len(n_rows)], use.names = FALSE)

keys_db <- keys(clariomshumanhttranscriptcluster.db, keytype = "PROBEID")
present <- ids %in% keys_db
mapped <- AnnotationDbi::select(
  clariomshumanhttranscriptcluster.db,
  keys = unique(ids[present]),
  columns = "SYMBOL",
  keytype = "PROBEID"
)
symbol_ok <- !is.na(mapped$SYMBOL) & mapped$SYMBOL != ""
probes_with_symbol <- unique(mapped$PROBEID[symbol_ok])
row_has_symbol <- ids %in% probes_with_symbol

summary <- data.frame(
  package = pkg,
  version = as.character(desc$Version),
  license = as.character(desc$License),
  description_path = desc_path,
  n_platform_rows = length(ids),
  n_db_keys = length(keys_db),
  n_rows_key_in_db = sum(present),
  n_rows_with_symbol = sum(row_has_symbol),
  n_unique_probesets_with_symbol = length(probes_with_symbol),
  n_mapping_rows = nrow(mapped),
  n_mapping_rows_with_symbol = sum(symbol_ok),
  n_unique_symbols = length(unique(mapped$SYMBOL[symbol_ok])),
  stringsAsFactors = FALSE
)
write.table(summary, file.path(out_dir, "gpl24324_map_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

symbols <- sort(unique(mapped$SYMBOL[symbol_ok]))
writeLines(symbols, file.path(out_dir, "gpl24324_symbols.txt"))

ns <- asNamespace("RNAAgeCalc")
obj_names <- ls(ns)
writeLines(obj_names, file.path(out_dir, "rnaagecalc_namespace_objects.txt"))
ga <- get("genelist_all", ns)
rows <- list()
for (sig in c("Peters", "GTExAge", "DESeq2")) {
  coefs <- ga$blood[[sig]]
  genes <- names(coefs)
  genes <- genes[genes != "(Intercept)"]
  rows[[sig]] <- data.frame(signature = sig, gene = genes, stringsAsFactors = FALSE)
}
write.table(do.call(rbind, rows), file.path(out_dir, "rnaage_coef_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("GPL24324 rows with symbol", sum(row_has_symbol), "of", length(ids), "\n")
cat("unique symbols", length(symbols), "\n")
