# Part B expression matrices. GSE113883 TPM to symbols. GSE196728 counts to FPKM if the species check passed.

suppressPackageStartupMessages({
  library(data.table)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

root <- Sys.getenv("TOD_ROOT", unset = normalizePath(".."))
if (!dir.exists(file.path(root, "results"))) {
  root <- normalizePath(".")
}
out <- file.path(root, "results", "stage3")
expr_dir <- file.path(out, "expr")
dir.create(expr_dir, recursive = TRUE, showWarnings = FALSE)

map_symbols <- function(ensembl) {
  keys <- unique(ensembl)
  mapped <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = keys,
    columns = "SYMBOL",
    keytype = "ENSEMBL"
  )
  n_symbol <- tapply(mapped$SYMBOL, mapped$ENSEMBL, function(x) length(unique(x[!is.na(x)])))
  multi <- names(n_symbol)[n_symbol > 1]
  one <- mapped[!is.na(mapped$SYMBOL) & !duplicated(mapped$ENSEMBL) & !(mapped$ENSEMBL %in% multi), ]
  list(table = one, n_multi = length(multi), n_keys = length(keys))
}

prepare_113883 <- function() {
  samples <- fread(file.path(out, "GSE113883_samples.tsv"))
  tpm_path <- file.path(root, "data/raw/GSE113883/GSE113883_TPM_processed_data.txt.gz")
  raw <- fread(tpm_path)
  gene_col <- names(raw)[1]
  setnames(raw, gene_col, "gene_raw")
  raw[, ensembl_versioned := sub(".*__", "", gene_raw)]
  raw[, ensembl := sub("\\..*$", "", ensembl_versioned)]
  raw[, prefix_symbol := sub("__.*$", "", gene_raw)]
  id_summary <- data.table(
    n_rows = nrow(raw),
    n_with_ensembl = sum(grepl("^ENSG", raw$ensembl)),
    n_versioned = sum(grepl("\\.", raw$ensembl_versioned)),
    id_example = paste(head(raw$gene_raw, 5), collapse = "|")
  )
  fwrite(id_summary, file.path(out, "GSE113883_id_format.tsv"), sep = "\t")
  mapped <- map_symbols(raw$ensembl)
  map_dt <- as.data.table(mapped$table)
  raw[map_dt, symbol := i.SYMBOL, on = c(ensembl = "ENSEMBL")]
  sample_cols <- setdiff(names(raw), c("gene_raw", "ensembl_versioned", "ensembl", "prefix_symbol", "symbol"))
  raw[, mean_tpm := rowMeans(.SD, na.rm = TRUE), .SDcols = sample_cols]
  n_unmapped <- sum(is.na(raw$symbol))
  keep <- raw[!is.na(symbol)]
  setorder(keep, symbol, -mean_tpm, ensembl)
  selected <- keep[, .SD[1], by = symbol]
  title_to_gsm <- setNames(samples$gsm, samples$sample_title)
  missing_titles <- setdiff(sample_cols, names(title_to_gsm))
  if (length(missing_titles)) stop("TPM columns not in sample table: ", paste(head(missing_titles, 5), collapse = ","))
  mat <- as.matrix(selected[, ..sample_cols])
  storage.mode(mat) <- "double"
  colnames(mat) <- unname(title_to_gsm[sample_cols])
  rownames(mat) <- selected$symbol
  mat <- mat[, samples$gsm, drop = FALSE]
  log2_mat <- log2(mat + 1)
  fwrite(cbind(data.table(gene = rownames(mat)), as.data.table(mat)), file.path(expr_dir, "GSE113883_tpm.tsv.gz"), sep = "\t")
  fwrite(cbind(data.table(gene = rownames(log2_mat)), as.data.table(log2_mat)), file.path(expr_dir, "GSE113883_log2tpm1.tsv.gz"), sep = "\t")
  prefix_match <- selected[, .(
    n_selected = .N,
    n_prefix_equals_db = sum(prefix_symbol == symbol)
  )]
  summary <- data.table(
    n_rows = nrow(raw),
    n_ensembl_keys = mapped$n_keys,
    n_ensembl_multi_symbol = mapped$n_multi,
    n_unmapped = n_unmapped,
    n_mapped_rows = nrow(keep),
    n_symbols = nrow(selected),
    n_prefix_equals_db_among_selected = prefix_match$n_prefix_equals_db,
    org_hs_eg_db = as.character(packageVersion("org.Hs.eg.db"))
  )
  fwrite(summary, file.path(out, "GSE113883_symbol_map_summary.tsv"), sep = "\t")
  print(summary)
}

prepare_196728 <- function() {
  decision <- fread(file.path(out, "gse196728_species_decision.tsv"))
  count_path <- file.path(expr_dir, "GSE196728_counts.tsv.gz")
  if (!isTRUE(decision$use_for_analysis) || !file.exists(count_path)) {
    message("GSE196728 expression skipped")
    return(invisible(NULL))
  }
  count_dt <- fread(count_path)
  count_genes <- count_dt[[1]]
  counts <- as.matrix(count_dt[, -1])
  storage.mode(counts) <- "double"
  rownames(counts) <- count_genes
  warns <- character()
  fpkm <- withCallingHandlers(
    RNAAgeCalc::count2FPKM(counts, idtype = "SYMBOL"),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  fpkm <- as.matrix(fpkm)
  log2_fpkm <- log2(fpkm + 1)
  fwrite(cbind(data.table(gene = rownames(fpkm)), as.data.table(fpkm)), file.path(expr_dir, "GSE196728_fpkm.tsv.gz"), sep = "\t")
  fwrite(cbind(data.table(gene = rownames(log2_fpkm)), as.data.table(log2_fpkm)), file.path(expr_dir, "GSE196728_log2fpkm1.tsv.gz"), sep = "\t")
  na_len <- sum(is.na(fpkm[, 1]) & is.na(fpkm[, ncol(fpkm)]))
  # A gene with missing length becomes NA across samples in getRPKM.
  n_na_gene <- sum(apply(fpkm, 1, function(x) all(is.na(x))))
  info <- data.table(
    n_genes = nrow(counts),
    n_samples = ncol(counts),
    n_genes_all_na_fpkm = n_na_gene,
    length_source = "RNAAgeCalc internal lengthDB, via org.Hs.eg.db SYMBOL to ENSEMBL, column bp_length",
    function_name = "RNAAgeCalc::count2FPKM",
    idtype = "SYMBOL",
    warnings = paste(unique(warns), collapse = " | ")
  )
  fwrite(info, file.path(out, "GSE196728_fpkm_method.tsv"), sep = "\t")
  print(info)
}

prepare_113883()
prepare_196728()

ver_path <- file.path(root, "config", "env_versions.txt")
note <- paste0(
  "\nstage3 ", format(Sys.time(), "%Y-%m-%d"),
  " org.Hs.eg.db ", as.character(packageVersion("org.Hs.eg.db")),
  " ", packageDescription("org.Hs.eg.db")$License, "\n"
)
cat(note, file = ver_path, append = TRUE)
capture.output(sessionInfo(), file = ver_path, append = TRUE)
cat("\n", file = ver_path, append = TRUE)
