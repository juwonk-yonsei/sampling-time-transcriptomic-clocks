#!/usr/bin/env Rscript
# Part E. Drop cross-dataset rhythmic genes from lin clocks. No bootstrap.

suppressPackageStartupMessages({
  library(limma)
  library(data.table)
  library(RNAAgeCalc)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
root <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(root, "scripts", "lib_cosinor.R"))
out_dir <- file.path(root, "results", "stage4", "h5")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

forms <- main_formula()
clocks <- c("RA_GTEx", "RA_Peters", "RA_default")
sig_of <- c(RA_GTEx = "GTExAge", RA_Peters = "Peters", RA_default = "DESeq2")
genelist <- get("genelist_all", asNamespace("RNAAgeCalc"))$blood
datasets <- c("GSE39445", "GSE56931", "GSE113883")

read_samples <- function(dataset) {
  path <- if (dataset == "GSE113883") {
    file.path(root, "results", "stage3", "GSE113883_samples.tsv")
  } else {
    file.path(root, "data", "processed", paste0(dataset, "_samples.tsv"))
  }
  samples <- fread(path)
  keep <- which(samples$condition == "control" & is.finite(samples$clock_time_h))
  as.data.frame(samples[keep, .(gsm, subject_id, clock_time_h)])
}

read_matrix <- function(path) {
  dt <- fread(path)
  mat <- as.matrix(dt[, -1])
  rownames(mat) <- dt[[1]]
  storage.mode(mat) <- "double"
  mat
}

rhythmic <- list()
count_rows <- list()
for (dataset in datasets) {
  cat("limma", dataset, "\n")
  flush.console()
  samples <- read_samples(dataset)
  if (dataset == "GSE113883") {
    tpm <- read_matrix(file.path(root, "results", "stage3", "expr", "GSE113883_tpm.tsv.gz"))
    log2 <- read_matrix(file.path(root, "results", "stage3", "expr", "GSE113883_log2tpm1.tsv.gz"))
    tpm <- tpm[, samples$gsm, drop = FALSE]
    log2 <- log2[, samples$gsm, drop = FALSE]
    keep <- rownames(tpm)[is.finite(rowMeans(tpm)) & rowMeans(tpm) >= 1]
    mat <- log2[keep, , drop = FALSE]
  } else {
    mat <- read_matrix(file.path(root, "data", "processed", paste0(dataset, "_expr_gene_log2.tsv.gz")))
    mat <- mat[, samples$gsm, drop = FALSE]
  }
  finite <- apply(mat, 1, function(x) all(is.finite(x)))
  mat <- mat[finite, , drop = FALSE]
  samples$subject_id <- factor(samples$subject_id)
  samples$cos_t <- cos(2 * pi * samples$clock_time_h / 24)
  samples$sin_t <- sin(2 * pi * samples$clock_time_h / 24)
  design <- model.matrix(~ subject_id + cos_t + sin_t, data = samples)
  fit <- eBayes(lmFit(mat, design))
  tt <- topTable(fit, coef = c("cos_t", "sin_t"), number = Inf, sort.by = "none")
  hits <- rownames(tt)[is.finite(tt$adj.P.Val) & tt$adj.P.Val < 0.05]
  rhythmic[[dataset]] <- hits
  writeLines(hits, file.path(out_dir, paste0("R_", sub("^GSE", "", dataset), ".txt")))
  count_rows[[dataset]] <- data.table(
    dataset = dataset,
    n_genes_tested = nrow(tt),
    n_rhythmic_q_lt_0_05 = length(hits),
    gene_filter = if (dataset == "GSE113883") "mean TPM >= 1 and finite log2" else "all finite genes"
  )
  cat(dataset, "tested", nrow(tt), "rhythmic", length(hits), "\n")
}
fwrite(rbindlist(count_rows), file.path(out_dir, "rhythmic_counts.tsv"), sep = "\t")

filter_of <- list(
  GSE39445 = union(rhythmic[["GSE56931"]], rhythmic[["GSE113883"]]),
  GSE56931 = union(rhythmic[["GSE39445"]], rhythmic[["GSE113883"]]),
  GSE113883 = union(rhythmic[["GSE39445"]], rhythmic[["GSE56931"]])
)

score_from <- function(mat, coefs, genes) {
  if (length(genes) == 0L) {
    return(rep(unname(coefs[["(Intercept)"]]), ncol(mat)))
  }
  unname(coefs[["(Intercept)"]]) + as.numeric(crossprod(coefs[genes], mat[genes, , drop = FALSE]))
}

rows <- list()
for (dataset in datasets) {
  cat("clocks", dataset, "\n")
  flush.console()
  samples <- read_samples(dataset)
  if (dataset == "GSE113883") {
    mat <- read_matrix(file.path(root, "results", "stage3", "expr", "GSE113883_log2tpm1.tsv.gz"))
  } else {
    mat <- read_matrix(file.path(root, "data", "processed", paste0(dataset, "_expr_gene_log2.tsv.gz")))
  }
  mat <- mat[, samples$gsm, drop = FALSE]
  mat <- mat[apply(mat, 1, function(x) all(is.finite(x))), , drop = FALSE]
  drop_set <- filter_of[[dataset]]
  for (clock in clocks) {
    coefs <- genelist[[sig_of[[clock]]]]
    present <- intersect(setdiff(names(coefs), "(Intercept)"), rownames(mat))
    removed <- intersect(present, drop_set)
    kept <- setdiff(present, removed)
    original <- score_from(mat, coefs, present)
    filtered <- score_from(mat, coefs, kept)
    names(original) <- colnames(mat)
    names(filtered) <- colnames(mat)
    means_o <- tapply(original, samples$subject_id, mean)
    means_f <- tapply(filtered, samples$subject_id, mean)
    for (version in c("original", "filtered")) {
      values <- if (version == "original") original else filtered
      dat <- samples
      dat$score <- unname(values)
      dat <- add_time_terms(dat)
      fit <- fit_cosinor(dat, forms$full, forms$null)
      rows[[paste(dataset, clock, version)]] <- data.table(
        dataset = dataset, clock = clock, version = version,
        n_genes_original = length(present),
        n_genes_removed = length(removed),
        n_genes_remaining = length(kept),
        n_subjects = fit$n_subjects, n_samples = fit$n_samples,
        r2_within = fit$r2_within, lrt_p = fit$lrt_p,
        two_a = fit$two_a, two_a_over_sd_between = fit$two_a_over_sd_between,
        sd_between = fit$sd_between, peak_h = fit$peak_h,
        pearson_subject_mean = cor(means_o, means_f),
        reml_converged = fit$reml_converged, ml_converged = fit$ml_converged,
        null_converged = fit$null_converged
      )
    }
  }
}
fwrite(rbindlist(rows), file.path(out_dir, "filtered_clocks.tsv"), sep = "\t")
cat("Part E done\n")
