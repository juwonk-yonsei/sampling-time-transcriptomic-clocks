#!/usr/bin/env Rscript
# Part E. RNAAgeCalc inspection, linear scores, and official predict_age.

suppressPackageStartupMessages({
  library(RNAAgeCalc)
  library(impute)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
root <- normalizePath(file.path(dirname(file_arg), ".."))
out_dir <- file.path(root, "results", "stage1")
clock_dir <- file.path(out_dir, "clocks")
gene_dir <- file.path(root, "data", "processed", "rnaage_genes")
dir.create(clock_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

ns <- asNamespace("RNAAgeCalc")
genelist_all <- get("genelist_all", envir = ns)
pkg <- system.file(package = "RNAAgeCalc")
sha256_file <- function(path) {
  line <- system2("sha256sum", path, stdout = TRUE)
  strsplit(line, " ", fixed = TRUE)[[1]][1]
}

hash_paths <- c(
  file.path(pkg, "DESCRIPTION"),
  file.path(pkg, "R", "sysdata.rdb"),
  file.path(pkg, "R", "sysdata.rdx"),
  file.path(pkg, "data", "fpkmExample.RData"),
  file.path(pkg, "data", "countExample.RData")
)
hash_tab <- data.frame(
  path = hash_paths,
  bytes = file.info(hash_paths)$size,
  sha256 = vapply(hash_paths, sha256_file, ""),
  stringsAsFactors = FALSE
)
write.table(hash_tab, file.path(out_dir, "rnaagecalc_file_hashes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

clocks <- data.frame(
  clock = c("RA_GTEx", "RA_Peters", "RA_default"),
  signature = c("GTExAge", "Peters", "DESeq2"),
  documented_n_genes = c(NA_integer_, 1497L, NA_integer_),
  documented_source = c(
    "vignette does not state a gene count",
    "RNAAge-vignette.Rmd: Peters signature contains the 1,497 genes",
    "vignette does not state a gene count; blood default when signature is omitted is DESeq2"
  ),
  stringsAsFactors = FALSE
)

audit_rows <- list()
gene_rows <- list()
for (i in seq_len(nrow(clocks))) {
  sig <- clocks$signature[i]
  coefs <- genelist_all[["blood"]][[sig]]
  n_genes <- length(coefs) - 1L
  documented <- clocks$documented_n_genes[i]
  match_doc <- if (is.na(documented)) NA else n_genes == documented
  pass_next <- if (is.na(match_doc)) TRUE else isTRUE(match_doc)
  audit_rows[[i]] <- data.frame(
    clock = clocks$clock[i],
    tissue = "blood",
    stype = "all",
    object = "genelist_all$blood",
    signature = sig,
    intercept_name = names(coefs)[1],
    intercept = unname(coefs[1]),
    n_genes = n_genes,
    documented_n_genes = documented,
    count_matches_document = match_doc,
    pass_to_next_stage = pass_next,
    stringsAsFactors = FALSE
  )
  gene_rows[[i]] <- data.frame(
    clock = clocks$clock[i],
    signature = sig,
    gene = names(coefs)[-1],
    coefficient = unname(coefs[-1]),
    stringsAsFactors = FALSE
  )
}
audit <- do.call(rbind, audit_rows)
genes <- do.call(rbind, gene_rows)
write.table(audit, file.path(out_dir, "rnaage_signature_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(genes, file.path(gene_dir, "blood_signature_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# What impute.knn does to an all-NA row, which is how predict_age adds absent genes.
demo <- rbind(MISSING = rep(NA_real_, 5), matrix(c(1, 2, 3, 4, 5, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE))
rownames(demo) <- c("MISSING", "G1", "G2")
imp <- impute.knn(demo, maxp = nrow(demo))$data
demo_out <- data.frame(
  sample = paste0("s", seq_len(ncol(demo))),
  imputed_missing_gene = imp["MISSING", ],
  column_mean_of_observed_genes = colMeans(demo[-1, , drop = FALSE]),
  stringsAsFactors = FALSE
)
demo_out$abs_diff <- abs(demo_out$imputed_missing_gene - demo_out$column_mean_of_observed_genes)
write.table(demo_out, file.path(out_dir, "impute_knn_allna_check.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

lin_predict <- function(log2_mat, coefs) {
  present <- intersect(names(coefs)[-1], rownames(log2_mat))
  used <- coefs[present]
  score <- as.numeric(crossprod(used, log2_mat[present, , drop = FALSE])) + unname(coefs[1])
  names(score) <- colnames(log2_mat)
  list(score = score, n_used = length(present), n_signature = length(coefs) - 1L)
}

read_log2 <- function(gse) {
  path <- file.path(root, "data", "processed", paste0(gse, "_expr_gene_log2.tsv.gz"))
  expr <- read.delim(gzfile(path), check.names = FALSE, row.names = 1)
  as.matrix(expr)
}

score_rows <- list()
message_rows <- list()
datasets <- c("GSE39445", "GSE56931")
for (gse in datasets) {
  log2_mat <- read_log2(gse)
  linear <- 2^log2_mat
  for (i in seq_len(nrow(clocks))) {
    sig <- clocks$signature[i]
    clock <- clocks$clock[i]
    coefs <- genelist_all[["blood"]][[sig]]
    lin <- lin_predict(log2_mat, coefs)
    score_rows[[length(score_rows) + 1]] <- data.frame(
      dataset = gse,
      gsm = names(lin$score),
      clock = clock,
      method = "lin",
      value = unname(lin$score),
      n_genes_used = lin$n_used,
      n_genes_signature = lin$n_signature,
      n_genes_imputed = NA_integer_,
      stringsAsFactors = FALSE
    )
    msgs <- character()
    timed <- system.time({
      fit <- withCallingHandlers(
        if (clock == "RA_default") {
          predict_age(linear, tissue = "blood", exprtype = "FPKM",
                      idtype = "SYMBOL", stype = "all")
        } else {
          predict_age(linear, tissue = "blood", exprtype = "FPKM",
                      idtype = "SYMBOL", stype = "all", signature = sig)
        },
        message = function(m) {
          msgs <<- c(msgs, conditionMessage(m))
          invokeRestart("muffleMessage")
        }
      )
    })
    n_missing <- sum(!names(coefs)[-1] %in% rownames(log2_mat))
    score_rows[[length(score_rows) + 1]] <- data.frame(
      dataset = gse,
      gsm = rownames(fit),
      clock = clock,
      method = "official",
      value = fit$RNAAge,
      n_genes_used = length(coefs) - 1L,
      n_genes_signature = length(coefs) - 1L,
      n_genes_imputed = n_missing,
      stringsAsFactors = FALSE
    )
    message_rows[[length(message_rows) + 1]] <- data.frame(
      dataset = gse,
      clock = clock,
      signature = sig,
      elapsed_sec = unname(timed[["elapsed"]]),
      user_sec = unname(timed[["user.self"]]),
      sys_sec = unname(timed[["sys.self"]]),
      n_genes_imputed = n_missing,
      messages = paste(msgs, collapse = " || "),
      stringsAsFactors = FALSE
    )
    message(gse, " ", clock, " official elapsed_sec=", unname(timed[["elapsed"]]))
  }
  rm(log2_mat, linear)
  gc()
}

scores <- do.call(rbind, score_rows)
messages <- do.call(rbind, message_rows)
write.table(messages, file.path(out_dir, "rnaage_official_runs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

dist_rows <- list()
for (gse in datasets) {
  sub <- scores[scores$dataset == gse, c("gsm", "clock", "method", "value")]
  write.table(sub, file.path(clock_dir, paste0(gse, "_rnaage.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  for (clock in unique(sub$clock)) {
    for (method in unique(sub$method)) {
      vals <- sub$value[sub$clock == clock & sub$method == method]
      dist_rows[[length(dist_rows) + 1]] <- data.frame(
        dataset = gse, clock = clock, method = method,
        n = length(vals), mean = mean(vals), sd = sd(vals),
        min = min(vals), max = max(vals),
        stringsAsFactors = FALSE
      )
    }
  }
}
write.table(do.call(rbind, dist_rows), file.path(out_dir, "rnaage_value_distribution.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("RNAAGE_DONE")
