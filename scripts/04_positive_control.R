#!/usr/bin/env Rscript
# Part D. Circadian positive control on development control samples.

suppressPackageStartupMessages({
  library(lme4)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
root <- normalizePath(file.path(dirname(file_arg), ".."))

genes <- c("PER1", "PER2", "PER3", "NR1D1", "NR1D2", "ARNTL", "CRY1", "CRY2", "DBP", "TEF")
datasets <- c("GSE39445", "GSE56931")
out_dir <- file.path(root, "results", "stage1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

peak_0_24 <- function(beta_sin, beta_cos) {
  peak <- (24 / (2 * pi)) * atan2(beta_sin, beta_cos)
  peak %% 24
}

circ_dist_0_12 <- function(a, b) {
  d <- abs(a - b) %% 24
  pmin(d, 24 - d)
}

fit_gene <- function(df) {
  df$cos_t <- cos(2 * pi * df$clock_time_h / 24)
  df$sin_t <- sin(2 * pi * df$clock_time_h / 24)
  df$subject_id <- factor(df$subject_id)
  warn <- character()
  capture <- function(expr) {
    withCallingHandlers(
      tryCatch(expr, error = function(e) e),
      warning = function(w) {
        warn <<- c(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  }
  reml <- capture(lmer(expr ~ cos_t + sin_t + (1 | subject_id), data = df, REML = TRUE))
  ml <- capture(lmer(expr ~ cos_t + sin_t + (1 | subject_id), data = df, REML = FALSE))
  null <- capture(lmer(expr ~ 1 + (1 | subject_id), data = df, REML = FALSE))
  failed <- function(x) inherits(x, "error")
  row <- data.frame(
    n_subjects = length(unique(df$subject_id)),
    n_samples = nrow(df),
    beta_cos = NA_real_,
    beta_sin = NA_real_,
    amplitude = NA_real_,
    peak_h = NA_real_,
    lrt_stat = NA_real_,
    lrt_df = 2,
    lrt_p = NA_real_,
    reml_ok = !failed(reml),
    ml_ok = !failed(ml),
    null_ok = !failed(null),
    reml_singular = if (!failed(reml)) isTRUE(isSingular(reml, tol = 1e-4)) else NA,
    ml_singular = if (!failed(ml)) isTRUE(isSingular(ml, tol = 1e-4)) else NA,
    null_singular = if (!failed(null)) isTRUE(isSingular(null, tol = 1e-4)) else NA,
    warnings = paste(unique(warn), collapse = " | "),
    stringsAsFactors = FALSE
  )
  if (!failed(reml)) {
    beta <- fixef(reml)
    row$beta_cos <- unname(beta["cos_t"])
    row$beta_sin <- unname(beta["sin_t"])
    row$amplitude <- sqrt(row$beta_cos^2 + row$beta_sin^2)
    row$peak_h <- peak_0_24(row$beta_sin, row$beta_cos)
  }
  if (!failed(ml) && !failed(null)) {
    stat <- as.numeric(2 * (logLik(ml) - logLik(null)))
    row$lrt_stat <- stat
    row$lrt_p <- pchisq(stat, df = 2, lower.tail = FALSE)
  }
  if (failed(reml)) row$warnings <- paste(row$warnings, conditionMessage(reml), sep = " | ")
  row
}

results <- list()
for (gse in datasets) {
  expr <- read.delim(gzfile(file.path(root, "data", "processed", paste0(gse, "_expr_gene_log2.tsv.gz"))),
                      check.names = FALSE, row.names = 1)
  meta <- read.delim(file.path(root, "data", "processed", paste0(gse, "_samples.tsv")),
                      check.names = FALSE, stringsAsFactors = FALSE)
  meta <- meta[meta$condition == "control" & meta$time_source == "given" & meta$clock_time_h != "NA", ]
  meta$clock_time_h <- as.numeric(meta$clock_time_h)
  keep <- intersect(meta$gsm, colnames(expr))
  meta <- meta[match(keep, meta$gsm), ]
  expr <- expr[, keep, drop = FALSE]
  for (gene in genes) {
    if (!(gene %in% rownames(expr))) {
      results[[length(results) + 1]] <- data.frame(
        dataset = gse, gene = gene, present = FALSE,
        n_subjects = NA, n_samples = NA,
        beta_cos = NA, beta_sin = NA, amplitude = NA, peak_h = NA,
        lrt_stat = NA, lrt_df = 2, lrt_p = NA,
        reml_ok = FALSE, ml_ok = FALSE, null_ok = FALSE,
        warnings = "absent",
        stringsAsFactors = FALSE
      )
      next
    }
    df <- data.frame(
      expr = as.numeric(expr[gene, ]),
      clock_time_h = meta$clock_time_h,
      subject_id = meta$subject_id,
      stringsAsFactors = FALSE
    )
    df <- df[is.finite(df$expr) & is.finite(df$clock_time_h), ]
    fit <- fit_gene(df)
    fit$dataset <- gse
    fit$gene <- gene
    fit$present <- TRUE
    results[[length(results) + 1]] <- fit
    message(gse, " ", gene, " peak=", fit$peak_h, " p=", fit$lrt_p)
  }
}

tab <- do.call(rbind, results)
tab <- tab[, c("dataset", "gene", "present", "n_subjects", "n_samples",
               "beta_cos", "beta_sin", "amplitude", "peak_h",
               "lrt_stat", "lrt_df", "lrt_p", "reml_ok", "ml_ok", "null_ok",
               "reml_singular", "ml_singular", "null_singular", "warnings")]
write.table(tab, file.path(out_dir, "positive_control.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

decision_rows <- list()
for (gse in datasets) {
  sub <- tab[tab$dataset == gse, ]
  peaks <- setNames(sub$peak_h, sub$gene)
  need_a <- c("PER1", "PER3", "NR1D1", "NR1D2")
  a_ok <- all(need_a %in% names(peaks)) && all(is.finite(peaks[need_a])) &&
    all(peaks[need_a] >= 4 & peaks[need_a] <= 12)
  b_ok <- all(c("PER3", "ARNTL") %in% names(peaks)) &&
    all(is.finite(peaks[c("PER3", "ARNTL")])) &&
    circ_dist_0_12(peaks["PER3"], peaks["ARNTL"]) >= 6
  dist <- if (all(c("PER3", "ARNTL") %in% names(peaks)) && all(is.finite(peaks[c("PER3", "ARNTL")]))) {
    circ_dist_0_12(peaks["PER3"], peaks["ARNTL"])
  } else NA_real_
  decision_rows[[gse]] <- data.frame(
    dataset = gse,
    criterion_a = a_ok,
    criterion_b = b_ok,
    per3_arntl_circ_dist_h = as.numeric(dist),
    decision = if (isTRUE(a_ok) && isTRUE(b_ok)) "pass" else "review",
    stringsAsFactors = FALSE
  )
}
dec <- do.call(rbind, decision_rows)
write.table(dec, file.path(out_dir, "positive_control_decision.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(dec)
