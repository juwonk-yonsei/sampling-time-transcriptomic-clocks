#!/usr/bin/env Rscript
# D1 check on GSE39445 control. Development data only. Does not open locked expression.

suppressPackageStartupMessages(library(data.table))
root <- "/HDD8T1/WORK/GeroSci/TOD"
setwd(root)
source(file.path(root, "scripts", "confirm", "run_pipeline.R"))
source(file.path(root, "scripts", "confirm", "run_locked_analysis.R"))

samples <- fread(file.path(root, "data", "processed", "GSE39445_samples.tsv"))
ctrl <- samples[condition == "control"]
ctrl[, clock_time_h := as.numeric(clock_time_h)]
expr <- read_gene_matrix(file.path(root, "data", "processed", "GSE39445_expr_gene_log2.tsv.gz"))
missing <- setdiff(ctrl$gsm, colnames(expr))
if (length(missing)) stop(paste("expression is missing", length(missing), "control samples"))
base <- expr[, ctrl$gsm, drop = FALSE]
rm(expr)
gc()

mask_random <- function(mat, fraction = 0.20, seed = 20261015) {
  set.seed(seed)
  finite <- which(is.finite(mat))
  n_mask <- as.integer(round(fraction * length(finite)))
  out <- mat
  out[sample(finite, n_mask)] <- NA_real_
  list(matrix = out, n_masked = n_mask, n_finite = length(finite))
}

mask_low_expression <- function(mat, fraction = 0.20, seed = 20261015) {
  set.seed(seed)
  mu <- rowMeans(mat, na.rm = TRUE)
  rk <- rank(mu, ties.method = "average", na.last = "keep")
  weight <- max(rk, na.rm = TRUE) + 1 - rk
  finite <- which(is.finite(mat))
  rows <- ((finite - 1L) %% nrow(mat)) + 1L
  n_mask <- as.integer(round(fraction * length(finite)))
  out <- mat
  out[sample(finite, n_mask, prob = weight[rows])] <- NA_real_
  list(matrix = out, n_masked = n_mask, n_finite = length(finite))
}

score_ra <- function(mat) {
  sc <- ra_linear_scores(mat)
  sc[, gsm := as.character(gsm)]
  sc
}

fit_one <- function(samples_i, scores, clock_id) {
  dat <- prepare_score(samples_i[, .(gsm, subject_id, clock_time_h)], scores, clock_id)
  fit <- model_m(dat)
  data.table(clock = clock_id, n_subjects = fit$n_subjects, n_samples = fit$n_samples, two_a = fit$two_a, lrt_p = fit$lrt_p, lrt_stat = fit$lrt_stat)
}

within_cor <- function(samples_i, original, treated, clock_id) {
  a <- original[clock == clock_id, .(gsm = as.character(gsm), original = value)]
  b <- treated[clock == clock_id, .(gsm = as.character(gsm), treated = value)]
  m <- merge(merge(a, b, by = "gsm"), samples_i[, .(gsm = as.character(gsm), subject_id)], by = "gsm")
  subjects <- unique(m$subject_id)
  rbindlist(lapply(subjects, function(sid) {
    block <- m[subject_id == sid]
    ok <- nrow(block) >= 3 && stats::sd(block$original) > 0 && stats::sd(block$treated) > 0
    r <- if (ok) unname(stats::cor(block$original, block$treated, method = "pearson")) else NA_real_
    data.table(subject_id = sid, n = nrow(block), pearson = r)
  }))
}

clocks <- c("RA_default", "RA_Peters", "RA_GTEx")
original <- score_ra(base)
base_fit <- rbindlist(lapply(clocks, function(clock_id) cbind(data.table(scheme = "unmasked"), fit_one(ctrl, original, clock_id))))
schemes <- list(random = mask_random(base), low_expression = mask_low_expression(base))
cor_rows <- list()
fit_rows <- list(base_fit)
audit <- list()
for (name in names(schemes)) {
  masked <- schemes[[name]]
  prepped <- d1_impute_matrix(masked$matrix)
  treated <- score_ra(prepped$matrix)
  if (any(!is.finite(treated$value))) stop(paste("D1 scores are not finite for", name))
  fit_rows[[name]] <- rbindlist(lapply(clocks, function(clock_id) cbind(data.table(scheme = name), fit_one(ctrl, treated, clock_id))))
  cor_rows[[name]] <- rbindlist(lapply(clocks, function(clock_id) {
    x <- within_cor(ctrl, original, treated, clock_id)
    x$clock <- clock_id
    x$scheme <- name
    x
  }))
  audit[[name]] <- data.table(
    scheme = name, n_masked = masked$n_masked, n_finite_before_mask = masked$n_finite,
    n_genes_dropped = prepped$n_genes_dropped, n_genes_kept = prepped$n_genes_kept,
    n_values_imputed = prepped$n_values_imputed
  )
}
out <- file.path(root, "results", "stage7b")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fwrite(rbindlist(fit_rows), file.path(out, "d1_development_model_m.tsv"), sep = "\t")
fwrite(rbindlist(cor_rows), file.path(out, "d1_development_within_pearson.tsv"), sep = "\t")
fwrite(rbindlist(audit), file.path(out, "d1_development_audit.tsv"), sep = "\t")
cat("d1 development check wrote", out, "\n")
