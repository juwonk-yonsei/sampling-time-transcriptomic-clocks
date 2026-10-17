#!/usr/bin/env Rscript
# Part B checks. Feature counts open locked L1. The D11 mask check uses GSE39445 only.

suppressPackageStartupMessages(library(data.table))
root <- "/HDD8T1/WORK/GeroSci/TOD"
setwd(root)
source(file.path(root, "scripts", "confirm", "run_pipeline.R"))
source(file.path(root, "scripts", "confirm", "run_locked_analysis.R"))

demo <- matrix(c(1, 5, NA, 6, 3, NA, NA, 8), nrow = 2, byrow = FALSE)
colnames(demo) <- c("a1", "a2", "b1", "b2")
rownames(demo) <- c("g1", "g2")
got <- d11_participant_mean_matrix(demo, c("A", "A", "B", "B"))
stopifnot(abs(got$matrix["g1", "a2"] - 1) < 1e-8)
stopifnot(abs(got$matrix["g1", "b2"] - 3) < 1e-8)
stopifnot(abs(got$matrix["g2", "b1"] - 8) < 1e-8)
empty <- matrix(c(NA, 4, NA, 6), nrow = 1)
colnames(empty) <- c("a1", "b1", "a2", "b2")
got_empty <- d11_participant_mean_matrix(empty, c("A", "B", "A", "B"))
stopifnot(abs(got_empty$matrix[1, "a1"] - 5) < 1e-8)
cat("d11 unit check passed\n")

out <- file.path(root, "results", "stage7c")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
work <- file.path(out, "l1_build")
dir.create(work, recursive = TRUE, showWarnings = FALSE)

l1_meta <- fread(file.path(root, "data", "processed", "GSE48113_samples.tsv"))
l1_meta <- l1_meta[!gsm %in% c("GSM1168837", "GSM1168735")]
day1 <- l1_meta[condition == "control" & period == "D1"]
day4 <- l1_meta[period == "D4"]
expect_n("L1 day1", nrow(day1), uniqueN(day1$subject_id), 147L, 22L)
expect_n("L1 day4", nrow(day4), uniqueN(day4$subject_id), 138L, 21L)
day1[, analysis_set := "day1"]
day4[, analysis_set := "day4"]
l1_samples <- rbind(day1, day4, fill = TRUE)
fwrite(l1_samples, file.path(work, "samples_scoring.tsv"), sep = "\t")
l1_raw <- series_table(load_locked_expression("L1", root))
l1_series_path <- file.path(work, "L1_series_table.tsv.gz")
fwrite(l1_raw, l1_series_path, sep = "\t")
rm(l1_raw)
gc()
status <- system2(python_bin, c(
  file.path(root, "scripts", "confirm", "preprocess_array.py"),
  "--matrix", l1_series_path,
  "--platform", file.path(root, "data", "raw", "gpl", "GPL15331_platform.txt"),
  "--symbol-column", "GENE_SYMBOL",
  "--samples", file.path(work, "samples_scoring.tsv"),
  "--out-expr", file.path(work, "log2_genes.tsv.gz"),
  "--out-audit", file.path(work, "preprocess_audit.tsv")
))
if (!identical(status, 0L)) stop("L1 preprocessing failed")
unlink(l1_series_path)
l1_log <- read_gene_matrix(file.path(work, "log2_genes.tsv.gz"))
l1_for <- l1_log[, l1_samples$gsm, drop = FALSE]
l1_d1 <- d1_impute_matrix(l1_for)
write_gene_matrix(l1_d1$matrix, file.path(work, "log2_genes_d1_imputed.tsv"))
d11_l1 <- d11_participant_mean_matrix(l1_for, l1_samples$subject_id)
write_gene_matrix(d11_l1$matrix, file.path(work, "log2_genes_d11.tsv"))
fwrite(data.table(
  n_missing = d11_l1$n_missing, n_filled_participant = d11_l1$n_filled_participant,
  n_filled_gene_mean = d11_l1$n_filled_gene_mean, n_genes_all_missing = d11_l1$n_genes_all_missing,
  n_nonfinite_after = d11_l1$n_nonfinite_after
), file.path(out, "d11_l1_audit.tsv"), sep = "\t")
rm(l1_log, l1_for, l1_d1, d11_l1)
gc()
count_script <- file.path(root, "scripts", "48_part7c_feature_count.py")
specs <- list(
  list(version = "frozen", expr = file.path(work, "log2_genes.tsv.gz")),
  list(version = "D2_knn", expr = file.path(work, "log2_genes_d1_imputed.tsv")),
  list(version = "D11_participant_mean", expr = file.path(work, "log2_genes_d11.tsv"))
)
for (spec in specs) {
  status <- system2(python_bin, c(
    count_script, "--expr", spec$expr, "--samples", file.path(work, "samples_scoring.tsv"),
    "--version", spec$version, "--out", file.path(out, paste0("tams_features_", spec$version, ".tsv"))
  ))
  if (!identical(status, 0L)) stop(paste("feature count failed", spec$version))
}
fwrite(rbindlist(lapply(specs, function(spec) fread(file.path(out, paste0("tams_features_", spec$version, ".tsv"))))), file.path(out, "tams_feature_counts.tsv"), sep = "\t")
unlink(file.path(work, c("log2_genes.tsv.gz", "log2_genes_d1_imputed.tsv", "log2_genes_d11.tsv")))
cat("feature counts wrote", file.path(out, "tams_feature_counts.tsv"), "\n")

mask_random <- function(mat, fraction = 0.20, seed = 20261015) {
  set.seed(seed)
  finite <- which(is.finite(mat))
  n_mask <- as.integer(round(fraction * length(finite)))
  out_m <- mat
  out_m[sample(finite, n_mask)] <- NA_real_
  list(matrix = out_m, n_masked = n_mask, n_finite = length(finite))
}

mask_low_expression <- function(mat, fraction = 0.20, seed = 20261015) {
  set.seed(seed)
  mu <- rowMeans(mat, na.rm = TRUE)
  rk <- rank(mu, ties.method = "average", na.last = "keep")
  weight <- max(rk, na.rm = TRUE) + 1 - rk
  finite <- which(is.finite(mat))
  rows <- ((finite - 1L) %% nrow(mat)) + 1L
  n_mask <- as.integer(round(fraction * length(finite)))
  out_m <- mat
  out_m[sample(finite, n_mask, prob = weight[rows])] <- NA_real_
  list(matrix = out_m, n_masked = n_mask, n_finite = length(finite))
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

samples <- fread(file.path(root, "data", "processed", "GSE39445_samples.tsv"))
ctrl <- samples[condition == "control"]
ctrl[, clock_time_h := as.numeric(clock_time_h)]
expr <- read_gene_matrix(file.path(root, "data", "processed", "GSE39445_expr_gene_log2.tsv.gz"))
missing <- setdiff(ctrl$gsm, colnames(expr))
if (length(missing)) stop(paste("expression is missing", length(missing), "control samples"))
base <- expr[, ctrl$gsm, drop = FALSE]
rm(expr)
gc()
clocks <- c("RA_default", "RA_Peters", "RA_GTEx")
original <- score_ra(base)
base_fit <- rbindlist(lapply(clocks, function(clock_id) cbind(data.table(scheme = "unmasked"), fit_one(ctrl, original, clock_id))))
schemes <- list(random = mask_random(base), low_expression = mask_low_expression(base))
cor_rows <- list()
fit_rows <- list(base_fit)
audit <- list()
for (name in names(schemes)) {
  masked <- schemes[[name]]
  prepped <- d11_participant_mean_matrix(masked$matrix, ctrl$subject_id)
  treated <- score_ra(prepped$matrix)
  if (any(!is.finite(treated$value))) stop(paste("D11 scores are not finite for", name))
  fit_rows[[name]] <- rbindlist(lapply(clocks, function(clock_id) cbind(data.table(scheme = name), fit_one(ctrl, treated, clock_id))))
  cor_rows[[name]] <- rbindlist(lapply(clocks, function(clock_id) {
    x <- within_cor(ctrl, original, treated, clock_id)
    x$clock <- clock_id
    x$scheme <- name
    x
  }))
  audit[[name]] <- data.table(
    scheme = name, n_masked = masked$n_masked, n_finite_before_mask = masked$n_finite,
    n_missing = prepped$n_missing, n_filled_participant = prepped$n_filled_participant,
    n_filled_gene_mean = prepped$n_filled_gene_mean, n_nonfinite_after = prepped$n_nonfinite_after
  )
}
fwrite(rbindlist(fit_rows), file.path(out, "d11_development_model_m.tsv"), sep = "\t")
fwrite(rbindlist(cor_rows), file.path(out, "d11_development_within_pearson.tsv"), sep = "\t")
fwrite(rbindlist(audit), file.path(out, "d11_development_audit.tsv"), sep = "\t")
cat("d11 development check wrote", out, "\n")
