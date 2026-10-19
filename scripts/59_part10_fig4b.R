#!/usr/bin/env Rscript
# Fig. 4B post hoc composition peaks. Does not write to results/confirm/.

suppressPackageStartupMessages({
  library(lme4)
  library(data.table)
  library(MCPcounter)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
root <- normalizePath(file.path(script_dir, ".."))
source(file.path(root, "scripts", "lib_cosinor.R"))
source(file.path(root, "scripts", "confirm", "lib_tests.R"))

out_dir <- file.path(root, "results", "stage10")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
say <- function(...) message(format(Sys.time(), "%H:%M:%S"), " ", paste(...))

read_gene_matrix <- function(path) {
  dt <- fread(path)
  genes <- as.character(dt[[1]])
  mat <- as.matrix(dt[, -1])
  storage.mode(mat) <- "double"
  rownames(mat) <- genes
  colnames(mat) <- as.character(colnames(mat))
  mat
}

mcp_scores <- function(log_mat) {
  genes <- read.table(
    file.path(root, "data", "processed", "ref", "MCPcounter_genes.txt"),
    sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "\"", comment.char = ""
  )
  est <- MCPcounter.estimate(log_mat, featuresType = "HUGO_symbols", genes = genes)
  as.data.frame(t(est))
}

boot_control <- lmerControl(calc.derivs = FALSE, check.conv.singular = "ignore")

fit_peak <- function(dat) {
  fit <- suppressWarnings(tryCatch(
    lmer(score ~ cos_t + sin_t + (1 | subject_id), data = dat, REML = TRUE, control = boot_control),
    error = function(e) e
  ))
  if (inherits(fit, "error") || !optimizer_converged(fit)) return(NA_real_)
  beta <- fixef(fit)
  peak_0_24(unname(beta["sin_t"]), unname(beta["cos_t"]))
}

boot_peaks <- function(dat, n_boot) {
  parts <- split(dat, dat$subject_id)
  ids <- names(parts)
  n <- length(ids)
  peaks <- rep(NA_real_, n_boot)
  for (i in seq_len(n_boot)) {
    draw <- sample.int(n, n, replace = TRUE)
    pieces <- vector("list", n)
    for (j in seq_len(n)) {
      part <- parts[[draw[j]]]
      part$subject_id <- j
      pieces[[j]] <- part
    }
    peaks[i] <- fit_peak(rbindlist(pieces))
  }
  peaks
}

load_development <- function(dataset, sample_path, mcp_path) {
  samples <- fread(sample_path)
  samples <- samples[condition == "control"]
  samples[, gsm := as.character(gsm)]
  samples[, subject_id := as.character(subject_id)]
  samples[, clock_time_h := as.numeric(clock_time_h)]
  mcp <- fread(mcp_path)
  mcp[, gsm := as.character(gsm)]
  merged <- merge(samples[, .(gsm, subject_id, clock_time_h)], mcp, by = "gsm")
  merged
}

say("load composition")
blocks <- list(
  GSE39445 = load_development(
    "GSE39445",
    file.path(root, "data", "processed", "GSE39445_samples.tsv"),
    file.path(root, "results", "stage2", "h2_mcp_scores_GSE39445.tsv")
  ),
  GSE56931 = load_development(
    "GSE56931",
    file.path(root, "data", "processed", "GSE56931_samples.tsv"),
    file.path(root, "results", "stage2", "h2_mcp_scores_GSE56931.tsv")
  ),
  GSE113883 = load_development(
    "GSE113883",
    file.path(root, "results", "stage3", "GSE113883_samples.tsv"),
    file.path(root, "results", "stage4", "rnaseq", "h2_mcp_scores_GSE113883.tsv")
  )
)

say("L1 MCPcounter")
l1_samples <- fread(file.path(root, "results", "confirm", "L1", "samples_scoring.tsv"))
l1_samples <- l1_samples[period == "D1" & analysis_set == "day1"]
l1_samples[, gsm := as.character(gsm)]
l1_samples[, subject_id := as.character(subject_id)]
l1_samples[, clock_time_h := as.numeric(clock_time_h)]
l1_log <- read_gene_matrix(file.path(root, "results", "confirm", "L1", "log2_genes.tsv.gz"))
missing <- setdiff(l1_samples$gsm, colnames(l1_log))
if (length(missing)) stop(paste("L1 expression missing", length(missing), "samples"))
est <- mcp_scores(l1_log[, l1_samples$gsm, drop = FALSE])
rm(l1_log)
gc()
est$gsm <- rownames(est)
l1 <- merge(l1_samples[, .(gsm, subject_id, clock_time_h)], as.data.table(est), by = "gsm")
blocks$GSE48113 <- l1

expected <- c(GSE39445 = 221, GSE56931 = 130, GSE113883 = 153, GSE48113 = 147)
for (dataset in names(expected)) {
  n <- nrow(blocks[[dataset]])
  if (n != expected[[dataset]]) stop(paste(dataset, "n_samples", n, "expected", expected[[dataset]]))
}

set.seed(20261019)
rows <- list()
for (dataset in names(expected)) {
  block <- blocks[[dataset]]
  keep <- rep(TRUE, nrow(block))
  for (population in blood_labels) keep <- keep & is.finite(block[[population]])
  n_dropped <- nrow(block) - sum(keep)
  used <- block[keep]
  say(dataset, "n", nrow(used), "dropped", n_dropped)
  for (population in blood_labels) {
    say(dataset, population)
    dat <- data.table(
      subject_id = used$subject_id,
      clock_time_h = used$clock_time_h,
      score = used[[population]]
    )
    dat <- add_time_terms(dat)
    point <- fit_cosinor(dat, forms$full, forms$null)
    peaks <- boot_peaks(dat, 2000L)
    interval <- peak_interval(peaks)
    finite <- peaks[is.finite(peaks)]
    pct <- if (length(finite)) stats::quantile(finite, c(0.025, 0.975), names = FALSE, type = 7) else c(NA_real_, NA_real_)
    rows[[paste(dataset, population)]] <- data.table(
      dataset = dataset,
      population = population,
      n_samples_before_drop = nrow(block),
      n_dropped = n_dropped,
      n_samples = point$n_samples,
      n_subjects = point$n_subjects,
      two_a = point$two_a,
      peak_h = point$peak_h,
      lrt_p = point$lrt_p,
      lrt_stat = point$lrt_stat,
      reml_converged = point$reml_converged,
      ml_converged = point$ml_converged,
      n_boot = 2000L,
      n_boot_finite = length(finite),
      peak_center = unname(interval["center"]),
      peak_halfwidth_95 = unname(interval["halfwidth_95"]),
      peak_p025 = pct[1],
      peak_p975 = pct[2],
      seed = 20261019L,
      bootstrap = "participant resample with replacement; peak_interval circular mean and 95% circular distance",
      bootstrap_control = "lmerControl(calc.derivs=FALSE, check.conv.singular='ignore'); convergence uses conv$opt",
      point_fit = "fit_cosinor model M, default lmer control",
      title = "post hoc"
    )
  }
}
fwrite(rbindlist(rows), file.path(out_dir, "fig4b_composition_peaks.tsv"), sep = "\t")
say("fig4b finished")
