#!/usr/bin/env Rscript
# Part F. Gate candidates G1-G3 on development control samples.
# Criteria are not changed after seeing results.

suppressPackageStartupMessages(library(data.table))

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
root <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(root, "scripts", "lib_cosinor.R"))
out_dir <- file.path(root, "results", "stage3", "gates")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clocks <- c("RA_GTEx", "RA_Peters", "RA_default")
sig_of <- c(RA_GTEx = "GTExAge", RA_Peters = "Peters", RA_default = "DESeq2")
forms <- main_formula()
ns <- asNamespace("RNAAgeCalc")
genelist <- get("genelist_all", ns)$blood

read_expr <- function(path) {
  expr <- fread(path)
  mat <- as.matrix(expr[, -1])
  rownames(mat) <- expr[[1]]
  storage.mode(mat) <- "double"
  mat
}

lin_named <- function(mat, clock) {
  coefs <- genelist[[sig_of[[clock]]]]
  intercept <- unname(coefs[["(Intercept)"]])
  genes <- setdiff(names(coefs), "(Intercept)")
  have <- intersect(genes, rownames(mat))
  used <- coefs[have]
  scores <- as.numeric(intercept + drop(crossprod(used, mat[have, , drop = FALSE])))
  names(scores) <- colnames(mat)
  scores
}

gene_named <- function(mat, gene) {
  symbol <- gene
  if (!gene %in% rownames(mat) && gene == "ARNTL" && "BMAL1" %in% rownames(mat)) symbol <- "BMAL1"
  if (!symbol %in% rownames(mat)) return(NULL)
  vals <- as.numeric(mat[symbol, ])
  names(vals) <- colnames(mat)
  list(values = vals, matrix_symbol = symbol, n_finite = sum(is.finite(vals)))
}

g1_rows <- list()
g2_rows <- list()
g2_clock_rows <- list()

run_g1 <- function(dataset, dat, values) {
  peaks <- list()
  fits <- list()
  for (gene in c("PER1", "PER3", "ARNTL")) {
    got <- values[[gene]]
    if (is.null(got)) {
      peaks[[gene]] <- NA_real_
      g1_rows[[paste(dataset, gene, "missing")]] <<- data.table(
        dataset = dataset, gene = gene, matrix_symbol = NA_character_, present = FALSE,
        n_finite = 0, peak_h = NA_real_, lrt_p = NA_real_, two_a = NA_real_,
        n_samples = 0, n_subjects = 0, reml_converged = FALSE
      )
      next
    }
    d <- dat
    d$score <- unname(got$values[d$gsm])
    d <- d[is.finite(d$score), ]
    if (nrow(d) == 0L) {
      peaks[[gene]] <- NA_real_
      g1_rows[[paste(dataset, gene)]] <<- data.table(
        dataset = dataset, gene = gene, matrix_symbol = got$matrix_symbol, present = TRUE,
        n_finite = 0, peak_h = NA_real_, lrt_p = NA_real_, two_a = NA_real_,
        n_samples = 0, n_subjects = 0, reml_converged = FALSE
      )
      next
    }
    d <- add_time_terms(d)
    fit <- fit_cosinor(d, forms$full, forms$null)
    peaks[[gene]] <- fit$peak_h
    fits[[gene]] <- fit
    g1_rows[[paste(dataset, gene)]] <<- data.table(
      dataset = dataset, gene = gene, matrix_symbol = got$matrix_symbol, present = TRUE,
      n_finite = sum(is.finite(got$values[dat$gsm])),
      peak_h = fit$peak_h, lrt_p = fit$lrt_p, two_a = fit$two_a,
      n_samples = fit$n_samples, n_subjects = fit$n_subjects,
      reml_converged = fit$reml_converged, ml_converged = fit$ml_converged
    )
  }
  dist <- if (isTRUE(is.finite(peaks[["PER3"]])) && isTRUE(is.finite(peaks[["ARNTL"]]))) {
    circ_dist_0_12(peaks[["PER3"]], peaks[["ARNTL"]])
  } else {
    NA_real_
  }
  p1 <- if (is.null(fits[["PER1"]])) NA_real_ else fits[["PER1"]]$lrt_p
  p3 <- if (is.null(fits[["PER3"]])) NA_real_ else fits[["PER3"]]$lrt_p
  pass <- isTRUE(is.finite(p1) && p1 < 0.01 && is.finite(p3) && p3 < 0.01 && is.finite(dist) && dist >= 6)
  data.table(
    dataset = dataset,
    gate = "G1",
    per1_lrt_p = if (is.null(p1)) NA_real_ else p1,
    per3_lrt_p = if (is.null(p3)) NA_real_ else p3,
    per3_arntl_circ_dist_h = dist,
    pass = pass
  )
}

window_means <- function(dat, values) {
  morning <- dat$clock_time_h >= 7 & dat$clock_time_h <= 9
  evening <- dat$clock_time_h >= 19 & dat$clock_time_h <= 21
  rows <- list()
  subjects <- sort(unique(dat$subject_id))
  for (sid in subjects) {
    m <- dat$subject_id == sid & morning & is.finite(values)
    e <- dat$subject_id == sid & evening & is.finite(values)
    if (!any(m) || !any(e)) next
    rows[[sid]] <- data.table(
      subject_id = sid,
      morning = mean(values[m]),
      evening = mean(values[e]),
      n_morning = sum(m),
      n_evening = sum(e)
    )
  }
  if (length(rows) == 0L) return(data.table())
  rbindlist(rows)
}

paired_summary <- function(wm) {
  if (nrow(wm) < 2L) {
    return(list(n = nrow(wm), mean_diff = if (nrow(wm) == 1L) wm$morning - wm$evening else NA_real_,
                sd_diff = NA_real_, p = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }
  diff <- wm$morning - wm$evening
  tt <- t.test(wm$morning, wm$evening, paired = TRUE)
  list(
    n = nrow(wm),
    mean_diff = unname(tt$estimate),
    sd_diff = stats::sd(diff),
    p = unname(tt$p.value),
    ci_low = unname(tt$conf.int[1]),
    ci_high = unname(tt$conf.int[2])
  )
}

analyze_dataset <- function(dataset, sample_path, expr_path, control_only, score_clocks) {
  cat(dataset, "\n")
  flush.console()
  samples <- fread(sample_path)
  if (control_only) {
    dat <- samples[which(samples$condition == "control" & is.finite(samples$clock_time_h))]
  } else {
    dat <- samples[which(is.finite(samples$clock_time_h))]
  }
  dat <- as.data.frame(dat[, .(gsm, subject_id, clock_time_h, hours_awake)])
  mat <- read_expr(expr_path)
  values <- list()
  for (gene in c("PER1", "PER3", "ARNTL")) values[[gene]] <- gene_named(mat, gene)
  decision <- run_g1(dataset, dat, values)
  for (gene in c("PER1", "PER3")) {
    got <- values[[gene]]
    if (is.null(got)) {
      g2_rows[[paste(dataset, gene)]] <<- data.table(
        dataset = dataset, outcome = gene, outcome_type = "gene", n_subjects = 0,
        mean_diff = NA_real_, sd_diff = NA_real_, p = NA_real_, ci95_low = NA_real_, ci95_high = NA_real_,
        pass_gene = FALSE
      )
      next
    }
    wm <- window_means(dat, unname(got$values[dat$gsm]))
    sm <- paired_summary(wm)
    pass_gene <- isTRUE(is.finite(sm$mean_diff) && sm$mean_diff > 0 && is.finite(sm$p) && sm$p < 0.01)
    g2_rows[[paste(dataset, gene)]] <<- data.table(
      dataset = dataset, outcome = gene, matrix_symbol = got$matrix_symbol, outcome_type = "gene",
      n_subjects = sm$n, mean_diff = sm$mean_diff, sd_diff = sm$sd_diff, p = sm$p,
      ci95_low = sm$ci_low, ci95_high = sm$ci_high, pass_gene = pass_gene
    )
  }
  if (score_clocks) {
    for (clock in clocks) {
      named <- lin_named(mat, clock)
      wm <- window_means(dat, unname(named[dat$gsm]))
      sm <- paired_summary(wm)
      g2_clock_rows[[paste(dataset, clock)]] <<- data.table(
        dataset = dataset, outcome = clock, outcome_type = "clock_lin",
        n_subjects = sm$n, mean_diff = sm$mean_diff, sd_diff = sm$sd_diff,
        ci95_low = sm$ci_low, ci95_high = sm$ci_high, p = sm$p
      )
    }
  } else {
    g2_clock_rows[[dataset]] <<- data.table(
      dataset = dataset, outcome = NA_character_, outcome_type = "clock_lin",
      n_subjects = NA_integer_, mean_diff = NA_real_, sd_diff = NA_real_,
      ci95_low = NA_real_, ci95_high = NA_real_, p = NA_real_,
      status = "not_run",
      reason = "positive control failed; lin clock scores were not produced in section 3.2"
    )
  }
  rm(mat)
  gc()
  decision
}

decisions <- list()
decisions[["GSE39445"]] <- analyze_dataset(
  "GSE39445",
  file.path(root, "data", "processed", "GSE39445_samples.tsv"),
  file.path(root, "data", "processed", "GSE39445_expr_gene_log2.tsv.gz"),
  TRUE, TRUE
)
decisions[["GSE56931"]] <- analyze_dataset(
  "GSE56931",
  file.path(root, "data", "processed", "GSE56931_samples.tsv"),
  file.path(root, "data", "processed", "GSE56931_expr_gene_log2.tsv.gz"),
  TRUE, TRUE
)
decisions[["GSE113883"]] <- analyze_dataset(
  "GSE113883",
  file.path(root, "results", "stage3", "GSE113883_samples.tsv"),
  file.path(root, "results", "stage3", "expr", "GSE113883_log2tpm1.tsv.gz"),
  TRUE, FALSE
)
decisions[["GSE196728"]] <- analyze_dataset(
  "GSE196728",
  file.path(root, "results", "stage3", "GSE196728_samples.tsv"),
  file.path(root, "results", "stage3", "expr", "GSE196728_log2fpkm1.tsv.gz"),
  TRUE, FALSE
)

g1 <- rbindlist(decisions)
g2 <- rbindlist(g2_rows, fill = TRUE)
pass_map <- g2[, .(n_pass = sum(pass_gene), n_gene = .N), by = dataset]
pass_map[, g2_pass := n_pass == 2L]

cat("G3\n")
flush.console()
samples <- fread(file.path(root, "data", "processed", "GSE39445_samples.tsv"))
dat <- samples[which(samples$condition == "control" & is.finite(samples$hours_awake))]
dat <- as.data.frame(dat[, .(gsm, subject_id, hours_awake)])
mat <- read_expr(file.path(root, "data", "processed", "GSE39445_expr_gene_log2.tsv.gz"))
g3_rows <- list()
for (gene in c("PER1", "PER3")) {
  got <- gene_named(mat, gene)
  d <- dat
  d$clock_time_h <- d$hours_awake
  d$score <- unname(got$values[d$gsm])
  d <- d[is.finite(d$score), ]
  d <- add_time_terms(d)
  fit <- fit_cosinor(d, forms$full, forms$null)
  g3_rows[[gene]] <- data.table(
    dataset = "GSE39445",
    gate = "G3",
    gene = gene,
    matrix_symbol = got$matrix_symbol,
    time = "hours_awake",
    n_subjects = fit$n_subjects,
    n_samples = fit$n_samples,
    peak_h = fit$peak_h,
    two_a = fit$two_a,
    lrt_p = fit$lrt_p,
    reml_converged = fit$reml_converged
  )
}

fwrite(rbindlist(g1_rows, fill = TRUE), file.path(out_dir, "g1_gene_fits.tsv"), sep = "\t")
fwrite(g1, file.path(out_dir, "g1_decision.tsv"), sep = "\t")
fwrite(g2, file.path(out_dir, "g2_genes.tsv"), sep = "\t")
fwrite(pass_map, file.path(out_dir, "g2_decision.tsv"), sep = "\t")
fwrite(rbindlist(g2_clock_rows, fill = TRUE), file.path(out_dir, "g2_clocks.tsv"), sep = "\t")
fwrite(rbindlist(g3_rows), file.path(out_dir, "g3_peaks.tsv"), sep = "\t")
cat("done\n")
