#!/usr/bin/env Rscript
# Part E. Preregistration calculations. Locked expression is not opened.

suppressPackageStartupMessages({
  library(lme4)
  library(data.table)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
root <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(root, "scripts", "lib_cosinor.R"))
out_dir <- file.path(root, "results", "stage5")
forms <- main_formula()
datasets <- c("GSE39445", "GSE56931", "GSE113883")
ra_clocks <- c("RA_GTEx", "RA_Peters", "RA_default")
tams_clocks <- c("TAms_chrono", "TAms_mort")
summary_clocks <- c(ra_clocks, tams_clocks)
n_rep <- 2000
seed <- 20261015

read_samples <- function(dataset) {
  path <- if (dataset == "GSE113883") {
    file.path(root, "results", "stage3", "GSE113883_samples.tsv")
  } else {
    file.path(root, "data", "processed", paste0(dataset, "_samples.tsv"))
  }
  samples <- fread(path)
  keep <- which(samples$condition == "control" & is.finite(samples$clock_time_h))
  samples[keep, .(gsm = as.character(gsm), subject_id = as.character(subject_id), clock_time_h)]
}

variance_terms <- function(data) {
  reml <- quiet_lmer(forms$full, data, TRUE)
  null_reml <- quiet_lmer(forms$null, data, TRUE)
  sigma_reml <- if (!inherits(reml, "error") && optimizer_converged(reml)) sigma(reml) else NA_real_
  sigma_u <- NA_real_
  if (!inherits(null_reml, "error") && optimizer_converged(null_reml)) {
    vc <- lme4::VarCorr(null_reml)
    sds <- if (length(vc) == 1L) attr(vc[[1]], "stddev") else numeric()
    if (length(sds) == 1L) sigma_u <- as.numeric(sds)
  }
  c(sigma_reml = sigma_reml, sigma_u = sigma_u)
}

fit_named <- function(samples, named) {
  dat <- as.data.frame(samples)
  dat$score <- unname(named[dat$gsm])
  if (anyNA(dat$score)) stop("missing RA score")
  dat <- add_time_terms(dat)
  fit <- fit_cosinor(dat, forms$full, forms$null)
  extra <- variance_terms(dat)
  fit$sigma_reml <- unname(extra["sigma_reml"])
  fit$sigma_u <- unname(extra["sigma_u"])
  fit$two_a_over_sigma <- fit$two_a / fit$sigma_reml
  fit
}

ra_rows <- list()
for (dataset in datasets) {
  samples <- read_samples(dataset)
  path <- if (dataset == "GSE113883") {
    file.path(root, "results", "stage4", "rnaseq", "clock_scores.tsv")
  } else {
    file.path(root, "results", "stage1", "clocks", paste0(dataset, "_rnaage.tsv"))
  }
  stored_scores <- fread(path)
  stored_scores[, gsm := as.character(gsm)]
  stored_fit <- if (dataset == "GSE113883") {
    fread(file.path(root, "results", "stage4", "rnaseq", "h1_point_estimates.tsv"))
  } else {
    fread(file.path(root, "results", "stage2", "h1_point_estimates.tsv"))
  }
  for (clock in ra_clocks) {
    idx <- which(stored_scores$clock == clock & stored_scores$method == "lin")
    sc <- stored_scores[idx]
    named <- sc$value
    names(named) <- sc$gsm
    fit <- fit_named(samples, named)
    old_idx <- which(stored_fit$dataset == dataset & stored_fit$clock == clock & stored_fit$method == "lin")
    old <- stored_fit[old_idx]
    if (nrow(old) != 1L) stop(paste("stored fit missing", dataset, clock))
    ra_rows[[paste(dataset, clock)]] <- data.table(
      dataset = dataset, clock = clock,
      lrt_p = fit$lrt_p, peak_h = fit$peak_h, two_a = fit$two_a,
      sigma_reml = fit$sigma_reml, two_a_over_sigma = fit$two_a_over_sigma,
      sigma_u = fit$sigma_u,
      stored_lrt_p = old$lrt_p, stored_peak_h = old$peak_h, stored_two_a = old$two_a,
      abs_diff_lrt_p = abs(fit$lrt_p - old$lrt_p),
      abs_diff_peak_h = circ_dist_0_12(fit$peak_h, old$peak_h),
      abs_diff_two_a = abs(fit$two_a - old$two_a)
    )
  }
}
ra <- rbindlist(ra_rows)
if (max(ra$abs_diff_lrt_p) > 1e-8 || max(ra$abs_diff_two_a) > 1e-8) {
  stop("refit of stored lin scores does not match the stored H1 estimates")
}
fwrite(ra, file.path(out_dir, "ra_sigma_refit.tsv"), sep = "\t")

points <- fread(file.path(out_dir, "h1_point_estimates.tsv"))
tams <- points[clock %in% tams_clocks, .(
  dataset, clock, lrt_p, peak_h, two_a, sigma_reml, two_a_over_sigma, sigma_u
)]
clock_long <- rbindlist(list(
  ra[, .(dataset, clock, lrt_p, peak_h, two_a, sigma_reml, two_a_over_sigma, sigma_u)],
  tams
), use.names = TRUE)
clock_summ <- list()
for (clock in summary_clocks) {
  idx <- which(clock_long$clock == clock)
  part <- clock_long[idx]
  part <- part[match(datasets, part$dataset)]
  phi <- circ_mean_h(part$peak_h)
  dist <- circ_dist_0_12(part$peak_h, phi)
  clock_summ[[clock]] <- data.table(
    clock = clock, phi_bar = phi, max_circ_dist = max(dist),
    dist_GSE39445 = dist[1], dist_GSE56931 = dist[2], dist_GSE113883 = dist[3]
  )
}
fwrite(clock_long, file.path(out_dir, "clock_summary_long.tsv"), sep = "\t")
fwrite(rbindlist(clock_summ), file.path(out_dir, "clock_summary_phase.tsv"), sep = "\t")

module_ids <- unique(points$clock[points$family %in% c("TM", "CMP")])
rhythm_rows <- list()
for (clock in module_ids) {
  idx <- which(points$clock == clock)
  part <- points[idx]
  part <- part[match(datasets, part$dataset)]
  phi <- circ_mean_h(part$peak_h)
  dist <- circ_dist_0_12(part$peak_h, phi)
  q <- part$q_bh
  rhythm_rows[[clock]] <- data.table(
    clock = clock,
    family = part$family[1],
    outcome = part$outcome[1],
    module = part$module[1],
    parent_clock = part$parent_clock[1],
    phi_bar = phi,
    max_circ_dist = max(dist),
    q_GSE39445 = q[1], q_GSE56931 = q[2], q_GSE113883 = q[3],
    dist_GSE39445 = dist[1], dist_GSE56931 = dist[2], dist_GSE113883 = dist[3],
    two_a_over_sigma_GSE39445 = part$two_a_over_sigma[1],
    two_a_over_sigma_GSE56931 = part$two_a_over_sigma[2],
    two_a_over_sigma_GSE113883 = part$two_a_over_sigma[3],
    mean_two_a_over_sigma = mean(part$two_a_over_sigma),
    q_all_below_0.05 = all(q < 0.05),
    peaks_within_3h = all(dist <= 3),
    reproducible = all(q < 0.05) && all(dist <= 3)
  )
}
rhythm <- rbindlist(rhythm_rows)
fwrite(rhythm, file.path(out_dir, "reproducible_modules.tsv"), sep = "\t")

g2 <- fread(file.path(root, "results", "stage4", "rnaseq", "g2_clocks.tsv"))
g2 <- g2[, .(clock = outcome, g2_mean_diff = mean_diff, g2_sd_diff = sd_diff, g2_observed_d = mean_diff / sd_diff)]
g113 <- fread(file.path(root, "results", "stage3", "GSE113883_samples.tsv"))
g113_keep <- which(g113$condition == "control" & is.finite(g113$clock_time_h))
g113_dat <- as.data.frame(g113[g113_keep, .(gsm = as.character(gsm), subject_id = as.character(subject_id), clock_time_h)])
tams_scores <- fread(file.path(out_dir, "clock_scores.tsv"))
tams_scores[, gsm := as.character(gsm)]
tams_g2 <- list()
for (clock in tams_clocks) {
  score_idx <- which(tams_scores$dataset == "GSE113883" & tams_scores$clock == clock)
  sc <- tams_scores[score_idx]
  named <- sc$value
  names(named) <- sc$gsm
  values <- unname(named[g113_dat$gsm])
  morning <- g113_dat$clock_time_h >= 7 & g113_dat$clock_time_h <= 9
  evening <- g113_dat$clock_time_h >= 19 & g113_dat$clock_time_h <= 21
  rows <- list()
  for (sid in sort(unique(g113_dat$subject_id))) {
    m <- g113_dat$subject_id == sid & morning & is.finite(values)
    e <- g113_dat$subject_id == sid & evening & is.finite(values)
    if (!any(m) || !any(e)) next
    rows[[sid]] <- data.table(morning = mean(values[m]), evening = mean(values[e]))
  }
  wm <- rbindlist(rows)
  diff <- wm$morning - wm$evening
  tt <- stats::t.test(wm$morning, wm$evening, paired = TRUE)
  tams_g2[[clock]] <- data.table(
    clock = clock, n_subjects = nrow(wm),
    g2_mean_diff = unname(tt$estimate), g2_sd_diff = stats::sd(diff),
    g2_observed_d = unname(tt$estimate) / stats::sd(diff)
  )
}
fwrite(rbindlist(tams_g2), file.path(out_dir, "g2_tams.tsv"), sep = "\t")
g2 <- rbindlist(list(g2, rbindlist(tams_g2)[, .(clock, g2_mean_diff, g2_sd_diff, g2_observed_d)]), use.names = TRUE)
targets <- c(summary_clocks, rhythm$clock[rhythm$reproducible])
d_rows <- list()
for (clock in targets) {
  idx <- which(clock_long$clock == clock)
  if (length(idx) == 0L) {
    idx <- which(points$clock == clock)
    part <- points[idx, .(dataset, clock, two_a, peak_h, sigma_reml)]
  } else {
    part <- clock_long[idx, .(dataset, clock, two_a, peak_h, sigma_reml)]
  }
  part <- part[match(datasets, part$dataset)]
  expected <- part$two_a * cos(2 * pi * (8 - part$peak_h) / 24)
  d_exp <- expected / (sqrt(2) * part$sigma_reml)
  d_rows[[clock]] <- data.table(
    clock = clock,
    d_GSE39445 = d_exp[1], d_GSE56931 = d_exp[2], d_GSE113883 = d_exp[3],
    d_exp_median = median(d_exp)
  )
}
d_exp <- rbindlist(d_rows)
g2_idx <- match(d_exp$clock, g2$clock)
d_exp[, g2_mean_diff := g2$g2_mean_diff[g2_idx]]
d_exp[, g2_sd_diff := g2$g2_sd_diff[g2_idx]]
d_exp[, g2_observed_d := g2$g2_observed_d[g2_idx]]
paired_t_power <- function(d, n = 58, alpha) {
  df <- n - 1
  crit <- stats::qt(1 - alpha / 2, df)
  ncp <- d * sqrt(n)
  stats::pt(-crit, df, ncp, lower.tail = TRUE) + (1 - stats::pt(crit, df, ncp, lower.tail = TRUE))
}
d_exp[, power_alpha_0.05 := paired_t_power(d_exp_median, 58, 0.05)]
d_exp[, power_alpha_0.01 := paired_t_power(d_exp_median, 58, 0.01)]
d_exp[, n_pairs := 58L]
fwrite(d_exp, file.path(out_dir, "gse295758_expected_d.tsv"), sep = "\t")

meta <- fread(file.path(root, "results", "stage3", "GSE295758_samples.tsv"))
counts <- meta[, .N, by = subject_id]
four <- counts$subject_id[counts$N == 4]
four_idx <- which(meta$subject_id %in% four)
verbatim <- meta[four_idx, .(
  gsm, subject_id, sample_title, raw_characteristics, age, sex, char_opioids, char_timepoint
)]
fwrite(verbatim, file.path(out_dir, "gse295758_four_sample_subjects.tsv"), sep = "\t")
subject_tbl <- unique(meta[, .(
  subject_id, age, sex, opioid_group = char_opioids
)])
subject_tbl <- merge(subject_tbl, counts, by = "subject_id", sort = FALSE)
setnames(subject_tbl, "N", "n_samples")
fwrite(subject_tbl, file.path(out_dir, "gse295758_subject_covariates.tsv"), sep = "\t")

design_48113 <- function() {
  samples <- fread(file.path(root, "data", "processed", "GSE48113_samples.tsv"))
  keep <- which(
    samples$condition == "control" & samples$period == "D1" &
      !(samples$gsm %in% c("GSM1168837", "GSM1168735")) &
      is.finite(samples$clock_time_h)
  )
  kept <- samples[keep, .(subject_id = as.character(subject_id), clock_time_h)]
  kept
}

design_122725 <- function() {
  samples <- fread(file.path(root, "data", "processed", "GSE122725_samples.tsv"))
  wake <- fread(file.path(root, "results", "stage2", "GSE122725_wake_ref.tsv"))
  keep <- which(samples$condition == "control" & samples$period == "period1")
  kept <- samples[keep, .(gsm = as.character(gsm), subject_id = as.character(subject_id))]
  wake[, gsm := as.character(gsm)]
  merged <- merge(kept, wake[, .(gsm, wake_ref_h)], by = "gsm", sort = FALSE)
  if (nrow(merged) != nrow(kept) || any(!is.finite(merged$wake_ref_h))) stop("wake time missing")
  merged[, .(subject_id, clock_time_h = wake_ref_h)]
}

designs <- list(GSE48113 = design_48113(), GSE122725 = design_122725())
params <- clock_long[clock %in% summary_clocks]
set.seed(seed)
sim_rows <- list()
for (design_name in names(designs)) {
  design <- as.data.frame(designs[[design_name]])
  design$subject_id <- factor(design$subject_id)
  n_subject <- nlevels(design$subject_id)
  for (clock in summary_clocks) {
    for (source in datasets) {
      hit <- which(params$clock == clock & params$dataset == source)
      par <- params[hit]
      cat("sim", design_name, clock, source, "\n")
      flush.console()
      p_lt_05 <- 0L
      p_lt_01 <- 0L
      n_ok <- 0L
      amp <- par$two_a / 2
      for (b in seq_len(n_rep)) {
        u <- stats::rnorm(n_subject, 0, par$sigma_u)
        names(u) <- levels(design$subject_id)
        dat <- design
        dat$score <- unname(u[dat$subject_id]) + amp * cos(2 * pi * (dat$clock_time_h - par$peak_h) / 24) +
          stats::rnorm(nrow(dat), 0, par$sigma_reml)
        dat <- add_time_terms(dat)
        ml <- quiet_lmer(forms$full, dat, FALSE)
        null <- quiet_lmer(forms$null, dat, FALSE)
        ok <- !inherits(ml, "error") && !inherits(null, "error") &&
          optimizer_converged(ml) && optimizer_converged(null)
        if (ok) {
          n_ok <- n_ok + 1L
          stat <- as.numeric(2 * (logLik(ml) - logLik(null)))
          p <- stats::pchisq(stat, df = 2, lower.tail = FALSE)
          if (p < 0.05) p_lt_05 <- p_lt_05 + 1L
          if (p < 0.01) p_lt_01 <- p_lt_01 + 1L
        }
      }
      one_row <- data.table(
        design = design_name, clock = clock, param_dataset = source,
        n_subjects = n_subject, n_samples = nrow(design),
        A = amp, phi = par$peak_h, sigma_u = par$sigma_u, sigma = par$sigma_reml,
        n_rep = n_rep, n_converged = n_ok,
        frac_p_lt_0.05 = p_lt_05 / n_rep, frac_p_lt_0.01 = p_lt_01 / n_rep,
        seed = seed
      )
      sim_rows[[paste(design_name, clock, source)]] <- one_row
      fwrite(rbindlist(sim_rows), file.path(out_dir, "locked_design_power.tsv"), sep = "\t")
    }
  }
}
cat("Part E done\n")
