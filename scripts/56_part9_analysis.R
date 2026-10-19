#!/usr/bin/env Rscript
# Stage 9 post-hoc analyses. Does not write to results/confirm/.

suppressPackageStartupMessages(library(data.table))

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
root <- normalizePath(file.path(script_dir, ".."))
source(file.path(root, "scripts", "confirm", "run_pipeline.R"))
source(file.path(root, "scripts", "confirm", "run_locked_analysis.R"))

out_dir <- file.path(root, "results", "stage9")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
title <- "post hoc"
clocks <- c(primary_clocks, secondary_clocks)
genes <- c("PER1", "PER2", "PER3", "NR1D1", "NR1D2", "ARNTL", "CRY1", "CRY2", "DBP", "TEF")
say <- function(...) message(format(Sys.time(), "%H:%M:%S"), " ", paste(...))

circ_sd_h <- function(hours) {
  ang <- hours * 2 * pi / 24
  rbar <- sqrt(mean(cos(ang))^2 + mean(sin(ang))^2)
  if (!is.finite(rbar) || rbar <= 0) return(NA_real_)
  sqrt(-2 * log(min(rbar, 1))) * 24 / (2 * pi)
}
signed_diff <- function(a, b) {
  d <- (a - b) %% 24
  ifelse(d > 12, d - 24, d)
}
fpr_t <- function(ncp, df) {
  crit <- qt(0.975, df)
  1 - pt(crit, df, ncp) + pt(-crit, df, ncp)
}
sigma_u_of <- function(dat) {
  null_reml <- quiet_lmer(forms$null, dat, TRUE)
  if (inherits(null_reml, "error") || !optimizer_converged(null_reml)) return(NA_real_)
  vc <- lme4::VarCorr(null_reml)
  sds <- if (length(vc) == 1L) attr(vc[[1]], "stddev") else numeric()
  if (length(sds) == 1L) as.numeric(sds) else NA_real_
}
fit_named <- function(samples, values) {
  used <- samples[is.finite(values) & is.finite(samples$clock_time_h)]
  score <- data.table(gsm = used$gsm, clock = "x", value = values[is.finite(values) & is.finite(samples$clock_time_h)])
  model_m(prepare_score(used, score, "x"))
}

say("Part A L2 wake_corr")
l2_samples <- fread(file.path(root, "results", "confirm", "L2", "samples_scoring.tsv"))
l2_samples[, gsm := as.character(gsm)]
period1 <- l2_samples[period == "period1"]
period1[, relclocktime := wake_ref_h]
period1[, wake_corr := (wake_ref_h - 8) %% 24]
l2_expr <- as.matrix(read_gene_matrix(file.path(root, "results", "confirm", "L2", "log2_genes.tsv")))

fit_symbol <- function(time_values, symbol) {
  one <- copy(period1)
  one[, clock_time_h := time_values]
  if (!(symbol %in% rownames(l2_expr))) {
    return(data.table(n_samples = 0L, n_subjects = 0L, peak_h = NA_real_, lrt_p = NA_real_, lrt_stat = NA_real_, status = "absent"))
  }
  values <- as.numeric(l2_expr[symbol, one$gsm])
  fit <- fit_named(one, values)
  data.table(
    n_samples = fit$n_samples, n_subjects = fit$n_subjects,
    peak_h = fit$peak_h, lrt_p = fit$lrt_p, lrt_stat = fit$lrt_stat, status = "fit"
  )
}

gate_rows <- list()
for (gene in c("PER1", "PER3")) {
  rel <- fit_symbol(period1$relclocktime %% 24, gene)
  corr <- fit_symbol(period1$wake_corr, gene)
  gate_rows[[gene]] <- data.table(
    gene = gene, n_samples = corr$n_samples, n_subjects = corr$n_subjects,
    peak_relclock = rel$peak_h, peak_wake_corr = corr$peak_h,
    in_prereg_window = is.finite(corr$peak_h) && peak_in_window(corr$peak_h, 18, 2),
    window_rule = "peak >= 18 or peak <= 2, the preregistered 18-26 h window modulo 24",
    lrt_p_relclock = rel$lrt_p, lrt_p_wake_corr = corr$lrt_p,
    lrt_stat_relclock = rel$lrt_stat, lrt_stat_wake_corr = corr$lrt_stat,
    lrt_p_difference = corr$lrt_p - rel$lrt_p,
    lrt_p_below_0.01 = is.finite(corr$lrt_p) && corr$lrt_p < 0.01,
    prereg_gate_pass = is.finite(corr$lrt_p) && corr$lrt_p < 0.01 && is.finite(corr$peak_h) && peak_in_window(corr$peak_h, 18, 2),
    title = title
  )
}
fwrite(rbindlist(gate_rows), file.path(out_dir, "l2_wake_corr_gate.tsv"), sep = "\t")

gse <- fread(file.path(root, "results", "stage8", "gene_peak_differences_arntl_alias.tsv"))
gene_rows <- list()
for (gene in genes) {
  symbol <- if (identical(gene, "ARNTL")) "BMAL1" else gene
  corr <- fit_symbol(period1$wake_corr, symbol)
  gene_id <- gene
  ref <- gse[gene == gene_id]
  gene_rows[[gene]] <- data.table(
    gene = gene, symbol_used = symbol, n_samples = corr$n_samples, n_subjects = corr$n_subjects,
    peak_wake_corr = corr$peak_h, lrt_p_wake_corr = corr$lrt_p, status = corr$status,
    peak_gse39445_hours_awake = if (nrow(ref) == 1L) ref$peak_gse39445 else NA_real_,
    signed_diff_h = signed_diff(corr$peak_h, if (nrow(ref) == 1L) ref$peak_gse39445 else NA_real_),
    title = title
  )
}
gene_tab <- rbindlist(gene_rows)
gene_tab[, circ_dist_0_12 := mapply(circ_dist_0_12, peak_wake_corr, peak_gse39445_hours_awake)]
fwrite(gene_tab, file.path(out_dir, "l2_wake_corr_genes.tsv"), sep = "\t")
finite_diff <- gene_tab[is.finite(signed_diff_h), signed_diff_h]
fwrite(data.table(
  n_genes = length(finite_diff),
  circ_mean_h = if (length(finite_diff)) circ_mean_h(finite_diff) else NA_real_,
  circ_sd_h = if (length(finite_diff)) circ_sd_h(finite_diff) else NA_real_,
  signed_diff_definition = "(peak_wake_corr - peak_gse39445_hours_awake) mod 24, then subtract 24 when the result is greater than 12; range (-12, 12]",
  circ_sd_definition = "sqrt(-2*log(R)) * 24 / (2*pi), R = mean resultant length of the signed differences",
  arntl_source = "L2 symbol BMAL1; GSE39445 symbol ARNTL, hours_awake axis, stage 8 table",
  title = title
), file.path(out_dir, "l2_wake_corr_gene_circular.tsv"), sep = "\t")

l2_scores <- rbind(
  fread(file.path(root, "results", "posthoc_stage7c", "l2_tage_scores.tsv")),
  ra_linear_scores(l2_expr[, period1$gsm, drop = FALSE])[, .(gsm, clock, value)],
  fill = TRUE
)
l2_scores[, gsm := as.character(gsm)]
clock_rows <- list()
for (clock_id in clocks) {
  rel_s <- copy(period1)
  rel_s[, clock_time_h := relclocktime %% 24]
  corr_s <- copy(period1)
  corr_s[, clock_time_h := wake_corr]
  one <- l2_scores[clock == clock_id]
  rel_fit <- model_m(prepare_score(rel_s, one, clock_id))
  corr_fit <- model_m(prepare_score(corr_s, one, clock_id))
  clock_rows[[clock_id]] <- data.table(
    clock = clock_id, n_subjects = corr_fit$n_subjects, n_samples = corr_fit$n_samples,
    peak_relclock = rel_fit$peak_h, peak_wake_corr = corr_fit$peak_h,
    lrt_p_relclock = rel_fit$lrt_p, lrt_p_wake_corr = corr_fit$lrt_p,
    lrt_p_difference = corr_fit$lrt_p - rel_fit$lrt_p,
    role = "reference", title = title
  )
}
fwrite(rbindlist(clock_rows), file.path(out_dir, "l2_wake_corr_clocks.tsv"), sep = "\t")
rm(l2_expr)
gc()

say("Part B parameters")
day1 <- fread(file.path(root, "results", "confirm", "L1", "samples_scoring.tsv"))
day1[, gsm := as.character(gsm)]
day1 <- day1[period == "D1"]
l1_scores <- rbind(
  fread(file.path(root, "results", "confirm", "L1", "tage_scores.tsv")),
  fread(file.path(root, "results", "confirm", "L1", "ra_scores.tsv")),
  fill = TRUE
)
l1_scores[, gsm := as.character(gsm)]
fam <- rbind(
  fread(file.path(root, "results", "confirm", "family_P.tsv")),
  fread(file.path(root, "results", "confirm", "family_S1.tsv")),
  fill = TRUE
)
f1 <- fread(file.path(root, "prereg", "frozen", "F1.tsv"))
pairs <- fread(file.path(root, "results", "confirm", "L3", "paired_deltas.tsv"))
param_rows <- list()
lab_fun <- list()
for (clock_id in clocks) {
  dat <- prepare_score(day1, l1_scores, clock_id)
  fit <- model_m(dat)
  sig_u <- sigma_u_of(dat)
  ref <- fam[dataset == "L1" & clock == clock_id]
  ref_l3 <- fam[dataset == "L3" & clock == clock_id]
  if (nrow(ref) != 1L || !is.finite(fit$two_a) || abs(fit$two_a - ref$two_a) / abs(ref$two_a) > 1e-6) {
    stop(paste("L1 two_a does not match family table", clock_id, fit$two_a, if (nrow(ref)) ref$two_a else NA))
  }
  if (!is.finite(fit$sigma_reml) || abs(fit$sigma_reml - ref$sigma_reml) / abs(ref$sigma_reml) > 1e-6) {
    stop(paste("L1 sigma does not match family table", clock_id))
  }
  phi_vals <- unique(f1[clock == clock_id, phi_bar])
  if (length(phi_vals) != 1L) stop(paste("phi_bar is not unique", clock_id))
  block <- pairs[clock == clock_id]
  morning_sd <- stats::sd(block$morning)
  if (nrow(ref_l3) != 1L || !is.finite(ref_l3$sd_delta) || abs(stats::sd(block$delta) - ref_l3$sd_delta) / abs(ref_l3$sd_delta) > 1e-6) {
    stop(paste("L3 sd_delta does not match", clock_id))
  }
  dcos <- cos(2 * pi * (8 - phi_vals) / 24) - cos(2 * pi * (20 - phi_vals) / 24)
  param_rows[[paste(clock_id, "lab")]] <- data.table(
    scenario = "lab", clock = clock_id, phi_bar = NA_real_, A = NA_real_, cosine_denominator = NA_real_,
    beta_cos = fit$beta_cos, beta_sin = fit$beta_sin, sigma = fit$sigma_reml, sigma_u = sig_u,
    total_sd = sqrt(fit$sigma_reml^2 + sig_u^2), sd_diff = sqrt(2) * fit$sigma_reml,
    mean_delta = NA_real_, ci_high = NA_real_, morning_sd = NA_real_,
    assumption = "f(t) = beta_cos*cos(2*pi*t/24) + beta_sin*sin(2*pi*t/24) from L1 Day 1 model M",
    unit = if (startsWith(clock_id, "RA_")) "year" else "score", title = title
  )
  for (scenario in c("free-living", "free-living upper")) {
    target <- if (scenario == "free-living") ref_l3$mean_delta else ref_l3$ci_high
    param_rows[[paste(clock_id, scenario)]] <- data.table(
      scenario = scenario, clock = clock_id, phi_bar = phi_vals, A = target / dcos, cosine_denominator = dcos,
      beta_cos = NA_real_, beta_sin = NA_real_, sigma = NA_real_, sigma_u = NA_real_,
      total_sd = morning_sd, sd_diff = ref_l3$sd_delta,
      mean_delta = ref_l3$mean_delta, ci_high = ref_l3$ci_high, morning_sd = morning_sd,
      assumption = "f(t) = A*cos(2*pi*(t-phi_bar)/24). phi_bar is F1. A makes f(8)-f(20) equal the L3 mean delta, or its 95 percent upper limit in the upper scenario. The cosine shape and phi_bar are assumptions.",
      unit = if (startsWith(clock_id, "RA_")) "year" else "score", title = title
    )
  }
  lab_fun[[clock_id]] <- c(beta_cos = fit$beta_cos, beta_sin = fit$beta_sin, sigma = fit$sigma_reml, sigma_u = sig_u)
}
params <- rbindlist(param_rows, fill = TRUE)
fwrite(params, file.path(out_dir, "sim_parameters.tsv"), sep = "\t")

rhythm <- function(scenario, clock_id, time) {
  scen <- scenario
  row <- params[scenario == scen & clock == clock_id]
  if (scenario == "lab") {
    par <- lab_fun[[clock_id]]
    unname(par["beta_cos"]) * cos(2 * pi * time / 24) + unname(par["beta_sin"]) * sin(2 * pi * time / 24)
  } else {
    row$A * cos(2 * pi * (time - row$phi_bar) / 24)
  }
}

say("Part B designs 1 and 2")
s_values <- c(0, 1, 2, 3, 4, 6, 8, 12)
n_paired <- c(20, 50, 100, 200)
n_group <- c(25, 50, 100, 200)
bias_rows <- list()
fpr_rows <- list()
small_rows <- list()
d2_rows <- list()
for (scenario in c("lab", "free-living", "free-living upper")) {
  for (clock_id in clocks) {
    scen <- scenario
    row <- params[scenario == scen & clock == clock_id]
    for (s in s_values) {
      b <- rhythm(scenario, clock_id, 8 + s) - rhythm(scenario, clock_id, 8)
      bias_rows[[paste(scenario, clock_id, s)]] <- data.table(
        scenario = scenario, clock = clock_id, s = s, b = b,
        b_over_total_sd = b / row$total_sd, total_sd = row$total_sd, unit = row$unit, title = title
      )
      for (n in n_paired) {
        ncp <- b / row$sd_diff * sqrt(n)
        fpr_rows[[paste(scenario, clock_id, s, n)]] <- data.table(
          scenario = scenario, clock = clock_id, s = s, n = n, df = n - 1, ncp = ncp,
          sd_diff = row$sd_diff, fpr = fpr_t(ncp, n - 1), alpha = 0.05, test = "paired t", title = title
        )
      }
      for (n in n_group) {
        ncp <- b / row$total_sd * sqrt(n / 2)
        d2_rows[[paste(scenario, clock_id, s, n)]] <- data.table(
          scenario = scenario, clock = clock_id, s = s, n_per_group = n, df = 2 * n - 2, ncp = ncp,
          total_sd = row$total_sd, b = b, fpr = fpr_t(ncp, 2 * n - 2),
          alpha = 0.05, test = "two-sample t, equal variance", title = title
        )
      }
    }
    fpr100 <- vapply(s_values, function(s) {
      b <- rhythm(scenario, clock_id, 8 + s) - rhythm(scenario, clock_id, 8)
      fpr_t(b / row$sd_diff * 10, 99)
    }, numeric(1))
    hit <- s_values[fpr100 > 0.10]
    small_rows[[paste(scenario, clock_id)]] <- data.table(
      scenario = scenario, clock = clock_id, n = 100, threshold = 0.10,
      smallest_s = if (length(hit)) min(hit) else NA_real_,
      status = if (length(hit)) "exceeds" else "none of the listed s",
      title = title
    )
  }
}
fwrite(rbindlist(bias_rows), file.path(out_dir, "sim_design1_bias.tsv"), sep = "\t")
fwrite(rbindlist(fpr_rows), file.path(out_dir, "sim_design1_fpr.tsv"), sep = "\t")
fwrite(rbindlist(small_rows), file.path(out_dir, "sim_design1_smallest_s.tsv"), sep = "\t")
fwrite(rbindlist(d2_rows), file.path(out_dir, "sim_design2_fpr.tsv"), sep = "\t")

say("Part B Monte Carlo design 1")
set.seed(20261019)
mc_rows <- list()
for (clock_id in clocks) {
  par <- lab_fun[[clock_id]]
  n <- 50
  s <- 3
  reject <- 0L
  failed <- 0L
  for (b in seq_len(2000L)) {
    u <- rnorm(n, 0, unname(par["sigma_u"]))
    e1 <- rnorm(n, 0, unname(par["sigma"]))
    e2 <- rnorm(n, 0, unname(par["sigma"]))
    d <- (rhythm("lab", clock_id, 8 + s) + u + e2) - (rhythm("lab", clock_id, 8) + u + e1)
    sd_d <- stats::sd(d)
    if (!is.finite(sd_d) || sd_d == 0) {
      failed <- failed + 1L
      next
    }
    tstat <- mean(d) / (sd_d / sqrt(n))
    p <- 2 * pt(-abs(tstat), n - 1)
    if (is.finite(p) && p < 0.05) reject <- reject + 1L
  }
  b_s <- rhythm("lab", clock_id, 11) - rhythm("lab", clock_id, 8)
  analytical <- fpr_t(b_s / (sqrt(2) * unname(par["sigma"])) * sqrt(n), n - 1)
  mc_rows[[clock_id]] <- data.table(
    scenario = "lab", clock = clock_id, n = n, s = s, n_rep = 2000L, n_reject = reject, n_failed = failed,
    mc_fpr = reject / 2000, analytical_fpr = analytical, seed = 20261019,
    generator = "model M: f(t) + random intercept + residual; paired difference", title = title
  )
}
fwrite(rbindlist(mc_rows), file.path(out_dir, "sim_design1_montecarlo.tsv"), sep = "\t")

say("Part B design 3")
set.seed(20261019)
d3_rows <- list()
for (W in c(0, 2, 4, 8)) {
  for (clock_id in clocks) {
    par <- lab_fun[[clock_id]]
    sigma <- unname(par["sigma"])
    delta <- 0.3 * sqrt(2) * sigma
    n <- 50
    estimates <- numeric(2000)
    reject <- 0L
    failed <- 0L
    for (b in seq_len(2000L)) {
      t1 <- runif(n, 8, 8 + W)
      t2 <- runif(n, 8, 8 + W)
      u <- rnorm(n, 0, unname(par["sigma_u"]))
      d <- rhythm("lab", clock_id, t2) + u + rnorm(n, 0, sigma) + delta - (rhythm("lab", clock_id, t1) + u + rnorm(n, 0, sigma))
      estimates[b] <- mean(d)
      sd_d <- stats::sd(d)
      if (!is.finite(sd_d) || sd_d == 0) {
        failed <- failed + 1L
        next
      }
      tstat <- mean(d) / (sd_d / sqrt(n))
      p <- 2 * pt(-abs(tstat), n - 1)
      if (is.finite(p) && p < 0.05) reject <- reject + 1L
    }
    d3_rows[[paste(W, clock_id)]] <- data.table(
      scenario = "lab", clock = clock_id, W = W, n = n, delta = delta, n_rep = 2000L,
      n_reject = reject, n_failed = failed, power = reject / 2000,
      mean_estimate = mean(estimates), bias = mean(estimates) - delta, seed = 20261019,
      time_draw = "independent Uniform(8, 8+W)", effect = "delta added to the post sample", title = title
    )
  }
}
fwrite(rbindlist(d3_rows), file.path(out_dir, "sim_design3_power.tsv"), sep = "\t")

say("Part C cosinor composition")
center_comp <- function(comp, subject) {
  out <- as.data.frame(comp)
  subject <- as.character(subject)
  for (j in seq_along(out)) out[[j]] <- out[[j]] - ave(out[[j]], subject, FUN = mean)
  out
}
fit_composition <- function(dat, comp) {
  unadj <- fit_cosinor(dat, forms$full, forms$null)
  safe <- safe_name(names(comp))
  names(comp) <- safe
  use <- cbind(dat, comp)
  full <- as.formula(paste("score ~ cos_t + sin_t +", paste(safe, collapse = " + "), "+ (1 | subject_id)"))
  null <- as.formula(paste("score ~", paste(safe, collapse = " + "), "+ (1 | subject_id)"))
  adj <- fit_cosinor(use, full, null)
  data.table(
    n_subjects = unadj$n_subjects, n_samples = unadj$n_samples,
    two_a_before = unadj$two_a, r2_within_before = unadj$r2_within,
    lrt_p_before = unadj$lrt_p, peak_before = unadj$peak_h,
    two_a_after = adj$two_a, r2_within_after = adj$r2_within,
    lrt_p_after = adj$lrt_p, peak_after = adj$peak_h,
    two_a_ratio = if (is.finite(unadj$two_a) && unadj$two_a != 0) adj$two_a / unadj$two_a else NA_real_,
    reml_converged_after = adj$reml_converged, ml_converged_after = adj$ml_converged,
    null_converged_after = adj$null_converged
  )
}
dataset_samples <- function(dataset) {
  if (dataset == "GSE113883") {
    samples <- fread(file.path(root, "results", "stage3", "GSE113883_samples.tsv"))
  } else if (dataset == "L1") {
    samples <- copy(day1)
  } else {
    samples <- fread(file.path(root, "data", "processed", paste0(dataset, "_samples.tsv")))
  }
  samples[, gsm := as.character(gsm)]
  if ("condition" %in% names(samples)) samples <- samples[condition == "control"]
  samples[is.finite(clock_time_h)]
}
score_table <- function(dataset) {
  if (dataset == "L1") return(l1_scores)
  if (dataset == "GSE113883") {
    tage <- fread(file.path(root, "results", "stage5", "clock_scores.tsv"))
    tage <- tage[dataset == "GSE113883" & clock %in% c("TAms_chrono", "TAms_mort"), .(gsm = as.character(gsm), clock, value)]
    ra <- fread(file.path(root, "results", "stage4", "rnaseq", "clock_scores.tsv"))
    ra <- ra[method == "lin", .(gsm = as.character(gsm), clock, value)]
    return(rbind(tage, ra, fill = TRUE))
  }
  ds <- dataset
  tage <- fread(file.path(root, "results", "stage5", "clock_scores.tsv"))
  tage <- tage[dataset == ds & condition == "control" & clock %in% c("TAms_chrono", "TAms_mort"), .(gsm = as.character(gsm), clock, value)]
  expr <- as.matrix(read_gene_matrix(file.path(root, "data", "processed", paste0(dataset, "_expr_gene_log2.tsv.gz"))))
  samples <- dataset_samples(dataset)
  ra <- ra_linear_scores(expr[, samples$gsm, drop = FALSE])[, .(gsm, clock, value)]
  rm(expr)
  gc()
  rbind(tage, ra, fill = TRUE)
}
read_mcp <- function(path) {
  tab <- fread(path)
  tab[, gsm := as.character(gsm)]
  missing <- setdiff(blood_labels, names(tab))
  if (length(missing)) stop(paste("MCP file missing", paste(missing, collapse = ", ")))
  tab
}

comp_rows <- list()
dev_specs <- list(
  GSE39445 = file.path(root, "results", "stage2", "h2_mcp_scores_GSE39445.tsv"),
  GSE56931 = file.path(root, "results", "stage2", "h2_mcp_scores_GSE56931.tsv"),
  GSE113883 = file.path(root, "results", "stage4", "rnaseq", "h2_mcp_scores_GSE113883.tsv")
)
for (dataset in names(dev_specs)) {
  say(dataset)
  samples <- dataset_samples(dataset)
  scores <- score_table(dataset)
  mcp <- read_mcp(dev_specs[[dataset]])
  for (clock_id in clocks) {
    dat <- prepare_score(samples, scores, clock_id)
    hit <- match(dat$gsm, mcp$gsm)
    if (anyNA(hit)) stop(paste("MCP gsm missing", dataset, clock_id))
    comp <- as.data.frame(mcp[hit, ..blood_labels])
    bad <- !apply(comp, 1, function(x) all(is.finite(x)))
    dat <- dat[!bad, , drop = FALSE]
    comp <- comp[!bad, , drop = FALSE]
    centered <- center_comp(comp, dat$subject_id)
    row <- fit_composition(dat, centered)
    row[, dataset := dataset]
    row[, clock := clock_id]
    row[, n_dropped := sum(bad)]
    row[, title := title]
    comp_rows[[paste(dataset, clock_id)]] <- row
    if (dataset != "L1") {
      ref <- fread(file.path(root, "results", "stage8", "cosinor_effects_boot.tsv"))
      ds <- dataset
      prev <- ref[dataset == ds & clock == clock_id]
      if (sum(bad) == 0L && nrow(prev) == 1L && is.finite(row$two_a_before) && abs(row$two_a_before - prev$two_a) / abs(prev$two_a) > 1e-4) {
        stop(paste("unadjusted two_a does not match stage 8", dataset, clock_id, row$two_a_before, prev$two_a))
      }
    }
  }
}

say("L1 composition")
l1_log <- as.matrix(read_gene_matrix(file.path(root, "results", "confirm", "L1", "log2_genes.tsv.gz")))
l1_mcp <- composition_matrix(l1_log, day1)
rm(l1_log)
gc()
for (clock_id in clocks) {
  dat <- prepare_score(day1, l1_scores, clock_id)
  hit <- match(dat$gsm, rownames(l1_mcp))
  if (anyNA(hit)) stop("L1 composition gsm missing")
  comp <- l1_mcp[hit, blood_labels, drop = FALSE]
  bad <- apply(comp, 1, function(x) any(!is.finite(x)))
  dat <- dat[!bad, , drop = FALSE]
  comp <- comp[!bad, , drop = FALSE]
  centered <- center_comp(comp, dat$subject_id)
  row <- fit_composition(dat, centered)
  row[, dataset := "L1"]
  row[, clock := clock_id]
  row[, n_dropped := sum(bad)]
  row[, title := title]
  comp_rows[[paste("L1", clock_id)]] <- row
}
fwrite(rbindlist(comp_rows, fill = TRUE), file.path(out_dir, "composition_cosinor.tsv"), sep = "\t")

say("L3 composition")
l3_samples <- fread(file.path(root, "results", "confirm", "L3", "samples_scoring.tsv"))
l3_samples[, gsm := as.character(gsm)]
l3_fpkm <- as.matrix(read_gene_matrix(file.path(root, "results", "confirm", "L3", "log2_fpkm.tsv")))
l3_mcp <- composition_matrix(l3_fpkm, l3_samples)
rm(l3_fpkm)
gc()
comp_delta <- paired_composition_delta(l3_samples, l3_mcp)
h5 <- fread(file.path(root, "results", "confirm", "h5.tsv"))
l3_rows <- list()
beta_rows <- list()
age_rows <- list()
cov <- fread(file.path(root, "results", "stage5", "gse295758_subject_covariates.tsv"))
for (clock_id in clocks) {
  block <- pairs[clock == clock_id]
  ord <- match(comp_delta$subjects, block$subject_id)
  if (anyNA(ord)) stop("L3 subjects do not match")
  delta <- block$delta[ord]
  morning <- block$morning[ord]
  evening <- block$evening[ord]
  subjects <- comp_delta$subjects
  x <- comp_delta$delta
  safe <- safe_name(names(x))
  names(x) <- safe
  frame <- data.frame(delta = delta, x, check.names = FALSE)
  fit <- stats::lm(delta ~ ., data = frame)
  coefs <- summary(fit)$coefficients
  ci <- tryCatch(stats::confint(fit)["(Intercept)", ], error = function(e) c(NA_real_, NA_real_))
  ref <- fam[dataset == "L3" & clock == clock_id]
  h5_row <- h5[clock == clock_id]
  l3_rows[[clock_id]] <- data.table(
    clock = clock_id, n = nrow(frame),
    mean_delta_before = ref$mean_delta, ci_low_before = ref$ci_low, ci_high_before = ref$ci_high, p_before = ref$p,
    mean_delta_composition = unname(stats::coef(fit)[["(Intercept)"]]),
    ci_low_composition = unname(ci[1]), ci_high_composition = unname(ci[2]),
    p_composition = unname(coefs["(Intercept)", "Pr(>|t|)"]),
    df_residual = fit$df.residual,
    mean_delta_time = if (nrow(h5_row) == 1L) h5_row$mean_delta_minus_c else NA_real_,
    ci_low_time = if (nrow(h5_row) == 1L) h5_row$ci_low else NA_real_,
    ci_high_time = if (nrow(h5_row) == 1L) h5_row$ci_high else NA_real_,
    time_source = if (nrow(h5_row) == 1L) "H5 mean_delta_minus_c" else "H5 c is not defined for this clock",
    title = title
  )
  beta <- stats::coef(fit)
  beta_rows[[clock_id]] <- data.table(clock = clock_id, term = names(beta), beta = as.numeric(beta), title = title)
  age <- cov$age[match(subjects, cov$subject_id)]
  half <- as.numeric(as.matrix(x) %*% beta[safe]) / 2
  morning_adj <- morning - half
  evening_adj <- evening + half
  person <- (morning + evening) / 2
  person_adj <- (morning_adj + evening_adj) / 2
  one_r <- function(x) {
    ok <- is.finite(x) & is.finite(age)
    if (sum(ok) < 3L) return(c(r = NA_real_, n = sum(ok)))
    c(r = unname(stats::cor(x[ok], age[ok])), n = sum(ok))
  }
  vals <- list(
    morning_before = one_r(morning), morning_after = one_r(morning_adj),
    evening_before = one_r(evening), evening_after = one_r(evening_adj),
    person_mean = one_r(person)
  )
  age_rows[[clock_id]] <- data.table(
    clock = clock_id,
    morning_r_before = vals$morning_before["r"], morning_r_after = vals$morning_after["r"],
    evening_r_before = vals$evening_before["r"], evening_r_after = vals$evening_after["r"],
    person_mean_r = vals$person_mean["r"], n = vals$person_mean["n"],
    max_abs_person_mean_change = max(abs(person_adj - person)),
    title = title
  )
}
fwrite(rbindlist(l3_rows), file.path(out_dir, "composition_l3.tsv"), sep = "\t")
fwrite(rbindlist(beta_rows), file.path(out_dir, "composition_l3_betas.tsv"), sep = "\t")
fwrite(rbindlist(age_rows), file.path(out_dir, "age_correlation_l3.tsv"), sep = "\t")

cos_out <- rbindlist(comp_rows, fill = TRUE)
sum_rows <- list()
for (i in seq_len(nrow(cos_out))) {
  sum_rows[[length(sum_rows) + 1L]] <- data.table(
    dataset = cos_out$dataset[i], clock = cos_out$clock[i], effect = "two_a",
    effect_before = cos_out$two_a_before[i], effect_composition = cos_out$two_a_after[i],
    effect_time = NA_real_, two_a_ratio = cos_out$two_a_ratio[i], title = title
  )
}
for (clock_id in clocks) {
  row <- l3_rows[[clock_id]]
  sum_rows[[length(sum_rows) + 1L]] <- data.table(
    dataset = "L3", clock = clock_id, effect = "mean_delta",
    effect_before = row$mean_delta_before, effect_composition = row$mean_delta_composition,
    effect_time = row$mean_delta_time, two_a_ratio = NA_real_, title = title
  )
}
fwrite(rbindlist(sum_rows), file.path(out_dir, "mitigation_summary.tsv"), sep = "\t")
say("stage9 analysis finished")
