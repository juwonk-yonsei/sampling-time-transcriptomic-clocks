#!/usr/bin/env Rscript
# Score-level rehearsal. Uses the confirmation test functions. No bootstrap.

suppressPackageStartupMessages({
  library(lme4)
  library(data.table)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
root <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(root, "scripts", "lib_cosinor.R"))
source(file.path(root, "scripts", "confirm", "lib_tests.R"))
out_dir <- file.path(root, "results", "stage6", "rehearsal", "scores")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

f1 <- fread(file.path(root, "prereg", "frozen", "F1.tsv"))
f4 <- fread(file.path(root, "prereg", "frozen", "F4.tsv"))
g2_ra <- fread(file.path(root, "results", "stage4", "rnaseq", "g2_clocks.tsv"))
g2_ta <- fread(file.path(root, "results", "stage5", "g2_tams.tsv"))
primary <- c("TAms_chrono", "TAms_mort", "RA_default")

meta_l1 <- fread(file.path(root, "data", "processed", "GSE48113_samples.tsv"))
d1_idx <- which(meta_l1$condition == "control" & meta_l1$period == "D1" & !(meta_l1$gsm %in% c("GSM1168837", "GSM1168735")) & is.finite(meta_l1$clock_time_h))
d4_idx <- which(meta_l1$period == "D4" & is.finite(meta_l1$clock_time_h))
l1_d1 <- meta_l1[d1_idx, .(subject_id = as.character(subject_id), clock_time_h)]
l1_d4 <- meta_l1[d4_idx, .(subject_id = as.character(subject_id), clock_time_h, day = 1)]
l1_d1[, day := 0]
l1_both <- rbind(l1_d1, l1_d4, use.names = TRUE)

meta_l2 <- fread(file.path(root, "data", "processed", "GSE122725_samples.tsv"))
wake <- fread(file.path(root, "results", "stage2", "GSE122725_wake_ref.tsv"))
wake[, gsm := as.character(gsm)]
meta_l2[, gsm := as.character(gsm)]
l2 <- merge(meta_l2, wake[, .(gsm, wake_ref_h)], by = "gsm")
l2[, t := wake_ref_h %% 24]
p1_idx <- which(l2$condition == "control" & l2$period == "period1")
p2_idx <- which(l2$condition == "night_oriented" & l2$period == "period2")
l2_p1 <- l2[p1_idx, .(subject_id = as.character(subject_id), clock_time_h = t, day = 0)]
l2_p2 <- l2[p2_idx, .(subject_id = as.character(subject_id), clock_time_h = t, day = 1)]
l2_both <- rbind(l2_p1, l2_p2, use.names = TRUE)
cov <- fread(file.path(root, "results", "stage5", "gse295758_subject_covariates.tsv"))
cov <- cov[, .(subject_id, age = as.numeric(age), sex, opioid_group)]
fwrite(data.table(
  design = c("L1_day1", "L1_day4", "L1_both", "L2_period1", "L2_period2", "L2_both", "L3"),
  n_samples = c(nrow(l1_d1), nrow(l1_d4), nrow(l1_both), nrow(l2_p1), nrow(l2_p2), nrow(l2_both), nrow(cov)),
  n_subjects = c(uniqueN(l1_d1$subject_id), uniqueN(l1_d4$subject_id), uniqueN(l1_both$subject_id),
                 uniqueN(l2_p1$subject_id), uniqueN(l2_p2$subject_id), uniqueN(l2_both$subject_id), nrow(cov))
), file.path(out_dir, "design_counts.tsv"), sep = "\t")

params <- function(clock) {
  hit <- which(f1$clock == clock & f1$dataset == "GSE113883")
  row <- f1[hit]
  c(A = row$two_a / 2, phi = row$peak_h, sigma = row$sigma, sigma_u = row$sigma_u)
}
sd_delta <- function(clock) {
  if (clock %in% g2_ta$clock) return(g2_ta$g2_sd_diff[match(clock, g2_ta$clock)])
  g2_ra$sd_diff[match(clock, g2_ra$outcome)]
}

binom_row <- function(k, n) {
  ci <- stats::binom.test(k, n)$conf.int
  c(frac = k / n, low = ci[1], high = ci[2])
}

sim_cosinor <- function(design, A, phi, sigma_u, sigma) {
  design <- as.data.frame(design)
  design$subject_id <- factor(design$subject_id)
  u <- stats::rnorm(nlevels(design$subject_id), 0, sigma_u)
  names(u) <- levels(design$subject_id)
  design$score <- unname(u[design$subject_id]) + A * cos(2 * pi * (design$clock_time_h - phi) / 24) +
    stats::rnorm(nrow(design), 0, sigma)
  fit_cosinor(add_time_terms(design), forms$full, forms$null)
}

signed_circ <- function(a, b) {
  d <- (a - b + 12) %% 24 - 12
  d
}

rows <- list()
add <- function(scenario, design, clock, kind, n_rep, p, bias) {
  finite <- is.finite(p)
  if (!any(finite)) {
    br <- c(frac = NA_real_, low = NA_real_, high = NA_real_)
  } else {
    br <- binom_row(sum(p < 0.05, na.rm = TRUE), n_rep)
  }
  rows[[length(rows) + 1]] <<- data.table(
    scenario = scenario, design = design, clock = clock, kind = kind, n_rep = n_rep,
    n_finite_p = sum(is.finite(p)), frac_p_lt_0.05 = unname(br["frac"]),
    frac_low = unname(br["low"]), frac_high = unname(br["high"]),
    bias_2a = bias["bias_2a"], bias_peak = bias["bias_peak"],
    bias_p_med = bias["bias_p_med"], bias_age = bias["bias_age"], bias_mean_delta_minus_c = bias["bias_mean_delta_minus_c"]
  )
}

run_h1 <- function(id, design_name, design, kind, n_rep) {
  set.seed(20261015 + as.integer(id))
  for (clock in primary) {
    par <- params(clock)
    A <- if (kind == "null") 0 else unname(par["A"])
    p <- numeric(n_rep)
    two_a <- numeric(n_rep)
    peak <- numeric(n_rep)
    for (b in seq_len(n_rep)) {
      fit <- sim_cosinor(design, A, unname(par["phi"]), unname(par["sigma_u"]), unname(par["sigma"]))
      p[b] <- fit$lrt_p
      two_a[b] <- fit$two_a
      peak[b] <- fit$peak_h
    }
    bias <- c(bias_2a = mean(two_a) - 2 * A, bias_peak = if (A == 0) NA_real_ else signed_circ(circ_mean_h(peak), unname(par["phi"])),
              bias_p_med = NA_real_, bias_age = NA_real_, bias_mean_delta_minus_c = NA_real_)
    add(paste0("H1_", design_name), design_name, clock, kind, n_rep, p, bias)
    cat("H1", design_name, clock, kind, "\n"); flush.console()
  }
}

run_h1("01", "L1", l1_d1, "null", 500)
run_h1("02", "L1", l1_d1, "planted", 200)
run_h1("03", "L2", l2_p1, "null", 500)
run_h1("04", "L2", l2_p1, "planted", 200)

run_h1_l3 <- function(id, kind, n_rep) {
  set.seed(20261015 + as.integer(id))
  for (clock in primary) {
    c_value <- f4$c[match(clock, f4$clock)]
    sd <- sd_delta(clock)
    mu <- if (kind == "null") 0 else c_value
    p <- numeric(n_rep)
    means <- numeric(n_rep)
    for (b in seq_len(n_rep)) {
      delta <- stats::rnorm(nrow(cov), mu, sd)
      tt <- t.test(delta)
      p[b] <- tt$p.value
      means[b] <- mean(delta)
    }
    bias <- c(bias_2a = NA_real_, bias_peak = NA_real_, bias_p_med = NA_real_, bias_age = NA_real_,
              bias_mean_delta_minus_c = mean(means - c_value))
    add("H1_L3", "L3", clock, kind, n_rep, p, bias)
    cat("H1 L3", clock, kind, "\n"); flush.console()
  }
}
run_h1_l3("05", "null", 500)
run_h1_l3("06", "planted", 200)

run_h3 <- function(id, kind, n_rep) {
  set.seed(20261015 + as.integer(id))
  age10 <- cov$age / 10
  for (clock in primary) {
    c_value <- f4$c[match(clock, f4$clock)]
    slope <- if (kind == "null") 0 else -0.25 * c_value
    sd <- sd_delta(clock)
    p <- numeric(n_rep)
    est <- numeric(n_rep)
    for (b in seq_len(n_rep)) {
      delta <- slope * age10 + stats::rnorm(nrow(cov), 0, sd)
      fit <- h3_age(delta, cov$age, cov$sex, cov$opioid_group)
      p[b] <- fit$age_p
      est[b] <- fit$age_coef
    }
    bias <- c(bias_2a = NA, bias_peak = NA, bias_p_med = NA, bias_age = mean(est) - slope, bias_mean_delta_minus_c = NA)
    add("H3", "L3", clock, kind, n_rep, p, bias)
    cat("H3", clock, kind, "\n"); flush.console()
  }
}
run_h3("07", "null", 500)
run_h3("08", "planted", 200)

run_h4 <- function(id, design_name, design, kind, n_rep) {
  set.seed(20261015 + as.integer(id))
  for (clock in primary) {
    par <- params(clock)
    A <- unname(par["A"])
    p <- numeric(n_rep)
    for (b in seq_len(n_rep)) {
      dat <- as.data.frame(design)
      dat$subject_id <- factor(dat$subject_id)
      u <- stats::rnorm(nlevels(dat$subject_id), 0, unname(par["sigma_u"]))
      names(u) <- levels(dat$subject_id)
      amp <- ifelse(dat$day == 1 & kind == "planted", A / 2, A)
      dat$score <- unname(u[dat$subject_id]) + amp * cos(2 * pi * (dat$clock_time_h - unname(par["phi"])) / 24) +
        stats::rnorm(nrow(dat), 0, unname(par["sigma"]))
      dat <- add_time_terms(dat)
      fit <- h4_interaction(dat)
      p[b] <- fit$lrt_p
    }
    bias <- c(bias_2a = NA, bias_peak = NA, bias_p_med = NA, bias_age = NA, bias_mean_delta_minus_c = NA)
    add(paste0("H4_", design_name), design_name, clock, kind, n_rep, p, bias)
    cat("H4", design_name, clock, kind, "\n"); flush.console()
  }
}
run_h4("09", "L1", l1_both, "null", 500)
run_h4("10", "L1", l1_both, "planted", 200)
run_h4("11", "L2", l2_both, "null", 500)
run_h4("12", "L2", l2_both, "planted", 200)

run_h2 <- function(id, kind, n_rep) {
  set.seed(20261015 + as.integer(id))
  for (clock in primary) {
    par <- params(clock)
    A <- unname(par["A"])
    phi <- unname(par["phi"])
    p <- rep(NA_real_, n_rep)
    est <- numeric(n_rep)
    target <- if (kind == "null") 0 else 0.6
    for (b in seq_len(n_rep)) {
      dat <- as.data.frame(l1_d1)
      dat$subject_id <- factor(dat$subject_id)
      theta <- cos(2 * pi * (dat$clock_time_h - phi) / 24)
      if (kind == "null") {
        dat[["T cells"]] <- stats::rnorm(nrow(dat))
        direct <- A
        alpha <- 0
      } else {
        dat[["T cells"]] <- theta + stats::rnorm(nrow(dat), 0, 0.5)
        alpha <- 0.6 * A
        direct <- 0.4 * A
      }
      u <- stats::rnorm(nlevels(dat$subject_id), 0, unname(par["sigma_u"]))
      names(u) <- levels(dat$subject_id)
      dat$score <- alpha * dat[["T cells"]] + direct * theta + unname(u[dat$subject_id]) +
        stats::rnorm(nrow(dat), 0, unname(par["sigma"]))
      dat <- add_time_terms(dat)
      fit <- h2_l1l2(dat, "T cells", NULL)
      est[b] <- fit$point$p_med
    }
    bias <- c(bias_2a = NA, bias_peak = NA, bias_p_med = mean(est) - target, bias_age = NA, bias_mean_delta_minus_c = NA)
    add("H2", "L1", clock, kind, n_rep, p, bias)
    cat("H2", clock, kind, "\n"); flush.console()
  }
}
run_h2("13", "null", 500)
run_h2("14", "planted", 200)

run_h5 <- function(id, scale, n_rep) {
  set.seed(20261015 + as.integer(id))
  for (clock in primary) {
    c_value <- f4$c[match(clock, f4$clock)]
    sd <- sd_delta(clock)
    mu <- scale * c_value
    p <- numeric(n_rep)
    bias_vals <- numeric(n_rep)
    for (b in seq_len(n_rep)) {
      delta <- stats::rnorm(nrow(cov), mu, sd)
      fit <- h5_transfer(delta, c_value)
      tt <- t.test(delta - c_value)
      p[b] <- tt$p.value
      bias_vals[b] <- fit$mean_delta_minus_c
    }
    bias <- c(bias_2a = NA, bias_peak = NA, bias_p_med = NA, bias_age = NA, bias_mean_delta_minus_c = mean(bias_vals))
    add("H5", "L3", clock, paste0("planted_", scale, "c"), n_rep, p, bias)
    cat("H5", clock, scale, "\n"); flush.console()
  }
}
run_h5("15", 1, 200)
run_h5("16", 0.5, 200)

out <- rbindlist(rows)
fwrite(out, file.path(out_dir, "score_rehearsal.tsv"), sep = "\t")
seeds <- data.table(
  rule = "set.seed(20261015 + scenario_id) at the start of each scenario; scenario_id is 1 to 16 in the order H1 L1 null, H1 L1 planted, H1 L2 null, H1 L2 planted, H1 L3 null, H1 L3 planted, H3 null, H3 planted, H4 L1 null, H4 L1 planted, H4 L2 null, H4 L2 planted, H2 null, H2 planted, H5 1c, H5 0.5c",
  base_seed = 20261015
)
fwrite(seeds, file.path(out_dir, "seed_rule.tsv"), sep = "\t")
cat("score rehearsal done", nrow(out), "\n")
