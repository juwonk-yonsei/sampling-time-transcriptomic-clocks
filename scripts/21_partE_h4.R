#!/usr/bin/env Rscript
# Part E. H4 time axes, sleep restriction, sleep state.

suppressPackageStartupMessages(library(data.table))

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
root <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(root, "scripts", "lib_cosinor.R"))
out_dir <- file.path(root, "results", "stage3", "h4")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clocks <- c("RA_GTEx", "RA_Peters", "RA_default")
sig_of <- c(RA_GTEx = "GTExAge", RA_Peters = "Peters", RA_default = "DESeq2")
blood_labels <- c(
  "T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells",
  "B lineage", "Monocytic lineage", "Myeloid dendritic cells", "Neutrophils"
)
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
  list(values = vals, matrix_symbol = symbol)
}

ml_aic <- function(data, formula) {
  fit <- quiet_lmer(formula, data, FALSE)
  if (inherits(fit, "error")) return(c(aic = NA_real_, converged = 0))
  c(aic = AIC(fit), converged = as.numeric(optimizer_converged(fit)))
}

lrt_ml <- function(data, full, reduced, df) {
  f <- quiet_lmer(full, data, FALSE)
  r <- quiet_lmer(reduced, data, FALSE)
  out <- data.table(
    lrt_df = df,
    lrt_stat = NA_real_,
    lrt_p = NA_real_,
    full_converged = !inherits(f, "error") && optimizer_converged(f),
    reduced_converged = !inherits(r, "error") && optimizer_converged(r),
    full_singular = if (!inherits(f, "error")) isTRUE(lme4::isSingular(f, tol = 1e-4)) else NA,
    reduced_singular = if (!inherits(r, "error")) isTRUE(lme4::isSingular(r, tol = 1e-4)) else NA
  )
  if (!inherits(f, "error") && !inherits(r, "error")) {
    stat <- as.numeric(2 * (logLik(f) - logLik(r)))
    out$lrt_stat <- stat
    out$lrt_p <- pchisq(stat, df = df, lower.tail = FALSE)
  }
  out
}

coef_condition <- function(data, formula, coef_name) {
  fit <- quiet_lmer(formula, data, FALSE)
  if (inherits(fit, "error") || !optimizer_converged(fit)) return(NA_real_)
  beta <- lme4::fixef(fit)
  if (!coef_name %in% names(beta)) return(NA_real_)
  unname(beta[[coef_name]])
}

sd_between <- function(score, subject_id) {
  stats::sd(tapply(score, subject_id, mean))
}

make_draws <- function(subjects, n_boot, seed) {
  set.seed(seed)
  rows <- vector("list", n_boot)
  for (b in seq_len(n_boot)) {
    picked <- sample(subjects, length(subjects), replace = TRUE)
    times <- ave(seq_along(picked), picked, FUN = seq_along)
    rows[[b]] <- data.table(
      draw = b,
      pick = seq_along(picked),
      subject_id = picked,
      new_subject_id = paste(picked, times, sep = "#")
    )
  }
  rbindlist(rows)
}

axis_rows <- list()

fit_axis <- function(dat, time_col, dataset, outcome, axis) {
  d <- dat
  d$clock_time_h <- d[[time_col]]
  d <- add_time_terms(d)
  d$score <- d$outcome_value
  fit <- fit_cosinor(d, forms$full, forms$null)
  aic <- ml_aic(d, forms$full)
  data.table(
    dataset = dataset,
    outcome = outcome,
    matrix_symbol = dat$matrix_symbol[1],
    axis = axis,
    time_col = time_col,
    n_subjects = fit$n_subjects,
    n_samples = fit$n_samples,
    two_a = fit$two_a,
    peak_h = fit$peak_h,
    lrt_p = fit$lrt_p,
    r2_within = fit$r2_within,
    aic_ml = unname(aic["aic"]),
    ml_converged = unname(aic["converged"]) == 1,
    reml_converged = fit$reml_converged,
    reml_singular = fit$reml_singular,
    n_unique_time_within_max = max(tapply(d$clock_time_h, d$subject_id, function(x) length(unique(x))))
  )
}

cat("5.1 GSE39445\n")
flush.console()
samples <- fread(file.path(root, "data", "processed", "GSE39445_samples.tsv"))
mat <- read_expr(file.path(root, "data", "processed", "GSE39445_expr_gene_log2.tsv.gz"))
base <- samples[which(samples$condition == "control" & is.finite(samples$clock_time_h) & is.finite(samples$circadian_ref_h) & is.finite(samples$hours_awake))]
base <- as.data.frame(base[, .(gsm, subject_id, clock_time_h, circadian_ref_h, hours_awake)])
outcomes <- list()
for (clock in clocks) outcomes[[clock]] <- lin_named(mat, clock)
for (gene in c("PER1", "PER3", "ARNTL")) {
  got <- gene_named(mat, gene)
  if (is.null(got)) stop(paste("missing", gene))
  outcomes[[gene]] <- got$values
  base[[paste0("symbol_", gene)]] <- got$matrix_symbol
}
rm(mat)
gc()
for (outcome in names(outcomes)) {
  dat <- base
  dat$outcome_value <- unname(outcomes[[outcome]][dat$gsm])
  dat$matrix_symbol <- if (outcome %in% c("PER1", "PER3", "ARNTL")) dat[[paste0("symbol_", outcome)]][1] else outcome
  if (anyNA(dat$outcome_value)) stop(paste("missing outcome", outcome))
  for (spec in list(
    c(axis = "T_clock", col = "clock_time_h"),
    c(axis = "T_circ", col = "circadian_ref_h"),
    c(axis = "T_wake", col = "hours_awake")
  )) {
    axis_rows[[paste("GSE39445", outcome, spec["axis"])]] <- fit_axis(
      dat, spec["col"], "GSE39445", outcome, spec["axis"]
    )
  }
}

cat("5.1 GSE113883 genes only\n")
flush.console()
samples_r <- fread(file.path(root, "results", "stage3", "GSE113883_samples.tsv"))
mat_r <- read_expr(file.path(root, "results", "stage3", "expr", "GSE113883_log2tpm1.tsv.gz"))
base_r <- samples_r[which(is.finite(samples_r$circadian_ref_h) & is.finite(samples_r$clock_time_h))]
base_r <- as.data.frame(base_r[, .(gsm, subject_id, clock_time_h, circadian_ref_h)])
for (gene in c("PER1", "PER3", "ARNTL")) {
  got <- gene_named(mat_r, gene)
  if (is.null(got)) stop(paste("missing", gene))
  dat <- base_r
  dat$outcome_value <- unname(got$values[dat$gsm])
  dat$matrix_symbol <- got$matrix_symbol
  if (anyNA(dat$outcome_value)) stop(paste("missing", gene))
  for (spec in list(
    c(axis = "T_clock", col = "clock_time_h"),
    c(axis = "T_dlmo", col = "circadian_ref_h")
  )) {
    axis_rows[[paste("GSE113883", gene, spec["axis"])]] <- fit_axis(
      dat, spec["col"], "GSE113883", gene, spec["axis"]
    )
  }
}
rm(mat_r)
gc()
fwrite(data.table(
  dataset = "GSE113883",
  clock_axes = "not_run",
  reason = "positive control failed; lin clock scores were not produced in section 3.2"
), file.path(out_dir, "gse113883_clock_axis_status.tsv"), sep = "\t")

axes <- rbindlist(axis_rows)
ref <- axes[axes$axis == "T_clock", .(dataset, outcome, aic_clock = aic_ml)]
axes <- merge(axes, ref, by = c("dataset", "outcome"), sort = FALSE)
axes[, delta_aic_vs_T_clock := aic_ml - aic_clock]
fwrite(axes, file.path(out_dir, "time_axis_comparison.tsv"), sep = "\t")

cat("5.2 sleep restriction\n")
flush.console()
samples <- fread(file.path(root, "data", "processed", "GSE39445_samples.tsv"))
mat <- read_expr(file.path(root, "data", "processed", "GSE39445_expr_gene_log2.tsv.gz"))
use <- samples[which(samples$condition %in% c("control", "sleep_restriction") & is.finite(samples$clock_time_h))]
use <- as.data.frame(use[, .(gsm, subject_id, condition, clock_time_h, hours_awake)])
use$condition <- factor(use$condition, levels = c("control", "sleep_restriction"))
use <- add_time_terms(use)
ctrl_ids <- sort(unique(use$subject_id[use$condition == "control"]))
rest_ids <- sort(unique(use$subject_id[use$condition == "sleep_restriction"]))
both <- sort(intersect(ctrl_ids, rest_ids))
draws <- make_draws(both, 2000L, 20261015L)
fwrite(draws, file.path(out_dir, "sleep_restriction_boot_draws.tsv"), sep = "\t")

form_mean <- score ~ cos_t + sin_t + condition + (1 | subject_id)
form_base <- score ~ cos_t + sin_t + (1 | subject_id)
form_full <- score ~ cos_t + sin_t + condition + condition:cos_t + condition:sin_t + (1 | subject_id)
form_mean_w <- score ~ cos_t + sin_t + condition + hours_awake + (1 | subject_id)
form_base_w <- score ~ cos_t + sin_t + hours_awake + (1 | subject_id)
form_full_w <- score ~ cos_t + sin_t + condition + condition:cos_t + condition:sin_t + hours_awake + (1 | subject_id)
coef_name <- "conditionsleep_restriction"

restrict_rows <- list()
cond_rows <- list()
boot_sum <- list()
for (clock in clocks) {
  cat(" clock", clock, "\n")
  flush.console()
  named <- lin_named(mat, clock)
  dat <- use
  dat$score <- unname(named[dat$gsm])
  if (anyNA(dat$score)) stop("missing lin score")
  ctrl <- dat[dat$condition == "control", ]
  sdb <- sd_between(ctrl$score, ctrl$subject_id)
  diff0 <- coef_condition(dat, form_mean, coef_name)
  if (!is.finite(diff0)) stop(paste("mean difference not finite", clock))
  lrt_main <- lrt_ml(dat, form_mean, form_base, 1)
  lrt_int <- lrt_ml(dat, form_full, form_mean, 2)
  dat_w <- dat[is.finite(dat$hours_awake), ]
  diff_w <- coef_condition(dat_w, form_mean_w, coef_name)
  lrt_main_w <- lrt_ml(dat_w, form_mean_w, form_base_w, 1)
  lrt_int_w <- lrt_ml(dat_w, form_full_w, form_mean_w, 2)
  restrict_rows[[clock]] <- data.table(
    dataset = "GSE39445",
    clock = clock,
    method = "lin",
    n_samples = nrow(dat),
    n_subjects_control = length(ctrl_ids),
    n_subjects_sleep_restriction = length(rest_ids),
    n_subjects_both = length(both),
    mean_diff_restriction_minus_control = diff0,
    sd_between_control = sdb,
    diff_over_sd_between = diff0 / sdb,
    lrt_condition_df = lrt_main$lrt_df,
    lrt_condition_stat = lrt_main$lrt_stat,
    lrt_condition_p = lrt_main$lrt_p,
    lrt_condition_full_converged = lrt_main$full_converged,
    lrt_condition_reduced_converged = lrt_main$reduced_converged,
    lrt_interaction_df = lrt_int$lrt_df,
    lrt_interaction_stat = lrt_int$lrt_stat,
    lrt_interaction_p = lrt_int$lrt_p,
    lrt_interaction_full_converged = lrt_int$full_converged,
    lrt_interaction_reduced_converged = lrt_int$reduced_converged,
    sensitivity = "hours_awake",
    n_samples_sensitivity = nrow(dat_w),
    mean_diff_sensitivity = diff_w,
    diff_over_sd_between_sensitivity = diff_w / sdb,
    lrt_condition_p_sensitivity = lrt_main_w$lrt_p,
    lrt_interaction_p_sensitivity = lrt_int_w$lrt_p,
    lrt_condition_converged_sensitivity = lrt_main_w$full_converged && lrt_main_w$reduced_converged,
    lrt_interaction_converged_sensitivity = lrt_int_w$full_converged && lrt_int_w$reduced_converged
  )
  for (cond in c("control", "sleep_restriction")) {
    part <- dat[dat$condition == cond, ]
    fit <- fit_cosinor(part, forms$full, forms$null)
    cond_rows[[paste(clock, cond)]] <- data.table(
      dataset = "GSE39445",
      clock = clock,
      method = "lin",
      condition = cond,
      n_subjects = fit$n_subjects,
      n_samples = fit$n_samples,
      two_a = fit$two_a,
      peak_h = fit$peak_h,
      lrt_p = fit$lrt_p,
      reml_converged = fit$reml_converged
    )
  }
  reps <- numeric(2000L)
  reps_w <- numeric(2000L)
  ok <- logical(2000L)
  ok_w <- logical(2000L)
  for (b in seq_len(2000L)) {
    one <- draws[which(draws$draw == b)]
    boot <- apply_draw(dat, one$subject_id, one$new_subject_id)
    val <- coef_condition(boot, form_mean, coef_name)
    reps[b] <- val
    ok[b] <- is.finite(val)
    boot_w <- boot[is.finite(boot$hours_awake), ]
    val_w <- coef_condition(boot_w, form_mean_w, coef_name)
    reps_w[b] <- val_w
    ok_w[b] <- is.finite(val_w)
    if (b %% 200L == 0L) {
      cat("  ", b, "\n")
      flush.console()
    }
  }
  ci <- percentile_ci(reps)
  ci_w <- percentile_ci(reps_w)
  fwrite(data.table(draw = seq_len(2000L), mean_diff = reps, converged = ok, mean_diff_sensitivity = reps_w, converged_sensitivity = ok_w),
         file.path(out_dir, paste0("sleep_restriction_reps_", clock, ".tsv")), sep = "\t")
  boot_sum[[clock]] <- data.table(
    dataset = "GSE39445",
    clock = clock,
    method = "lin",
    n_boot = 2000L,
    n_finite = sum(ok),
    mean_diff = diff0,
    mean_diff_p025 = ci[1],
    mean_diff_p975 = ci[2],
    n_finite_sensitivity = sum(ok_w),
    mean_diff_sensitivity = diff_w,
    mean_diff_sensitivity_p025 = ci_w[1],
    mean_diff_sensitivity_p975 = ci_w[2]
  )
}
rm(mat)
gc()
fwrite(rbindlist(restrict_rows), file.path(out_dir, "sleep_restriction_lrt.tsv"), sep = "\t")
fwrite(rbindlist(cond_rows), file.path(out_dir, "sleep_restriction_by_condition.tsv"), sep = "\t")
fwrite(rbindlist(boot_sum), file.path(out_dir, "sleep_restriction_bootstrap.tsv"), sep = "\t")

cat("5.3 sleep state\n")
flush.console()
samples <- fread(file.path(root, "data", "processed", "GSE56931_samples.tsv"))
mat <- read_expr(file.path(root, "data", "processed", "GSE56931_expr_gene_log2.tsv.gz"))
mcp <- fread(file.path(root, "results", "stage2", "h2_mcp_scores_GSE56931.tsv"))
use <- samples[which(samples$period %in% c("baseline", "deprivation") & is.finite(samples$clock_time_h))]
use <- as.data.frame(use[, .(gsm, subject_id, period, clock_time_h)])
use <- merge(use, as.data.frame(mcp), by = "gsm", sort = FALSE)
use$period <- factor(use$period, levels = c("baseline", "deprivation"))
use <- add_time_terms(use)
form_mean <- score ~ cos_t + sin_t + period + (1 | subject_id)
form_base <- score ~ cos_t + sin_t + (1 | subject_id)
form_full <- score ~ cos_t + sin_t + period + period:cos_t + period:sin_t + (1 | subject_id)
coef_name <- "perioddeprivation"

state_rows <- list()
state_cond <- list()
for (clock in clocks) {
  named <- lin_named(mat, clock)
  dat <- use
  dat$score <- unname(named[dat$gsm])
  if (anyNA(dat$score)) stop("missing lin score")
  lrt_main <- lrt_ml(dat, form_mean, form_base, 1)
  lrt_int <- lrt_ml(dat, form_full, form_mean, 2)
  diff0 <- coef_condition(dat, form_mean, coef_name)
  if (!is.finite(diff0)) stop(paste("sleep-state difference not finite", clock))
  state_rows[[clock]] <- data.table(
    dataset = "GSE56931",
    outcome = clock,
    outcome_type = "clock_lin",
    model = "model1",
    n_samples = nrow(dat),
    n_subjects = length(unique(dat$subject_id)),
    mean_diff_deprivation_minus_baseline = diff0,
    lrt_period_df = lrt_main$lrt_df,
    lrt_period_stat = lrt_main$lrt_stat,
    lrt_period_p = lrt_main$lrt_p,
    lrt_period_full_converged = lrt_main$full_converged,
    lrt_period_reduced_converged = lrt_main$reduced_converged,
    lrt_interaction_df = lrt_int$lrt_df,
    lrt_interaction_stat = lrt_int$lrt_stat,
    lrt_interaction_p = lrt_int$lrt_p,
    lrt_interaction_full_converged = lrt_int$full_converged,
    lrt_interaction_reduced_converged = lrt_int$reduced_converged
  )
  for (per in c("baseline", "deprivation")) {
    part <- dat[dat$period == per, ]
    fit <- fit_cosinor(part, forms$full, forms$null)
    state_cond[[paste(clock, per)]] <- data.table(
      dataset = "GSE56931",
      outcome = clock,
      period = per,
      n_subjects = fit$n_subjects,
      n_samples = fit$n_samples,
      two_a = fit$two_a,
      peak_h = fit$peak_h,
      lrt_p = fit$lrt_p,
      reml_converged = fit$reml_converged
    )
  }
}

wald_rows <- list()
wald_one <- function(dat, outcome, outcome_type) {
  dat$tod <- factor(dat$clock_time_h, levels = c(0, 4, 8, 12, 16, 20))
  if (anyNA(dat$tod)) stop("clock time outside 0,4,8,12,16,20")
  fit <- quiet_lmer(score ~ tod * period + (1 | subject_id), dat, TRUE)
  converged <- !inherits(fit, "error") && optimizer_converged(fit)
  beta <- if (converged) lme4::fixef(fit) else NULL
  vv <- if (converged) as.matrix(vcov(fit)) else NULL
  for (hour in c(0, 4, 8, 12, 16, 20)) {
    n_base <- sum(dat$period == "baseline" & dat$clock_time_h == hour)
    n_dep <- sum(dat$period == "deprivation" & dat$clock_time_h == hour)
    est <- NA_real_
    se <- NA_real_
    if (converged) {
      L <- setNames(rep(0, length(beta)), names(beta))
      L["perioddeprivation"] <- 1
      if (hour != 0) {
        key <- paste0("tod", hour, ":perioddeprivation")
        if (!key %in% names(beta)) stop(paste("missing contrast", key))
        L[key] <- 1
      }
      est <- sum(L * beta)
      se <- sqrt(as.numeric(t(L) %*% vv %*% L))
    }
    wald_rows[[paste(outcome, hour)]] <<- data.table(
      dataset = "GSE56931",
      outcome = outcome,
      outcome_type = outcome_type,
      model = "model2",
      clock_time_h = hour,
      n_baseline = n_base,
      n_deprivation = n_dep,
      diff_deprivation_minus_baseline = est,
      se = se,
      ci95_low = est - qnorm(0.975) * se,
      ci95_high = est + qnorm(0.975) * se,
      reml_converged = converged,
      reml_singular = if (!inherits(fit, "error")) isTRUE(lme4::isSingular(fit, tol = 1e-4)) else NA
    )
  }
}
for (clock in clocks) {
  dat <- use
  dat$score <- unname(lin_named(mat, clock)[dat$gsm])
  wald_one(dat, clock, "clock_lin")
}
for (pop in blood_labels) {
  dat <- use
  dat$score <- dat[[pop]]
  if (anyNA(dat$score)) stop(paste("missing", pop))
  wald_one(dat, pop, "composition")
}
rm(mat)
gc()
fwrite(rbindlist(state_rows), file.path(out_dir, "sleep_state_model1.tsv"), sep = "\t")
fwrite(rbindlist(state_cond), file.path(out_dir, "sleep_state_by_period.tsv"), sep = "\t")
fwrite(rbindlist(wald_rows), file.path(out_dir, "sleep_state_model2_wald.tsv"), sep = "\t")
cat("done\n")
