# Hypothesis tests for the confirmation pipeline. Source scripts/lib_cosinor.R before this file.

forms <- main_formula()
blood_labels <- c(
  "T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells",
  "B lineage", "Monocytic lineage", "Myeloid dendritic cells", "Neutrophils"
)
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

prepare_score <- function(samples, scores, clock_id) {
  idx <- which(scores$clock == clock_id)
  sc <- scores[idx]
  named <- sc$value
  names(named) <- sc$gsm
  dat <- as.data.frame(samples)
  dat$score <- unname(named[dat$gsm])
  if (anyNA(dat$score)) stop(paste("missing score", clock_id))
  dat$subject_id <- factor(dat$subject_id)
  add_time_terms(dat)
}

model_m <- function(dat) {
  fit <- fit_cosinor(dat, forms$full, forms$null)
  extra <- quiet_lmer(forms$full, dat, TRUE)
  fit$sigma_reml <- if (!inherits(extra, "error") && optimizer_converged(extra)) sigma(extra) else NA_real_
  fit
}

peak_in_window <- function(peak, start, end) {
  # Window may cross midnight when end < start, e.g. 18 to 2.
  p <- peak %% 24
  if (start <= end) return(p >= start && p <= end)
  p >= start || p <= end
}

gate_cosinor <- function(expr_gene, samples, start, end) {
  rows <- list()
  for (gene in c("PER1", "PER3")) {
    if (!(gene %in% rownames(expr_gene))) {
      rows[[gene]] <- data.table(gene = gene, lrt_p = NA_real_, peak_h = NA_real_, pass = FALSE, reason = "gene_absent")
      next
    }
    dat <- as.data.frame(samples)
    dat$score <- as.numeric(expr_gene[gene, dat$gsm])
    dat <- add_time_terms(dat)
    fit <- model_m(dat)
    in_window <- is.finite(fit$peak_h) && peak_in_window(fit$peak_h, start, end)
    rows[[gene]] <- data.table(
      gene = gene, lrt_p = fit$lrt_p, peak_h = fit$peak_h,
      pass = is.finite(fit$lrt_p) && fit$lrt_p < 0.01 && in_window,
      reason = ""
    )
  }
  out <- rbindlist(rows)
  out[, dataset_pass := all(pass)]
  out
}

gate_l3 <- function(gene_log, samples) {
  # samples need columns gsm, subject_id, phase in {morning, evening}
  rows <- list()
  for (gene in "PER1") {
    if (!(gene %in% rownames(gene_log))) {
      return(data.table(gene = gene, mean_delta = NA_real_, p = NA_real_, pass = FALSE, dataset_pass = FALSE))
    }
    vals <- gene_log[gene, samples$gsm]
    wide <- data.table(subject_id = samples$subject_id, phase = samples$phase, value = as.numeric(vals))
    morning <- wide[phase == "morning", .(morning = mean(value)), by = subject_id]
    evening <- wide[phase == "evening", .(evening = mean(value)), by = subject_id]
    both <- merge(morning, evening, by = "subject_id")
    diff <- both$morning - both$evening
    tt <- t.test(diff)
    rows[[gene]] <- data.table(
      gene = gene, mean_delta = unname(tt$estimate), p = unname(tt$p.value),
      pass = unname(tt$estimate) > 0 && unname(tt$p.value) < 0.01,
      dataset_pass = unname(tt$estimate) > 0 && unname(tt$p.value) < 0.01
    )
  }
  rbindlist(rows)
}

paired_delta <- function(samples, scores, clock_id) {
  idx <- which(scores$clock == clock_id)
  sc <- scores[idx]
  named <- sc$value
  names(named) <- sc$gsm
  wide <- data.table(subject_id = samples$subject_id, phase = samples$phase, gsm = samples$gsm)
  wide[, value := unname(named[gsm])]
  morning <- wide[phase == "morning", .(morning = mean(value)), by = subject_id]
  evening <- wide[phase == "evening", .(evening = mean(value)), by = subject_id]
  both <- merge(morning, evening, by = "subject_id")
  diff <- both$morning - both$evening
  tt <- t.test(diff)
  data.table(
    clock = clock_id, n = length(diff), mean_delta = unname(tt$estimate),
    sd_delta = stats::sd(diff), t_stat = unname(tt$statistic),
    p = unname(tt$p.value), ci_low = unname(tt$conf.int[1]), ci_high = unname(tt$conf.int[2]),
    subjects = list(both$subject_id), morning = list(both$morning), evening = list(both$evening),
    delta = list(diff)
  )
}

participant_draws <- function(subjects, n_boot = 2000, seed = 20261015) {
  set.seed(seed)
  subjects <- as.character(subjects)
  rows <- vector("list", n_boot)
  for (b in seq_len(n_boot)) {
    drawn <- sample(subjects, length(subjects), replace = TRUE)
    rows[[b]] <- data.table(draw = b, subject_id = drawn, new_subject_id = paste0(drawn, "#", seq_along(drawn)))
  }
  rbindlist(rows)
}

bootstrap_model_m <- function(dat, draws) {
  reps <- vector("list", uniqueN(draws$draw))
  for (b in seq_len(length(reps))) {
    one <- draws[draw == b]
    boot_dat <- apply_draw(dat, one$subject_id, one$new_subject_id)
    reps[[b]] <- fit_cosinor(boot_dat, forms$full, forms$null)
  }
  rbindlist(lapply(reps, as.data.table))
}

h2_l1l2 <- function(dat, composition, n_boot_draws = NULL) {
  safe <- safe_name(composition)
  names(safe) <- composition
  for (label in composition) dat[[safe[[label]]]] <- dat[[label]]
  adj <- as.formula(paste("score ~ cos_t + sin_t +", paste(unname(safe), collapse = " + "), "+ (1 | subject_id)"))
  unadj <- model_m(dat)
  adj_fit <- quiet_lmer(adj, dat, TRUE)
  if (inherits(adj_fit, "error") || !optimizer_converged(adj_fit)) {
    return(list(point = data.table(p_med = NA_real_, lrt_p = unadj$lrt_p), reps = NULL))
  }
  beta <- lme4::fixef(adj_fit)
  ru <- c(unadj$beta_cos, unadj$beta_sin)
  ra <- c(unname(beta["cos_t"]), unname(beta["sin_t"]))
  ru2 <- sum(ru * ru)
  p_med <- 1 - sum(ra * ru) / ru2
  point <- data.table(p_med = p_med, lrt_p = unadj$lrt_p, two_a = unadj$two_a)
  reps <- NULL
  if (!is.null(n_boot_draws)) {
    vals <- numeric(uniqueN(n_boot_draws$draw))
    for (b in seq_along(vals)) {
      one <- n_boot_draws[draw == b]
      boot_dat <- apply_draw(dat, one$subject_id, one$new_subject_id)
      u <- quiet_lmer(forms$full, boot_dat, TRUE)
      a <- quiet_lmer(adj, boot_dat, TRUE)
      if (inherits(u, "error") || inherits(a, "error") || !optimizer_converged(u) || !optimizer_converged(a)) {
        vals[b] <- NA_real_
      } else {
        ub <- lme4::fixef(u)
        ab <- lme4::fixef(a)
        ru_b <- c(unname(ub["cos_t"]), unname(ub["sin_t"]))
        ra_b <- c(unname(ab["cos_t"]), unname(ab["sin_t"]))
        vals[b] <- 1 - sum(ra_b * ru_b) / sum(ru_b * ru_b)
      }
    }
    reps <- vals
    ci <- percentile_ci(vals)
    point[, p_med_p025 := ci[1]]
    point[, p_med_p975 := ci[2]]
  }
  list(point = point, reps = reps)
}

h3_age <- function(delta, age, sex, opioid) {
  frame <- data.frame(delta = delta, age10 = age / 10, sex = factor(sex), opioid = factor(opioid))
  fit <- lm(delta ~ age10 + sex + opioid, data = frame)
  coef <- summary(fit)$coefficients
  data.table(
    age_coef = unname(coef["age10", "Estimate"]),
    age_p = unname(coef["age10", "Pr(>|t|)"]),
    n = nrow(frame)
  )
}

h4_interaction <- function(dat) {
  full <- score ~ cos_t + sin_t + day + day:cos_t + day:sin_t + (1 | subject_id)
  null <- score ~ cos_t + sin_t + day + (1 | subject_id)
  ml <- quiet_lmer(full, dat, FALSE)
  base <- quiet_lmer(null, dat, FALSE)
  reml <- quiet_lmer(full, dat, TRUE)
  ok <- !inherits(ml, "error") && !inherits(base, "error") && optimizer_converged(ml) && optimizer_converged(base)
  stat <- if (ok) as.numeric(2 * (logLik(ml) - logLik(base))) else NA_real_
  p <- if (ok) pchisq(stat, df = 2, lower.tail = FALSE) else NA_real_
  data.table(lrt_stat = stat, lrt_p = p, converged = ok, reml_converged = !inherits(reml, "error") && optimizer_converged(reml))
}

h5_transfer <- function(delta, c_value) {
  tt <- t.test(delta - c_value)
  data.table(
    mean_delta = mean(delta),
    mean_delta_minus_c = unname(tt$estimate),
    ci_low = unname(tt$conf.int[1]),
    ci_high = unname(tt$conf.int[2]),
    ratio = mean(delta) / c_value
  )
}

flag_3mad <- function(log_mat) {
  med <- apply(log_mat, 2, median, na.rm = TRUE)
  center <- median(med, na.rm = TRUE)
  mad <- median(abs(med - center), na.rm = TRUE)
  flagged <- if (is.finite(mad) && mad > 0) abs(med - center) > 3 * mad else rep(FALSE, length(med))
  data.table(gsm = colnames(log_mat), median_log = med, flag_3mad = flagged)
}

holm_family <- function(pvals) {
  p.adjust(pvals, method = "holm")
}

sensitivity_intercept <- function(delta, sex, opioid) {
  frame <- data.frame(delta = delta, sex = factor(sex), opioid = factor(opioid))
  fit <- lm(delta ~ sex + opioid, data = frame)
  coef <- summary(fit)$coefficients
  data.table(
    intercept = unname(coef["(Intercept)", "Estimate"]),
    intercept_p = unname(coef["(Intercept)", "Pr(>|t|)"]),
    n = nrow(frame)
  )
}

h2_l3 <- function(delta_clock, delta_comp) {
  frame <- data.frame(delta = delta_clock, delta_comp, check.names = FALSE)
  fit <- stats::lm(delta ~ ., data = frame)
  alpha <- unname(stats::coef(fit)[["(Intercept)"]])
  data.table(alpha = alpha, mean_delta = mean(delta_clock), p_med = 1 - alpha / mean(delta_clock))
}

h2_l3_bootstrap <- function(delta_clock, delta_comp, n_boot = 2000, seed = 20261015) {
  set.seed(seed)
  n <- length(delta_clock)
  vals <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    comp_b <- delta_comp[idx, , drop = FALSE]
    vals[b] <- h2_l3(delta_clock[idx], comp_b)$p_med
  }
  point <- h2_l3(delta_clock, delta_comp)
  ci <- percentile_ci(vals)
  point[, p_med_p025 := ci[1]]
  point[, p_med_p975 := ci[2]]
  point
}
