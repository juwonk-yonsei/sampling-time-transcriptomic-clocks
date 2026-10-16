# Shared cosinor helpers for stage 2. Estimates use REML. Likelihood metrics use ML.

peak_0_24 <- function(beta_sin, beta_cos) {
  ((24 / (2 * pi)) * atan2(beta_sin, beta_cos)) %% 24
}

circ_dist_0_12 <- function(a, b) {
  d <- abs(a - b) %% 24
  pmin(d, 24 - d)
}

circ_mean_h <- function(hours) {
  hours <- hours[is.finite(hours)]
  ang <- hours / 24 * 2 * pi
  ((atan2(mean(sin(ang)), mean(cos(ang)))) %% (2 * pi)) * 24 / (2 * pi)
}

quiet_lmer <- function(formula, data, reml) {
  warn <- character()
  fit <- withCallingHandlers(
    tryCatch(
      lme4::lmer(formula, data = data, REML = reml, na.action = na.exclude),
      error = function(e) e
    ),
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (inherits(fit, "error")) warn <- c(warn, conditionMessage(fit))
  attr(fit, "warn_text") <- paste(unique(warn), collapse = " | ")
  fit
}

optimizer_converged <- function(fit) {
  if (inherits(fit, "error")) return(FALSE)
  opt <- fit@optinfo$conv$opt
  if (is.null(opt)) return(TRUE)
  isTRUE(all(unlist(opt) == 0))
}

fit_cosinor <- function(data, formula_full, formula_null) {
  data$subject_id <- factor(data$subject_id)
  reml <- quiet_lmer(formula_full, data, TRUE)
  ml <- quiet_lmer(formula_full, data, FALSE)
  null <- quiet_lmer(formula_null, data, FALSE)
  out <- data.frame(
    n_subjects = length(unique(data$subject_id)),
    n_samples = nrow(data),
    beta_cos = NA_real_,
    beta_sin = NA_real_,
    beta_intercept = NA_real_,
    amplitude = NA_real_,
    two_a = NA_real_,
    peak_h = NA_real_,
    r2_within = NA_real_,
    sd_between = NA_real_,
    sd_resid_null = NA_real_,
    two_a_over_sd_between = NA_real_,
    two_a_over_sd_resid = NA_real_,
    lrt_stat = NA_real_,
    lrt_df = 2,
    lrt_p = NA_real_,
    icc_null = NA_real_,
    reml_ok = !inherits(reml, "error"),
    ml_ok = !inherits(ml, "error"),
    null_ok = !inherits(null, "error"),
    reml_converged = optimizer_converged(reml),
    ml_converged = optimizer_converged(ml),
    null_converged = optimizer_converged(null),
    reml_singular = if (!inherits(reml, "error")) isTRUE(lme4::isSingular(reml, tol = 1e-4)) else NA,
    ml_singular = if (!inherits(ml, "error")) isTRUE(lme4::isSingular(ml, tol = 1e-4)) else NA,
    null_singular = if (!inherits(null, "error")) isTRUE(lme4::isSingular(null, tol = 1e-4)) else NA,
    warnings = paste(
      attr(reml, "warn_text"),
      attr(ml, "warn_text"),
      attr(null, "warn_text"),
      sep = " || "
    ),
    stringsAsFactors = FALSE
  )
  if (!inherits(reml, "error")) {
    beta <- lme4::fixef(reml)
    out$beta_intercept <- unname(beta["(Intercept)"])
    out$beta_cos <- unname(beta["cos_t"])
    out$beta_sin <- unname(beta["sin_t"])
    out$amplitude <- sqrt(out$beta_cos^2 + out$beta_sin^2)
    out$two_a <- 2 * out$amplitude
    out$peak_h <- peak_0_24(out$beta_sin, out$beta_cos)
  }
  if (!inherits(ml, "error") && !inherits(null, "error")) {
    stat <- as.numeric(2 * (logLik(ml) - logLik(null)))
    out$lrt_stat <- stat
    out$lrt_p <- pchisq(stat, df = 2, lower.tail = FALSE)
    v_full <- sigma(ml)^2
    v_null <- sigma(null)^2
    out$r2_within <- if (is.finite(v_null) && v_null > 0) 1 - v_full / v_null else NA_real_
    out$sd_resid_null <- sqrt(v_null)
    vc <- lme4::VarCorr(null)
    sds <- attr(vc[[1]], "stddev")
    if (length(sds) == 1L) {
      sb <- as.numeric(sds)
      out$icc_null <- sb^2 / (sb^2 + v_null)
    }
  }
  means <- tapply(data$score, data$subject_id, mean)
  out$sd_between <- stats::sd(means)
  if (is.finite(out$two_a) && is.finite(out$sd_between) && out$sd_between > 0) {
    out$two_a_over_sd_between <- out$two_a / out$sd_between
  }
  if (is.finite(out$two_a) && is.finite(out$sd_resid_null) && out$sd_resid_null > 0) {
    out$two_a_over_sd_resid <- out$two_a / out$sd_resid_null
  }
  out
}

main_formula <- function() {
  list(
    full = score ~ cos_t + sin_t + (1 | subject_id),
    null = score ~ 1 + (1 | subject_id)
  )
}

add_time_terms <- function(data) {
  data$cos_t <- cos(2 * pi * data$clock_time_h / 24)
  data$sin_t <- sin(2 * pi * data$clock_time_h / 24)
  data
}

apply_draw <- function(data, subject_ids, new_ids) {
  pieces <- vector("list", length(subject_ids))
  for (i in seq_along(subject_ids)) {
    part <- data[data$subject_id == subject_ids[[i]], , drop = FALSE]
    part$subject_id <- new_ids[[i]]
    pieces[[i]] <- part
  }
  do.call(rbind, pieces)
}

percentile_ci <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(NA_real_, NA_real_))
  stats::quantile(x, c(0.025, 0.975), names = FALSE, type = 7)
}

peak_interval <- function(hours) {
  hours <- hours[is.finite(hours)]
  if (length(hours) == 0) {
    return(c(center = NA_real_, halfwidth_95 = NA_real_))
  }
  center <- circ_mean_h(hours)
  half <- as.numeric(stats::quantile(circ_dist_0_12(hours, center), 0.95, names = FALSE, type = 7))
  c(center = center, halfwidth_95 = half)
}
