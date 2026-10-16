#!/usr/bin/env Rscript
# Part C. H1 development exploration. Not a confirmation test.

suppressPackageStartupMessages({
  library(lme4)
  library(data.table)
  library(circular)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
root <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(root, "scripts", "lib_cosinor.R"))
out_dir <- file.path(root, "results", "stage2")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 20261015
n_boot <- 2000
n_random <- 1000
clocks <- c("RA_GTEx", "RA_Peters", "RA_default")
sig_of <- c(RA_GTEx = "GTExAge", RA_Peters = "Peters", RA_default = "DESeq2")
datasets <- c("GSE39445", "GSE56931")
forms <- main_formula()

ns <- asNamespace("RNAAgeCalc")
genelist <- get("genelist_all", ns)$blood

read_samples <- function(dataset) {
  fread(file.path(root, "data", "processed", paste0(dataset, "_samples.tsv")))
}

read_scores <- function(dataset) {
  fread(file.path(root, "results", "stage1", "clocks", paste0(dataset, "_rnaage.tsv")))
}

field_map <- function(text, key) {
  parts <- strsplit(text, "||", fixed = TRUE)[[1]]
  prefix <- paste0(key, ":")
  for (part in parts) {
    part <- trimws(part)
    if (startsWith(part, prefix)) return(trimws(substring(part, nchar(prefix) + 1)))
  }
  NA_character_
}

make_draws <- function(subjects) {
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

within_compare <- function(wide) {
  wide[, dev_lin := value_lin - mean(value_lin), by = subject_id]
  wide[, dev_off := value_official - mean(value_official), by = subject_id]
  wide[, diff := value_official - value_lin]
  wide[, diff_dev := diff - mean(diff), by = subject_id]
  data.table(
    pearson_within = cor(wide$dev_lin, wide$dev_off, use = "complete.obs"),
    sd_within_official_minus_lin = sd(wide$diff_dev),
    n = nrow(wide)
  )
}

full_expr <- new.env(parent = emptyenv())

read_full_expr <- function(dataset) {
  if (exists(dataset, envir = full_expr, inherits = FALSE)) {
    return(get(dataset, envir = full_expr, inherits = FALSE))
  }
  expr <- fread(file.path(root, "data", "processed", paste0(dataset, "_expr_gene_log2.tsv.gz")))
  mat <- as.matrix(expr[, -1])
  rownames(mat) <- expr$gene
  storage.mode(mat) <- "double"
  assign(dataset, mat, envir = full_expr)
  mat
}

lin_from_expr <- function(dataset, clock) {
  sig <- sig_of[[clock]]
  coefs <- genelist[[sig]]
  intercept <- unname(coefs[["(Intercept)"]])
  genes <- names(coefs)
  genes <- genes[genes != "(Intercept)"]
  beta <- coefs[genes]
  mat_all <- read_full_expr(dataset)
  have <- intersect(genes, rownames(mat_all))
  mat <- mat_all[have, , drop = FALSE]
  used <- beta[have]
  scores <- as.numeric(intercept + drop(crossprod(used, mat)))
  names(scores) <- colnames(mat)
  list(scores = scores, coef = used, matrix = mat)
}

control_frame <- function(samples, score_named) {
  dat <- samples[condition == "control" & !is.na(clock_time_h)]
  dat <- dat[, .(gsm, subject_id, clock_time_h, hours_awake, raw_characteristics)]
  dat[, score := score_named[gsm]]
  dat <- dat[!is.na(score)]
  add_time_terms(as.data.frame(dat))
}

cat("building bootstrap draws\n")
for (dataset in datasets) {
  samples <- read_samples(dataset)
  subjects <- sort(unique(samples[condition == "control" & !is.na(clock_time_h), subject_id]))
  draws <- make_draws(subjects)
  draws[, dataset := dataset]
  fwrite(draws, file.path(out_dir, paste0("boot_draws_", dataset, ".tsv")), sep = "\t")
  cat(dataset, "subjects", length(subjects), "draw rows", nrow(draws), "\n")
}

cat("official vs lin\n")
cmp_rows <- list()
for (dataset in datasets) {
  samples <- read_samples(dataset)
  scores <- read_scores(dataset)
  setnames(scores, "clock", "clock_name")
  ctrl <- samples[condition == "control" & !is.na(clock_time_h), .(gsm, subject_id)]
  for (clock in clocks) {
    lin <- scores[clock_name == clock & method == "lin", .(gsm, value_lin = value)]
    off <- scores[clock_name == clock & method == "official", .(gsm, value_official = value)]
    wide <- merge(merge(ctrl, lin, by = "gsm"), off, by = "gsm")
    stats <- within_compare(wide)
    stats[, dataset := dataset]
    stats[, clock := clock]
    cmp_rows[[paste(dataset, clock)]] <- stats
  }
}
cmp <- rbindlist(cmp_rows)
cmp[, use_lin_only := pearson_within >= 0.99]
fwrite(cmp, file.path(out_dir, "h1_official_lin_within.tsv"), sep = "\t")
print(cmp)

cat("reconstruct lin scores\n")
recon <- list()
expr_cache <- list()
for (dataset in datasets) {
  for (clock in clocks) {
    built <- lin_from_expr(dataset, clock)
    stored <- read_scores(dataset)
    setnames(stored, "clock", "clock_name")
    stored <- stored[clock_name == clock & method == "lin"]
    common <- intersect(names(built$scores), stored$gsm)
    diff <- built$scores[common] - stored$value[match(common, stored$gsm)]
    recon[[paste(dataset, clock)]] <- data.table(
      dataset = dataset,
      clock = clock,
      n_genes_used = length(built$coef),
      n_compared = length(common),
      max_abs_diff = max(abs(diff)),
      cor = cor(built$scores[common], stored$value[match(common, stored$gsm)])
    )
    expr_cache[[paste(dataset, clock)]] <- built
  }
}
recon <- rbindlist(recon)
fwrite(recon, file.path(out_dir, "h1_lin_reconstruction.tsv"), sep = "\t")
print(recon)
if (any(recon$max_abs_diff > 1e-6)) stop("lin reconstruction does not match stored scores")

use_lin_only <- function(ds, ck) {
  hit <- cmp$use_lin_only[cmp$dataset == ds & cmp$clock == ck]
  isTRUE(length(hit) == 1L && isTRUE(hit[[1]]))
}

summarize_boot <- function(rep, dataset, clock, method) {
  ok <- rep$reml_converged & rep$ml_converged & rep$null_converged
  summ <- data.table(
    dataset = dataset,
    clock = clock,
    method = method,
    n_boot = nrow(rep),
    n_converged = sum(ok),
    n_not_converged = sum(!ok)
  )
  for (metric in c("two_a", "r2_within", "two_a_over_sd_between", "two_a_over_sd_resid")) {
    ci <- percentile_ci(rep[[metric]][ok])
    summ[, paste0(metric, "_p025") := ci[1]]
    summ[, paste0(metric, "_p975") := ci[2]]
  }
  pk <- peak_interval(rep$peak_h[ok])
  summ[, peak_center := unname(pk["center"])]
  summ[, peak_halfwidth_95 := unname(pk["halfwidth_95"])]
  summ
}

point_rows <- list()
boot_sum_rows <- list()
for (dataset in datasets) {
  samples <- read_samples(dataset)
  scores <- read_scores(dataset)
  setnames(scores, "clock", "clock_name")
  draws <- fread(file.path(out_dir, paste0("boot_draws_", dataset, ".tsv")))
  for (clock in clocks) {
    methods <- if (isTRUE(use_lin_only(dataset, clock))) c("lin") else c("lin", "official")
    full_methods <- if (isTRUE(use_lin_only(dataset, clock))) "lin" else c("lin", "official")
    # point estimates for both methods always
    for (method in c("lin", "official")) {
      meth <- method
      sc <- scores[clock_name == clock & method == meth]
      named <- sc$value
      names(named) <- sc$gsm
      dat <- control_frame(samples, named)
      fit <- fit_cosinor(dat, forms$full, forms$null)
      fit$dataset <- dataset
      fit$clock <- clock
      fit$method <- method
      fit$role <- if (method == "lin" || !isTRUE(use_lin_only(dataset, clock))) "primary" else "appendix_point"
      point_rows[[paste(dataset, clock, method)]] <- fit
      do_boot <- method == "lin" || !isTRUE(use_lin_only(dataset, clock))
      if (!do_boot) next
      rep_path <- file.path(out_dir, paste0("h1_boot_reps_", dataset, "_", clock, "_", method, ".tsv"))
      if (file.exists(rep_path) && nrow(fread(rep_path)) == n_boot) {
        cat("reuse bootstrap", dataset, clock, method, "\n")
        flush.console()
        boot_sum_rows[[paste(dataset, clock, method)]] <- summarize_boot(fread(rep_path), dataset, clock, method)
        next
      }
      cat("bootstrap", dataset, clock, method, "\n")
      flush.console()
      reps <- vector("list", n_boot)
      for (b in seq_len(n_boot)) {
        one <- draws[draw == b]
        boot_dat <- apply_draw(dat, one$subject_id, one$new_subject_id)
        reps[[b]] <- fit_cosinor(boot_dat, forms$full, forms$null)
        if (b %% 200 == 0) {
          cat("  ", b, "\n")
          flush.console()
        }
      }
      rep <- rbindlist(lapply(reps, as.data.table))
      rep[, draw := seq_len(n_boot)]
      boot_sum_rows[[paste(dataset, clock, method)]] <- summarize_boot(rep, dataset, clock, method)
      fwrite(
        rep[, .(draw, two_a, r2_within, two_a_over_sd_between, two_a_over_sd_resid, peak_h,
                reml_converged, ml_converged, null_converged)],
        rep_path,
        sep = "\t"
      )
    }
  }
}
points <- rbindlist(point_rows, fill = TRUE)
fwrite(points, file.path(out_dir, "h1_point_estimates.tsv"), sep = "\t")
boot_sum <- rbindlist(boot_sum_rows)
fwrite(boot_sum, file.path(out_dir, "h1_bootstrap_ci.tsv"), sep = "\t")

cat("random signatures\n")
rand_sum <- list()
rand_done <- file.path(out_dir, "h1_random_signature.tsv")
if (file.exists(rand_done) && nrow(fread(rand_done)) == length(datasets) * length(clocks)) {
  cat("reuse random signatures\n")
} else {
for (dataset in datasets) {
  samples_dt <- read_samples(dataset)
  for (clock in clocks) {
    built <- expr_cache[[paste(dataset, clock)]]
    mat_all <- read_full_expr(dataset)
    mu <- rowMeans(mat_all, na.rm = TRUE)
    ok_gene <- is.finite(mu)
    mu <- mu[ok_gene]
    dec <- as.integer(ceiling(rank(mu, ties.method = "average") / length(mu) * 10))
    dec[dec < 1] <- 1L
    dec[dec > 10] <- 10L
    names(dec) <- names(mu)
    clock_genes <- intersect(names(built$coef), names(dec))
    coef <- built$coef[clock_genes]
    non_clock <- setdiff(names(dec), clock_genes)
    pools <- split(non_clock, dec[non_clock])
    need <- table(dec[clock_genes])
    pool_n <- vapply(names(need), function(d) length(pools[[d]]), integer(1))
    if (any(pool_n < as.integer(need))) {
      stop(paste(dataset, clock, "decile pool smaller than clock-gene count",
                 paste(names(need), as.integer(need), pool_n, sep = ":", collapse = " ")))
    }
    set.seed(seed)
    r2 <- rep(NA_real_, n_random)
    n_fail <- 0L
    gsm <- colnames(mat_all)
    for (i in seq_len(n_random)) {
      chosen <- character(length(clock_genes))
      for (d in names(need)) {
        pool <- pools[[d]]
        k <- as.integer(need[[d]])
        pick <- sample(pool, k)
        idx <- which(as.character(dec[clock_genes]) == d)
        chosen[idx] <- pick
      }
      shuffled <- sample(unname(coef))
      score_all <- drop(crossprod(shuffled, mat_all[chosen, , drop = FALSE]))
      names(score_all) <- gsm
      dat <- control_frame(samples_dt, score_all)
      fit <- fit_cosinor(dat, forms$full, forms$null)
      if (!isTRUE(fit$ml_converged) || !isTRUE(fit$null_converged)) n_fail <- n_fail + 1L
      r2[i] <- fit$r2_within
    }
    hit <- points$dataset == dataset & points$clock == clock & points$method == "lin"
    real <- points$r2_within[hit]
    if (length(real) != 1L) stop(paste("expected one lin R2", dataset, clock))
    finite <- r2[is.finite(r2)]
    rand_sum[[paste(dataset, clock)]] <- data.table(
      dataset = dataset,
      clock = clock,
      method = "lin",
      n_random = n_random,
      n_finite = length(finite),
      n_not_converged = n_fail,
      random_median = median(finite),
      random_p05 = as.numeric(quantile(finite, 0.05, type = 7)),
      random_p95 = as.numeric(quantile(finite, 0.95, type = 7)),
      clock_r2_within = real,
      percentile_rank = 100 * mean(finite <= real)
    )
    fwrite(
      data.table(draw = seq_len(n_random), r2_within = r2),
      file.path(out_dir, paste0("h1_random_reps_", dataset, "_", clock, ".tsv")),
      sep = "\t"
    )
    cat("random done", dataset, clock, "\n")
  }
}
fwrite(rbindlist(rand_sum), file.path(out_dir, "h1_random_signature.tsv"), sep = "\t")
}

cat("sensitivity\n")
qc <- fread(file.path(root, "results", "stage1", "expression_qc_summary.tsv"))
sens_rows <- list()
for (dataset in datasets) {
  samples <- read_samples(dataset)
  scores <- read_scores(dataset)
  setnames(scores, "clock", "clock_name")
  flagged <- strsplit(qc$flagged_gsm[qc$dataset == dataset], ",", fixed = TRUE)[[1]]
  for (clock in clocks) {
    sc <- scores[clock_name == clock & method == "lin"]
    named <- sc$value
    names(named) <- sc$gsm
    base <- control_frame(samples, named)
    # S1
    if (dataset == "GSE39445") {
      dat <- base[!is.na(base$hours_awake), ]
      fit <- fit_cosinor(
        dat,
        score ~ cos_t + sin_t + hours_awake + (1 | subject_id),
        score ~ hours_awake + (1 | subject_id)
      )
      fit$dataset <- dataset
      fit$clock <- clock
      fit$sensitivity <- "S1"
      fit$null_formula <- "score ~ hours_awake + (1 | subject_id)"
      sens_rows[[paste(dataset, clock, "S1")]] <- fit
    }
    # S2
    fit <- fit_cosinor(
      base,
      score ~ cos_t + sin_t + (1 + cos_t + sin_t | subject_id),
      score ~ 1 + (1 + cos_t + sin_t | subject_id)
    )
    fit$dataset <- dataset
    fit$clock <- clock
    fit$sensitivity <- "S2"
    fit$null_formula <- "score ~ 1 + (1 + cos_t + sin_t | subject_id)"
    sens_rows[[paste(dataset, clock, "S2")]] <- fit
    # S5
    dat <- base[!base$gsm %in% flagged, ]
    fit <- fit_cosinor(dat, forms$full, forms$null)
    fit$dataset <- dataset
    fit$clock <- clock
    fit$sensitivity <- "S5"
    fit$null_formula <- "score ~ 1 + (1 | subject_id)"
    fit$n_removed <- sum(base$gsm %in% flagged)
    sens_rows[[paste(dataset, clock, "S5")]] <- fit
    # S3
    pieces <- split(base, base$subject_id)
    s3 <- lapply(pieces, function(part) {
      rng <- max(part$clock_time_h) - min(part$clock_time_h)
      if (nrow(part) < 5 || rng < 16) return(NULL)
      ols <- lm(score ~ cos_t + sin_t, data = part)
      b <- coef(ols)
      a <- sqrt(b[["cos_t"]]^2 + b[["sin_t"]]^2)
      data.frame(
        subject_id = part$subject_id[[1]],
        n = nrow(part),
        clock_range = rng,
        two_a = 2 * a,
        peak_h = peak_0_24(b[["sin_t"]], b[["cos_t"]]),
        r2 = summary(ols)$r.squared
      )
    })
    s3 <- do.call(rbind, s3)
    ray <- rayleigh.test(circular(s3$peak_h, units = "hours", template = "clock24"))
    sens_rows[[paste(dataset, clock, "S3")]] <- data.frame(
      dataset = dataset,
      clock = clock,
      sensitivity = "S3",
      n_subjects = nrow(s3),
      n_samples = sum(s3$n),
      two_a = median(s3$two_a),
      two_a_q25 = as.numeric(quantile(s3$two_a, 0.25)),
      two_a_q75 = as.numeric(quantile(s3$two_a, 0.75)),
      peak_h = circ_mean_h(s3$peak_h),
      r2_within = median(s3$r2),
      lrt_p = NA_real_,
      rayleigh_stat = unname(ray$statistic),
      rayleigh_p = ray$p.value,
      null_formula = "per-subject OLS; mixed-model LRT not fit",
      stringsAsFactors = FALSE
    )
  }
}
# S4 GSE56931 design window
samples <- read_samples("GSE56931")
samples[, day := vapply(raw_characteristics, field_map, character(1), key = "day")]
samples[, hour := vapply(raw_characteristics, field_map, character(1), key = "hour")]
keep <- samples$condition == "control" & (
  (samples$day == "1" & samples$hour %in% c("8", "12", "16", "20")) |
    (samples$day == "2" & samples$hour %in% c("0", "4"))
)
scores <- read_scores("GSE56931")
setnames(scores, "clock", "clock_name")
for (clock in clocks) {
  sc <- scores[clock_name == clock & method == "lin"]
  named <- sc$value
  names(named) <- sc$gsm
  dat <- samples[keep & !is.na(clock_time_h)]
  dat <- as.data.frame(dat[, .(gsm, subject_id, clock_time_h, hours_awake)])
  dat$score <- named[dat$gsm]
  dat <- add_time_terms(dat[!is.na(dat$score), ])
  fit <- fit_cosinor(dat, forms$full, forms$null)
  fit$dataset <- "GSE56931"
  fit$clock <- clock
  fit$sensitivity <- "S4"
  fit$null_formula <- "score ~ 1 + (1 | subject_id)"
  sens_rows[[paste("GSE56931", clock, "S4")]] <- fit
}
sens <- rbindlist(sens_rows, fill = TRUE)
fwrite(sens, file.path(out_dir, "h1_sensitivity.tsv"), sep = "\t")

cat("figures\n")
for (dataset in datasets) {
  samples <- read_samples(dataset)
  scores <- read_scores(dataset)
  setnames(scores, "clock", "clock_name")
  for (clock in clocks) {
    sc <- scores[clock_name == clock & method == "lin"]
    named <- sc$value
    names(named) <- sc$gsm
    dat <- control_frame(samples, named)
    hit <- which(points$dataset == dataset & points$clock == clock & points$method == "lin")
    fit <- points[hit, ]
    grid <- seq(0, 24, length.out = 241)
    curve <- fit$beta_intercept + fit$beta_cos * cos(2 * pi * grid / 24) + fit$beta_sin * sin(2 * pi * grid / 24)
    png(file.path(fig_dir, paste0(dataset, "_", clock, "_lin.png")), width = 900, height = 640)
    plot(dat$clock_time_h, dat$score, xlab = "clock_time_h", ylab = "score",
         main = paste(dataset, clock, "lin"), pch = 16, col = "#9aa0a6")
    parts <- split(dat, dat$subject_id)
    for (part in parts) {
      part <- part[order(part$clock_time_h), ]
      lines(part$clock_time_h, part$score, col = "#4c78a8")
    }
    lines(grid, curve, col = "#e45756", lwd = 2)
    dev.off()
  }
}
cat("Part C done\n")
