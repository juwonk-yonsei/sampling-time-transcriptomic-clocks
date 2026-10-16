#!/usr/bin/env Rscript
# Part C. GSE113883 RNA-seq development. Criterion (b) is not applied.
# GSE196728 is not read.

suppressPackageStartupMessages({
  library(lme4)
  library(data.table)
  library(circular)
  library(MCPcounter)
  library(RNAAgeCalc)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
root <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(root, "scripts", "lib_cosinor.R"))
out_dir <- file.path(root, "results", "stage4", "rnaseq")
fig_dir <- file.path(root, "results", "stage4", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 20261015
n_boot <- 2000
n_random <- 1000
dataset <- "GSE113883"
clocks <- c("RA_GTEx", "RA_Peters", "RA_default")
sig_of <- c(RA_GTEx = "GTExAge", RA_Peters = "Peters", RA_default = "DESeq2")
forms <- main_formula()
genelist <- get("genelist_all", asNamespace("RNAAgeCalc"))$blood
blood_labels <- c(
  "T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells",
  "B lineage", "Monocytic lineage", "Myeloid dendritic cells", "Neutrophils"
)

read_matrix <- function(path) {
  dt <- fread(path)
  mat <- as.matrix(dt[, -1])
  rownames(mat) <- dt[[1]]
  storage.mode(mat) <- "double"
  mat
}

samples <- fread(file.path(root, "results", "stage3", "GSE113883_samples.tsv"))
wake <- fread(file.path(root, "results", "stage4", "gse113883_hours_since_wake.tsv"))
samples <- merge(samples, wake[, .(gsm, hours_since_wake)], by = "gsm", sort = FALSE)
if (nrow(samples) != 153L) stop("sample count")
if (any(!is.finite(samples$clock_time_h)) || any(!is.finite(samples$hours_since_wake))) stop("missing time")
if (!all(samples$condition == "control")) stop("condition")

log2_mat <- read_matrix(file.path(root, "results", "stage3", "expr", "GSE113883_log2tpm1.tsv.gz"))
tpm <- read_matrix(file.path(root, "results", "stage3", "expr", "GSE113883_tpm.tsv.gz"))
log2_mat <- log2_mat[, samples$gsm, drop = FALSE]
tpm <- tpm[, samples$gsm, drop = FALSE]
if (!identical(rownames(log2_mat), rownames(tpm))) stop("gene order")

cat("QC\n")
flush.console()
mean_tpm <- rowMeans(tpm)
keep_qc <- is.finite(mean_tpm) & mean_tpm >= 1
qc_mat <- log2_mat[keep_qc, , drop = FALSE]
sample_median <- apply(qc_mat, 2, median, na.rm = TRUE)
center <- median(sample_median)
mad <- median(abs(sample_median - center))
flagged <- if (mad == 0) rep(FALSE, length(sample_median)) else abs(sample_median - center) > (3 * mad)
fwrite(data.table(
  gsm = colnames(qc_mat),
  median_log2tpm1 = sample_median,
  flag_3mad = flagged
), file.path(out_dir, "qc_sample_medians.tsv"), sep = "\t")
fwrite(data.table(
  dataset = dataset,
  value_scale = "log2(TPM+1)",
  gene_filter = "mean TPM >= 1",
  n_genes = sum(keep_qc),
  n_samples = ncol(qc_mat),
  median_of_sample_medians = center,
  mad_of_sample_medians = mad,
  mad_constant = "none",
  n_flagged_3mad = sum(flagged),
  flagged_gsm = paste(colnames(qc_mat)[flagged], collapse = ",")
), file.path(out_dir, "qc_summary.tsv"), sep = "\t")
cat("QC genes", sum(keep_qc), "flagged", sum(flagged), "\n")

cat("clock scores\n")
flush.console()
finite_lin <- apply(log2_mat, 1, function(x) all(is.finite(x)))
finite_tpm <- apply(tpm, 1, function(x) all(is.finite(x)))
lin_mat <- log2_mat[finite_lin, , drop = FALSE]
linear <- tpm[finite_tpm, , drop = FALSE]
score_rows <- list()
for (clock in clocks) {
  coefs <- genelist[[sig_of[[clock]]]]
  present <- intersect(setdiff(names(coefs), "(Intercept)"), rownames(lin_mat))
  used <- coefs[present]
  lin_score <- unname(coefs[["(Intercept)"]]) + as.numeric(crossprod(used, lin_mat[present, , drop = FALSE]))
  score_rows[[paste(clock, "lin")]] <- data.table(
    dataset = dataset, gsm = colnames(lin_mat), clock = clock, method = "lin",
    value = lin_score, n_genes_used = length(present),
    n_genes_signature = length(coefs) - 1L, n_genes_imputed = NA_integer_,
    exprtype_argument = NA_character_, value_kind = "log2(TPM+1)"
  )
  msgs <- character()
  fit <- withCallingHandlers(
    if (clock == "RA_default") {
      predict_age(linear, tissue = "blood", exprtype = "FPKM", idtype = "SYMBOL", stype = "all")
    } else {
      predict_age(linear, tissue = "blood", exprtype = "FPKM", idtype = "SYMBOL", stype = "all", signature = sig_of[[clock]])
    },
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  n_missing <- sum(!setdiff(names(coefs), "(Intercept)") %in% rownames(linear))
  score_rows[[paste(clock, "official")]] <- data.table(
    dataset = dataset, gsm = rownames(fit), clock = clock, method = "official",
    value = fit$RNAAge, n_genes_used = length(coefs) - 1L - n_missing,
    n_genes_signature = length(coefs) - 1L, n_genes_imputed = n_missing,
    exprtype_argument = "FPKM", value_kind = "TPM",
    messages = paste(unique(msgs), collapse = " | ")
  )
  cat(clock, "lin genes", length(present), "imputed", n_missing, "\n")
}
scores <- rbindlist(score_rows, fill = TRUE)
fwrite(scores, file.path(out_dir, "clock_scores.tsv"), sep = "\t")
fwrite(scores[, .(
  n = .N, mean = mean(value), sd = stats::sd(value), min = min(value), max = max(value),
  n_genes_used = n_genes_used[1], n_genes_signature = n_genes_signature[1],
  n_genes_imputed = n_genes_imputed[1]
), by = .(dataset, clock, method)], file.path(out_dir, "clock_score_distribution.tsv"), sep = "\t")

within_compare <- function(wide) {
  wide[, dev_lin := value_lin - mean(value_lin), by = subject_id]
  wide[, dev_off := value_official - mean(value_official), by = subject_id]
  wide[, diff := value_official - value_lin]
  wide[, diff_dev := diff - mean(diff), by = subject_id]
  data.table(
    pearson_within = cor(wide$dev_lin, wide$dev_off, use = "complete.obs"),
    sd_within_official_minus_lin = stats::sd(wide$diff_dev),
    n = nrow(wide)
  )
}
cmp_rows <- list()
for (clock in clocks) {
  lin_idx <- which(scores$clock == clock & scores$method == "lin")
  off_idx <- which(scores$clock == clock & scores$method == "official")
  lin <- scores[lin_idx, .(gsm, value_lin = value)]
  off <- scores[off_idx, .(gsm, value_official = value)]
  wide <- merge(merge(samples[, .(gsm, subject_id)], lin, by = "gsm"), off, by = "gsm")
  stats <- within_compare(wide)
  stats[, dataset := dataset]
  stats[, clock := clock]
  cmp_rows[[clock]] <- stats
}
cmp <- rbindlist(cmp_rows)
cmp[, use_lin_only := pearson_within >= 0.99]
fwrite(cmp, file.path(out_dir, "official_lin_within.tsv"), sep = "\t")
print(cmp)

named_score <- function(clock, method) {
  sc_idx <- which(scores$clock == clock & scores$method == method)
  sc <- scores[sc_idx]
  named <- sc$value
  names(named) <- sc$gsm
  named
}

control_frame <- function(score_named, time_col = "clock_time_h") {
  dat <- samples[, .(gsm, subject_id, clock_time_h = get(time_col))]
  dat[, score := score_named[gsm]]
  dat <- dat[is.finite(score) & is.finite(clock_time_h)]
  add_time_terms(as.data.frame(dat))
}

make_draws <- function(subjects) {
  set.seed(seed)
  rows <- vector("list", n_boot)
  for (b in seq_len(n_boot)) {
    picked <- sample(subjects, length(subjects), replace = TRUE)
    times <- ave(seq_along(picked), picked, FUN = seq_along)
    rows[[b]] <- data.table(
      draw = b, pick = seq_along(picked), subject_id = picked,
      new_subject_id = paste(picked, times, sep = "#")
    )
  }
  rbindlist(rows)
}

subjects <- sort(unique(samples$subject_id))
draw_path <- file.path(out_dir, "boot_draws_GSE113883.tsv")
if (!file.exists(draw_path)) {
  fwrite(make_draws(subjects), draw_path, sep = "\t")
}
draws <- fread(draw_path)

summarize_boot <- function(rep, clock, method) {
  ok <- rep$reml_converged & rep$ml_converged & rep$null_converged
  summ <- data.table(
    dataset = dataset, clock = clock, method = method,
    n_boot = nrow(rep), n_converged = sum(ok), n_not_converged = sum(!ok)
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
for (clock in clocks) {
  for (method in c("lin", "official")) {
    dat <- control_frame(named_score(clock, method))
    fit <- fit_cosinor(dat, forms$full, forms$null)
    fit$dataset <- dataset
    fit$clock <- clock
    fit$method <- method
    lin_only <- isTRUE(cmp$use_lin_only[cmp$clock == clock][[1]])
    fit$role <- if (method == "lin" || !lin_only) "primary" else "appendix_point"
    point_rows[[paste(clock, method)]] <- as.data.table(fit)
    if (method == "official" && lin_only) next
    rep_path <- file.path(out_dir, paste0("h1_boot_reps_", clock, "_", method, ".tsv"))
    if (file.exists(rep_path) && nrow(fread(rep_path, nrows = 1)) > 0 && nrow(fread(rep_path)) == n_boot) {
      boot_sum_rows[[paste(clock, method)]] <- summarize_boot(fread(rep_path), clock, method)
      next
    }
    cat("bootstrap", clock, method, "\n")
    flush.console()
    reps <- vector("list", n_boot)
    for (b in seq_len(n_boot)) {
      draw_idx <- which(draws$draw == b)
      one <- draws[draw_idx]
      boot_dat <- apply_draw(dat, one$subject_id, one$new_subject_id)
      reps[[b]] <- fit_cosinor(boot_dat, forms$full, forms$null)
      if (b %% 200L == 0L) {
        cat(" ", b, "\n")
        flush.console()
      }
    }
    rep <- rbindlist(lapply(reps, as.data.table))
    rep[, draw := seq_len(n_boot)]
    fwrite(rep[, .(draw, two_a, r2_within, two_a_over_sd_between, two_a_over_sd_resid, peak_h,
                   reml_converged, ml_converged, null_converged)], rep_path, sep = "\t")
    boot_sum_rows[[paste(clock, method)]] <- summarize_boot(rep, clock, method)
  }
}
points <- rbindlist(point_rows, fill = TRUE)
fwrite(points, file.path(out_dir, "h1_point_estimates.tsv"), sep = "\t")
fwrite(rbindlist(boot_sum_rows), file.path(out_dir, "h1_bootstrap_ci.tsv"), sep = "\t")

cat("random signatures\n")
flush.console()
mu <- rowMeans(lin_mat)
dec <- as.integer(ceiling(rank(mu, ties.method = "average") / length(mu) * 10))
dec[dec < 1] <- 1L
dec[dec > 10] <- 10L
names(dec) <- names(mu)
rand_sum <- list()
for (clock in clocks) {
  rep_path <- file.path(out_dir, paste0("h1_random_reps_", clock, ".tsv"))
  coefs <- genelist[[sig_of[[clock]]]]
  clock_genes <- intersect(setdiff(names(coefs), "(Intercept)"), names(dec))
  coef <- coefs[clock_genes]
  if (file.exists(rep_path) && nrow(fread(rep_path)) == n_random) {
    r2 <- fread(rep_path)$r2_within
    n_fail <- sum(!is.finite(r2))
  } else {
    non_clock <- setdiff(names(dec), clock_genes)
    pools <- split(non_clock, dec[non_clock])
    need <- table(dec[clock_genes])
    set.seed(seed)
    r2 <- rep(NA_real_, n_random)
    n_fail <- 0L
    for (i in seq_len(n_random)) {
      chosen <- character(length(clock_genes))
      for (d in names(need)) {
        pick <- sample(pools[[d]], as.integer(need[[d]]))
        idx <- which(as.character(dec[clock_genes]) == d)
        chosen[idx] <- pick
      }
      shuffled <- sample(unname(coef))
      score_all <- drop(crossprod(shuffled, lin_mat[chosen, , drop = FALSE]))
      names(score_all) <- colnames(lin_mat)
      dat <- control_frame(score_all)
      fit <- fit_cosinor(dat, forms$full, forms$null)
      if (!isTRUE(fit$ml_converged) || !isTRUE(fit$null_converged)) n_fail <- n_fail + 1L
      r2[i] <- fit$r2_within
      if (i %% 100L == 0L) {
        cat(" random", clock, i, "\n")
        flush.console()
      }
    }
    fwrite(data.table(draw = seq_len(n_random), r2_within = r2), rep_path, sep = "\t")
  }
  real <- points$r2_within[points$clock == clock & points$method == "lin"]
  finite <- r2[is.finite(r2)]
  rand_sum[[clock]] <- data.table(
    dataset = dataset, clock = clock, method = "lin",
    n_random = n_random, n_finite = length(finite), n_not_converged = n_fail,
    random_median = median(finite),
    random_p05 = as.numeric(quantile(finite, 0.05, type = 7)),
    random_p95 = as.numeric(quantile(finite, 0.95, type = 7)),
    clock_r2_within = real,
    percentile_rank = 100 * mean(finite <= real)
  )
}
fwrite(rbindlist(rand_sum), file.path(out_dir, "h1_random_signature.tsv"), sep = "\t")

cat("sensitivity\n")
flush.console()
sens_rows <- list()
for (clock in clocks) {
  base <- control_frame(named_score(clock, "lin"))
  fit <- fit_cosinor(
    base,
    score ~ cos_t + sin_t + (1 + cos_t + sin_t | subject_id),
    score ~ 1 + (1 + cos_t + sin_t | subject_id)
  )
  fit$dataset <- dataset
  fit$clock <- clock
  fit$sensitivity <- "S2"
  fit$null_formula <- "score ~ 1 + (1 + cos_t + sin_t | subject_id)"
  sens_rows[[paste(clock, "S2")]] <- as.data.table(fit)
  pieces <- split(base, base$subject_id)
  s3 <- lapply(pieces, function(part) {
    rng <- max(part$clock_time_h) - min(part$clock_time_h)
    if (nrow(part) < 5 || rng < 16) return(NULL)
    ols <- lm(score ~ cos_t + sin_t, data = part)
    b <- coef(ols)
    a <- sqrt(b[["cos_t"]]^2 + b[["sin_t"]]^2)
    data.frame(
      subject_id = part$subject_id[[1]], n = nrow(part), clock_range = rng,
      two_a = 2 * a, peak_h = peak_0_24(b[["sin_t"]], b[["cos_t"]]),
      r2 = summary(ols)$r.squared
    )
  })
  s3 <- do.call(rbind, s3)
  ray <- rayleigh.test(circular(s3$peak_h, units = "hours", template = "clock24"))
  sens_rows[[paste(clock, "S3")]] <- data.table(
    dataset = dataset, clock = clock, sensitivity = "S3",
    n_subjects = nrow(s3), n_samples = sum(s3$n),
    two_a = median(s3$two_a),
    two_a_q25 = as.numeric(quantile(s3$two_a, 0.25)),
    two_a_q75 = as.numeric(quantile(s3$two_a, 0.75)),
    peak_h = circ_mean_h(s3$peak_h),
    r2_within = median(s3$r2),
    lrt_p = NA_real_,
    rayleigh_stat = unname(ray$statistic),
    rayleigh_p = ray$p.value,
    null_formula = "per-subject OLS; mixed-model LRT not fit",
    clock_range_definition = "max(clock_time_h) - min(clock_time_h)"
  )
}
fwrite(rbindlist(sens_rows, fill = TRUE), file.path(out_dir, "h1_sensitivity.tsv"), sep = "\t")

cat("figures\n")
flush.console()
for (clock in clocks) {
  dat <- control_frame(named_score(clock, "lin"))
  fit_idx <- which(points$clock == clock & points$method == "lin")
  fit <- points[fit_idx]
  grid <- seq(0, 24, length.out = 241)
  curve <- fit$beta_intercept + fit$beta_cos * cos(2 * pi * grid / 24) + fit$beta_sin * sin(2 * pi * grid / 24)
  png(file.path(fig_dir, paste0("GSE113883_", clock, "_lin.png")), width = 900, height = 640)
  plot(dat$clock_time_h, dat$score, xlab = "clock_time_h", ylab = "score",
       main = paste("GSE113883", clock, "lin"), pch = 16, col = "#9aa0a6")
  parts <- split(dat, dat$subject_id)
  for (part in parts) {
    part <- part[order(part$clock_time_h), ]
    lines(part$clock_time_h, part$score, col = "#4c78a8")
  }
  lines(grid, curve, col = "#e45756", lwd = 2)
  dev.off()
}

cat("MCPcounter\n")
flush.console()
genes <- read.table(
  file.path(root, "data/processed/ref/MCPcounter_genes.txt"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,
  quote = "\"", comment.char = ""
)
marker_symbol <- genes[["HUGO symbols"]]
marker_pop <- genes[["Cell population"]]
est <- MCPcounter.estimate(lin_mat, featuresType = "HUGO_symbols", genes = genes)
score_dt <- as.data.table(t(est), keep.rownames = "gsm")
fwrite(score_dt, file.path(out_dir, "h2_mcp_scores_GSE113883.tsv"), sep = "\t")
marker_rows <- list()
for (pop in unique(marker_pop)) {
  marks <- marker_symbol[marker_pop == pop]
  marker_rows[[pop]] <- data.table(
    dataset = dataset, population = pop, n_marker = length(marks),
    n_marker_in_matrix = sum(marks %in% rownames(lin_mat))
  )
}
fwrite(rbindlist(marker_rows), file.path(out_dir, "h2_marker_coverage.tsv"), sep = "\t")
base <- as.data.frame(samples[, .(gsm, subject_id, clock_time_h)])
base <- merge(base, as.data.frame(score_dt), by = "gsm", sort = FALSE)
base <- add_time_terms(base)
cm <- cor(base[, blood_labels], use = "complete.obs")
corr_rows <- list()
for (i in seq_along(blood_labels)) {
  for (j in seq_along(blood_labels)) {
    if (j <= i) next
    corr_rows[[paste(i, j)]] <- data.table(
      dataset = dataset, population_a = blood_labels[i], population_b = blood_labels[j],
      pearson_r = cm[i, j], abs_r_gt_0_9 = abs(cm[i, j]) > 0.9
    )
  }
}
fwrite(rbindlist(corr_rows), file.path(out_dir, "h2_population_correlations.tsv"), sep = "\t")
pop_rows <- list()
for (pop in blood_labels) {
  dat <- base
  dat$score <- dat[[pop]]
  fit <- fit_cosinor(dat, forms$full, forms$null)
  fit$dataset <- dataset
  fit$population <- pop
  pop_rows[[pop]] <- as.data.table(fit)
}
fwrite(rbindlist(pop_rows, fill = TRUE), file.path(out_dir, "h2_population_cosinor.tsv"), sep = "\t")

cat("projection\n")
flush.console()
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
safe <- safe_name(blood_labels)
names(safe) <- blood_labels
for (label in blood_labels) base[[safe[[label]]]] <- base[[label]]
adj_full <- as.formula(paste(
  "score ~ cos_t + sin_t +", paste(safe, collapse = " + "), "+ (1 | subject_id)"
))
reml_beta <- function(data, formula) {
  fit <- quiet_lmer(formula, data, TRUE)
  if (inherits(fit, "error") || !optimizer_converged(fit)) {
    return(c(cos = NA_real_, sin = NA_real_, converged = 0))
  }
  beta <- lme4::fixef(fit)
  c(cos = unname(beta["cos_t"]), sin = unname(beta["sin_t"]), converged = 1)
}
projection <- function(ru, ra) {
  ru2 <- sum(ru * ru)
  if (!is.finite(ru2) || ru2 == 0) {
    return(c(p_med = NA_real_, two_a_orth = NA_real_, ru_norm = NA_real_))
  }
  scale <- sum(ra * ru) / ru2
  resid <- ra - scale * ru
  c(p_med = 1 - scale, two_a_orth = 2 * sqrt(sum(resid * resid)), ru_norm = sqrt(ru2))
}
proj_points <- list()
proj_ci <- list()
within_rows <- list()
for (clock in clocks) {
  dat <- base
  dat$score <- unname(named_score(clock, "lin")[dat$gsm])
  if (anyNA(dat$score)) stop("missing lin score")
  unadj <- fit_cosinor(dat, forms$full, forms$null)
  adj_beta <- reml_beta(dat, adj_full)
  ru <- c(unadj$beta_cos, unadj$beta_sin)
  ra <- c(unname(adj_beta["cos"]), unname(adj_beta["sin"]))
  proj <- projection(ru, ra)
  proj_points[[clock]] <- data.table(
    dataset = dataset, clock = clock, method = "lin",
    n_subjects = unadj$n_subjects, n_samples = unadj$n_samples,
    lrt_p_unadj = unadj$lrt_p,
    ru_cos = ru[1], ru_sin = ru[2], ru_norm = unname(proj["ru_norm"]),
    ra_cos = ra[1], ra_sin = ra[2],
    p_med = unname(proj["p_med"]), two_a_orth = unname(proj["two_a_orth"]),
    reml_converged_unadj = unadj$reml_converged,
    reml_converged_adj = isTRUE(adj_beta["converged"] == 1),
    two_a_unadj = unadj$two_a
  )
  dev_clock <- ave(dat$score, dat$subject_id, FUN = function(x) x - mean(x))
  for (pop in blood_labels) {
    dev_pop <- ave(dat[[pop]], dat$subject_id, FUN = function(x) x - mean(x))
    within_rows[[paste(clock, pop)]] <- data.table(
      dataset = dataset, clock = clock, population = pop,
      n_samples = nrow(dat), n_subjects = length(unique(dat$subject_id)),
      pearson_within = cor(dev_clock, dev_pop, use = "complete.obs")
    )
  }
  rep_path <- file.path(out_dir, paste0("projection_reps_", clock, ".tsv"))
  if (file.exists(rep_path) && nrow(fread(rep_path)) == n_boot) {
    rep <- fread(rep_path)
  } else {
    cat(" projection", clock, "\n")
    flush.console()
    reps <- vector("list", n_boot)
    for (b in seq_len(n_boot)) {
      draw_idx <- which(draws$draw == b)
      one <- draws[draw_idx]
      boot_dat <- apply_draw(dat, one$subject_id, one$new_subject_id)
      u <- reml_beta(boot_dat, forms$full)
      a <- reml_beta(boot_dat, adj_full)
      ok <- isTRUE(u["converged"] == 1) && isTRUE(a["converged"] == 1)
      pr <- if (ok) projection(c(u["cos"], u["sin"]), c(a["cos"], a["sin"])) else c(p_med = NA_real_, two_a_orth = NA_real_, ru_norm = NA_real_)
      reps[[b]] <- data.table(
        draw = b, reml_converged_both = ok,
        p_med = unname(pr["p_med"]), two_a_orth = unname(pr["two_a_orth"])
      )
      if (b %% 200L == 0L) {
        cat("  ", b, "\n")
        flush.console()
      }
    }
    rep <- rbindlist(reps)
    fwrite(rep, rep_path, sep = "\t")
  }
  p_ci <- percentile_ci(rep$p_med)
  a_ci <- percentile_ci(rep$two_a_orth)
  proj_ci[[clock]] <- data.table(
    dataset = dataset, clock = clock, method = "lin",
    n_boot = n_boot,
    n_converged_both_reml = sum(rep$reml_converged_both),
    n_finite_p_med = sum(is.finite(rep$p_med)),
    p_med = unname(proj["p_med"]),
    p_med_p025 = p_ci[1], p_med_p975 = p_ci[2],
    two_a_orth = unname(proj["two_a_orth"]),
    two_a_orth_p025 = a_ci[1], two_a_orth_p975 = a_ci[2],
    lrt_p_unadj = unadj$lrt_p, ru_norm = unname(proj["ru_norm"])
  )
}
fwrite(rbindlist(proj_points), file.path(out_dir, "projection_point.tsv"), sep = "\t")
fwrite(rbindlist(proj_ci), file.path(out_dir, "projection_ci.tsv"), sep = "\t")
fwrite(rbindlist(within_rows), file.path(out_dir, "clock_composition_within.tsv"), sep = "\t")

cat("G2 clocks\n")
flush.console()
window_means <- function(dat, values) {
  morning <- dat$clock_time_h >= 7 & dat$clock_time_h <= 9
  evening <- dat$clock_time_h >= 19 & dat$clock_time_h <= 21
  rows <- list()
  for (sid in sort(unique(dat$subject_id))) {
    m <- dat$subject_id == sid & morning & is.finite(values)
    e <- dat$subject_id == sid & evening & is.finite(values)
    if (!any(m) || !any(e)) next
    rows[[sid]] <- data.table(
      subject_id = sid, morning = mean(values[m]), evening = mean(values[e]),
      n_morning = sum(m), n_evening = sum(e)
    )
  }
  rbindlist(rows)
}
paired_summary <- function(wm) {
  diff <- wm$morning - wm$evening
  tt <- t.test(wm$morning, wm$evening, paired = TRUE)
  list(
    n = nrow(wm), mean_diff = unname(tt$estimate), sd_diff = stats::sd(diff),
    p = unname(tt$p.value), ci_low = unname(tt$conf.int[1]), ci_high = unname(tt$conf.int[2])
  )
}
g2_rows <- list()
dat_g2 <- as.data.frame(samples[, .(gsm, subject_id, clock_time_h)])
for (clock in clocks) {
  named <- named_score(clock, "lin")
  wm <- window_means(dat_g2, unname(named[dat_g2$gsm]))
  sm <- paired_summary(wm)
  g2_rows[[clock]] <- data.table(
    dataset = dataset, outcome = clock, outcome_type = "clock_lin",
    window_morning = "07:00-09:00 inclusive", window_evening = "19:00-21:00 inclusive",
    contrast = "morning - evening",
    n_subjects = sm$n, mean_diff = sm$mean_diff, sd_diff = sm$sd_diff,
    ci95_low = sm$ci_low, ci95_high = sm$ci_high, p = sm$p
  )
}
fwrite(rbindlist(g2_rows), file.path(out_dir, "g2_clocks.tsv"), sep = "\t")

cat("time axes\n")
flush.console()
ml_aic <- function(data, formula) {
  fit <- quiet_lmer(formula, data, FALSE)
  if (inherits(fit, "error")) return(c(aic = NA_real_, converged = 0))
  c(aic = AIC(fit), converged = as.numeric(optimizer_converged(fit)))
}
gene_named <- function(gene) {
  symbol <- gene
  if (!gene %in% rownames(log2_mat) && gene == "BMAL1" && "BMAL1" %in% rownames(log2_mat)) symbol <- "BMAL1"
  if (!symbol %in% rownames(log2_mat)) stop(paste("missing", gene))
  vals <- as.numeric(log2_mat[symbol, ])
  names(vals) <- colnames(log2_mat)
  list(values = vals, matrix_symbol = symbol)
}
outcomes <- list()
for (clock in clocks) {
  outcomes[[clock]] <- list(values = named_score(clock, "lin"), matrix_symbol = clock)
}
for (gene in c("PER1", "PER3", "BMAL1")) outcomes[[gene]] <- gene_named(gene)
axis_rows <- list()
for (outcome in names(outcomes)) {
  got <- outcomes[[outcome]]
  for (spec in list(c(axis = "T_clock", col = "clock_time_h"), c(axis = "T_wake", col = "hours_since_wake"))) {
    d <- as.data.frame(samples[, .(gsm, subject_id)])
    d$clock_time_h <- samples[[spec[["col"]]]]
    d$score <- unname(got$values[d$gsm])
    d <- add_time_terms(d)
    fit <- fit_cosinor(d, forms$full, forms$null)
    aic <- ml_aic(d, forms$full)
    axis_rows[[paste(outcome, spec[["axis"]])]] <- data.table(
      dataset = dataset, outcome = outcome, matrix_symbol = got$matrix_symbol,
      axis = spec[["axis"]], time_col = spec[["col"]],
      n_subjects = fit$n_subjects, n_samples = fit$n_samples,
      two_a = fit$two_a, peak_h = fit$peak_h, lrt_p = fit$lrt_p, r2_within = fit$r2_within,
      aic_ml = unname(aic["aic"]), ml_converged = unname(aic["converged"]) == 1,
      reml_converged = fit$reml_converged
    )
  }
}
axes <- rbindlist(axis_rows)
ref_idx <- which(axes$axis == "T_clock")
ref <- axes[ref_idx, .(outcome, aic_clock = aic_ml)]
axes <- merge(axes, ref, by = "outcome", sort = FALSE)
axes[, delta_aic_vs_T_clock := aic_ml - aic_clock]
fwrite(axes, file.path(out_dir, "time_axis_comparison.tsv"), sep = "\t")
fwrite(data.table(
  dataset = dataset, axis = "T_dlmo", performed = FALSE,
  reason = "salivary_dlmo_circadian_time values are outside 18-26",
  n_samples_other_axes = 153L
), file.path(out_dir, "t_dlmo_status.tsv"), sep = "\t")
cat("Part C done\n")
