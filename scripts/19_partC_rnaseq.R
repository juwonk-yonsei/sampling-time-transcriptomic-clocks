# Part C. RNA-seq development H1 and H2. Positive control gates later sections.

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
out_dir <- file.path(root, "results", "stage3", "rnaseq")
fig_dir <- file.path(root, "results", "stage3", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 20261015
n_boot <- 2000
n_random <- 1000
clocks <- c("RA_GTEx", "RA_Peters", "RA_default")
sig_of <- c(RA_GTEx = "GTExAge", RA_Peters = "Peters", RA_default = "DESeq2")
pc_genes <- c("PER1", "PER2", "PER3", "NR1D1", "NR1D2", "ARNTL", "CRY1", "CRY2", "DBP", "TEF")
candidate <- c("GSE113883", "GSE196728")
forms <- main_formula()
genelist <- get("genelist_all", asNamespace("RNAAgeCalc"))$blood
blood_labels <- c(
  "T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells",
  "B lineage", "Monocytic lineage", "Myeloid dendritic cells", "Neutrophils"
)

read_samples <- function(dataset) {
  fread(file.path(root, "results", "stage3", paste0(dataset, "_samples.tsv")))
}

read_log2 <- function(dataset) {
  path <- if (dataset == "GSE113883") {
    file.path(root, "results/stage3/expr/GSE113883_log2tpm1.tsv.gz")
  } else {
    file.path(root, "results/stage3/expr/GSE196728_log2fpkm1.tsv.gz")
  }
  dt <- fread(path)
  mat <- as.matrix(dt[, -1])
  rownames(mat) <- dt[[1]]
  storage.mode(mat) <- "double"
  mat
}

read_linear <- function(dataset) {
  if (dataset == "GSE113883") {
    dt <- fread(file.path(root, "results/stage3/expr/GSE113883_tpm.tsv.gz"))
  } else {
    dt <- fread(file.path(root, "results/stage3/expr/GSE196728_counts.tsv.gz"))
  }
  mat <- as.matrix(dt[, -1])
  rownames(mat) <- dt[[1]]
  storage.mode(mat) <- "double"
  mat
}

control_samples <- function(samples) {
  samples[condition == "control" & !is.na(clock_time_h)]
}

exprtype_of <- function(dataset) if (dataset == "GSE113883") "FPKM" else "counts"

cat("positive control\n")
pc_rows <- list()
pass <- character()
for (dataset in candidate) {
  samples <- control_samples(read_samples(dataset))
  mat <- read_log2(dataset)
  mat <- mat[, samples$gsm, drop = FALSE]
  peaks <- list()
  for (gene in pc_genes) {
    matrix_symbol <- if (gene %in% rownames(mat)) {
      gene
    } else if (gene == "ARNTL" && "BMAL1" %in% rownames(mat)) {
      "BMAL1"
    } else {
      NA_character_
    }
    if (is.na(matrix_symbol)) {
      peaks[[gene]] <- NA_real_
      pc_rows[[paste(dataset, gene)]] <- data.table(
        dataset = dataset, gene = gene, matrix_symbol = NA_character_, present = FALSE,
        peak_h = NA_real_, lrt_p = NA_real_
      )
      next
    }
    dat <- data.frame(
      score = as.numeric(mat[matrix_symbol, ]),
      subject_id = samples$subject_id,
      clock_time_h = samples$clock_time_h,
      stringsAsFactors = FALSE
    )
    dat <- dat[is.finite(dat$score), ]
    dat <- add_time_terms(dat)
    fit <- fit_cosinor(dat, forms$full, forms$null)
    peaks[[gene]] <- fit$peak_h
    pc_rows[[paste(dataset, gene)]] <- data.table(
      dataset = dataset, gene = gene, matrix_symbol = matrix_symbol, present = TRUE, peak_h = fit$peak_h,
      lrt_p = fit$lrt_p, two_a = fit$two_a, n_samples = fit$n_samples,
      n_subjects = fit$n_subjects, reml_converged = fit$reml_converged
    )
  }
  per1 <- peaks[["PER1"]]
  per3 <- peaks[["PER3"]]
  arntl <- peaks[["ARNTL"]]
  dist <- if (isTRUE(is.finite(per3)) && isTRUE(is.finite(arntl))) circ_dist_0_12(per3, arntl) else NA_real_
  a_prime <- is.finite(per1) && is.finite(per3) && per1 >= 4 && per1 <= 12 && per3 >= 4 && per3 <= 12
  b <- is.finite(dist) && dist >= 6
  decision <- if (isTRUE(a_prime) && isTRUE(b)) "pass" else "fail"
  if (decision == "pass") pass <- c(pass, dataset)
  pc_rows[[paste(dataset, "decision")]] <- data.table(
    dataset = dataset, gene = "DECISION", present = TRUE,
    peak_h = NA_real_, lrt_p = NA_real_, two_a = NA_real_,
    n_samples = nrow(samples), n_subjects = uniqueN(samples$subject_id),
    reml_converged = NA, criterion_a_prime = a_prime, criterion_b = b,
    per3_arntl_circ_dist_h = dist, decision = decision
  )
  cat(dataset, decision, "PER1", per1, "PER3", per3, "dist", dist, "\n")
}
fwrite(rbindlist(pc_rows, fill = TRUE), file.path(out_dir, "positive_control.tsv"), sep = "\t")
if (length(pass) == 0) {
  cat("no dataset passed; later sections skipped\n")
  quit(save = "no", status = 0)
}
datasets <- pass

cat("clock scores\n")
score_rows <- list()
for (dataset in datasets) {
  log2_mat <- read_log2(dataset)
  linear <- read_linear(dataset)
  keep <- apply(linear, 1, function(x) all(is.finite(x)))
  linear <- linear[keep, , drop = FALSE]
  log2_keep <- apply(log2_mat, 1, function(x) all(is.finite(x)))
  log2_mat <- log2_mat[log2_keep, , drop = FALSE]
  for (clock in clocks) {
    coefs <- genelist[[sig_of[[clock]]]]
    present <- intersect(setdiff(names(coefs), "(Intercept)"), rownames(log2_mat))
    used <- coefs[present]
    lin_score <- unname(coefs[["(Intercept)"]]) + as.numeric(crossprod(used, log2_mat[present, , drop = FALSE]))
    score_rows[[paste(dataset, clock, "lin")]] <- data.table(
      dataset = dataset, gsm = colnames(log2_mat), clock = clock, method = "lin",
      value = lin_score, n_genes_used = length(present),
      n_genes_signature = length(coefs) - 1L, n_genes_imputed = NA_integer_
    )
    msgs <- character()
    fit <- withCallingHandlers(
      if (clock == "RA_default") {
        predict_age(linear, tissue = "blood", exprtype = exprtype_of(dataset), idtype = "SYMBOL", stype = "all")
      } else {
        predict_age(linear, tissue = "blood", exprtype = exprtype_of(dataset), idtype = "SYMBOL", stype = "all", signature = sig_of[[clock]])
      },
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
    n_missing <- sum(!setdiff(names(coefs), "(Intercept)") %in% rownames(linear))
    score_rows[[paste(dataset, clock, "official")]] <- data.table(
      dataset = dataset, gsm = rownames(fit), clock = clock, method = "official",
      value = fit$RNAAge, n_genes_used = length(coefs) - 1L,
      n_genes_signature = length(coefs) - 1L, n_genes_imputed = n_missing,
      messages = paste(unique(msgs), collapse = " | ")
    )
    cat(dataset, clock, "lin genes", length(present), "imputed", n_missing, "\n")
  }
}
scores <- rbindlist(score_rows, fill = TRUE)
fwrite(scores, file.path(out_dir, "clock_scores.tsv"), sep = "\t")
dist_rows <- scores[, .(
  n = .N, mean = mean(value), sd = stats::sd(value), min = min(value), max = max(value),
  n_genes_used = n_genes_used[1], n_genes_imputed = n_genes_imputed[1]
), by = .(dataset, clock, method)]
fwrite(dist_rows, file.path(out_dir, "clock_score_distribution.tsv"), sep = "\t")

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
for (dataset in datasets) {
  samples <- control_samples(read_samples(dataset))
  for (clock in clocks) {
    lin <- scores[scores$dataset == dataset & scores$clock == clock & scores$method == "lin", .(gsm, value_lin = value)]
    off <- scores[scores$dataset == dataset & scores$clock == clock & scores$method == "official", .(gsm, value_official = value)]
    wide <- merge(merge(samples[, .(gsm, subject_id)], lin, by = "gsm"), off, by = "gsm")
    stats <- within_compare(wide)
    stats[, dataset := dataset]
    stats[, clock := clock]
    cmp_rows[[paste(dataset, clock)]] <- stats
  }
}
cmp <- rbindlist(cmp_rows)
cmp[, use_lin_only := pearson_within >= 0.99]
fwrite(cmp, file.path(out_dir, "official_lin_within.tsv"), sep = "\t")
print(cmp)

use_lin_only <- function(ds, ck) {
  hit <- cmp$use_lin_only[cmp$dataset == ds & cmp$clock == ck]
  isTRUE(length(hit) == 1L && isTRUE(hit[[1]]))
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

control_frame <- function(samples, score_named) {
  dat <- control_samples(samples)
  dat <- dat[, .(gsm, subject_id, clock_time_h)]
  dat[, score := score_named[gsm]]
  dat <- dat[is.finite(score)]
  add_time_terms(as.data.frame(dat))
}

named_score <- function(dataset, clock, method) {
  sc <- scores[scores$dataset == dataset & scores$clock == clock & scores$method == method]
  named <- sc$value
  names(named) <- sc$gsm
  named
}

for (dataset in datasets) {
  samples <- read_samples(dataset)
  subjects <- sort(unique(control_samples(samples)$subject_id))
  path <- file.path(out_dir, paste0("boot_draws_", dataset, ".tsv"))
  if (!file.exists(path)) {
    draws <- make_draws(subjects)
    fwrite(draws, path, sep = "\t")
  }
}

summarize_boot <- function(rep, dataset, clock, method) {
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
for (dataset in datasets) {
  samples <- read_samples(dataset)
  draws <- fread(file.path(out_dir, paste0("boot_draws_", dataset, ".tsv")))
  for (clock in clocks) {
    for (method in c("lin", "official")) {
      dat <- control_frame(samples, named_score(dataset, clock, method))
      fit <- fit_cosinor(dat, forms$full, forms$null)
      fit$dataset <- dataset
      fit$clock <- clock
      fit$method <- method
      fit$role <- if (method == "lin" || !isTRUE(use_lin_only(dataset, clock))) "primary" else "appendix_point"
      point_rows[[paste(dataset, clock, method)]] <- as.data.table(fit)
      do_boot <- method == "lin" || !isTRUE(use_lin_only(dataset, clock))
      if (!do_boot) next
      rep_path <- file.path(out_dir, paste0("h1_boot_reps_", dataset, "_", clock, "_", method, ".tsv"))
      if (file.exists(rep_path) && nrow(fread(rep_path)) == n_boot) {
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
      }
      rep <- rbindlist(lapply(reps, as.data.table))
      rep[, draw := seq_len(n_boot)]
      fwrite(rep[, .(draw, two_a, r2_within, two_a_over_sd_between, two_a_over_sd_resid, peak_h,
                     reml_converged, ml_converged, null_converged)], rep_path, sep = "\t")
      boot_sum_rows[[paste(dataset, clock, method)]] <- summarize_boot(rep, dataset, clock, method)
    }
  }
}
points <- rbindlist(point_rows, fill = TRUE)
fwrite(points, file.path(out_dir, "h1_point_estimates.tsv"), sep = "\t")
fwrite(rbindlist(boot_sum_rows), file.path(out_dir, "h1_bootstrap_ci.tsv"), sep = "\t")

if ("GSE113883" %in% datasets) {
  cat("random signatures GSE113883\n")
  rand_sum <- list()
  samples <- read_samples("GSE113883")
  mat_all <- read_log2("GSE113883")
  finite_gene <- apply(mat_all, 1, function(x) all(is.finite(x)))
  mat_all <- mat_all[finite_gene, , drop = FALSE]
  mu <- rowMeans(mat_all)
  dec <- as.integer(ceiling(rank(mu, ties.method = "average") / length(mu) * 10))
  dec[dec < 1] <- 1L
  dec[dec > 10] <- 10L
  names(dec) <- names(mu)
  for (clock in clocks) {
    coefs <- genelist[[sig_of[[clock]]]]
    clock_genes <- intersect(setdiff(names(coefs), "(Intercept)"), names(dec))
    coef <- coefs[clock_genes]
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
      score_all <- drop(crossprod(shuffled, mat_all[chosen, , drop = FALSE]))
      names(score_all) <- colnames(mat_all)
      dat <- control_frame(samples, score_all)
      fit <- fit_cosinor(dat, forms$full, forms$null)
      if (!isTRUE(fit$ml_converged) || !isTRUE(fit$null_converged)) n_fail <- n_fail + 1L
      r2[i] <- fit$r2_within
    }
    hit <- which(points$dataset == "GSE113883" & points$clock == clock & points$method == "lin")
    real <- points$r2_within[hit]
    finite <- r2[is.finite(r2)]
    rand_sum[[clock]] <- data.table(
      dataset = "GSE113883", clock = clock, method = "lin",
      n_random = n_random, n_finite = length(finite), n_not_converged = n_fail,
      random_median = median(finite),
      random_p05 = as.numeric(quantile(finite, 0.05, type = 7)),
      random_p95 = as.numeric(quantile(finite, 0.95, type = 7)),
      clock_r2_within = real,
      percentile_rank = 100 * mean(finite <= real)
    )
    fwrite(data.table(draw = seq_len(n_random), r2_within = r2),
           file.path(out_dir, paste0("h1_random_reps_GSE113883_", clock, ".tsv")), sep = "\t")
    cat("random", clock, "\n")
  }
  fwrite(rbindlist(rand_sum), file.path(out_dir, "h1_random_signature.tsv"), sep = "\t")
}

cat("sensitivity\n")
sens_rows <- list()
for (dataset in datasets) {
  samples <- read_samples(dataset)
  for (clock in clocks) {
    base <- control_frame(samples, named_score(dataset, clock, "lin"))
    fit <- fit_cosinor(
      base,
      score ~ cos_t + sin_t + (1 + cos_t + sin_t | subject_id),
      score ~ 1 + (1 + cos_t + sin_t | subject_id)
    )
    fit$dataset <- dataset
    fit$clock <- clock
    fit$sensitivity <- "S2"
    fit$null_formula <- "score ~ 1 + (1 + cos_t + sin_t | subject_id)"
    sens_rows[[paste(dataset, clock, "S2")]] <- as.data.table(fit)
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
    sens_rows[[paste(dataset, clock, "S3")]] <- data.table(
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
      null_formula = "per-subject OLS; mixed-model LRT not fit"
    )
  }
}
fwrite(rbindlist(sens_rows, fill = TRUE), file.path(out_dir, "h1_sensitivity.tsv"), sep = "\t")

cat("figures\n")
for (dataset in datasets) {
  samples <- read_samples(dataset)
  for (clock in clocks) {
    dat <- control_frame(samples, named_score(dataset, clock, "lin"))
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

cat("MCPcounter\n")
genes <- read.table(
  file.path(root, "data/processed/ref/MCPcounter_genes.txt"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,
  quote = "\"", comment.char = ""
)
marker_symbol <- genes[["HUGO symbols"]]
marker_pop <- genes[["Cell population"]]
score_rows <- list()
pop_rows <- list()
corr_rows <- list()
for (dataset in datasets) {
  mat <- read_log2(dataset)
  mat <- mat[apply(mat, 1, function(x) all(is.finite(x))), , drop = FALSE]
  est <- MCPcounter.estimate(mat, featuresType = "HUGO_symbols", genes = genes)
  score_dt <- as.data.table(t(est), keep.rownames = "gsm")
  fwrite(score_dt, file.path(out_dir, paste0("h2_mcp_scores_", dataset, ".tsv")), sep = "\t")
  for (pop in unique(marker_pop)) {
    marks <- marker_symbol[marker_pop == pop]
    score_rows[[paste(dataset, pop)]] <- data.table(
      dataset = dataset, population = pop, n_marker = length(marks),
      n_marker_in_matrix = sum(marks %in% rownames(mat))
    )
  }
  samples <- read_samples(dataset)
  base <- as.data.frame(control_samples(samples)[, .(gsm, subject_id, clock_time_h)])
  base <- merge(base, as.data.frame(score_dt), by = "gsm", sort = FALSE)
  base <- add_time_terms(base)
  cm <- cor(base[, blood_labels], use = "complete.obs")
  for (i in seq_along(blood_labels)) {
    for (j in seq_along(blood_labels)) {
      if (j <= i) next
      corr_rows[[paste(dataset, i, j)]] <- data.table(
        dataset = dataset, population_a = blood_labels[i], population_b = blood_labels[j],
        pearson_r = cm[i, j], abs_r_gt_0_9 = abs(cm[i, j]) > 0.9
      )
    }
  }
  for (pop in blood_labels) {
    dat <- base
    dat$score <- dat[[pop]]
    fit <- fit_cosinor(dat, forms$full, forms$null)
    fit$dataset <- dataset
    fit$population <- pop
    pop_rows[[paste(dataset, pop)]] <- as.data.table(fit)
  }
}
fwrite(rbindlist(score_rows), file.path(out_dir, "h2_marker_coverage.tsv"), sep = "\t")
fwrite(rbindlist(corr_rows), file.path(out_dir, "h2_population_correlations.tsv"), sep = "\t")
fwrite(rbindlist(pop_rows, fill = TRUE), file.path(out_dir, "h2_population_cosinor.tsv"), sep = "\t")
cat("Part C done\n")
