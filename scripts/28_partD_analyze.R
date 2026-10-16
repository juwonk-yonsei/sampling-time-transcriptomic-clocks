#!/usr/bin/env Rscript
# Part D. Module-clock cosinor on development control samples.
# Bootstrap and projection are limited to the two composite clocks.

suppressPackageStartupMessages({
  library(lme4)
  library(data.table)
  library(RNAAgeCalc)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
root <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(root, "scripts", "lib_cosinor.R"))
out_dir <- file.path(root, "results", "stage4", "module_clocks")
scores <- fread(file.path(out_dir, "clock_scores.tsv"))
forms <- main_formula()
n_boot <- 2000
datasets <- c("GSE39445", "GSE56931", "GSE113883")
composites <- c("TA_chrono", "TA_mort")
blood_labels <- c(
  "T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells",
  "B lineage", "Monocytic lineage", "Myeloid dendritic cells", "Neutrophils"
)
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
genelist <- get("genelist_all", asNamespace("RNAAgeCalc"))$blood

read_samples <- function(dataset) {
  path <- if (dataset == "GSE113883") {
    file.path(root, "results", "stage3", "GSE113883_samples.tsv")
  } else {
    file.path(root, "data", "processed", paste0(dataset, "_samples.tsv"))
  }
  samples <- fread(path)
  keep <- which(samples$condition == "control" & is.finite(samples$clock_time_h))
  samples[keep, .(gsm, subject_id, clock_time_h)]
}

read_log2 <- function(dataset) {
  path <- if (dataset == "GSE113883") {
    file.path(root, "results", "stage3", "expr", "GSE113883_log2tpm1.tsv.gz")
  } else {
    file.path(root, "data", "processed", paste0(dataset, "_expr_gene_log2.tsv.gz"))
  }
  dt <- fread(path)
  mat <- as.matrix(dt[, -1])
  rownames(mat) <- dt[[1]]
  storage.mode(mat) <- "double"
  mat
}

peters_lin <- function(dataset, gsms) {
  mat <- read_log2(dataset)
  mat <- mat[, gsms, drop = FALSE]
  keep <- apply(mat, 1, function(x) all(is.finite(x)))
  mat <- mat[keep, , drop = FALSE]
  coefs <- genelist[["Peters"]]
  present <- intersect(setdiff(names(coefs), "(Intercept)"), rownames(mat))
  used <- coefs[present]
  score <- unname(coefs[["(Intercept)"]]) + as.numeric(crossprod(used, mat[present, , drop = FALSE]))
  names(score) <- colnames(mat)
  score
}

draw_path <- function(dataset) {
  if (dataset == "GSE113883") {
    file.path(root, "results", "stage4", "rnaseq", "boot_draws_GSE113883.tsv")
  } else {
    file.path(root, "results", "stage2", paste0("boot_draws_", dataset, ".tsv"))
  }
}

mcp_path <- function(dataset) {
  if (dataset == "GSE113883") {
    file.path(root, "results", "stage4", "rnaseq", "h2_mcp_scores_GSE113883.tsv")
  } else {
    file.path(root, "results", "stage2", paste0("h2_mcp_scores_", dataset, ".tsv"))
  }
}

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

point_rows <- list()
boot_rows <- list()
within_rows <- list()
for (dataset in datasets) {
  cat(dataset, "\n")
  flush.console()
  samples <- read_samples(dataset)
  peters <- peters_lin(dataset, samples$gsm)
  draws <- fread(draw_path(dataset))
  subjects <- sort(unique(samples$subject_id))
  draw_subjects <- sort(unique(draws$subject_id))
  if (!identical(subjects, draw_subjects)) stop(paste("draw subjects differ", dataset))
  clock_names <- unique(scores$clock)
  for (clock in clock_names) {
    idx <- which(scores$dataset == dataset & scores$clock == clock)
    sc <- scores[idx]
    named <- sc$value
    names(named) <- sc$gsm
    dat <- as.data.frame(samples)
    dat$score <- unname(named[dat$gsm])
    if (anyNA(dat$score)) stop(paste("missing score", dataset, clock))
    dat <- add_time_terms(dat)
    fit <- fit_cosinor(dat, forms$full, forms$null)
    fit$dataset <- dataset
    fit$clock <- clock
    fit$outcome <- sc$outcome[1]
    fit$module <- sc$module[1]
    point_rows[[paste(dataset, clock)]] <- as.data.table(fit)
    dev_clock <- ave(dat$score, dat$subject_id, FUN = function(x) x - mean(x))
    dev_peters <- ave(unname(peters[dat$gsm]), dat$subject_id, FUN = function(x) x - mean(x))
    if (clock %in% composites) {
      within_rows[[paste(dataset, clock)]] <- data.table(
        dataset = dataset, clock = clock, other = "RA_Peters_lin",
        n_samples = nrow(dat), n_subjects = length(unique(dat$subject_id)),
        pearson_within = cor(dev_clock, dev_peters, use = "complete.obs")
      )
      rep_path <- file.path(out_dir, paste0("h1_boot_reps_", dataset, "_", clock, ".tsv"))
      if (file.exists(rep_path) && nrow(fread(rep_path)) == n_boot) {
        rep <- fread(rep_path)
      } else {
        cat(" bootstrap", clock, "\n")
        flush.console()
        reps <- vector("list", n_boot)
        for (b in seq_len(n_boot)) {
          draw_idx <- which(draws$draw == b)
          one <- draws[draw_idx]
          boot_dat <- apply_draw(dat, one$subject_id, one$new_subject_id)
          reps[[b]] <- fit_cosinor(boot_dat, forms$full, forms$null)
        }
        rep <- rbindlist(lapply(reps, as.data.table))
        rep[, draw := seq_len(n_boot)]
        fwrite(rep[, .(draw, two_a, r2_within, two_a_over_sd_between, two_a_over_sd_resid, peak_h,
                       reml_converged, ml_converged, null_converged)], rep_path, sep = "\t")
      }
      ok <- rep$reml_converged & rep$ml_converged & rep$null_converged
      summ <- data.table(
        dataset = dataset, clock = clock, n_boot = nrow(rep),
        n_converged = sum(ok), n_not_converged = sum(!ok)
      )
      for (metric in c("two_a", "r2_within", "two_a_over_sd_between", "two_a_over_sd_resid")) {
        ci <- percentile_ci(rep[[metric]][ok])
        summ[, paste0(metric, "_p025") := ci[1]]
        summ[, paste0(metric, "_p975") := ci[2]]
      }
      pk <- peak_interval(rep$peak_h[ok])
      summ[, peak_center := unname(pk["center"])]
      summ[, peak_halfwidth_95 := unname(pk["halfwidth_95"])]
      boot_rows[[paste(dataset, clock)]] <- summ
    }
  }
}
points <- rbindlist(point_rows, fill = TRUE)
mod_idx <- which(points$module != "composite")
modules <- points[mod_idx]
modules[, q_bh := p.adjust(lrt_p, method = "BH"), by = .(dataset, outcome)]
points <- merge(
  points,
  modules[, .(dataset, clock, q_bh)],
  by = c("dataset", "clock"),
  all.x = TRUE,
  sort = FALSE
)
fwrite(points, file.path(out_dir, "h1_point_estimates.tsv"), sep = "\t")
fwrite(rbindlist(boot_rows), file.path(out_dir, "h1_bootstrap_ci.tsv"), sep = "\t")
fwrite(rbindlist(within_rows), file.path(out_dir, "composite_peters_within.tsv"), sep = "\t")

cat("projection\n")
flush.console()
proj_points <- list()
proj_ci <- list()
for (dataset in datasets) {
  samples <- as.data.frame(read_samples(dataset))
  mcp <- fread(mcp_path(dataset))
  base <- merge(samples, as.data.frame(mcp), by = "gsm", sort = FALSE)
  if (anyNA(base[, blood_labels])) stop(paste("missing composition", dataset))
  base <- add_time_terms(base)
  safe <- safe_name(blood_labels)
  names(safe) <- blood_labels
  for (label in blood_labels) base[[safe[[label]]]] <- base[[label]]
  adj_full <- as.formula(paste(
    "score ~ cos_t + sin_t +", paste(unname(safe), collapse = " + "), "+ (1 | subject_id)"
  ))
  draws <- fread(draw_path(dataset))
  for (clock in composites) {
    idx <- which(scores$dataset == dataset & scores$clock == clock)
    sc <- scores[idx]
    named <- sc$value
    names(named) <- sc$gsm
    dat <- base
    dat$score <- unname(named[dat$gsm])
    unadj <- fit_cosinor(dat, forms$full, forms$null)
    adj_beta <- reml_beta(dat, adj_full)
    ru <- c(unadj$beta_cos, unadj$beta_sin)
    ra <- c(unname(adj_beta["cos"]), unname(adj_beta["sin"]))
    proj <- projection(ru, ra)
    proj_points[[paste(dataset, clock)]] <- data.table(
      dataset = dataset, clock = clock,
      n_subjects = unadj$n_subjects, n_samples = unadj$n_samples,
      lrt_p_unadj = unadj$lrt_p,
      p_med = unname(proj["p_med"]), two_a_orth = unname(proj["two_a_orth"]),
      ru_norm = unname(proj["ru_norm"]),
      reml_converged_unadj = unadj$reml_converged,
      reml_converged_adj = isTRUE(adj_beta["converged"] == 1),
      two_a_unadj = unadj$two_a
    )
    rep_path <- file.path(out_dir, paste0("projection_reps_", dataset, "_", clock, ".tsv"))
    if (file.exists(rep_path) && nrow(fread(rep_path)) == n_boot) {
      rep <- fread(rep_path)
    } else {
      cat(" projection", dataset, clock, "\n")
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
      }
      rep <- rbindlist(reps)
      fwrite(rep, rep_path, sep = "\t")
    }
    p_ci <- percentile_ci(rep$p_med)
    a_ci <- percentile_ci(rep$two_a_orth)
    proj_ci[[paste(dataset, clock)]] <- data.table(
      dataset = dataset, clock = clock, n_boot = n_boot,
      n_converged_both_reml = sum(rep$reml_converged_both),
      p_med = unname(proj["p_med"]),
      p_med_p025 = p_ci[1], p_med_p975 = p_ci[2],
      two_a_orth = unname(proj["two_a_orth"]),
      two_a_orth_p025 = a_ci[1], two_a_orth_p975 = a_ci[2],
      lrt_p_unadj = unadj$lrt_p
    )
  }
}
fwrite(rbindlist(proj_points), file.path(out_dir, "projection_point.tsv"), sep = "\t")
fwrite(rbindlist(proj_ci), file.path(out_dir, "projection_ci.tsv"), sep = "\t")
cat("Part D done\n")
