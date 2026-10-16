#!/usr/bin/env Rscript
# Confirmation steps. Locked expression is opened only by load_locked_expression().

suppressPackageStartupMessages({
  library(lme4)
  library(data.table)
  library(MCPcounter)
})

find_tod_root <- function() {
  starts <- getwd()
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if (length(file_arg)) starts <- c(dirname(sub("^--file=", "", file_arg[1])), starts)
  for (start in starts) {
    cur <- normalizePath(start)
    for (step in seq_len(8)) {
      if (file.exists(file.path(cur, "scripts", "lib_cosinor.R"))) return(cur)
      parent <- dirname(cur)
      if (parent == cur) break
      cur <- parent
    }
  }
  stop("TOD root not found")
}

args_all <- commandArgs(trailingOnly = FALSE)
root <- find_tod_root()
source(file.path(root, "scripts", "lib_cosinor.R"))
source(file.path(root, "scripts", "confirm", "lib_lock.R"))
source(file.path(root, "scripts", "confirm", "lib_tests.R"))

python_bin <- "/HDD8T1/WORK/GeroSci/tod4_clock_env/bin/python"
primary_clocks <- c("TAms_chrono", "TAms_mort", "RA_default")
secondary_clocks <- c("RA_Peters", "RA_GTEx")
ra_map <- c(RA_default = "DESeq2", RA_Peters = "Peters", RA_GTEx = "GTExAge")

ra_linear_scores <- function(log_mat) {
  log_mat <- as.matrix(log_mat)
  storage.mode(log_mat) <- "double"
  genelist <- get("genelist_all", asNamespace("RNAAgeCalc"))$blood
  rows <- list()
  for (clock in names(ra_map)) {
    coefs <- genelist[[ra_map[[clock]]]]
    present <- intersect(setdiff(names(coefs), "(Intercept)"), rownames(log_mat))
    used <- coefs[present]
    score <- unname(coefs[["(Intercept)"]]) + as.numeric(crossprod(used, log_mat[present, , drop = FALSE]))
    rows[[clock]] <- data.table(gsm = colnames(log_mat), clock = clock, value = score, subject_id = NA_character_)
  }
  rbindlist(rows)
}

counts_to_log2_fpkm <- function(counts) {
  fpkm_fun <- get("count2FPKM", asNamespace("RNAAgeCalc"))
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  fpkm <- as.matrix(fpkm_fun(counts, idtype = "SYMBOL"))
  storage.mode(fpkm) <- "double"
  log2(fpkm + 1)
}

log2_cpm <- function(counts) {
  lib <- colSums(counts)
  cpm <- sweep(counts, 2, lib, "/") * 1e6
  log2(cpm + 1)
}

mcp_scores <- function(log_mat) {
  genes <- read.table(
    file.path(root, "data", "processed", "ref", "MCPcounter_genes.txt"),
    sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "\"", comment.char = ""
  )
  est <- MCPcounter.estimate(log_mat, featuresType = "HUGO_symbols", genes = genes)
  as.data.frame(t(est))
}

score_tage <- function(expr_path, samples_path, mode, out_scores, out_components) {
  status <- system2(python_bin, c(
    file.path(root, "scripts", "confirm", "score_clocks.py"),
    "--expr", expr_path, "--samples", samples_path, "--mode", mode,
    "--out", out_scores, "--components", out_components
  ))
  if (!identical(status, 0L)) stop("tAge scoring failed")
  list(scores = fread(out_scores), components = fread(out_components))
}

family_tables <- function(samples_m, score_dt, clocks, draws, out_dir, family) {
  point <- list()
  boot <- list()
  frames <- list()
  for (clock in clocks) {
    dat <- prepare_score(samples_m, score_dt, clock)
    fit <- model_m(dat)
    fit$clock <- clock
    fit$family <- family
    point[[clock]] <- as.data.table(fit)
    fwrite(as.data.table(dat)[, .(subject_id, clock_time_h, cos_t, sin_t, score)], file.path(out_dir, paste0("frame_", family, "_", clock, ".tsv")), sep = "\t")
    if (!is.null(draws)) {
      rep <- bootstrap_model_m(dat, draws)
      ok <- rep$reml_converged & rep$ml_converged & rep$null_converged
      ci <- percentile_ci(rep$two_a[ok])
      pk <- peak_interval(rep$peak_h[ok])
      boot[[clock]] <- data.table(
        clock = clock, n_boot = nrow(rep), n_converged = sum(ok),
        two_a_p025 = ci[1], two_a_p975 = ci[2],
        peak_center = unname(pk["center"]), peak_halfwidth_95 = unname(pk["halfwidth_95"])
      )
    }
  }
  points <- rbindlist(point, fill = TRUE)
  points[, q_holm := holm_family(lrt_p)]
  fwrite(points, file.path(out_dir, paste0(family, "_point.tsv")), sep = "\t")
  if (length(boot)) fwrite(rbindlist(boot), file.path(out_dir, paste0(family, "_bootstrap.tsv")), sep = "\t")
  points
}

run_gate_or_stop <- function(gate, out_dir) {
  fwrite(gate, file.path(out_dir, "gate.tsv"), sep = "\t")
  passed <- isTRUE(gate$dataset_pass[1])
  if (!passed) {
    fwrite(data.table(status = "gate fail"), file.path(out_dir, "stopped.tsv"), sep = "\t")
  }
  passed
}

analyse_cosinor_dataset <- function(spec) {
  dir.create(spec$out_dir, recursive = TRUE, showWarnings = FALSE)
  samples <- spec$samples
  fwrite(samples, file.path(spec$out_dir, "samples_used.tsv"), sep = "\t")
  if (!run_gate_or_stop(spec$gate, spec$out_dir)) return(invisible(FALSE))
  score_dt <- spec$scores
  score_dt[, gsm := as.character(gsm)]
  primary_samples <- samples[samples$analysis_set == spec$primary_set, ]
  draws <- participant_draws(unique(primary_samples$subject_id), spec$n_boot, 20261015)
  fwrite(draws, file.path(spec$out_dir, "boot_draws.tsv"), sep = "\t")
  p_tab <- family_tables(primary_samples, score_dt, primary_clocks, draws, spec$out_dir, "P")
  family_tables(primary_samples, score_dt, secondary_clocks, draws, spec$out_dir, "S1")
  comp_cols <- intersect(blood_labels, colnames(spec$composition))
  if (isTRUE(spec$drop_neutrophils)) comp_cols <- setdiff(comp_cols, "Neutrophils")
  h2_rows <- list()
  for (clock_id in primary_clocks) {
    idx <- which(p_tab$clock == clock_id)
    if (!isTRUE(p_tab$q_holm[idx] < 0.05)) next
    dat <- prepare_score(primary_samples, score_dt, clock_id)
    comp <- spec$composition[dat$gsm, comp_cols, drop = FALSE]
    for (col in comp_cols) dat[[col]] <- comp[[col]]
    h2 <- h2_l1l2(dat, comp_cols, draws)
    h2$point$clock <- clock_id
    h2_rows[[clock_id]] <- h2$point
  }
  if (length(h2_rows)) fwrite(rbindlist(h2_rows, fill = TRUE), file.path(spec$out_dir, "h2.tsv"), sep = "\t")
  if (!is.null(spec$h4_samples)) {
    h4_rows <- list()
    for (clock_id in primary_clocks) {
      dat <- prepare_score(spec$h4_samples, score_dt, clock_id)
      dat$day <- spec$h4_samples$day
      fit <- h4_interaction(dat)
      fit$clock <- clock_id
      h4_rows[[clock_id]] <- fit
    }
    h4 <- rbindlist(h4_rows)
    h4[, q_holm := holm_family(lrt_p)]
    fwrite(h4, file.path(spec$out_dir, "h4.tsv"), sep = "\t")
  }
  if (!is.null(spec$module_scores)) {
    mod <- spec$module_scores
    f3 <- spec$f3_clocks
    mod_rows <- list()
    for (clock_id in unique(mod$clock)) {
      dat <- prepare_score(primary_samples, mod, clock_id)
      fit <- model_m(dat)
      fit$clock <- clock_id
      fit$in_f3 <- clock_id %in% f3
      mod_rows[[clock_id]] <- as.data.table(fit)
    }
    mod_tab <- rbindlist(mod_rows, fill = TRUE)
    idx <- which(!mod_tab$in_f3)
    mod_tab[, q_bh := NA_real_]
    if (length(idx)) set(mod_tab, i = idx, j = "q_bh", value = p.adjust(mod_tab$lrt_p[idx], method = "BH"))
    fwrite(mod_tab, file.path(spec$out_dir, "s4_point.tsv"), sep = "\t")
  }
  sens <- list()
  flags <- spec$flags
  keep <- primary_samples$gsm %in% flags$gsm[!flags$flag_3mad]
  if (sum(keep) > 5) {
    for (clock_id in primary_clocks) {
      dat <- prepare_score(primary_samples[keep, ], score_dt, clock_id)
      fit <- model_m(dat)
      fit$clock <- clock_id
      sens[[clock_id]] <- as.data.table(fit)
    }
    fwrite(rbindlist(sens, fill = TRUE), file.path(spec$out_dir, "sensitivity_drop_flags.tsv"), sep = "\t")
  }
  slope_rows <- list()
  for (clock_id in primary_clocks) {
    dat <- prepare_score(primary_samples, score_dt, clock_id)
    fit <- quiet_lmer(score ~ cos_t + sin_t + (1 + cos_t + sin_t | subject_id), dat, TRUE)
    slope_rows[[clock_id]] <- data.table(
      clock = clock_id,
      converged = !inherits(fit, "error") && optimizer_converged(fit)
    )
  }
  fwrite(rbindlist(slope_rows), file.path(spec$out_dir, "sensitivity_random_slopes.tsv"), sep = "\t")
  fwrite(data.table(status = "complete"), file.path(spec$out_dir, "stopped.tsv"), sep = "\t")
  invisible(TRUE)
}

analyse_paired_dataset <- function(spec) {
  dir.create(spec$out_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(spec$samples, file.path(spec$out_dir, "samples_used.tsv"), sep = "\t")
  if (!run_gate_or_stop(spec$gate, spec$out_dir)) return(invisible(FALSE))
  rows <- list()
  delta_rows <- list()
  for (clock_id in c(primary_clocks, secondary_clocks)) {
    res <- paired_delta(spec$samples, spec$scores, clock_id)
    delta <- res$delta[[1]]
    res_subjects <- res$subjects[[1]]
    res_morning <- res$morning[[1]]
    res_evening <- res$evening[[1]]
    res$delta <- NULL
    res$subjects <- NULL
    res$morning <- NULL
    res$evening <- NULL
    res$family <- if (clock_id %in% primary_clocks) "P" else "S1"
    rows[[clock_id]] <- res
    delta_rows[[clock_id]] <- data.table(
      clock = clock_id, subject_id = res_subjects, morning = res_morning,
      evening = res_evening, delta = delta
    )
  }
  tab <- rbindlist(rows)
  tab[, q_holm := NA_real_]
  for (fam in unique(tab$family)) {
    idx <- which(tab$family == fam)
    set(tab, i = idx, j = "q_holm", value = holm_family(tab$p[idx]))
  }
  fwrite(tab, file.path(spec$out_dir, "paired_tests.tsv"), sep = "\t")
  fwrite(rbindlist(delta_rows), file.path(spec$out_dir, "paired_deltas.tsv"), sep = "\t")
  if (!is.null(spec$covariates)) {
    h3_rows <- list()
    h5_rows <- list()
    intercept_rows <- list()
    cov <- spec$covariates
    for (clock_id in primary_clocks) {
      ord <- match(delta_rows[[clock_id]]$subject_id, cov$subject_id)
      if (anyNA(ord)) stop("covariates do not cover the paired subjects")
      cov_i <- cov[ord]
      delta <- delta_rows[[clock_id]]$delta
      h3_rows[[clock_id]] <- cbind(data.table(clock = clock_id), h3_age(delta, cov_i$age, cov_i$sex, cov_i$opioid_group))
      c_value <- spec$c_values[[clock_id]]
      h5_rows[[clock_id]] <- cbind(data.table(clock = clock_id, c = c_value), h5_transfer(delta, c_value))
      intercept_rows[[clock_id]] <- cbind(data.table(clock = clock_id), sensitivity_intercept(delta, cov_i$sex, cov_i$opioid_group))
    }
    h3 <- rbindlist(h3_rows)
    h3[, q_holm := holm_family(age_p)]
    fwrite(h3, file.path(spec$out_dir, "h3.tsv"), sep = "\t")
    fwrite(rbindlist(h5_rows), file.path(spec$out_dir, "h5.tsv"), sep = "\t")
    fwrite(rbindlist(intercept_rows), file.path(spec$out_dir, "sensitivity_intercept.tsv"), sep = "\t")
    p_primary <- tab[family == "P"]
    if (!is.null(spec$delta_comp)) {
      h2_rows <- list()
      for (clock_id in primary_clocks) {
        idx <- which(p_primary$clock == clock_id)
        if (!isTRUE(p_primary$q_holm[idx] < 0.05)) next
        delta <- delta_rows[[clock_id]]$delta
        fit <- h2_l3_bootstrap(delta, spec$delta_comp, spec$n_boot, 20261015)
        fit$clock <- clock_id
        h2_rows[[clock_id]] <- fit
      }
      if (length(h2_rows)) fwrite(rbindlist(h2_rows, fill = TRUE), file.path(spec$out_dir, "h2.tsv"), sep = "\t")
    }
  }
  if (!is.null(spec$module_scores)) {
    mod_rows <- list()
    for (clock_id in unique(spec$module_scores$clock)) {
      res <- paired_delta(spec$samples, spec$module_scores, clock_id)
      res$delta <- NULL
      res$subjects <- NULL
      res$morning <- NULL
      res$evening <- NULL
      res$in_f3 <- clock_id %in% spec$f3_clocks
      mod_rows[[clock_id]] <- res
    }
    mod_tab <- rbindlist(mod_rows, fill = TRUE)
    idx <- which(!mod_tab$in_f3)
    mod_tab[, q_bh := NA_real_]
    if (length(idx)) set(mod_tab, i = idx, j = "q_bh", value = p.adjust(mod_tab$p[idx], method = "BH"))
    fwrite(mod_tab, file.path(spec$out_dir, "s4_paired.tsv"), sep = "\t")
  }
  batch_n <- if (is.null(spec$batch2_n)) 0L else spec$batch2_n
  fwrite(data.table(batch2_pairs = batch_n, note = "batch-2 pairs are absent in this run"), file.path(spec$out_dir, "sensitivity_batch2.tsv"), sep = "\t")
  combat_note <- if (is.null(spec$combat_path)) "combat-vst file absent; not downloaded" else spec$combat_path
  fwrite(data.table(combat = combat_note), file.path(spec$out_dir, "sensitivity_combat.tsv"), sep = "\t")
  fwrite(data.table(status = "complete"), file.path(spec$out_dir, "stopped.tsv"), sep = "\t")
  invisible(TRUE)
}

run_locked <- function() {
  check_manifest(root)
  # Step 2. The load below returns only after the unlock file and TOD_UNLOCK match the manifest hash.
  load_locked_expression("L1", root)
  stop("locked confirmation steps 3-12 are instruction 7")
}

check_manifest <- function(root) {
  manifest <- file.path(root, "prereg", "manifest.tsv")
  if (!file.exists(manifest)) stop("manifest is missing")
  tab <- fread(manifest)
  bad <- character()
  for (i in seq_len(nrow(tab))) {
    path <- file.path(root, tab$path[i])
    if (!file.exists(path)) {
      bad <- c(bad, paste("missing", tab$path[i]))
      next
    }
    size <- file.info(path)$size
    if (!identical(as.numeric(size), as.numeric(tab$bytes[i]))) bad <- c(bad, paste("bytes", tab$path[i]))
    hash <- strsplit(system2("sha256sum", path, stdout = TRUE), " ", fixed = TRUE)[[1]][1]
    if (!identical(hash, tab$sha256[i])) bad <- c(bad, paste("hash", tab$path[i]))
  }
  if (length(bad)) stop(paste(bad, collapse = "; "))
  invisible(TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
file_arg_now <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg_now) && grepl("run_pipeline\\.R$", file_arg_now[1]) && "--mode" %in% args) {
  mode <- args[match("--mode", args) + 1]
  if (mode == "locked") run_locked()
}
