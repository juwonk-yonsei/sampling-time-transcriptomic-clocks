#!/usr/bin/env Rscript
# Part D. Rhythm-vector projection on microarray development sets.
# RNA-seq lin scores were not produced because the positive control failed.

suppressPackageStartupMessages({
  library(data.table)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
root <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(root, "scripts", "lib_cosinor.R"))
stage2 <- file.path(root, "results", "stage2")
out_dir <- file.path(root, "results", "stage3", "h2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- c("GSE39445", "GSE56931")
clocks <- c("RA_GTEx", "RA_Peters", "RA_default")
sig_of <- c(RA_GTEx = "GTExAge", RA_Peters = "Peters", RA_default = "DESeq2")
blood_labels <- c(
  "T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells",
  "B lineage", "Monocytic lineage", "Myeloid dendritic cells", "Neutrophils"
)
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
forms <- main_formula()

ns <- asNamespace("RNAAgeCalc")
genelist <- get("genelist_all", ns)$blood

read_expr <- function(dataset) {
  expr <- fread(file.path(root, "data", "processed", paste0(dataset, "_expr_gene_log2.tsv.gz")))
  mat <- as.matrix(expr[, -1])
  rownames(mat) <- expr$gene
  storage.mode(mat) <- "double"
  mat
}

lin_scores <- function(mat, clock) {
  coefs <- genelist[[sig_of[[clock]]]]
  intercept <- unname(coefs[["(Intercept)"]])
  genes <- setdiff(names(coefs), "(Intercept)")
  have <- intersect(genes, rownames(mat))
  used <- coefs[have]
  scores <- as.numeric(intercept + drop(crossprod(used, mat[have, , drop = FALSE])))
  names(scores) <- colnames(mat)
  scores
}

reml_beta <- function(data, formula) {
  fit <- quiet_lmer(formula, data, TRUE)
  if (inherits(fit, "error") || !optimizer_converged(fit)) {
    return(c(cos = NA_real_, sin = NA_real_, converged = FALSE))
  }
  beta <- lme4::fixef(fit)
  c(cos = unname(beta["cos_t"]), sin = unname(beta["sin_t"]), converged = TRUE)
}

projection <- function(ru, ra) {
  ru2 <- sum(ru * ru)
  if (!is.finite(ru2) || ru2 == 0) {
    return(c(p_med = NA_real_, two_a_orth = NA_real_, ru_norm = NA_real_))
  }
  scale <- sum(ra * ru) / ru2
  resid <- ra - scale * ru
  c(
    p_med = 1 - scale,
    two_a_orth = 2 * sqrt(sum(resid * resid)),
    ru_norm = sqrt(ru2)
  )
}

stage2_points <- fread(file.path(stage2, "h1_point_estimates.tsv"))
stage2_points <- stage2_points[method == "lin"]

point_rows <- list()
ci_rows <- list()
within_rows <- list()
match_rows <- list()

fwrite(data.table(
  dataset = c("GSE113883", "GSE196728"),
  projection = "not_run",
  within_correlation = "not_run",
  reason = "positive control failed; lin clock scores and MCPcounter scores were not produced in section 3.2-3.4"
), file.path(out_dir, "rnaseq_projection_status.tsv"), sep = "\t")

for (dataset in datasets) {
  cat("dataset", dataset, "\n")
  flush.console()
  samples <- fread(file.path(root, "data", "processed", paste0(dataset, "_samples.tsv")))
  mcp <- fread(file.path(stage2, paste0("h2_mcp_scores_", dataset, ".tsv")))
  mat <- read_expr(dataset)
  base <- samples[condition == "control" & !is.na(clock_time_h), .(gsm, subject_id, clock_time_h)]
  base <- merge(base, mcp, by = "gsm", sort = FALSE)
  if (anyNA(base[, ..blood_labels])) stop("missing composition score")
  base <- as.data.frame(base)
  base <- add_time_terms(base)
  safe <- safe_name(blood_labels)
  names(safe) <- blood_labels
  for (label in blood_labels) base[[safe[[label]]]] <- base[[label]]
  adj_full <- as.formula(paste(
    "score ~ cos_t + sin_t +", paste(safe, collapse = " + "), "+ (1 | subject_id)"
  ))
  draws <- fread(file.path(stage2, paste0("boot_draws_", dataset, ".tsv")))
  subjects <- sort(unique(base$subject_id))
  draw_subjects <- sort(unique(draws$subject_id))
  if (!identical(subjects, draw_subjects)) stop("bootstrap subjects do not match control subjects")

  for (clock in clocks) {
    cat(" clock", clock, "\n")
    flush.console()
    named <- lin_scores(mat, clock)
    dat <- base
    dat$score <- unname(named[dat$gsm])
    if (anyNA(dat$score)) stop("missing lin score")
    unadj <- fit_cosinor(dat, forms$full, forms$null)
    adj_beta <- reml_beta(dat, adj_full)
    ru <- c(unadj$beta_cos, unadj$beta_sin)
    ra <- c(unname(adj_beta["cos"]), unname(adj_beta["sin"]))
    proj <- projection(ru, ra)
    keep <- which(stage2_points$dataset == dataset & stage2_points$clock == clock)
    stored <- stage2_points[keep]
    if (nrow(stored) != 1L) stop("stage2 point row count")
    match_rows[[paste(dataset, clock)]] <- data.table(
      dataset = dataset,
      clock = clock,
      two_a = unadj$two_a,
      two_a_stage2 = stored$two_a,
      abs_diff_two_a = abs(unadj$two_a - stored$two_a),
      lrt_p = unadj$lrt_p,
      lrt_p_stage2 = stored$lrt_p
    )
    if (abs(unadj$two_a - stored$two_a) > 1e-6) {
      stop(paste("unadjusted 2A does not match stage 2", dataset, clock))
    }
    point_rows[[paste(dataset, clock)]] <- data.table(
      dataset = dataset,
      clock = clock,
      method = "lin",
      n_subjects = unadj$n_subjects,
      n_samples = unadj$n_samples,
      lrt_p_unadj = unadj$lrt_p,
      ru_cos = ru[1],
      ru_sin = ru[2],
      ru_norm = unname(proj["ru_norm"]),
      ra_cos = ra[1],
      ra_sin = ra[2],
      p_med = unname(proj["p_med"]),
      two_a_orth = unname(proj["two_a_orth"]),
      reml_converged_unadj = unadj$reml_converged,
      reml_converged_adj = isTRUE(adj_beta["converged"] == 1),
      two_a_unadj = unadj$two_a
    )

    dev_clock <- ave(dat$score, dat$subject_id, FUN = function(x) x - mean(x))
    for (pop in blood_labels) {
      dev_pop <- ave(dat[[pop]], dat$subject_id, FUN = function(x) x - mean(x))
      within_rows[[paste(dataset, clock, pop)]] <- data.table(
        dataset = dataset,
        clock = clock,
        population = pop,
        n_samples = nrow(dat),
        n_subjects = length(unique(dat$subject_id)),
        pearson_within = cor(dev_clock, dev_pop, use = "complete.obs")
      )
    }

    reps <- vector("list", 2000L)
    for (b in seq_len(2000L)) {
      one <- draws[which(draws$draw == b)]
      boot_dat <- apply_draw(dat, one$subject_id, one$new_subject_id)
      u <- reml_beta(boot_dat, forms$full)
      a <- reml_beta(boot_dat, adj_full)
      ok <- isTRUE(u["converged"] == 1) && isTRUE(a["converged"] == 1)
      pr <- if (ok) projection(c(u["cos"], u["sin"]), c(a["cos"], a["sin"])) else c(p_med = NA_real_, two_a_orth = NA_real_, ru_norm = NA_real_)
      reps[[b]] <- data.table(
        draw = b,
        reml_converged_both = ok,
        p_med = unname(pr["p_med"]),
        two_a_orth = unname(pr["two_a_orth"])
      )
      if (b %% 200L == 0L) {
        cat("  ", b, "\n")
        flush.console()
      }
    }
    rep <- rbindlist(reps)
    fwrite(rep, file.path(out_dir, paste0("projection_reps_", dataset, "_", clock, ".tsv")), sep = "\t")
    p_ci <- percentile_ci(rep$p_med)
    a_ci <- percentile_ci(rep$two_a_orth)
    ci_rows[[paste(dataset, clock)]] <- data.table(
      dataset = dataset,
      clock = clock,
      method = "lin",
      n_boot = 2000L,
      n_converged_both_reml = sum(rep$reml_converged_both),
      n_finite_p_med = sum(is.finite(rep$p_med)),
      p_med = unname(proj["p_med"]),
      p_med_p025 = p_ci[1],
      p_med_p975 = p_ci[2],
      two_a_orth = unname(proj["two_a_orth"]),
      two_a_orth_p025 = a_ci[1],
      two_a_orth_p975 = a_ci[2],
      lrt_p_unadj = unadj$lrt_p,
      ru_norm = unname(proj["ru_norm"])
    )
  }
  rm(mat)
  gc()
}

fwrite(rbindlist(point_rows), file.path(out_dir, "projection_point.tsv"), sep = "\t")
fwrite(rbindlist(ci_rows), file.path(out_dir, "projection_ci.tsv"), sep = "\t")
fwrite(rbindlist(within_rows), file.path(out_dir, "clock_composition_within.tsv"), sep = "\t")
fwrite(rbindlist(match_rows), file.path(out_dir, "unadj_match_stage2.tsv"), sep = "\t")
cat("done\n")
