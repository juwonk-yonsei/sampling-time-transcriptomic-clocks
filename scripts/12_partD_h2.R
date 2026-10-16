#!/usr/bin/env Rscript
# Part D. H2 development exploration: MCPcounter composition.

suppressPackageStartupMessages({
  library(MCPcounter)
  library(data.table)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
root <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(root, "scripts", "lib_cosinor.R"))
out_dir <- file.path(root, "results", "stage2")

n_boot <- 2000
datasets <- c("GSE39445", "GSE56931")
clocks <- c("RA_GTEx", "RA_Peters", "RA_default")
sig_of <- c(RA_GTEx = "GTExAge", RA_Peters = "Peters", RA_default = "DESeq2")
blood_labels <- c(
  "T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells",
  "B lineage", "Monocytic lineage", "Myeloid dendritic cells", "Neutrophils"
)
score_only_labels <- c("Endothelial cells", "Fibroblasts")
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

genes_path <- file.path(root, "data", "processed", "ref", "MCPcounter_genes.txt")
genes <- read.table(
  genes_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
  check.names = FALSE, quote = "\"", comment.char = ""
)
marker_symbol <- genes[["HUGO symbols"]]
marker_pop <- genes[["Cell population"]]

ns <- asNamespace("RNAAgeCalc")
genelist <- get("genelist_all", ns)$blood

read_samples <- function(dataset) {
  fread(file.path(root, "data", "processed", paste0(dataset, "_samples.tsv")))
}

read_expr <- function(dataset) {
  expr <- fread(file.path(root, "data", "processed", paste0(dataset, "_expr_gene_log2.tsv.gz")))
  mat <- as.matrix(expr[, -1])
  rownames(mat) <- expr$gene
  storage.mode(mat) <- "double"
  mat
}

control_base <- function(samples) {
  dat <- samples[condition == "control" & !is.na(clock_time_h)]
  as.data.frame(dat[, .(gsm, subject_id, clock_time_h)])
}

overlap_rows <- list()
for (clock in clocks) {
  coefs <- genelist[[sig_of[[clock]]]]
  clock_genes <- setdiff(names(coefs), "(Intercept)")
  for (pop in unique(marker_pop)) {
    marks <- marker_symbol[marker_pop == pop]
    overlap_rows[[paste(clock, pop)]] <- data.frame(
      clock = clock,
      population = pop,
      n_marker = length(marks),
      n_clock_coef = length(clock_genes),
      n_overlap = length(intersect(clock_genes, marks)),
      stringsAsFactors = FALSE
    )
  }
}
fwrite(rbindlist(overlap_rows), file.path(out_dir, "h2_marker_clock_overlap.tsv"), sep = "\t")

forms_unadj <- main_formula()
pop_rows <- list()
clock_rows <- list()
corr_rows <- list()
score_rows <- list()
boot_rows <- list()

for (dataset in datasets) {
  cat("scores", dataset, "\n")
  flush.console()
  mat <- read_expr(dataset)
  est <- MCPcounter.estimate(mat, featuresType = "HUGO_symbols", genes = genes)
  if (!identical(colnames(est), colnames(mat))) stop("MCPcounter columns do not match expression")
  score_dt <- as.data.table(t(est), keep.rownames = "gsm")
  fwrite(score_dt, file.path(out_dir, paste0("h2_mcp_scores_", dataset, ".tsv")), sep = "\t")
  present <- rownames(mat)
  for (pop in unique(marker_pop)) {
    marks <- marker_symbol[marker_pop == pop]
    score_rows[[paste(dataset, pop, "coverage")]] <- data.frame(
      dataset = dataset,
      population = pop,
      n_marker = length(marks),
      n_marker_in_matrix = sum(marks %in% present),
      stringsAsFactors = FALSE
    )
  }
  samples <- read_samples(dataset)
  base <- control_base(samples)
  base <- merge(base, as.data.frame(score_dt), by = "gsm", sort = FALSE)
  if (anyNA(base[, blood_labels])) stop("missing composition score on a control sample")
  base <- add_time_terms(base)
  for (pop in c(blood_labels, score_only_labels)) {
    vals_all <- score_dt[[pop]]
    vals_ctrl <- base[[pop]]
    score_rows[[paste(dataset, pop, "summary")]] <- data.frame(
      dataset = dataset,
      population = pop,
      in_clock_model = pop %in% blood_labels,
      n_all_samples = length(vals_all),
      mean_all = mean(vals_all),
      sd_all = stats::sd(vals_all),
      min_all = min(vals_all),
      max_all = max(vals_all),
      n_control = length(vals_ctrl),
      mean_control = mean(vals_ctrl),
      sd_control = stats::sd(vals_ctrl),
      min_control = min(vals_ctrl),
      max_control = max(vals_ctrl),
      stringsAsFactors = FALSE
    )
  }
  cm <- cor(base[, blood_labels], use = "complete.obs")
  for (i in seq_along(blood_labels)) {
    for (j in seq_along(blood_labels)) {
      if (j <= i) next
      corr_rows[[paste(dataset, i, j)]] <- data.frame(
        dataset = dataset,
        sample_set = "control_nonmissing_clock_time_h",
        n = nrow(base),
        population_a = blood_labels[[i]],
        population_b = blood_labels[[j]],
        pearson_r = cm[i, j],
        abs_r_gt_0.9 = abs(cm[i, j]) > 0.9,
        stringsAsFactors = FALSE
      )
    }
  }
  for (pop in blood_labels) {
    dat <- base
    dat$score <- dat[[pop]]
    fit <- fit_cosinor(dat, forms_unadj$full, forms_unadj$null)
    if (!isTRUE(fit$reml_ok)) stop(paste("population fit failed", dataset, pop, fit$warnings))
    fit$dataset <- dataset
    fit$population <- pop
    pop_rows[[paste(dataset, pop)]] <- fit
  }
  draws <- fread(file.path(out_dir, paste0("boot_draws_", dataset, ".tsv")))
  scores <- fread(file.path(root, "results", "stage1", "clocks", paste0(dataset, "_rnaage.tsv")))
  safe <- safe_name(blood_labels)
  names(safe) <- blood_labels
  for (label in blood_labels) base[[safe[[label]]]] <- base[[label]]
  adj_full <- as.formula(paste(
    "score ~ cos_t + sin_t +", paste(safe, collapse = " + "), "+ (1 | subject_id)"
  ))
  adj_null <- as.formula(paste(
    "score ~", paste(safe, collapse = " + "), "+ (1 | subject_id)"
  ))
  for (clock in clocks) {
    cat("clock", dataset, clock, "\n")
    flush.console()
    keep <- which(scores$clock == clock & scores$method == "lin")
    sc <- scores[keep]
    named <- sc$value
    names(named) <- sc$gsm
    dat <- base
    dat$score <- named[dat$gsm]
    if (anyNA(dat$score)) stop("missing lin score")
    unadj <- fit_cosinor(dat, forms_unadj$full, forms_unadj$null)
    adj <- fit_cosinor(dat, adj_full, adj_null)
    ratio <- 1 - adj$two_a / unadj$two_a
    clock_rows[[paste(dataset, clock)]] <- data.frame(
      dataset = dataset,
      clock = clock,
      method = "lin",
      n_subjects = unadj$n_subjects,
      n_samples = unadj$n_samples,
      two_a_unadj = unadj$two_a,
      r2_within_unadj = unadj$r2_within,
      peak_h_unadj = unadj$peak_h,
      lrt_p_unadj = unadj$lrt_p,
      reml_converged_unadj = unadj$reml_converged,
      two_a_adj = adj$two_a,
      r2_within_adj = adj$r2_within,
      peak_h_adj = adj$peak_h,
      lrt_p_adj = adj$lrt_p,
      reml_converged_adj = adj$reml_converged,
      ml_converged_adj = adj$ml_converged,
      null_converged_adj = adj$null_converged,
      reduction_ratio = ratio,
      adj_full_formula = paste(deparse(adj_full), collapse = ""),
      adj_null_formula = paste(deparse(adj_null), collapse = ""),
      stringsAsFactors = FALSE
    )
    reps <- vector("list", n_boot)
    for (b in seq_len(n_boot)) {
      idx <- which(draws$draw == b)
      one <- draws[idx]
      boot_dat <- apply_draw(dat, one$subject_id, one$new_subject_id)
      u <- fit_cosinor(boot_dat, forms_unadj$full, forms_unadj$null)
      a <- fit_cosinor(boot_dat, adj_full, adj_null)
      ok <- isTRUE(u$reml_converged) && isTRUE(a$reml_converged)
      reps[[b]] <- data.frame(
        draw = b,
        two_a_unadj = u$two_a,
        two_a_adj = a$two_a,
        reduction_ratio = if (ok && is.finite(u$two_a) && u$two_a != 0) 1 - a$two_a / u$two_a else NA_real_,
        reml_converged_unadj = u$reml_converged,
        reml_converged_adj = a$reml_converged,
        stringsAsFactors = FALSE
      )
      if (b %% 200 == 0) {
        cat("  ", b, "\n")
        flush.console()
      }
    }
    rep <- rbindlist(reps)
    fwrite(rep, file.path(out_dir, paste0("h2_reduction_reps_", dataset, "_", clock, ".tsv")), sep = "\t")
    ci <- percentile_ci(rep$reduction_ratio)
    boot_rows[[paste(dataset, clock)]] <- data.frame(
      dataset = dataset,
      clock = clock,
      method = "lin",
      n_boot = n_boot,
      n_converged_both_reml = sum(is.finite(rep$reduction_ratio)),
      reduction_ratio = ratio,
      reduction_p025 = ci[1],
      reduction_p975 = ci[2],
      stringsAsFactors = FALSE
    )
  }
  rm(mat)
  gc()
}

fwrite(rbindlist(pop_rows, fill = TRUE), file.path(out_dir, "h2_population_cosinor.tsv"), sep = "\t")
fwrite(rbindlist(clock_rows, fill = TRUE), file.path(out_dir, "h2_clock_adjusted.tsv"), sep = "\t")
fwrite(rbindlist(corr_rows), file.path(out_dir, "h2_population_correlations.tsv"), sep = "\t")
fwrite(rbindlist(score_rows, fill = TRUE), file.path(out_dir, "h2_score_summary.tsv"), sep = "\t")
fwrite(rbindlist(boot_rows), file.path(out_dir, "h2_reduction_ci.tsv"), sep = "\t")
cat("Part D done\n")
