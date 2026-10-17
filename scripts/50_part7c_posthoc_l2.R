#!/usr/bin/env Rscript
# Post-hoc L2 checks. Not part of family P or replication.

suppressPackageStartupMessages(library(data.table))

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
root <- normalizePath(file.path(script_dir, ".."))
source(file.path(root, "scripts", "confirm", "run_pipeline.R"))
source(file.path(root, "scripts", "confirm", "run_locked_analysis.R"))

out_dir <- file.path(root, "results", "posthoc_stage7c")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suffix <- "(사후, 관문 실패 자료)"

quotes <- data.table(
  source = c("GEO sample description", "paper figure 1 legend"),
  citation = c(
    "GSE122725 series matrix header, !Sample_description, GSM3484078",
    "Kervezee et al., J Biol Rhythms, PMID 30621487, DOI 10.1177/0748730418821776"
  ),
  sentence = c(
    "PBMC sample taken at 0.8 h after habitual wake time during day-oriented schedule in Subject T10 in the BrightLight group",
    "relative clock time refers to the time of day relative to the start of each participant's habitual sleep period."
  ),
  title = paste("relclocktime source sentences", suffix)
)
fwrite(quotes, file.path(out_dir, "relclocktime_quotes.tsv"), sep = "\t")

genes <- c("PER1", "PER2", "PER3", "NR1D1", "NR1D2", "ARNTL", "CRY1", "CRY2", "DBP", "TEF")

peak_one <- function(samples, expr, gene) {
  empty <- data.table(gene = gene, n_samples = 0L, n_subjects = 0L, two_a = NA_real_, peak_h = NA_real_, lrt_stat = NA_real_, lrt_p = NA_real_, status = "absent")
  if (!(gene %in% rownames(expr))) return(empty)
  values <- as.numeric(expr[gene, samples$gsm])
  keep <- is.finite(values) & is.finite(samples$clock_time_h)
  used <- samples[keep]
  used_values <- values[keep]
  if (nrow(used) < 2L) {
    empty$n_samples <- nrow(used)
    empty$status <- "fewer than 2 finite samples"
    return(empty)
  }
  score <- data.table(gsm = used$gsm, clock = gene, value = used_values)
  dat <- prepare_score(used, score, gene)
  fit <- model_m(dat)
  data.table(
    gene = gene, n_samples = nrow(used), n_subjects = uniqueN(used$subject_id),
    two_a = fit$two_a, peak_h = fit$peak_h, lrt_stat = fit$lrt_stat, lrt_p = fit$lrt_p, status = "fit"
  )
}

l2_samples <- fread(file.path(root, "results", "confirm", "L2", "samples_scoring.tsv"))
l2_samples[, gsm := as.character(gsm)]
l2_expr <- as.matrix(read_gene_matrix(file.path(root, "results", "confirm", "L2", "log2_genes.tsv")))
period1 <- l2_samples[period == "period1"]
period1[, clock_time_h := wake_ref_h %% 24]
l2_peaks <- rbindlist(lapply(genes, function(gene) peak_one(period1, l2_expr, gene)))
l2_peaks[, dataset := "L2_period1"]
l2_peaks[, time_axis := "wake_ref_h %% 24"]

dev <- fread(file.path(root, "data", "processed", "GSE39445_samples.tsv"))
dev <- dev[condition == "control"]
dev[, gsm := as.character(gsm)]
dev[, clock_time_h := hours_awake]
dev_expr <- as.matrix(read_gene_matrix(file.path(root, "data", "processed", "GSE39445_expr_gene_log2.tsv.gz")))
dev_peaks <- rbindlist(lapply(genes, function(gene) peak_one(dev, dev_expr, gene)))
dev_peaks[, dataset := "GSE39445_control"]
dev_peaks[, time_axis := "hours_awake"]

signed_diff <- function(a, b) {
  d <- (a - b) %% 24
  ifelse(d > 12, d - 24, d)
}
cmp <- merge(
  l2_peaks[, .(gene, peak_l2 = peak_h, status_l2 = status)],
  dev_peaks[, .(gene, peak_gse39445 = peak_h, status_gse39445 = status)],
  by = "gene"
)
cmp[, signed_diff_h := signed_diff(peak_l2, peak_gse39445)]
cmp[, circ_dist_0_12 := circ_dist_0_12(peak_l2, peak_gse39445)]
cmp[, title := paste("clock-gene peak comparison", suffix)]
setcolorder(cmp, c("gene", "peak_l2", "peak_gse39445", "signed_diff_h", "circ_dist_0_12", "status_l2", "status_gse39445", "title"))
cmp <- cmp[match(genes, gene)]

finite_diff <- cmp[is.finite(signed_diff_h), signed_diff_h]
circ_sd_h <- function(hours) {
  ang <- hours * 2 * pi / 24
  rbar <- sqrt(mean(cos(ang))^2 + mean(sin(ang))^2)
  if (!is.finite(rbar) || rbar <= 0) return(NA_real_)
  sqrt(-2 * log(min(rbar, 1))) * 24 / (2 * pi)
}
summary_row <- data.table(
  n_genes = length(finite_diff),
  circ_mean_h = if (length(finite_diff)) circ_mean_h(finite_diff) else NA_real_,
  circ_sd_h = if (length(finite_diff)) circ_sd_h(finite_diff) else NA_real_,
  signed_diff_definition = "(peak_l2 - peak_gse39445) mod 24, then subtract 24 when the result is greater than 12; range (-12, 12]",
  circ_sd_definition = "sqrt(-2*log(R)) * 24 / (2*pi), R = mean resultant length of the signed differences",
  title = paste("circular mean and SD of signed peak differences", suffix)
)
l2_peaks[, title := paste("L2 period 1 clock-gene peaks", suffix)]
dev_peaks[, title := paste("GSE39445 control clock-gene peaks", suffix)]
fwrite(l2_peaks, file.path(out_dir, "l2_period1_gene_peaks.tsv"), sep = "\t")
fwrite(dev_peaks, file.path(out_dir, "gse39445_gene_peaks.tsv"), sep = "\t")
fwrite(cmp, file.path(out_dir, "gene_peak_differences.tsv"), sep = "\t")
fwrite(summary_row, file.path(out_dir, "gene_peak_circular_summary.tsv"), sep = "\t")
rm(l2_expr, dev_expr)
gc()

say <- function(...) message(format(Sys.time(), "%H:%M:%S"), " ", paste(...))
say("scoring L2 clocks")
scored <- score_tage(
  file.path(root, "results", "confirm", "L2", "log2_genes.tsv"),
  file.path(root, "results", "confirm", "L2", "samples_scoring.tsv"),
  "array",
  file.path(out_dir, "l2_tage_scores.tsv"),
  file.path(out_dir, "l2_tage_components.tsv")
)
scores <- rbind(scored$scores, scored$components, fill = TRUE)
scores[, gsm := as.character(gsm)]
log_mat <- as.matrix(read_gene_matrix(file.path(root, "results", "confirm", "L2", "log2_genes.tsv")))
ra <- ra_linear_scores(log_mat[, l2_samples$gsm, drop = FALSE])
scores <- rbind(scores, ra, fill = TRUE)
scores <- scores[!duplicated(scores, by = c("gsm", "clock"))]
rm(log_mat)
gc()

f3 <- fread(file.path(root, "prereg", "frozen", "F3.tsv"))
h1_clocks <- c(primary_clocks, secondary_clocks, f3$clock)
h1_samples <- l2_samples[period == "period1"]
h1_samples[, clock_time_h := wake_ref_h %% 24]
h1_rows <- lapply(h1_clocks, function(clock_id) {
  if (!(clock_id %in% scores$clock)) {
    return(data.table(dataset = "L2", period = "period1", clock = clock_id, n_samples = 0L, p = NA_real_, status = "score column absent"))
  }
  one <- scores[clock == clock_id]
  values <- one$value[match(h1_samples$gsm, one$gsm)]
  keep <- is.finite(values)
  used <- h1_samples[keep]
  if (nrow(used) < 2L) {
    return(data.table(dataset = "L2", period = "period1", clock = clock_id, n_samples = nrow(used), p = NA_real_, status = "fewer than 2 finite samples"))
  }
  dat <- prepare_score(used, one, clock_id)
  fit <- model_m(dat)
  data.table(
    dataset = "L2", period = "period1", clock = clock_id,
    n_subjects = uniqueN(used$subject_id), n_samples = nrow(used),
    two_a = fit$two_a, peak_h = fit$peak_h, lrt_stat = fit$lrt_stat, p = fit$lrt_p,
    status = "fit"
  )
})
h1 <- rbindlist(h1_rows, fill = TRUE)
h1[, q_bh := p.adjust(p, method = "BH")]
h1[, label := "exploratory"]
h1[, note := "An LRT is invariant to adding a constant to the whole time axis."]
h1[, title := paste("H1 exploratory model M on L2 period 1", suffix)]
fwrite(h1, file.path(out_dir, "h1_exploratory.tsv"), sep = "\t")

h4_rows <- lapply(primary_clocks, function(clock_id) {
  one <- scores[clock == clock_id]
  values <- one$value[match(l2_samples$gsm, one$gsm)]
  keep <- is.finite(values) & is.finite(l2_samples$clock_time_h)
  used <- l2_samples[keep]
  used[, day := fifelse(period == "period1", 0L, 1L)]
  dat <- as.data.table(prepare_score(used, one, clock_id))
  dat[, day := used$day]
  fit <- h4_interaction(dat)
  fit[, dataset := "L2"]
  fit[, clock := clock_id]
  fit[, n_samples := nrow(used)]
  fit[, n_period1 := sum(used$day == 0L)]
  fit[, n_period2 := sum(used$day == 1L)]
  fit[, day_coding := "period1=0, period2=1"]
  fit
})
h4 <- rbindlist(h4_rows, fill = TRUE)
h4[, q_bh := p.adjust(lrt_p, method = "BH")]
h4[, label := "exploratory"]
h4[, note := "An LRT is invariant to adding a constant to the whole time axis."]
h4[, title := paste("H4 exploratory interaction LRT on L2 period 1 and period 2", suffix)]
fwrite(h4, file.path(out_dir, "h4_exploratory.tsv"), sep = "\t")
say("posthoc finished")
