#!/usr/bin/env Rscript
# Expression-level pipeline rehearsals. Calls the confirmation functions.
# Does not open data/locked.

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(MCPcounter)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
root <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(root, "scripts", "confirm", "run_pipeline.R"))
source(file.path(root, "scripts", "confirm", "preprocess_l2.R"))
source(file.path(root, "scripts", "confirm", "preprocess_l3.R"))

python_bin <- "/HDD8T1/WORK/GeroSci/tod4_clock_env/bin/python"
f3 <- fread(file.path(root, "prereg", "frozen", "F3.tsv"))
f4 <- fread(file.path(root, "prereg", "frozen", "F4.tsv"))
f2 <- fread(file.path(root, "prereg", "frozen", "F2.tsv"))

read_gene_matrix <- function(path) {
  tab <- fread(path)
  gene <- tab[[1]]
  mat <- as.matrix(tab[, -1])
  storage.mode(mat) <- "double"
  rownames(mat) <- gene
  colnames(mat) <- names(tab)[-1]
  mat
}

finish_cosinor <- function(out_dir, samples, log_mat, gate, primary_set, drop_neutrophils, n_boot) {
  fwrite(gate, file.path(out_dir, "gate.tsv"), sep = "\t")
  if (!isTRUE(gate$dataset_pass[1])) {
    fwrite(data.table(status = "gate fail"), file.path(out_dir, "stopped.tsv"), sep = "\t")
    return(invisible(FALSE))
  }
  score_path <- file.path(out_dir, "scores_input.tsv.gz")
  fwrite(as.data.table(log_mat, keep.rownames = "gene"), score_path, sep = "\t")
  samples_path <- file.path(out_dir, "samples_input.tsv")
  fwrite(samples, samples_path, sep = "\t")
  scored <- score_tage(score_path, samples_path, "array", file.path(out_dir, "tage_scores.tsv"), file.path(out_dir, "components.tsv"))
  ra <- ra_linear_scores(log_mat)
  ra[, subject_id := samples$subject_id[match(gsm, samples$gsm)]]
  scores <- rbind(scored$scores, ra, fill = TRUE)
  flags <- flag_3mad(log_mat)
  fwrite(flags, file.path(out_dir, "flags_3mad.tsv"), sep = "\t")
  comp <- mcp_scores(log_mat)
  comp <- comp[samples$gsm, , drop = FALSE]
  h4 <- samples
  spec <- list(
    out_dir = out_dir,
    samples = samples,
    scores = scores,
    gate = gate,
    primary_set = primary_set,
    composition = comp,
    drop_neutrophils = drop_neutrophils,
    h4_samples = h4,
    module_scores = rbind(scored$scores[grepl("^TM_", clock)], scored$components),
    f3_clocks = f3$clock,
    flags = flags,
    n_boot = n_boot
  )
  analyse_cosinor_dataset(spec)
  sign <- merge(fread(file.path(out_dir, "P_point.tsv")), f2[, .(clock, predicted_sign)], by = "clock", all.x = TRUE)
  fwrite(sign, file.path(out_dir, "P_with_prediction.tsv"), sep = "\t")
}

run_r1 <- function() {
  out_dir <- file.path(root, "results", "stage6", "rehearsal", "R1")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  meta <- fread(file.path(root, "data", "processed", "GSE39445_samples.tsv"))
  keep <- which(is.finite(meta$clock_time_h) & meta$condition %in% c("control", "sleep_restriction"))
  meta <- meta[keep]
  meta[, analysis_set := fifelse(condition == "control", "Day1", "Day4")]
  meta[, day := fifelse(condition == "control", 0, 1)]
  samples <- meta[, .(gsm = as.character(gsm), subject_id = as.character(subject_id), clock_time_h, analysis_set, day, condition)]
  fwrite(samples, file.path(out_dir, "samples_input.tsv"), sep = "\t")
  status <- system2(python_bin, c(
    file.path(root, "scripts", "confirm", "preprocess_array.py"),
    "--matrix", file.path(root, "data", "raw", "GSE39445_series_matrix.txt.gz"),
    "--platform", file.path(root, "data", "raw", "gpl", "GPL15331_platform.txt"),
    "--symbol-column", "GENE_SYMBOL",
    "--samples", file.path(out_dir, "samples_input.tsv"),
    "--out-expr", file.path(out_dir, "gene_log2.tsv.gz"),
    "--out-audit", file.path(out_dir, "preprocess_audit.tsv")
  ))
  if (!identical(status, 0L)) stop("R1 preprocess failed")
  log_mat <- read_gene_matrix(file.path(out_dir, "gene_log2.tsv.gz"))
  day1 <- samples[analysis_set == "Day1"]
  gate <- gate_cosinor(log_mat, day1, 4, 12)
  finish_cosinor(out_dir, samples, log_mat, gate, "Day1", FALSE, 2000)
}

run_r2 <- function() {
  out_dir <- file.path(root, "results", "stage6", "rehearsal", "R2")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  suppressPackageStartupMessages(library(clariomshumanhttranscriptcluster.db))
  ids <- AnnotationDbi::keys(clariomshumanhttranscriptcluster.db, keytype = "PROBEID")
  meta <- fread(file.path(root, "data", "processed", "GSE122725_samples.tsv"))
  wake <- fread(file.path(root, "results", "stage2", "GSE122725_wake_ref.tsv"))
  meta[, gsm := as.character(gsm)]
  wake[, gsm := as.character(gsm)]
  meta <- merge(meta, wake[, .(gsm, wake_ref_h)], by = "gsm")
  meta[, clock_time_h := wake_ref_h %% 24]
  meta[, analysis_set := fifelse(period == "period1" & condition == "control", "period1",
                                 fifelse(period == "period2", "period2", "other"))]
  meta[, day := fifelse(analysis_set == "period2", 1, 0)]
  samples <- meta[, .(gsm, subject_id = as.character(subject_id), clock_time_h, analysis_set, day, condition, period, light_group)]
  fwrite(samples, file.path(out_dir, "samples_input.tsv"), sep = "\t")
  light <- unique(samples[, .(subject_id, light_group)])
  fwrite(light[, .N, by = light_group], file.path(out_dir, "light_groups.tsv"), sep = "\t")
  set.seed(20261015)
  raw <- matrix(stats::rnorm(length(ids) * nrow(samples)), nrow = length(ids), dimnames = list(ids, samples$gsm))
  scale <- geo2r_scale_matrix(raw)
  log_mat <- to_log2_matrix(raw, scale)
  collapsed <- collapse_clariom(log_mat)
  fwrite(data.table(
    scale_call = scale, n_probes = collapsed$n_probes, n_multi_symbol = collapsed$n_multi,
    n_genes = collapsed$n_genes, seed = 20261015, probe_matrix_stored = FALSE
  ), file.path(out_dir, "preprocess_audit.tsv"), sep = "\t")
  gene <- collapsed$matrix
  fwrite(as.data.table(gene, keep.rownames = "gene"), file.path(out_dir, "gene_log2.tsv.gz"), sep = "\t")
  period1 <- samples[analysis_set == "period1"]
  gate <- gate_cosinor(gene, period1, 18, 2)
  if (!isTRUE(gate$dataset_pass[1])) {
    fwrite(gate, file.path(out_dir, "gate.tsv"), sep = "\t")
    fwrite(data.table(status = "gate fail"), file.path(out_dir, "stopped.tsv"), sep = "\t")
    return(invisible(FALSE))
  }
  finish_cosinor(out_dir, samples[analysis_set %in% c("period1", "period2")], gene[, samples[analysis_set %in% c("period1", "period2")]$gsm, drop = FALSE], gate, "period1", TRUE, 2000)
}

run_r3_one <- function(label, counts_path) {
  out_dir <- file.path(root, "results", "stage6", "rehearsal", "R3", label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  parsed <- read_l3_counts(counts_path, file.path(root, "prereg", "frozen", "gencode_v38_basic_transcripts.tsv"))
  fwrite(data.table(
    branch = parsed$layout$branch, gene_column = parsed$layout$gene_column,
    transcript_column = parsed$layout$transcript_column, n_genes = parsed$n_genes
  ), file.path(out_dir, "layout.tsv"), sep = "\t")
  samples <- fread(file.path(root, "results", "stage6", "rehearsal", "R3", "inputs", "samples.tsv"))
  cov <- fread(file.path(root, "results", "stage6", "rehearsal", "R3", "inputs", "covariates.tsv"))
  samples[, gsm := as.character(gsm)]
  samples[, subject_id := as.character(subject_id)]
  cov[, subject_id := as.character(subject_id)]
  counts <- parsed$counts[, samples$gsm, drop = FALSE]
  cpm <- log2_cpm(counts)
  gate <- gate_l3(cpm, samples)
  fwrite(gate, file.path(out_dir, "gate.tsv"), sep = "\t")
  if (!isTRUE(gate$dataset_pass[1])) {
    fwrite(data.table(status = "gate fail"), file.path(out_dir, "stopped.tsv"), sep = "\t")
    return(invisible(FALSE))
  }
  fpkm <- counts_to_log2_fpkm(counts)
  count_path <- file.path(out_dir, "gene_counts.tsv.gz")
  fwrite(as.data.table(counts, keep.rownames = "gene"), count_path, sep = "\t")
  fwrite(samples, file.path(out_dir, "samples_input.tsv"), sep = "\t")
  scored <- score_tage(count_path, file.path(out_dir, "samples_input.tsv"), "counts", file.path(out_dir, "tage_scores.tsv"), file.path(out_dir, "components.tsv"))
  ra <- ra_linear_scores(fpkm)
  scores <- rbind(scored$scores, ra, fill = TRUE)
  flags <- flag_3mad(fpkm)
  fwrite(flags, file.path(out_dir, "flags_3mad.tsv"), sep = "\t")
  comp <- mcp_scores(fpkm)
  wide_subjects <- NULL
  delta_comp <- NULL
  labels <- intersect(blood_labels, colnames(comp))
  pieces <- list()
  for (label_i in labels) {
    wide <- data.table(subject_id = samples$subject_id, phase = samples$phase, value = comp[samples$gsm, label_i])
    morning <- wide[phase == "morning", .(morning = mean(value)), by = subject_id]
    evening <- wide[phase == "evening", .(evening = mean(value)), by = subject_id]
    both <- merge(morning, evening, by = "subject_id")
    pieces[[label_i]] <- both$morning - both$evening
    wide_subjects <- both$subject_id
  }
  delta_comp <- do.call(cbind, pieces)
  rownames(delta_comp) <- wide_subjects
  c_values <- as.list(f4$c)
  names(c_values) <- f4$clock
  spec <- list(
    out_dir = out_dir,
    samples = samples,
    scores = scores,
    gate = gate,
    covariates = cov,
    c_values = c_values,
    delta_comp = delta_comp,
    module_scores = rbind(scored$scores[grepl("^TM_", clock)], scored$components),
    f3_clocks = f3$clock,
    n_boot = 2000,
    batch2_n = 0L,
    combat_path = NULL
  )
  analyse_paired_dataset(spec)
  pred <- fread(file.path(root, "prereg", "frozen", "F2.tsv"))
  tab <- fread(file.path(out_dir, "paired_tests.tsv"))
  tab <- merge(tab, pred[, .(clock, predicted_sign)], by = "clock", all.x = TRUE)
  tab[, sign_agree := sign(mean_delta) == fifelse(predicted_sign == "+", 1, fifelse(predicted_sign == "-", -1, 0))]
  fwrite(tab, file.path(out_dir, "paired_with_prediction.tsv"), sep = "\t")
}

args <- commandArgs(trailingOnly = TRUE)
which_run <- if ("--which" %in% args) args[match("--which", args) + 1] else stop("--which is required")
if (which_run == "R1") run_r1()
if (which_run == "R2") run_r2()
if (which_run == "R3") {
  status <- system2(python_bin, file.path(root, "scripts", "42_partD_make_r3.py"))
  if (!identical(status, 0L)) stop("R3 input build failed")
  run_r3_one("gene_column", file.path(root, "results", "stage6", "rehearsal", "R3", "inputs", "counts_gene_column.tsv.gz"))
  run_r3_one("transcript_column", file.path(root, "results", "stage6", "rehearsal", "R3", "inputs", "counts_transcript_column.tsv.gz"))
}
cat(which_run, "pipeline rehearsal finished\n")
