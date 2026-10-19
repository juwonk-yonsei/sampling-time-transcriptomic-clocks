#!/usr/bin/env Rscript
# Stage 8 post-hoc tables. Does not write to results/confirm/.

suppressPackageStartupMessages(library(data.table))

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
root <- normalizePath(file.path(script_dir, ".."))
source(file.path(root, "scripts", "confirm", "run_pipeline.R"))
source(file.path(root, "scripts", "confirm", "run_locked_analysis.R"))

out_dir <- file.path(root, "results", "stage8")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
label <- "post hoc, exploratory"

say <- function(...) message(format(Sys.time(), "%H:%M:%S"), " ", paste(...))

quotes <- data.table(
  source = c(
    "Kervezee et al. 2019 methods",
    "Kervezee et al. 2019 methods",
    "Kervezee et al. 2019 figure 1 legend"
  ),
  citation = "Kervezee et al., J Biol Rhythms, PMID 30621487, DOI 10.1177/0748730418821776",
  sentence = c(
    "Briefly, participants maintained a regular, selfselected sleep/wake schedule with a nocturnal sleep period of 8 h for at least one week before admission to the laboratory, which was verified by daily sleep logs and actigraphy.",
    "Upon admission to the laboratory, participants were involved in a 24-h baseline measurement period that consisted of an 8-h sleep period during their habitual sleep times and a 16-h constant posture procedure, involving a constant semi-recumbent position, minimal activity, and hourly isocaloric snacks.",
    "relative clock time refers to the time of day relative to the start of each participant's habitual sleep period."
  ),
  habitual_sleep_start_mean_range = "확인 못 함",
  reason = "The paper states an 8 h habitual sleep period and defines relative clock time from its start. It does not state a clock time mean or range for habitual sleep onset. Protocol details are referred to Cuesta et al., 2017.",
  time_transform = "not performed",
  title = paste("L2 time-axis source", label)
)
fwrite(quotes, file.path(out_dir, "l2_sleep_start_search.tsv"), sep = "\t")

circ_sd_h <- function(hours) {
  ang <- hours * 2 * pi / 24
  rbar <- sqrt(mean(cos(ang))^2 + mean(sin(ang))^2)
  if (!is.finite(rbar) || rbar <= 0) return(NA_real_)
  sqrt(-2 * log(min(rbar, 1))) * 24 / (2 * pi)
}
signed_diff <- function(a, b) {
  d <- (a - b) %% 24
  ifelse(d > 12, d - 24, d)
}
peak_gene <- function(samples, expr, gene) {
  if (!(gene %in% rownames(expr))) {
    return(data.table(gene = gene, n_samples = 0L, n_subjects = 0L, peak_h = NA_real_, lrt_p = NA_real_, status = "absent"))
  }
  values <- as.numeric(expr[gene, samples$gsm])
  keep <- is.finite(values) & is.finite(samples$clock_time_h)
  used <- samples[keep]
  if (nrow(used) < 2L) {
    return(data.table(gene = gene, n_samples = nrow(used), n_subjects = uniqueN(used$subject_id), peak_h = NA_real_, lrt_p = NA_real_, status = "fewer than 2 finite samples"))
  }
  score <- data.table(gsm = used$gsm, clock = gene, value = values[keep])
  fit <- model_m(prepare_score(used, score, gene))
  data.table(gene = gene, n_samples = nrow(used), n_subjects = uniqueN(used$subject_id), peak_h = fit$peak_h, lrt_p = fit$lrt_p, status = "fit")
}

say("ARNTL alias")
l2_samples <- fread(file.path(root, "results", "confirm", "L2", "samples_scoring.tsv"))
l2_samples[, gsm := as.character(gsm)]
period1 <- l2_samples[period == "period1"]
period1[, clock_time_h := wake_ref_h %% 24]
l2_expr <- as.matrix(read_gene_matrix(file.path(root, "results", "confirm", "L2", "log2_genes.tsv")))
bmal <- peak_gene(period1, l2_expr, "BMAL1")
bmal[, alias_of := "ARNTL"]
bmal[, title := paste("L2 BMAL1 used as ARNTL", label)]
fwrite(bmal, file.path(out_dir, "l2_bmal1_as_arntl.tsv"), sep = "\t")
rm(l2_expr)

prev <- fread(file.path(root, "results", "posthoc_stage7c", "gene_peak_differences.tsv"))
prev[gene == "ARNTL", peak_l2 := bmal$peak_h]
prev[gene == "ARNTL", status_l2 := bmal$status]
prev[, signed_diff_h := signed_diff(peak_l2, peak_gse39445)]
prev[, circ_dist_0_12 := mapply(circ_dist_0_12, peak_l2, peak_gse39445)]
prev[, title := paste("clock-gene peaks with BMAL1 alias", label)]
finite_diff <- prev[is.finite(signed_diff_h), signed_diff_h]
summary_row <- data.table(
  n_genes = length(finite_diff),
  circ_mean_h = if (length(finite_diff)) circ_mean_h(finite_diff) else NA_real_,
  circ_sd_h = if (length(finite_diff)) circ_sd_h(finite_diff) else NA_real_,
  signed_diff_definition = "(peak_l2 - peak_gse39445) mod 24, then subtract 24 when the result is greater than 12; range (-12, 12]",
  circ_sd_definition = "sqrt(-2*log(R)) * 24 / (2*pi), R = mean resultant length of the signed differences",
  arntl_source = "L2 symbol BMAL1; GSE39445 symbol ARNTL",
  title = paste("circular mean and SD after ARNTL alias", label)
)
fwrite(prev, file.path(out_dir, "gene_peak_differences_arntl_alias.tsv"), sep = "\t")
fwrite(summary_row, file.path(out_dir, "gene_peak_circular_summary_arntl_alias.tsv"), sep = "\t")

say("L3 observed d")
pairs <- fread(file.path(root, "results", "confirm", "L3", "paired_deltas.tsv"))
f2 <- fread(file.path(root, "prereg", "frozen", "F2.tsv"))
set.seed(20261015)
d_rows <- lapply(c(primary_clocks, secondary_clocks), function(clock_id) {
  block <- pairs[clock == clock_id]
  delta <- block$morning - block$evening
  n <- length(delta)
  point_d <- mean(delta) / stats::sd(delta)
  boots <- replicate(2000L, {
    draw <- sample(delta, n, replace = TRUE)
    mean(draw) / stats::sd(draw)
  })
  qs <- quantile(boots, c(0.025, 0.975), type = 7, names = FALSE)
  info <- f2[clock == clock_id]
  df <- n - 1
  ncp <- point_d * sqrt(n)
  crit <- qt(0.975, df)
  power <- 1 - pt(crit, df, ncp) + pt(-crit, df, ncp)
  data.table(
    clock = clock_id, n_pairs = n, d = point_d, d_p025 = qs[1], d_p975 = qs[2],
    d_exp_median = info$d_exp_median, d_GSE113883 = info$d_GSE113883,
    predicted_sign = info$predicted_sign, power_alpha_0.05 = power,
    power_definition = "paired t, df=n-1, ncp=d*sqrt(n), two-sided alpha=0.05",
    title = paste("L3 observed d and development prediction", label)
  )
})
fwrite(rbindlist(d_rows), file.path(out_dir, "l3_d_vs_prediction.tsv"), sep = "\t")

say("D11 module model M")
f3 <- fread(file.path(root, "prereg", "frozen", "F3.tsv"))
d11_scores <- rbind(
  fread(file.path(root, "results", "confirm", "L1", "tage_scores_d11.tsv")),
  fread(file.path(root, "results", "confirm", "L1", "tage_components_d11.tsv")),
  fill = TRUE
)
d11_scores[, gsm := as.character(gsm)]
day1 <- fread(file.path(root, "results", "confirm", "L1", "samples_scoring.tsv"))
day1[, gsm := as.character(gsm)]
day1 <- day1[period == "D1"]
d11_rows <- lapply(f3$clock, function(clock_id) {
  one <- d11_scores[clock == clock_id]
  values <- one$value[match(day1$gsm, one$gsm)]
  keep <- is.finite(values)
  used <- day1[keep]
  if (nrow(used) < 2L) {
    return(data.table(clock = clock_id, n_samples = nrow(used), p = NA_real_, status = "fewer than 2 finite samples"))
  }
  fit <- model_m(prepare_score(used, one, clock_id))
  phi <- f3$phi_bar[match(clock_id, f3$clock)]
  dist <- circ_dist_0_12(fit$peak_h, phi)
  data.table(
    clock = clock_id, n_subjects = uniqueN(used$subject_id), n_samples = nrow(used),
    two_a = fit$two_a, peak_h = fit$peak_h, lrt_stat = fit$lrt_stat, p = fit$lrt_p,
    phi_bar = phi, circ_dist = dist,
    posthoc_confirmed = is.finite(fit$lrt_p) && fit$lrt_p < 0.05 && is.finite(dist) && dist <= 3,
    criterion = "LRT p < 0.05 and circular distance to phi_bar <= 3 h",
    status = "fit",
    title = paste("L1 F3 outcomes on D11 scores", label)
  )
})
fwrite(rbindlist(d11_rows, fill = TRUE), file.path(out_dir, "l1_f3_d11_model_m.tsv"), sep = "\t")

fit_clock_boot <- function(dataset, clock_id, samples, scores, role) {
  one <- scores[clock == clock_id]
  values <- one$value[match(samples$gsm, one$gsm)]
  keep <- is.finite(values) & is.finite(samples$clock_time_h)
  used <- samples[keep]
  dat <- prepare_score(used, one, clock_id)
  fit <- model_m(dat)
  if (identical(role, "development")) {
    frozen <- fread(file.path(root, "prereg", "frozen", "F1.tsv"))
    ds <- dataset
    ref <- frozen[dataset == ds & clock == clock_id]
    if (nrow(ref) != 1L || !is.finite(fit$two_a) || abs(fit$two_a - ref$two_a) / abs(ref$two_a) > 1e-4) {
      stop(paste("point two_a does not match F1", dataset, clock_id, fit$two_a, if (nrow(ref)) ref$two_a else NA))
    }
  }
  draws <- participant_draws(unique(as.character(used$subject_id)))
  boot <- boot_effects(dat, draws)
  row <- point_row(dataset, clock_id, role, fit)
  attach_boot(row, boot)
  row[, title := paste("cosinor effect", label)]
  row
}

say("development and L2 bootstrap")
dev_rows <- list()
spec <- list(
  list(dataset = "GSE39445", samples = file.path(root, "data", "processed", "GSE39445_samples.tsv"), expr = file.path(root, "data", "processed", "GSE39445_expr_gene_log2.tsv.gz")),
  list(dataset = "GSE56931", samples = file.path(root, "data", "processed", "GSE56931_samples.tsv"), expr = file.path(root, "data", "processed", "GSE56931_expr_gene_log2.tsv.gz"))
)
stage5 <- fread(file.path(root, "results", "stage5", "clock_scores.tsv"))
stage5[, gsm := as.character(gsm)]
for (item in spec) {
  say(item$dataset)
  samples <- fread(item$samples)
  samples <- samples[condition == "control"]
  samples[, gsm := as.character(gsm)]
  tage <- stage5[dataset == item$dataset & condition == "control" & clock %in% c("TAms_chrono", "TAms_mort"),
                 .(gsm, clock, value)]
  expr <- as.matrix(read_gene_matrix(item$expr))
  ra <- ra_linear_scores(expr[, samples$gsm, drop = FALSE])
  scores <- rbind(tage, ra[, .(gsm, clock, value)], fill = TRUE)
  rm(expr)
  gc()
  for (clock_id in c(primary_clocks, secondary_clocks)) {
    say(item$dataset, clock_id)
    dev_rows[[length(dev_rows) + 1L]] <- fit_clock_boot(item$dataset, clock_id, samples, scores, "development")
  }
}

say("GSE113883")
s113 <- fread(file.path(root, "results", "stage3", "GSE113883_samples.tsv"))
s113 <- s113[condition == "control"]
s113[, gsm := as.character(gsm)]
tage113 <- stage5[dataset == "GSE113883" & clock %in% c("TAms_chrono", "TAms_mort"), .(gsm, clock, value)]
ra113 <- fread(file.path(root, "results", "stage4", "rnaseq", "clock_scores.tsv"))
ra113 <- ra113[method == "lin", .(gsm = as.character(gsm), clock, value)]
scores113 <- rbind(tage113, ra113, fill = TRUE)
for (clock_id in c(primary_clocks, secondary_clocks)) {
  say("GSE113883", clock_id)
  dev_rows[[length(dev_rows) + 1L]] <- fit_clock_boot("GSE113883", clock_id, s113, scores113, "development")
}

say("L2 period 1")
l2_scores <- rbind(
  fread(file.path(root, "results", "posthoc_stage7c", "l2_tage_scores.tsv")),
  fread(file.path(root, "results", "posthoc_stage7c", "l2_tage_components.tsv")),
  fill = TRUE
)
l2_scores[, gsm := as.character(gsm)]
l2_expr <- as.matrix(read_gene_matrix(file.path(root, "results", "confirm", "L2", "log2_genes.tsv")))
l2_ra <- ra_linear_scores(l2_expr[, period1$gsm, drop = FALSE])
l2_scores <- rbind(l2_scores[, .(gsm, clock, value)], l2_ra[, .(gsm, clock, value)], fill = TRUE)
rm(l2_expr)
gc()
for (clock_id in c(primary_clocks, secondary_clocks)) {
  say("L2", clock_id)
  dev_rows[[length(dev_rows) + 1L]] <- fit_clock_boot("L2", clock_id, period1, l2_scores, "post hoc")
}
fwrite(rbindlist(dev_rows, fill = TRUE), file.path(out_dir, "cosinor_effects_boot.tsv"), sep = "\t")
say("stage8 posthoc finished")
