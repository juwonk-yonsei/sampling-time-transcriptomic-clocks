# Steps 3-12 of the locked confirmation. Sourced by run_locked().
# Test definitions stay in lib_tests.R and lib_cosinor.R.

source(file.path(root, "scripts", "confirm", "preprocess_l2.R"))
source(file.path(root, "scripts", "confirm", "preprocess_l3.R"))

read_geo_series_matrix <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (!length(line)) stop("series_matrix_table_begin was not found")
    if (identical(line, "!series_matrix_table_begin")) break
  }
  lines <- readLines(con, warn = FALSE)
  end <- match("!series_matrix_table_end", lines)
  if (!is.na(end)) lines <- lines[seq_len(end - 1L)]
  lines <- lines[nzchar(lines)]
  data.table::fread(paste(lines, collapse = "\n"), sep = "\t")
}

l3_layout_from_header <- function(header_line) {
  sep <- if (grepl("\t", header_line, fixed = TRUE)) "\t" else ","
  cols <- strsplit(header_line, sep, fixed = TRUE)[[1]]
  cols <- gsub('^"|"$', "", cols)
  gene_names <- c("gene_id", "gene_name", "gene_symbol", "GENE", "Gene", "gene", "symbol", "SYMBOL")
  transcript_names <- c("transcript_id", "transcript", "Transcript")
  gene_col <- gene_names[gene_names %in% cols]
  transcript_col <- transcript_names[transcript_names %in% cols]
  out <- list(
    columns = cols,
    gene_column = if (length(gene_col)) gene_col[1] else NA_character_,
    transcript_column = if (length(transcript_col)) transcript_col[1] else NA_character_,
    branch = if (length(gene_col)) "gene_column" else "transcript_to_gene",
    separator = if (identical(sep, "\t")) "tab" else "comma"
  )
  l3_layout_last <<- out
  out
}

say <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), paste(..., collapse = " "), "\n")
  flush.console()
}

verify_frozen_manifest <- function() {
  manifest <- fread(file.path(root, "prereg", "manifest.tsv"))
  changed <- fread(file.path(root, "results", "stage7b", "changed_files.tsv"))
  bad <- character()
  for (i in seq_len(nrow(manifest))) {
    rel <- manifest$path[i]
    path <- file.path(root, rel)
    size <- file.info(path)$size
    hash <- strsplit(system2("sha256sum", path, stdout = TRUE), " ", fixed = TRUE)[[1]][1]
    hit <- which(changed$path == rel)
    if (length(hit)) {
      if (!identical(hash, changed$sha256[hit[1]])) bad <- c(bad, paste("changed hash", rel))
      say("manifest changed file", rel, hash)
      next
    }
    if (!identical(as.numeric(size), as.numeric(manifest$bytes[i])) || !identical(hash, manifest$sha256[i])) {
      bad <- c(bad, rel)
    }
  }
  if (length(bad)) stop(paste("manifest mismatch:", paste(bad, collapse = ", ")))
  invisible(TRUE)
}

d1_impute_matrix <- function(log_mat) {
  frac_bad <- rowMeans(!is.finite(log_mat))
  keep <- frac_bad <= 0.5
  mat <- log_mat[keep, , drop = FALSE]
  n_nonfinite <- sum(!is.finite(mat))
  mat[!is.finite(mat)] <- NA_real_
  suppressPackageStartupMessages(library(impute))
  imputed <- impute::impute.knn(mat, k = 10)$data
  rownames(imputed) <- rownames(mat)
  colnames(imputed) <- colnames(mat)
  if (any(!is.finite(imputed))) stop("D1 imputation left a non-finite value")
  list(
    matrix = imputed, n_genes_in = nrow(log_mat), n_genes_dropped = sum(!keep),
    n_genes_kept = sum(keep), n_values_imputed = n_nonfinite
  )
}

d1_complete_matrix <- function(log_mat) {
  keep <- rowSums(!is.finite(log_mat)) == 0
  list(matrix = log_mat[keep, , drop = FALSE], n_genes_dropped = sum(!keep), n_genes_kept = sum(keep))
}

drop_all_nonfinite_rows <- function(log_mat) {
  keep <- rowSums(is.finite(log_mat)) > 0
  list(matrix = log_mat[keep, , drop = FALSE], n_genes_dropped = sum(!keep), n_genes_kept = sum(keep))
}

clean_token <- function(x) gsub('^"|"$', "", as.character(x))

write_gene_matrix <- function(mat, path) {
  dt <- as.data.table(mat, keep.rownames = "gene")
  fwrite(dt, path, sep = "\t")
}

read_gene_matrix <- function(path) {
  dt <- fread(path)
  genes <- as.character(dt[[1]])
  mat <- as.matrix(dt[, -1])
  storage.mode(mat) <- "double"
  rownames(mat) <- genes
  colnames(mat) <- as.character(colnames(mat))
  mat
}

series_table <- function(dt) {
  setnames(dt, clean_token(names(dt)))
  id <- names(dt)[1]
  set(dt, j = id, value = clean_token(dt[[id]]))
  dt
}

expect_n <- function(label, n_samples, n_subjects, want_samples, want_subjects) {
  if (!identical(n_samples, want_samples) || !identical(n_subjects, want_subjects)) {
    stop(paste(label, n_samples, n_subjects, "expected", want_samples, want_subjects))
  }
}

point_row <- function(dataset, clock, family, fit) {
  sigma <- fit$sigma_reml
  ratio <- if (is.finite(sigma) && sigma > 0 && is.finite(fit$two_a)) fit$two_a / sigma else NA_real_
  data.table(
    dataset = dataset, clock = clock, family = family,
    n_subjects = fit$n_subjects, n_samples = fit$n_samples,
    beta_cos = fit$beta_cos, beta_sin = fit$beta_sin,
    two_a = fit$two_a, r2_within = fit$r2_within,
    sigma_reml = sigma, two_a_over_sigma = ratio,
    peak_h = fit$peak_h, lrt_stat = fit$lrt_stat, p = fit$lrt_p
  )
}

boot_effects <- function(dat, draws) {
  n <- uniqueN(draws$draw)
  two_a <- r2 <- ratio <- peak <- rep(NA_real_, n)
  for (b in seq_len(n)) {
    one <- draws[draw == b]
    boot_dat <- apply_draw(dat, one$subject_id, one$new_subject_id)
    fit <- model_m(boot_dat)
    ok <- isTRUE(fit$reml_converged) && isTRUE(fit$ml_converged) && isTRUE(fit$null_converged) &&
      is.finite(fit$sigma_reml) && fit$sigma_reml > 0 && is.finite(fit$two_a)
    if (!ok) next
    two_a[b] <- fit$two_a
    r2[b] <- fit$r2_within
    ratio[b] <- fit$two_a / fit$sigma_reml
    peak[b] <- fit$peak_h
    if (b %% 500L == 0L) say("bootstrap", b, "of", n)
  }
  list(two_a = two_a, r2 = r2, ratio = ratio, peak = peak, n_boot = n, n_converged = sum(is.finite(two_a)))
}

attach_boot <- function(row, boot) {
  ci_a <- percentile_ci(boot$two_a)
  ci_r <- percentile_ci(boot$r2)
  ci_s <- percentile_ci(boot$ratio)
  pk <- peak_interval(boot$peak)
  row[, two_a_p025 := ci_a[1]]
  row[, two_a_p975 := ci_a[2]]
  row[, r2_p025 := ci_r[1]]
  row[, r2_p975 := ci_r[2]]
  row[, ratio_p025 := ci_s[1]]
  row[, ratio_p975 := ci_s[2]]
  row[, peak_center := unname(pk["center"])]
  row[, peak_halfwidth_95 := unname(pk["halfwidth_95"])]
  row[, n_boot := boot$n_boot]
  row[, n_converged := boot$n_converged]
  row
}

apply_holm <- function(tab) {
  tab <- copy(tab)
  tab[, q_holm := NA_real_]
  tab[, decision := "gate fail"]
  tested <- which(is.finite(tab$p))
  if (length(tested)) {
    q <- holm_family(tab$p[tested])
    set(tab, i = tested, j = "q_holm", value = q)
    decision <- ifelse(q < 0.05, "significant", "not significant")
    set(tab, i = tested, j = "decision", value = decision)
  }
  tab
}

aggregate_loaded_l3 <- function(raw, gencode_path, layout) {
  gencode <- fread(gencode_path)
  if (identical(layout$branch, "gene_column")) {
    ids <- as.character(raw[[layout$gene_column]])
    sample_cols <- setdiff(names(raw), c(layout$gene_column, layout$transcript_column))
    sample_cols <- sample_cols[!is.na(sample_cols)]
    kind <- if (all(grepl("^ENSG[0-9]", base_id(ids)))) "ensembl" else "symbol"
    if (identical(kind, "ensembl")) {
      map <- unique(gencode[, .(gene_id_base, gene_name)])
      names_out <- map$gene_name[match(base_id(ids), map$gene_id_base)]
    } else {
      names_out <- ids
    }
  } else {
    id_col <- layout$transcript_column
    if (is.na(id_col) || !nzchar(id_col) || !(id_col %in% names(raw))) id_col <- names(raw)[1]
    ids <- as.character(raw[[id_col]])
    sample_cols <- setdiff(names(raw), id_col)
    map <- unique(gencode[, .(transcript_id_base, gene_name)])
    names_out <- map$gene_name[match(base_id(ids), map$transcript_id_base)]
    kind <- "transcript"
  }
  keep <- !is.na(names_out) & nzchar(names_out)
  mat <- as.matrix(raw[keep, sample_cols, with = FALSE])
  storage.mode(mat) <- "double"
  counts <- sum_by_gene(mat, names_out[keep])
  list(
    counts = counts, kind = kind, n_rows = nrow(raw), n_mapped_rows = sum(keep),
    n_genes = nrow(counts), id_column = if (exists("id_col")) id_col else layout$gene_column
  )
}

parse_l3_columns <- function(columns) {
  pat <- "^([0-9]+)_(D|N)(?:_batch([12]))?$"
  hit <- regmatches(columns, regexec(pat, columns))
  if (any(lengths(hit) == 0L)) stop("an L3 column name does not match the subject-time pattern")
  rows <- vector("list", length(columns))
  for (i in seq_along(columns)) {
    g <- hit[[i]]
    batch <- if (length(g) >= 4L && nzchar(g[4])) g[4] else "1"
    rows[[i]] <- data.table(
      gsm = columns[i],
      subject_id = paste0("subject", g[2]),
      phase = if (g[3] == "D") "morning" else "evening",
      clock_time_h = if (g[3] == "D") 8 else 20,
      batch = batch
    )
  }
  rbindlist(rows)
}

condition_effects <- function(dat) {
  rows <- list()
  for (level in c(0, 1)) {
    part <- dat[dat$day == level, , drop = FALSE]
    fit <- model_m(part)
    rows[[as.character(level)]] <- data.table(
      day = level, n_subjects = fit$n_subjects, n_samples = fit$n_samples,
      two_a = fit$two_a, peak_h = fit$peak_h
    )
  }
  rbindlist(rows)
}

mean_difference <- function(dat) {
  mean(dat$score[dat$day == 1]) - mean(dat$score[dat$day == 0])
}

fit_clocks_cosinor <- function(dataset, samples, scores, clocks, family, draws, out_dir) {
  rows <- list()
  for (clock in clocks) {
    say(dataset, family, clock)
    dat <- prepare_score(samples, scores, clock)
    fwrite(
      as.data.table(dat)[, .(subject_id, clock_time_h, cos_t, sin_t, score)],
      file.path(out_dir, paste0("frame_", family, "_", clock, ".tsv")),
      sep = "\t"
    )
    fit <- model_m(dat)
    row <- point_row(dataset, clock, family, fit)
    boot <- boot_effects(dat, draws)
    rows[[clock]] <- attach_boot(row, boot)
  }
  rbindlist(rows)
}

module_cosinor <- function(samples, scores, clocks) {
  rows <- list()
  for (clock in clocks) {
    if (!any(scores$clock == clock)) next
    dat <- prepare_score(samples, scores, clock)
    fit <- model_m(dat)
    row <- point_row("L1orL2", clock, "S4", fit)
    rows[[clock]] <- row
  }
  rbindlist(rows)
}

h4_for <- function(dataset, samples, scores, clocks) {
  rows <- list()
  cond_rows <- list()
  for (clock in clocks) {
    dat <- prepare_score(samples, scores, clock)
    dat$day <- samples$day
    fit <- h4_interaction(dat)
    fit$dataset <- dataset
    fit$clock <- clock
    fit$p <- fit$lrt_p
    fit$mean_difference <- mean_difference(dat)
    rows[[clock]] <- fit
    cond <- condition_effects(dat)
    cond$dataset <- dataset
    cond$clock <- clock
    cond_rows[[clock]] <- cond
  }
  list(tests = rbindlist(rows), conditions = rbindlist(cond_rows))
}

random_slopes <- function(samples, scores, clocks) {
  rows <- list()
  for (clock in clocks) {
    dat <- prepare_score(samples, scores, clock)
    fit <- quiet_lmer(score ~ cos_t + sin_t + (1 + cos_t + sin_t | subject_id), dat, TRUE)
    ok <- !inherits(fit, "error") && optimizer_converged(fit)
    two_a <- NA_real_
    peak <- NA_real_
    if (ok) {
      beta <- lme4::fixef(fit)
      cos_b <- unname(beta["cos_t"])
      sin_b <- unname(beta["sin_t"])
      two_a <- 2 * sqrt(cos_b^2 + sin_b^2)
      peak <- peak_0_24(sin_b, cos_b)
    }
    rows[[clock]] <- data.table(clock = clock, converged = ok, two_a = two_a, peak_h = peak)
  }
  rbindlist(rows)
}

drop_flag_cosinor <- function(samples, scores, flags, clocks) {
  keep <- samples$gsm %in% flags$gsm[!flags$flag_3mad]
  used <- samples[keep, ]
  rows <- list()
  for (clock in clocks) {
    dat <- prepare_score(used, scores, clock)
    fit <- model_m(dat)
    rows[[clock]] <- point_row("sensitivity", clock, "drop_flags", fit)
  }
  rbindlist(rows)
}

composition_matrix <- function(log_mat, samples) {
  est <- mcp_scores(log_mat[, samples$gsm, drop = FALSE])
  if (!identical(rownames(est), samples$gsm)) {
    if (!setequal(rownames(est), samples$gsm)) stop("MCPcounter sample names do not match")
    est <- est[samples$gsm, , drop = FALSE]
  }
  missing <- setdiff(blood_labels, colnames(est))
  if (length(missing)) stop(paste("MCPcounter is missing", paste(missing, collapse = ", ")))
  est[, blood_labels, drop = FALSE]
}

h2_after_holm <- function(dataset, samples, scores, composition, clocks, q_tab, draws, drop_neutrophils) {
  comp_cols <- blood_labels
  if (drop_neutrophils) comp_cols <- setdiff(comp_cols, "Neutrophils")
  rows <- list()
  for (clock in clocks) {
    idx <- which(q_tab$dataset == dataset & q_tab$clock == clock)
    if (!length(idx) || !isTRUE(q_tab$decision[idx] == "significant")) next
    say("H2", dataset, clock)
    dat <- prepare_score(samples, scores, clock)
    comp <- composition[match(dat$gsm, rownames(composition)), comp_cols, drop = FALSE]
    if (anyNA(comp)) stop("composition rows do not match the score samples")
    for (col in comp_cols) dat[[col]] <- comp[[col]]
    h2 <- h2_l1l2(dat, comp_cols, draws)
    h2$point$dataset <- dataset
    h2$point$clock <- clock
    h2$point$decision <- if (is.finite(h2$point$p_med_p025) && h2$point$p_med_p025 > 0) "supported" else "not supported"
    rows[[clock]] <- h2$point
  }
  if (length(rows)) rbindlist(rows, fill = TRUE) else data.table()
}

paired_composition_delta <- function(samples, composition) {
  subjects <- NULL
  cols <- list()
  for (label in colnames(composition)) {
    wide <- data.table(
      subject_id = samples$subject_id, phase = samples$phase,
      value = composition[samples$gsm, label]
    )
    morning <- wide[phase == "morning", .(morning = mean(value)), by = subject_id]
    evening <- wide[phase == "evening", .(evening = mean(value)), by = subject_id]
    both <- merge(morning, evening, by = "subject_id")
    both <- both[order(subject_id)]
    if (is.null(subjects)) subjects <- both$subject_id
    if (!identical(subjects, both$subject_id)) stop("composition subject order changed")
    cols[[label]] <- both$morning - both$evening
  }
  list(subjects = subjects, delta = as.data.frame(cols, check.names = FALSE))
}

run_confirmation <- function() {
  verify_frozen_manifest()
  out_root <- file.path(root, "results", "confirm")
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  fwrite(data.table(run = "run4", started = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")), file.path(out_root, "run_id.tsv"), sep = "\t")
  f1 <- fread(file.path(root, "prereg", "frozen", "F1.tsv"))
  f2 <- fread(file.path(root, "prereg", "frozen", "F2.tsv"))
  f3 <- fread(file.path(root, "prereg", "frozen", "F3.tsv"))
  f4 <- fread(file.path(root, "prereg", "frozen", "F4.tsv"), colClasses = list(character = "c"))
  if (!identical(nrow(f3), 8L)) stop("F3 does not have 8 outcomes")
  phi <- f1[, .(phi_bar = unique(phi_bar)), by = clock]
  if (any(phi[, .N, by = clock]$N != 1L)) stop("a clock has more than one phi_bar")
  c_values <- setNames(as.numeric(f4$c), f4$clock)
  c_text <- setNames(f4$c, f4$clock)

  say("L1 load")
  l1_dir <- file.path(out_root, "L1")
  dir.create(l1_dir, recursive = TRUE, showWarnings = FALSE)
  l1_meta <- fread(file.path(root, "data", "processed", "GSE48113_samples.tsv"))
  l1_meta <- l1_meta[!gsm %in% c("GSM1168837", "GSM1168735")]
  day1 <- l1_meta[condition == "control" & period == "D1"]
  day4 <- l1_meta[period == "D4"]
  expect_n("L1 day1", nrow(day1), uniqueN(day1$subject_id), 147L, 22L)
  expect_n("L1 day4", nrow(day4), uniqueN(day4$subject_id), 138L, 21L)
  day1[, clock_time_h := as.numeric(clock_time_h)]
  day4[, clock_time_h := as.numeric(clock_time_h)]
  day1[, analysis_set := "day1"]
  day4[, analysis_set := "day4"]
  l1_samples <- rbind(day1, day4, fill = TRUE)
  fwrite(l1_samples, file.path(l1_dir, "samples_scoring.tsv"), sep = "\t")
  l1_raw <- series_table(load_locked_expression("L1", root))
  l1_cache <- file.path(out_root, "cache")
  dir.create(l1_cache, recursive = TRUE, showWarnings = FALSE)
  l1_series_path <- file.path(l1_cache, "L1_series_table.tsv.gz")
  fwrite(l1_raw, l1_series_path, sep = "\t")
  rm(l1_raw)
  gc()
  say("L1 preprocess")
  status <- system2(python_bin, c(
    file.path(root, "scripts", "confirm", "preprocess_array.py"),
    "--matrix", l1_series_path,
    "--platform", file.path(root, "data", "raw", "gpl", "GPL15331_platform.txt"),
    "--symbol-column", "GENE_SYMBOL",
    "--samples", file.path(l1_dir, "samples_scoring.tsv"),
    "--out-expr", file.path(l1_dir, "log2_genes.tsv.gz"),
    "--out-audit", file.path(l1_dir, "preprocess_audit.tsv")
  ))
  if (!identical(status, 0L)) stop("L1 preprocessing failed")
  unlink(l1_series_path)
  l1_log <- read_gene_matrix(file.path(l1_dir, "log2_genes.tsv.gz"))
  l1_flags <- flag_3mad(l1_log[, l1_samples$gsm, drop = FALSE])
  fwrite(l1_flags, file.path(l1_dir, "flags_3mad.tsv"), sep = "\t")
  l1_gate_samples <- day1[, .(gsm, subject_id, clock_time_h = as.numeric(clock_time_h))]
  l1_gate <- gate_cosinor(l1_log, l1_gate_samples, 4, 12)
  fwrite(l1_gate, file.path(l1_dir, "gate.tsv"), sep = "\t")
  say("L1 gate", l1_gate$dataset_pass[1])
  l1_pass <- isTRUE(l1_gate$dataset_pass[1])
  l1_scores <- NULL
  l1_components <- NULL
  l1_comp <- NULL
  l1_draws <- NULL
  if (l1_pass) {
    scored <- score_tage(
      file.path(l1_dir, "log2_genes.tsv.gz"),
      file.path(l1_dir, "samples_scoring.tsv"),
      "array", file.path(l1_dir, "tage_scores.tsv"), file.path(l1_dir, "tage_components.tsv")
    )
    say("L1 D1 impute")
    l1_for_ra <- l1_log[, l1_samples$gsm, drop = FALSE]
    l1_d1 <- d1_impute_matrix(l1_for_ra)
    fwrite(data.table(
      n_genes_in = l1_d1$n_genes_in, n_genes_dropped = l1_d1$n_genes_dropped,
      n_genes_kept = l1_d1$n_genes_kept, n_values_imputed = l1_d1$n_values_imputed, k = 10
    ), file.path(l1_dir, "d1_impute_audit.tsv"), sep = "\t")
    write_gene_matrix(l1_d1$matrix, file.path(l1_dir, "log2_genes_d1_imputed.tsv"))
    ra <- ra_linear_scores(l1_d1$matrix)
    if (any(!is.finite(ra$value))) stop("L1 RA primary scores are not finite after D1")
    fwrite(ra, file.path(l1_dir, "ra_scores.tsv"), sep = "\t")
    l1_cc <- d1_complete_matrix(l1_for_ra)
    fwrite(data.table(n_genes_dropped = l1_cc$n_genes_dropped, n_genes_kept = l1_cc$n_genes_kept), file.path(l1_dir, "d1_complete_audit.tsv"), sep = "\t")
    ra_cc <- ra_linear_scores(l1_cc$matrix)
    fwrite(ra_cc, file.path(l1_dir, "ra_scores_complete_case.tsv"), sep = "\t")
    l1_scores <- rbind(scored$scores, ra, fill = TRUE)
    l1_scores[, gsm := as.character(gsm)]
    l1_components <- scored$components
    l1_components[, gsm := as.character(gsm)]
    l1_comp <- composition_matrix(l1_log, l1_samples)
    l1_draws <- participant_draws(unique(day1$subject_id), 2000, 20261015)
    fwrite(l1_draws, file.path(l1_dir, "boot_draws.tsv"), sep = "\t")
  } else {
    fwrite(data.table(status = "gate fail"), file.path(l1_dir, "stopped.tsv"), sep = "\t")
  }

  say("L2 load")
  l2_dir <- file.path(out_root, "L2")
  dir.create(l2_dir, recursive = TRUE, showWarnings = FALSE)
  l2_meta <- fread(file.path(root, "data", "processed", "GSE122725_samples.tsv"))
  wake <- fread(file.path(root, "results", "stage2", "GSE122725_wake_ref.tsv"))
  l2_meta <- merge(l2_meta, wake[, .(gsm, wake_ref_h)], by = "gsm")
  if (anyNA(l2_meta$wake_ref_h)) stop("L2 wake reference has a missing value")
  l2_meta[, clock_time_h := as.numeric(wake_ref_h) %% 24]
  period1 <- l2_meta[condition == "control" & period == "period1"]
  period2 <- l2_meta[condition == "night_oriented" & period == "period2"]
  expect_n("L2 period1", nrow(period1), uniqueN(period1$subject_id), 70L, 11L)
  expect_n("L2 period2", nrow(period2), uniqueN(period2$subject_id), 71L, 11L)
  light <- period1[, .(n_participants = uniqueN(subject_id)), by = light_group]
  fwrite(light, file.path(l2_dir, "light_groups.tsv"), sep = "\t")
  l2_samples <- rbind(cbind(period1, analysis_set = "period1"), cbind(period2, analysis_set = "period2"), fill = TRUE)
  fwrite(l2_samples, file.path(l2_dir, "samples_scoring.tsv"), sep = "\t")
  l2_raw <- series_table(load_locked_expression("L2", root))
  id_col <- names(l2_raw)[1]
  l2_mat <- as.matrix(l2_raw[, -1, with = FALSE])
  storage.mode(l2_mat) <- "double"
  rownames(l2_mat) <- as.character(l2_raw[[id_col]])
  colnames(l2_mat) <- names(l2_raw)[-1]
  rm(l2_raw)
  gc()
  l2_mat <- l2_mat[, l2_samples$gsm, drop = FALSE]
  l2_scale <- geo2r_scale_matrix(l2_mat)
  l2_log_probes <- to_log2_matrix(l2_mat, l2_scale)
  l2_collapsed <- collapse_clariom(l2_log_probes)
  fwrite(data.table(
    scale_call = l2_scale, n_probes = l2_collapsed$n_probes,
    n_multi = l2_collapsed$n_multi, n_genes = l2_collapsed$n_genes, n_samples = ncol(l2_collapsed$matrix)
  ), file.path(l2_dir, "preprocess_audit.tsv"), sep = "\t")
  write_gene_matrix(l2_collapsed$matrix, file.path(l2_dir, "log2_genes.tsv"))
  l2_log <- l2_collapsed$matrix
  rm(l2_mat, l2_log_probes)
  gc()
  l2_flags <- flag_3mad(l2_log)
  fwrite(l2_flags, file.path(l2_dir, "flags_3mad.tsv"), sep = "\t")
  l2_gate_samples <- period1[, .(gsm, subject_id, clock_time_h)]
  l2_gate <- gate_cosinor(l2_log, l2_gate_samples, 18, 2)
  fwrite(l2_gate, file.path(l2_dir, "gate.tsv"), sep = "\t")
  say("L2 gate", l2_gate$dataset_pass[1])
  l2_pass <- isTRUE(l2_gate$dataset_pass[1])
  l2_scores <- NULL
  l2_components <- NULL
  l2_comp <- NULL
  l2_draws <- NULL
  if (l2_pass) {
    scored <- score_tage(
      file.path(l2_dir, "log2_genes.tsv"),
      file.path(l2_dir, "samples_scoring.tsv"),
      "array", file.path(l2_dir, "tage_scores.tsv"), file.path(l2_dir, "tage_components.tsv")
    )
    ra <- ra_linear_scores(l2_log[, l2_samples$gsm, drop = FALSE])
    fwrite(ra, file.path(l2_dir, "ra_scores.tsv"), sep = "\t")
    l2_scores <- rbind(scored$scores, ra, fill = TRUE)
    l2_scores[, gsm := as.character(gsm)]
    l2_components <- scored$components
    l2_components[, gsm := as.character(gsm)]
    l2_comp <- composition_matrix(l2_log, l2_samples)
    l2_draws <- participant_draws(unique(period1$subject_id), 2000, 20261015)
    fwrite(l2_draws, file.path(l2_dir, "boot_draws.tsv"), sep = "\t")
  } else {
    fwrite(data.table(status = "gate fail"), file.path(l2_dir, "stopped.tsv"), sep = "\t")
  }

  say("L3 load")
  l3_dir <- file.path(out_root, "L3")
  dir.create(l3_dir, recursive = TRUE, showWarnings = FALSE)
  l3_raw <- load_locked_expression("L3", root)
  l3_layout <- l3_layout_last
  fwrite(data.table(column = l3_layout$columns), file.path(l3_dir, "column_names.tsv"), sep = "\t")
  fwrite(data.table(
    branch = l3_layout$branch, separator = l3_layout$separator,
    gene_column = l3_layout$gene_column, transcript_column = l3_layout$transcript_column,
    n_columns = length(l3_layout$columns), fread_names = paste(names(l3_raw), collapse = ",")
  ), file.path(l3_dir, "layout.tsv"), sep = "\t")
  l3_agg <- aggregate_loaded_l3(l3_raw, file.path(root, "prereg", "frozen", "gencode_v38_basic_transcripts.tsv"), l3_layout)
  rm(l3_raw)
  gc()
  if (l3_agg$n_genes < 1000L) stop(paste("L3 gene map returned", l3_agg$n_genes, "genes from", l3_agg$n_mapped_rows, "rows"))
  l3_samples_all <- parse_l3_columns(colnames(l3_agg$counts))
  l3_primary <- l3_samples_all[batch == "1"]
  l3_batch2 <- l3_samples_all[batch == "2"]
  expect_n("L3 primary", nrow(l3_primary), uniqueN(l3_primary$subject_id), 116L, 58L)
  cov <- fread(file.path(root, "results", "stage5", "gse295758_subject_covariates.tsv"))
  if (!all(unique(l3_primary$subject_id) %in% cov$subject_id)) stop("L3 covariates do not cover the primary subjects")
  fwrite(l3_primary, file.path(l3_dir, "samples_scoring.tsv"), sep = "\t")
  l3_counts <- l3_agg$counts[, l3_primary$gsm, drop = FALSE]
  write_gene_matrix(l3_counts, file.path(l3_dir, "gene_counts.tsv"))
  fwrite(data.table(
    branch = l3_layout$branch, kind = l3_agg$kind, n_rows = l3_agg$n_rows,
    n_mapped_rows = l3_agg$n_mapped_rows, n_genes = l3_agg$n_genes,
    n_primary_samples = ncol(l3_counts), n_batch2_samples = nrow(l3_batch2)
  ), file.path(l3_dir, "preprocess_audit.tsv"), sep = "\t")
  l3_cpm <- log2_cpm(l3_counts)
  l3_fpkm <- counts_to_log2_fpkm(l3_counts)
  # count2FPKM returns a data.frame, which renames columns such as 1_D to X1_D.
  colnames(l3_fpkm) <- colnames(l3_counts)
  write_gene_matrix(l3_fpkm, file.path(l3_dir, "log2_fpkm.tsv"))
  l3_flags <- flag_3mad(l3_fpkm)
  fwrite(l3_flags, file.path(l3_dir, "flags_3mad.tsv"), sep = "\t")
  l3_gate <- gate_l3(l3_cpm, l3_primary)
  fwrite(l3_gate, file.path(l3_dir, "gate.tsv"), sep = "\t")
  say("L3 gate", l3_gate$dataset_pass[1])
  l3_pass <- isTRUE(l3_gate$dataset_pass[1])
  l3_scores <- NULL
  l3_components <- NULL
  l3_deltas <- NULL
  l3_comp <- NULL
  if (l3_pass) {
    scored <- score_tage(
      file.path(l3_dir, "gene_counts.tsv"),
      file.path(l3_dir, "samples_scoring.tsv"),
      "counts", file.path(l3_dir, "tage_scores.tsv"), file.path(l3_dir, "tage_components.tsv")
    )
    l3_ra_prep <- drop_all_nonfinite_rows(l3_fpkm[, l3_primary$gsm, drop = FALSE])
    fwrite(data.table(n_genes_dropped = l3_ra_prep$n_genes_dropped, n_genes_kept = l3_ra_prep$n_genes_kept), file.path(l3_dir, "d3_ra_audit.tsv"), sep = "\t")
    ra <- ra_linear_scores(l3_ra_prep$matrix)
    if (any(!is.finite(ra$value))) stop("L3 RA primary scores are not finite after D3")
    fwrite(ra, file.path(l3_dir, "ra_scores.tsv"), sep = "\t")
    l3_scores <- rbind(scored$scores, ra, fill = TRUE)
    l3_scores[, gsm := as.character(gsm)]
    l3_components <- scored$components
    l3_components[, gsm := as.character(gsm)]
    l3_comp <- composition_matrix(l3_fpkm, l3_primary)
  } else {
    fwrite(data.table(status = "gate fail"), file.path(l3_dir, "stopped.tsv"), sep = "\t")
  }

  day1 <- day1[, .(gsm = as.character(gsm), subject_id = as.character(subject_id), clock_time_h = as.numeric(clock_time_h))]
  day4 <- day4[, .(gsm = as.character(gsm), subject_id = as.character(subject_id), clock_time_h = as.numeric(clock_time_h))]
  period1 <- period1[, .(gsm = as.character(gsm), subject_id = as.character(subject_id), clock_time_h = as.numeric(clock_time_h))]
  period2 <- period2[, .(gsm = as.character(gsm), subject_id = as.character(subject_id), clock_time_h = as.numeric(clock_time_h))]
  l3_primary <- l3_primary[, .(gsm = as.character(gsm), subject_id = as.character(subject_id), phase = as.character(phase), clock_time_h = as.numeric(clock_time_h), batch = as.character(batch))]
  rm(l1_log, l2_log, l3_cpm, l3_fpkm, l3_counts)
  gc()

  blank_cosinor <- function(dataset, clocks, family) {
    rbindlist(lapply(clocks, function(clock) data.table(
      dataset = dataset, clock = clock, family = family, p = NA_real_,
      n_subjects = NA_integer_, n_samples = NA_integer_
    )))
  }
  p_parts <- list()
  s1_parts <- list()
  if (l1_pass) {
    p_parts[["L1"]] <- fit_clocks_cosinor("L1", day1, l1_scores, primary_clocks, "P", l1_draws, l1_dir)
    s1_parts[["L1"]] <- fit_clocks_cosinor("L1", day1, l1_scores, secondary_clocks, "S1", l1_draws, l1_dir)
  } else {
    p_parts[["L1"]] <- blank_cosinor("L1", primary_clocks, "P")
    s1_parts[["L1"]] <- blank_cosinor("L1", secondary_clocks, "S1")
  }
  if (l2_pass) {
    p_parts[["L2"]] <- fit_clocks_cosinor("L2", period1, l2_scores, primary_clocks, "P", l2_draws, l2_dir)
    s1_parts[["L2"]] <- fit_clocks_cosinor("L2", period1, l2_scores, secondary_clocks, "S1", l2_draws, l2_dir)
  } else {
    p_parts[["L2"]] <- blank_cosinor("L2", primary_clocks, "P")
    s1_parts[["L2"]] <- blank_cosinor("L2", secondary_clocks, "S1")
  }
  if (l3_pass) {
    l3_rows <- list()
    l3_delta_rows <- list()
    for (clock in c(primary_clocks, secondary_clocks)) {
      res <- paired_delta(l3_primary, l3_scores, clock)
      l3_delta_rows[[clock]] <- data.table(
        clock = clock, subject_id = res$subjects[[1]], morning = res$morning[[1]],
        evening = res$evening[[1]], delta = res$delta[[1]]
      )
      l3_rows[[clock]] <- data.table(
        dataset = "L3", clock = clock,
        family = if (clock %in% primary_clocks) "P" else "S1",
        n_subjects = res$n, n_samples = res$n * 2L,
        mean_delta = res$mean_delta, sd_delta = res$sd_delta,
        d = res$mean_delta / res$sd_delta,
        t_stat = res$t_stat, p = res$p, ci_low = res$ci_low, ci_high = res$ci_high
      )
    }
    l3_tab <- rbindlist(l3_rows)
    l3_deltas <- rbindlist(l3_delta_rows)
    fwrite(l3_deltas, file.path(l3_dir, "paired_deltas.tsv"), sep = "\t")
    p_parts[["L3"]] <- l3_tab[family == "P"]
    s1_parts[["L3"]] <- l3_tab[family == "S1"]
  } else {
    p_parts[["L3"]] <- blank_cosinor("L3", primary_clocks, "P")
    s1_parts[["L3"]] <- blank_cosinor("L3", secondary_clocks, "S1")
  }

  family_p <- apply_holm(rbindlist(p_parts, fill = TRUE))
  family_s1 <- apply_holm(rbindlist(s1_parts, fill = TRUE))
  fwrite(family_p, file.path(out_root, "family_P.tsv"), sep = "\t")
  fwrite(family_s1, file.path(out_root, "family_S1.tsv"), sep = "\t")

  rep_rows <- lapply(primary_clocks, function(clock_id) {
    sub <- family_p[clock == clock_id]
    n_sig <- sum(sub$decision == "significant")
    data.table(
      clock = clock_id,
      n_significant = n_sig,
      n_gate_fail = sum(sub$decision == "gate fail"),
      decision = if (n_sig >= 2L) "replicated" else "not replicated"
    )
  })
  replication <- rbindlist(rep_rows)
  fwrite(replication, file.path(out_root, "replication.tsv"), sep = "\t")

  sec_rows <- list()
  for (clock_id in primary_clocks) {
    phi_bar <- phi$phi_bar[match(clock_id, phi$clock)]
    pred <- f2$predicted_sign[match(clock_id, f2$clock)]
    if (l1_pass) {
      peak <- family_p[dataset == "L1" & clock == clock_id]$peak_h
      dist <- circ_dist_0_12(peak, phi_bar)
      sec_rows[[paste("L1", clock_id)]] <- data.table(
        dataset = "L1", clock = clock_id, phi_bar = phi_bar, peak_h = peak,
        circ_dist = dist, within_3h = is.finite(dist) && dist <= 3
      )
    } else {
      sec_rows[[paste("L1", clock_id)]] <- data.table(dataset = "L1", clock = clock_id, decision = "gate fail")
    }
    if (l3_pass) {
      delta <- family_p[dataset == "L3" & clock == clock_id]$mean_delta
      obs <- if (!is.finite(delta) || delta == 0) "0" else if (delta > 0) "+" else "-"
      sec_rows[[paste("L3", clock_id)]] <- data.table(
        dataset = "L3", clock = clock_id, predicted_sign = pred, observed_sign = obs,
        sign_agrees = identical(obs, pred)
      )
    } else {
      sec_rows[[paste("L3", clock_id)]] <- data.table(dataset = "L3", clock = clock_id, decision = "gate fail")
    }
  }
  fwrite(rbindlist(sec_rows, fill = TRUE), file.path(out_root, "h1_secondary.tsv"), sep = "\t")

  h2_rows <- list()
  if (l1_pass) h2_rows[["L1"]] <- h2_after_holm("L1", day1, l1_scores, l1_comp, primary_clocks, family_p, l1_draws, FALSE)
  if (l2_pass) h2_rows[["L2"]] <- h2_after_holm("L2", period1, l2_scores, l2_comp, primary_clocks, family_p, l2_draws, TRUE)
  if (l3_pass) {
    comp_delta <- paired_composition_delta(l3_primary, l3_comp)
    for (clock_id in primary_clocks) {
      idx <- which(family_p$dataset == "L3" & family_p$clock == clock_id)
      if (!isTRUE(family_p$decision[idx] == "significant")) next
      say("H2 L3", clock_id)
      block <- l3_deltas[clock == clock_id]
      delta <- block$delta
      ord <- match(comp_delta$subjects, block$subject_id)
      if (anyNA(ord)) stop("L3 composition subjects do not match the clock delta")
      delta <- delta[ord]
      fit <- h2_l3_bootstrap(delta, comp_delta$delta, 2000, 20261015)
      fit$dataset <- "L3"
      fit$clock <- clock_id
      fit$decision <- if (is.finite(fit$p_med_p025) && fit$p_med_p025 > 0) "supported" else "not supported"
      h2_rows[[paste("L3", clock_id)]] <- fit
    }
  }
  h2_tab <- rbindlist(h2_rows, fill = TRUE)
  if (nrow(h2_tab)) fwrite(h2_tab, file.path(out_root, "h2.tsv"), sep = "\t")

  h3_rows <- list()
  h5_rows <- list()
  intercept_rows <- list()
  if (l3_pass) {
    for (clock_id in primary_clocks) {
      block <- l3_deltas[clock == clock_id]
      ord <- match(block$subject_id, cov$subject_id)
      cov_i <- cov[ord]
      age_fit <- h3_age(block$delta, cov_i$age, cov_i$sex, cov_i$opioid_group)
      frame <- data.frame(delta = block$delta, age10 = cov_i$age / 10, sex = factor(cov_i$sex), opioid = factor(cov_i$opioid_group))
      lm_fit <- lm(delta ~ age10 + sex + opioid, data = frame)
      ci <- confint(lm_fit)["age10", ]
      age_fit$clock <- clock_id
      age_fit$ci_low <- unname(ci[1])
      age_fit$ci_high <- unname(ci[2])
      age_fit$direction <- if (age_fit$age_coef < 0) "-" else if (age_fit$age_coef > 0) "+" else "0"
      h3_rows[[clock_id]] <- age_fit
      h5_rows[[clock_id]] <- cbind(
        data.table(clock = clock_id, c = c_text[[clock_id]], c_numeric = c_values[[clock_id]]),
        h5_transfer(block$delta, c_values[[clock_id]])
      )
      intercept_rows[[clock_id]] <- cbind(
        data.table(clock = clock_id),
        sensitivity_intercept(block$delta, cov_i$sex, cov_i$opioid_group)
      )
    }
    h3 <- apply_holm(rbindlist(h3_rows)[, p := age_p])
    fwrite(h3, file.path(out_root, "h3.tsv"), sep = "\t")
    h5 <- rbindlist(h5_rows)
    h5[, l3_significant := family_p[dataset == "L3"]$decision[match(clock, family_p[dataset == "L3"]$clock)] == "significant"]
    h5[, decision := ifelse(!is.na(ci_low) & !is.na(ci_high) & ci_low <= 0 & ci_high >= 0 & l3_significant, "calibrated", "not calibrated")]
    fwrite(h5, file.path(out_root, "h5.tsv"), sep = "\t")
    fwrite(rbindlist(intercept_rows), file.path(out_root, "sensitivity_intercept.tsv"), sep = "\t")
    fwrite(data.table(
      sex_reference = levels(factor(cov$sex))[1],
      opioid_reference = levels(factor(cov$opioid_group))[1]
    ), file.path(out_root, "sensitivity_intercept_reference.tsv"), sep = "\t")
  } else {
    fwrite(blank_cosinor("L3", primary_clocks, "S2")[, decision := "gate fail"], file.path(out_root, "h3.tsv"), sep = "\t")
    fwrite(data.table(dataset = "L3", decision = "gate fail"), file.path(out_root, "h5.tsv"), sep = "\t")
  }

  h4_tests <- list()
  h4_cond <- list()
  if (l1_pass) {
    h4_day1 <- day1[, .(gsm, subject_id, clock_time_h)]
    h4_day1[, day := 0]
    h4_day4 <- day4[, .(gsm, subject_id, clock_time_h)]
    h4_day4[, day := 1]
    h4_samples <- rbind(h4_day1, h4_day4)
    h4 <- h4_for("L1", h4_samples, l1_scores, primary_clocks)
    h4_tests[["L1"]] <- h4$tests
    h4_cond[["L1"]] <- h4$conditions
  }
  if (l2_pass) {
    h4_p1 <- period1[, .(gsm, subject_id, clock_time_h)]
    h4_p1[, day := 0]
    h4_p2 <- period2[, .(gsm, subject_id, clock_time_h)]
    h4_p2[, day := 1]
    h4_samples <- rbind(h4_p1, h4_p2)
    h4 <- h4_for("L2", h4_samples, l2_scores, primary_clocks)
    h4_tests[["L2"]] <- h4$tests
    h4_cond[["L2"]] <- h4$conditions
  }
  if (!l1_pass) h4_tests[["L1"]] <- data.table(dataset = "L1", clock = primary_clocks, p = NA_real_)
  if (!l2_pass) h4_tests[["L2"]] <- data.table(dataset = "L2", clock = primary_clocks, p = NA_real_)
  h4_tab <- apply_holm(rbindlist(h4_tests, fill = TRUE))
  fwrite(h4_tab, file.path(out_root, "h4.tsv"), sep = "\t")
  if (length(h4_cond)) fwrite(rbindlist(h4_cond, fill = TRUE), file.path(out_root, "h4_conditions.tsv"), sep = "\t")

  s4_one <- function(dataset, samples, score_dt, paired) {
    clocks <- unique(score_dt$clock)
    f3_clocks <- intersect(clocks, f3$clock)
    other <- setdiff(clocks, f3$clock)
    rows <- list()
    if (paired) {
      for (clock_id in clocks) {
        res <- paired_delta(samples, score_dt, clock_id)
        info <- f3[clock == clock_id]
        rows[[clock_id]] <- data.table(
          dataset = dataset, clock = clock_id, in_f3 = clock_id %in% f3$clock,
          n_subjects = res$n, mean_delta = res$mean_delta, p = res$p,
          observed_sign = if (!is.finite(res$mean_delta) || res$mean_delta == 0) "0" else if (res$mean_delta > 0) "+" else "-",
          parent_clock = if (startsWith(clock_id, "CMP_TAms_chrono")) "TAms_chrono" else if (startsWith(clock_id, "CMP_TAms_mort")) "TAms_mort" else ""
        )
      }
    } else {
      for (clock_id in clocks) {
        dat <- prepare_score(samples, score_dt, clock_id)
        fit <- model_m(dat)
        rows[[clock_id]] <- data.table(
          dataset = dataset, clock = clock_id, in_f3 = clock_id %in% f3$clock,
          n_subjects = fit$n_subjects, n_samples = fit$n_samples,
          two_a = fit$two_a, peak_h = fit$peak_h, lrt_stat = fit$lrt_stat, p = fit$lrt_p,
          parent_clock = if (startsWith(clock_id, "CMP_TAms_chrono")) "TAms_chrono" else if (startsWith(clock_id, "CMP_TAms_mort")) "TAms_mort" else ""
        )
      }
    }
    tab <- rbindlist(rows, fill = TRUE)
    tab[, q_bh := NA_real_]
    tab[, confirmed := NA]
    if (!paired && dataset == "L1") {
      for (i in which(tab$in_f3)) {
        phi_bar <- f3$phi_bar[match(tab$clock[i], f3$clock)]
        dist <- circ_dist_0_12(tab$peak_h[i], phi_bar)
        set(tab, i = i, j = "confirmed", value = is.finite(tab$p[i]) && tab$p[i] < 0.05 && is.finite(dist) && dist <= 3)
        set(tab, i = i, j = "circ_dist", value = dist)
        set(tab, i = i, j = "phi_bar", value = phi_bar)
      }
    }
    if (paired) {
      for (i in which(tab$in_f3)) {
        info <- f3[f3$clock == tab$clock[i]]
        obs <- tab$observed_sign[i]
        eligible <- is.finite(info$d_exp_median) && abs(info$d_exp_median) >= 0.2
        set(tab, i = i, j = "predicted_sign", value = info$predicted_sign)
        set(tab, i = i, j = "d_exp_median", value = info$d_exp_median)
        set(tab, i = i, j = "sign_eligible", value = eligible)
        set(tab, i = i, j = "sign_agrees", value = if (eligible) identical(obs, info$predicted_sign) else NA)
      }
    }
    exploratory <- which(!tab$in_f3)
    parents <- unique(tab$parent_clock[exploratory])
    parents <- parents[nzchar(parents)]
    for (parent in parents) {
      idx <- which(!tab$in_f3 & tab$parent_clock == parent)
      if (length(idx)) set(tab, i = idx, j = "q_bh", value = p.adjust(tab$p[idx], method = "BH"))
    }
    mod_idx <- which(!tab$in_f3 & !nzchar(tab$parent_clock))
    if (length(mod_idx)) set(tab, i = mod_idx, j = "q_bh", value = p.adjust(tab$p[mod_idx], method = "BH"))
    tab
  }

  s4_parts <- list()
  if (l1_pass) {
    s4_parts[["L1"]] <- s4_one("L1", day1, rbind(l1_scores[grepl("^TM_", clock)], l1_components, fill = TRUE), FALSE)
    fwrite(drop_flag_cosinor(day1, l1_scores, l1_flags, primary_clocks), file.path(l1_dir, "sensitivity_drop_flags.tsv"), sep = "\t")
    fwrite(random_slopes(day1, l1_scores, primary_clocks), file.path(l1_dir, "sensitivity_random_slopes.tsv"), sep = "\t")
    fwrite(data.table(status = "complete"), file.path(l1_dir, "stopped.tsv"), sep = "\t")
  }
  if (l2_pass) {
    s4_parts[["L2"]] <- s4_one("L2", period1, rbind(l2_scores[grepl("^TM_", clock)], l2_components, fill = TRUE), FALSE)
    fwrite(drop_flag_cosinor(period1, l2_scores, l2_flags, primary_clocks), file.path(l2_dir, "sensitivity_drop_flags.tsv"), sep = "\t")
    fwrite(random_slopes(period1, l2_scores, primary_clocks), file.path(l2_dir, "sensitivity_random_slopes.tsv"), sep = "\t")
    fwrite(data.table(status = "complete"), file.path(l2_dir, "stopped.tsv"), sep = "\t")
  }
  if (l3_pass) {
    s4_parts[["L3"]] <- s4_one("L3", l3_primary, rbind(l3_scores[grepl("^TM_", clock)], l3_components, fill = TRUE), TRUE)
    keep_gsm <- l3_flags$gsm[!l3_flags$flag_3mad]
    flagged_subjects <- unique(l3_primary$subject_id[!l3_primary$gsm %in% keep_gsm])
    kept <- l3_primary[!subject_id %in% flagged_subjects]
    drop_rows <- list()
    for (clock in primary_clocks) {
      res <- paired_delta(kept, l3_scores, clock)
      drop_rows[[clock]] <- data.table(
        clock = clock, n_subjects = res$n, mean_delta = res$mean_delta,
        ci_low = res$ci_low, ci_high = res$ci_high, p = res$p
      )
    }
    fwrite(rbindlist(drop_rows), file.path(l3_dir, "sensitivity_drop_flags.tsv"), sep = "\t")
    fwrite(data.table(status = "complete"), file.path(l3_dir, "stopped.tsv"), sep = "\t")
  }
  if (length(s4_parts)) {
    s4 <- rbindlist(s4_parts, fill = TRUE)
    fwrite(s4, file.path(out_root, "s4.tsv"), sep = "\t")
    if (l1_pass) {
      n_confirmed <- s4[dataset == "L1" & in_f3 == TRUE, sum(confirmed)]
      fwrite(data.table(n_confirmed = n_confirmed, n_f3 = 8L), file.path(out_root, "s4_l1_confirmed.tsv"), sep = "\t")
    }
  }

  batch_names <- l3_layout$columns[grepl("batch", l3_layout$columns)]
  fwrite(data.table(
    n_batch2_columns = nrow(l3_batch2),
    batch_columns = paste(batch_names, collapse = ","),
    note = "batch-2 columns are absent from the supplement header"
  ), file.path(l3_dir, "sensitivity_batch2.tsv"), sep = "\t")

  say("L3 combat sensitivity")
  combat_status <- tryCatch({
    combat_raw <- load_locked_expression("L3_combat", root)
    combat_layout <- l3_layout_last
    combat_agg <- aggregate_loaded_l3(combat_raw, file.path(root, "prereg", "frozen", "gencode_v38_basic_transcripts.tsv"), combat_layout)
    rm(combat_raw)
    gc()
    combat_counts <- combat_agg$counts[, l3_primary$gsm, drop = FALSE]
    combat_path_counts <- file.path(l3_dir, "combat_gene_counts.tsv")
    write_gene_matrix(combat_counts, combat_path_counts)
    combat_fpkm <- counts_to_log2_fpkm(combat_counts)
    colnames(combat_fpkm) <- colnames(combat_counts)
    scored <- score_tage(
      combat_path_counts, file.path(l3_dir, "samples_scoring.tsv"),
      "counts", file.path(l3_dir, "combat_tage_scores.tsv"), file.path(l3_dir, "combat_tage_components.tsv")
    )
    ra <- ra_linear_scores(drop_all_nonfinite_rows(combat_fpkm[, l3_primary$gsm, drop = FALSE])$matrix)
    combat_scores <- rbind(scored$scores, ra, fill = TRUE)
    combat_scores[, gsm := as.character(gsm)]
    rows <- lapply(primary_clocks, function(clock) {
      res <- paired_delta(l3_primary, combat_scores, clock)
      data.table(
        clock = clock, n_subjects = res$n, mean_delta = res$mean_delta,
        sd_delta = res$sd_delta, ci_low = res$ci_low, ci_high = res$ci_high,
        t_stat = res$t_stat, p = res$p, status = "scored"
      )
    })
    fwrite(rbindlist(rows), file.path(l3_dir, "sensitivity_combat.tsv"), sep = "\t")
    "scored"
  }, error = function(e) {
    fwrite(data.table(status = "error", message = conditionMessage(e)), file.path(l3_dir, "sensitivity_combat.tsv"), sep = "\t")
    conditionMessage(e)
  })
  say("combat", combat_status)

  if (l1_pass) {
    say("D1 complete-case and D2 sensitivity")
    ra_cc <- fread(file.path(l1_dir, "ra_scores_complete_case.tsv"))
    cc_rows <- lapply(c(primary_clocks[3], secondary_clocks), function(clock_id) {
      dat <- prepare_score(day1, ra_cc, clock_id)
      fit <- model_m(dat)
      point_row("L1", clock_id, "D1_complete", fit)
    })
    fwrite(rbindlist(cc_rows), file.path(l1_dir, "sensitivity_d1_complete_case.tsv"), sep = "\t")
    d2 <- score_tage(
      file.path(l1_dir, "log2_genes_d1_imputed.tsv"),
      file.path(l1_dir, "samples_scoring.tsv"),
      "array", file.path(l1_dir, "tage_scores_d1_imputed.tsv"), file.path(l1_dir, "tage_components_d1_imputed.tsv")
    )
    d2_scores <- rbind(d2$scores, d2$components, fill = TRUE)
    d2_scores[, gsm := as.character(gsm)]
    d2_clocks <- c("TAms_chrono", "TAms_mort", f3$clock)
    d2_rows <- lapply(d2_clocks, function(clock_id) {
      dat <- prepare_score(day1, d2_scores, clock_id)
      fit <- model_m(dat)
      point_row("L1", clock_id, "D2_imputed", fit)
    })
    fwrite(rbindlist(d2_rows), file.path(l1_dir, "sensitivity_d2_imputed.tsv"), sep = "\t")
  }

  if (l3_pass) {
    say("D6 comprehensive sensitivity")
    comp_raw <- load_locked_expression("L3", root)
    comp_layout <- l3_layout_last
    comp_agg <- aggregate_loaded_l3(
      comp_raw,
      file.path(root, "results", "stage7b", "gencode_v38_comprehensive_transcripts.tsv"),
      comp_layout
    )
    rm(comp_raw)
    gc()
    comp_counts <- comp_agg$counts[, l3_primary$gsm, drop = FALSE]
    fwrite(data.table(
      n_rows = comp_agg$n_rows, n_mapped_rows = comp_agg$n_mapped_rows, n_genes = comp_agg$n_genes
    ), file.path(l3_dir, "d6_comprehensive_audit.tsv"), sep = "\t")
    comp_path <- file.path(l3_dir, "gene_counts_comprehensive.tsv")
    write_gene_matrix(comp_counts, comp_path)
    comp_scored <- score_tage(
      comp_path, file.path(l3_dir, "samples_scoring.tsv"), "counts",
      file.path(l3_dir, "tage_scores_comprehensive.tsv"), file.path(l3_dir, "tage_components_comprehensive.tsv")
    )
    comp_fpkm <- counts_to_log2_fpkm(comp_counts)
    colnames(comp_fpkm) <- colnames(comp_counts)
    comp_ra <- ra_linear_scores(drop_all_nonfinite_rows(comp_fpkm)$matrix)
    comp_scores <- rbind(comp_scored$scores, comp_ra, fill = TRUE)
    comp_scores[, gsm := as.character(gsm)]
    d6_rows <- lapply(primary_clocks, function(clock_id) {
      res <- paired_delta(l3_primary, comp_scores, clock_id)
      data.table(
        clock = clock_id, n_subjects = res$n, n_samples = res$n * 2L,
        mean_delta = res$mean_delta, sd_delta = res$sd_delta, d = res$mean_delta / res$sd_delta,
        t_stat = res$t_stat, p = res$p, ci_low = res$ci_low, ci_high = res$ci_high
      )
    })
    fwrite(rbindlist(d6_rows), file.path(l3_dir, "sensitivity_d6_comprehensive.tsv"), sep = "\t")
    unlink(comp_path)
  }
  say("confirmation finished")
  invisible(TRUE)
}
