#!/usr/bin/env Rscript
# Part D. Cosinor on the rescored clocks. Bootstrap is limited to TAms and TArd.

suppressPackageStartupMessages({
  library(lme4)
  library(data.table)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
root <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(root, "scripts", "lib_cosinor.R"))
out_dir <- file.path(root, "results", "stage5")
dir.create(file.path(out_dir, "boot"), showWarnings = FALSE, recursive = TRUE)

forms <- main_formula()
n_boot <- 2000
datasets <- c("GSE39445", "GSE56931", "GSE113883")
boot_clocks <- c("TAms_chrono", "TAms_mort", "TArd_chrono", "TArd_mort")
proj_clocks <- c("TAms_chrono", "TAms_mort")
blood_labels <- c(
  "T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells",
  "B lineage", "Monocytic lineage", "Myeloid dendritic cells", "Neutrophils"
)
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

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

read_samples <- function(dataset) {
  path <- if (dataset == "GSE113883") {
    file.path(root, "results", "stage3", "GSE113883_samples.tsv")
  } else {
    file.path(root, "data", "processed", paste0(dataset, "_samples.tsv"))
  }
  samples <- fread(path)
  keep <- which(samples$condition == "control" & is.finite(samples$clock_time_h))
  samples[keep, .(gsm = as.character(gsm), subject_id = as.character(subject_id), clock_time_h)]
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

variance_terms <- function(data) {
  reml <- quiet_lmer(forms$full, data, TRUE)
  null_reml <- quiet_lmer(forms$null, data, TRUE)
  sigma_reml <- if (!inherits(reml, "error") && optimizer_converged(reml)) sigma(reml) else NA_real_
  sigma_u <- NA_real_
  if (!inherits(null_reml, "error") && optimizer_converged(null_reml)) {
    vc <- lme4::VarCorr(null_reml)
    sds <- if (length(vc) == 1L) attr(vc[[1]], "stddev") else numeric()
    if (length(sds) == 1L) sigma_u <- as.numeric(sds)
  }
  c(sigma_reml = sigma_reml, sigma_u = sigma_u)
}

attach_score <- function(samples, score_dt, clock_id) {
  idx <- which(score_dt$clock_id == clock_id)
  sc <- score_dt[idx]
  named <- sc$value
  names(named) <- sc$gsm
  dat <- as.data.frame(samples)
  dat$score <- unname(named[dat$gsm])
  if (anyNA(dat$score)) stop(paste("missing score", clock_id))
  add_time_terms(dat)
}

scores <- fread(file.path(out_dir, "clock_scores.tsv"))
components <- fread(file.path(out_dir, "component_scores.tsv"))
scores[, gsm := as.character(gsm)]
components[, gsm := as.character(gsm)]
clock_meta <- unique(scores[, .(
  clock_id = clock, family = fifelse(catalog == "human_catalog" & module != "composite", "TM",
    fifelse(catalog == "rodent_module_rescore", "TMrd", "TA")),
  outcome, module, role, catalog, parent_clock = ""
)])
comp_meta <- unique(components[, .(
  clock_id = paste0("CMP_", clock, "_", component),
  family = "CMP",
  outcome = fifelse(grepl("mort", clock), "mortality", "chronological"),
  module = component,
  role = "primary",
  catalog = "component",
  parent_clock = clock
)])
analysis <- rbindlist(list(
  scores[, .(dataset, gsm, clock_id = clock, value)],
  components[, .(dataset, gsm, clock_id = paste0("CMP_", clock, "_", component), value)]
), use.names = TRUE)
meta <- rbindlist(list(clock_meta, comp_meta), use.names = TRUE)

point_rows <- list()
boot_rows <- list()
for (dataset in datasets) {
  cat(dataset, "\n")
  flush.console()
  samples <- read_samples(dataset)
  draws <- fread(draw_path(dataset))
  subjects <- sort(unique(samples$subject_id))
  draw_subjects <- sort(unique(as.character(draws$subject_id)))
  if (!identical(subjects, draw_subjects)) stop(paste("draw subjects differ", dataset))
  present <- which(analysis$dataset == dataset)
  one <- analysis[present]
  clock_ids <- meta$clock_id
  for (clock_id in clock_ids) {
    dat <- attach_score(samples, one, clock_id)
    fit <- fit_cosinor(dat, forms$full, forms$null)
    extra <- variance_terms(dat)
    fit$sigma_reml <- unname(extra["sigma_reml"])
    fit$sigma_u <- unname(extra["sigma_u"])
    fit$two_a_over_sigma <- if (is.finite(fit$two_a) && is.finite(fit$sigma_reml) && fit$sigma_reml > 0) {
      fit$two_a / fit$sigma_reml
    } else {
      NA_real_
    }
    info_idx <- which(meta$clock_id == clock_id)
    info <- meta[info_idx]
    fit$dataset <- dataset
    fit$clock <- clock_id
    fit$family <- info$family
    fit$outcome <- info$outcome
    fit$module <- info$module
    fit$role <- info$role
    fit$parent_clock <- info$parent_clock
    point_rows[[paste(dataset, clock_id)]] <- as.data.table(fit)
    if (clock_id %in% boot_clocks) {
      rep_path <- file.path(out_dir, "boot", paste0("h1_boot_reps_", dataset, "_", clock_id, ".tsv"))
      if (file.exists(rep_path) && nrow(fread(rep_path)) == n_boot) {
        rep <- fread(rep_path)
      } else {
        cat(" bootstrap", clock_id, "\n")
        flush.console()
        reps <- vector("list", n_boot)
        for (b in seq_len(n_boot)) {
          draw_idx <- which(draws$draw == b)
          drawn <- draws[draw_idx]
          boot_dat <- apply_draw(dat, drawn$subject_id, drawn$new_subject_id)
          reps[[b]] <- fit_cosinor(boot_dat, forms$full, forms$null)
        }
        rep <- rbindlist(lapply(reps, as.data.table))
        rep[, draw := seq_len(n_boot)]
        fwrite(rep[, .(draw, two_a, r2_within, two_a_over_sd_between, two_a_over_sd_resid, peak_h,
                       reml_converged, ml_converged, null_converged)], rep_path, sep = "\t")
      }
      ok <- rep$reml_converged & rep$ml_converged & rep$null_converged
      summ <- data.table(
        dataset = dataset, clock = clock_id, n_boot = nrow(rep),
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
      boot_rows[[paste(dataset, clock_id)]] <- summ
    }
  }
}

points <- rbindlist(point_rows, fill = TRUE)
tm_idx <- which(points$family == "TM")
tm <- points[tm_idx]
tm[, q_bh := p.adjust(lrt_p, method = "BH"), by = .(dataset, outcome)]
cmp_idx <- which(points$family == "CMP")
cmp <- points[cmp_idx]
cmp[, q_bh := p.adjust(lrt_p, method = "BH"), by = .(dataset, parent_clock)]
qtab <- rbindlist(list(
  tm[, .(dataset, clock, q_bh)],
  cmp[, .(dataset, clock, q_bh)]
), use.names = TRUE)
points <- merge(points, qtab, by = c("dataset", "clock"), all.x = TRUE, sort = FALSE)
fwrite(points, file.path(out_dir, "h1_point_estimates.tsv"), sep = "\t")
fwrite(rbindlist(boot_rows), file.path(out_dir, "h1_bootstrap_ci.tsv"), sep = "\t")

stage4 <- fread(file.path(root, "results", "stage4", "module_clocks", "h1_point_estimates.tsv"))
stage4_map <- data.table(
  clock = points[family %in% c("TA", "TMrd"), unique(clock)],
  stage4_clock = NA_character_
)
for (i in seq_len(nrow(stage4_map))) {
  clock_id <- stage4_map$clock[i]
  stage4_map$stage4_clock[i] <- if (clock_id == "TArd_chrono") {
    "TA_chrono"
  } else if (clock_id == "TArd_mort") {
    "TA_mort"
  } else {
    sub("^TMrd_", "TM_", clock_id)
  }
}
new_keep <- points[family %in% c("TA", "TMrd"), .(
  dataset, clock, stage4_r2 = r2_within, stage4_lrt_p = lrt_p, stage4_peak = peak_h
)]
# column names above are the new results; rename before the join
setnames(new_keep, c("stage4_r2", "stage4_lrt_p", "stage4_peak"), c("r2_within", "lrt_p", "peak_h"))
old_keep <- stage4[, .(dataset, stage4_clock = clock, stage4_r2 = r2_within, stage4_lrt_p = lrt_p, stage4_peak = peak_h)]
comparison <- merge(merge(new_keep, stage4_map, by = "clock"), old_keep, by = c("dataset", "stage4_clock"))
fwrite(comparison, file.path(out_dir, "h1_vs_stage4.tsv"), sep = "\t")

cat("projection\n")
flush.console()
proj_points <- list()
proj_ci <- list()
for (dataset in datasets) {
  samples <- as.data.frame(read_samples(dataset))
  mcp <- fread(mcp_path(dataset))
  mcp[, gsm := as.character(gsm)]
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
  present <- which(analysis$dataset == dataset)
  one <- analysis[present]
  for (clock_id in proj_clocks) {
    dat <- base
    named_idx <- which(one$clock_id == clock_id)
    sc <- one[named_idx]
    named <- sc$value
    names(named) <- sc$gsm
    dat$score <- unname(named[dat$gsm])
    if (anyNA(dat$score)) stop(paste("missing projection score", dataset, clock_id))
    unadj <- fit_cosinor(dat, forms$full, forms$null)
    adj_beta <- reml_beta(dat, adj_full)
    ru <- c(unadj$beta_cos, unadj$beta_sin)
    ra <- c(unname(adj_beta["cos"]), unname(adj_beta["sin"]))
    proj <- projection(ru, ra)
    proj_points[[paste(dataset, clock_id)]] <- data.table(
      dataset = dataset, clock = clock_id,
      n_subjects = unadj$n_subjects, n_samples = unadj$n_samples,
      lrt_p_unadj = unadj$lrt_p,
      p_med = unname(proj["p_med"]), two_a_orth = unname(proj["two_a_orth"]),
      ru_norm = unname(proj["ru_norm"]),
      reml_converged_unadj = unadj$reml_converged,
      reml_converged_adj = isTRUE(adj_beta["converged"] == 1),
      two_a_unadj = unadj$two_a
    )
    rep_path <- file.path(out_dir, "boot", paste0("projection_reps_", dataset, "_", clock_id, ".tsv"))
    if (file.exists(rep_path) && nrow(fread(rep_path)) == n_boot) {
      rep <- fread(rep_path)
    } else {
      cat(" projection", dataset, clock_id, "\n")
      flush.console()
      reps <- vector("list", n_boot)
      for (b in seq_len(n_boot)) {
        draw_idx <- which(draws$draw == b)
        drawn <- draws[draw_idx]
        boot_dat <- apply_draw(dat, drawn$subject_id, drawn$new_subject_id)
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
    proj_ci[[paste(dataset, clock_id)]] <- data.table(
      dataset = dataset, clock = clock_id, n_boot = n_boot,
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
