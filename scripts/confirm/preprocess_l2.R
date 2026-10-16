# L2 preprocessing. Does not open locked series matrices.

refuse_locked_path <- function(path) {
  norm <- tryCatch(normalizePath(path, mustWork = FALSE), error = function(e) as.character(path))
  if (grepl("/data/locked(/|$)", norm)) stop("refusing to open a locked expression path")
}

geo2r_scale_matrix <- function(values) {
  flat <- as.numeric(values)
  flat <- flat[is.finite(flat)]
  qs <- stats::quantile(flat, c(0, 0.25, 0.99, 1), names = FALSE, type = 7)
  if (qs[3] > 100 || ((qs[4] - qs[1]) > 50 && qs[2] > 0)) return("linear")
  "log"
}

to_log2_matrix <- function(values, scale) {
  out <- values
  if (identical(scale, "linear")) {
    out[out <= 0] <- NA_real_
    out <- log2(out)
  }
  out
}

collapse_clariom <- function(log_mat) {
  suppressPackageStartupMessages(library(clariomshumanhttranscriptcluster.db))
  ids <- rownames(log_mat)
  ann <- AnnotationDbi::select(
    clariomshumanhttranscriptcluster.db,
    keys = ids,
    columns = "SYMBOL",
    keytype = "PROBEID"
  )
  ann <- ann[!is.na(ann$SYMBOL) & nzchar(ann$SYMBOL), , drop = FALSE]
  n_symbol <- tapply(ann$SYMBOL, ann$PROBEID, function(x) length(unique(x)))
  n_multi <- sum(n_symbol > 1, na.rm = TRUE)
  means <- rowMeans(log_mat, na.rm = TRUE)
  ann$mean_log2 <- unname(means[ann$PROBEID])
  ann <- ann[order(ann$SYMBOL, -ann$mean_log2, ann$PROBEID), , drop = FALSE]
  ann <- ann[!duplicated(ann$SYMBOL), , drop = FALSE]
  out <- log_mat[ann$PROBEID, , drop = FALSE]
  rownames(out) <- ann$SYMBOL
  list(matrix = out, n_probes = nrow(log_mat), n_multi = n_multi, n_genes = nrow(out))
}
