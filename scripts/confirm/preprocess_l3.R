# L3 count preprocessing for files the caller names.
# Locked paths are refused. load_locked_expression is the only opener of locked files.

source_lock <- function() {
  if (!exists("l3_layout_from_header", mode = "function")) {
    args_all <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args_all, value = TRUE)
    start <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg[1])) else getwd()
    cur <- normalizePath(start)
    root <- NULL
    for (step in seq_len(8)) {
      if (file.exists(file.path(cur, "scripts", "confirm", "lib_lock.R"))) {
        root <- cur
        break
      }
      parent <- dirname(cur)
      if (parent == cur) break
      cur <- parent
    }
    if (is.null(root)) stop("lib_lock.R not found")
    sys.source(file.path(root, "scripts", "confirm", "lib_lock.R"), envir = parent.frame())
  }
}

base_id <- function(x) sub("\\..*$", "", x)

sum_by_gene <- function(counts, genes) {
  dt <- data.table::as.data.table(counts)
  dt[, gene := genes]
  out <- dt[, lapply(.SD, sum, na.rm = TRUE), by = gene]
  mat <- as.matrix(out[, -1])
  rownames(mat) <- out$gene
  colnames(mat) <- colnames(counts)
  mat
}

read_l3_counts <- function(path, gencode_path) {
  source_lock()
  norm <- tryCatch(normalizePath(path, mustWork = FALSE), error = function(e) as.character(path))
  if (grepl("/data/locked(/|$)", norm)) stop("refusing to open a locked expression path")
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  header <- readLines(con, n = 1, warn = FALSE)
  close(con)
  layout <- l3_layout_from_header(header)
  raw <- data.table::fread(path)
  sample_cols <- setdiff(names(raw), c(layout$gene_column, layout$transcript_column))
  sample_cols <- sample_cols[!is.na(sample_cols)]
  gencode <- data.table::fread(gencode_path)
  if (identical(layout$branch, "gene_column")) {
    ids <- as.character(raw[[layout$gene_column]])
    kind <- if (all(grepl("^ENSG[0-9]", base_id(ids)))) "ensembl" else "symbol"
    if (identical(kind, "ensembl")) {
      key <- base_id(ids)
      map <- unique(gencode[, .(gene_id_base, gene_name)])
      names_out <- map$gene_name[match(key, map$gene_id_base)]
    } else {
      names_out <- ids
    }
    keep <- !is.na(names_out) & nzchar(names_out)
    mat <- as.matrix(raw[keep, sample_cols, with = FALSE])
    storage.mode(mat) <- "double"
    counts <- sum_by_gene(mat, names_out[keep])
  } else {
    id_col <- layout$transcript_column
    if (is.na(id_col)) id_col <- names(raw)[1]
    ids <- as.character(raw[[id_col]])
    map <- unique(gencode[, .(transcript_id_base, gene_name)])
    names_out <- map$gene_name[match(base_id(ids), map$transcript_id_base)]
    keep <- !is.na(names_out) & nzchar(names_out)
    sample_cols <- setdiff(names(raw), id_col)
    mat <- as.matrix(raw[keep, sample_cols, with = FALSE])
    storage.mode(mat) <- "double"
    counts <- sum_by_gene(mat, names_out[keep])
  }
  list(counts = counts, layout = layout, n_genes = nrow(counts))
}
