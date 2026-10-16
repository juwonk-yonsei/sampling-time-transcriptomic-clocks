# The only function that opens locked expression.

confirm_root <- function() {
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

locked_expression_path <- function(which_data, root = confirm_root()) {
  if (which_data == "L1") return(file.path(root, "data", "locked", "GSE48113_series_matrix.txt.gz"))
  if (which_data == "L2") return(file.path(root, "data", "locked", "GSE122725_series_matrix.txt.gz"))
  if (which_data == "L3") return(file.path(root, "data", "locked", "GSE295758", "GSE295758_rawTranscriptCounts.csv.gz"))
  if (which_data == "L3_combat") return(file.path(root, "data", "locked", "GSE295758", "GSE295758_combatVstTranscriptCounts.csv.gz"))
  stop("unknown locked dataset")
}

freeze_hash <- function(root = confirm_root()) {
  manifest <- file.path(root, "prereg", "manifest.tsv")
  if (!file.exists(manifest)) stop("manifest is missing")
  out <- system2("sha256sum", manifest, stdout = TRUE)
  strsplit(out, " ", fixed = TRUE)[[1]][1]
}

load_locked_expression <- function(which_data, root = confirm_root()) {
  unlock <- file.path(root, "config", "unlock_approved.txt")
  if (!file.exists(unlock)) stop("config/unlock_approved.txt is missing; locked expression stays closed")
  approved <- trimws(readLines(unlock, n = 1, warn = FALSE))
  freeze <- freeze_hash(root)
  env_hash <- Sys.getenv("TOD_UNLOCK", unset = "")
  if (!identical(approved, freeze) || !identical(env_hash, freeze)) {
    stop("unlock hash does not match the manifest SHA-256; locked expression stays closed")
  }
  path <- locked_expression_path(which_data, root)
  if (!file.exists(path)) stop("locked expression file is not on disk")
  if (which_data %in% c("L3", "L3_combat")) {
    con <- gzfile(path, "rt")
    header <- readLines(con, n = 1, warn = FALSE)
    close(con)
    layout <- l3_layout_from_header(header)
    table <- data.table::fread(path)
    attr(table, "l3_branch") <- layout$branch
    return(table)
  }
  read_geo_series_matrix(path)
}

read_geo_series_matrix <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (!length(line)) stop("series_matrix_table_begin was not found")
    if (identical(line, "!series_matrix_table_begin")) break
  }
  data.table::fread(con, sep = "\t")
}

l3_layout_from_header <- function(header_line) {
  cols <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
  cols <- gsub('^"|"$', "", cols)
  gene_names <- c("gene_id", "gene_name", "gene_symbol", "GENE", "Gene", "gene", "symbol", "SYMBOL")
  transcript_names <- c("transcript_id", "transcript", "Transcript")
  gene_col <- gene_names[gene_names %in% cols]
  transcript_col <- transcript_names[transcript_names %in% cols]
  list(
    columns = cols,
    gene_column = if (length(gene_col)) gene_col[1] else NA_character_,
    transcript_column = if (length(transcript_col)) transcript_col[1] else NA_character_,
    branch = if (length(gene_col)) "gene_column" else "transcript_to_gene"
  )
}
