# Download functions for the L3 supplements. Not called by run_confirm.sh.
# Instruction 7 runs them after unlock.

l3_supplement_urls <- function() {
  c(
    GSE295758_rawTranscriptCounts.csv.gz = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE295nnn/GSE295758/suppl/GSE295758_rawTranscriptCounts.csv.gz",
    GSE295758_combatVstTranscriptCounts.csv.gz = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE295nnn/GSE295758/suppl/GSE295758_combatVstTranscriptCounts.csv.gz"
  )
}

download_l3_supplements <- function(dest_dir) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  urls <- l3_supplement_urls()
  for (name in names(urls)) {
    dest <- file.path(dest_dir, name)
    utils::download.file(urls[[name]], dest, mode = "wb", quiet = FALSE)
  }
  invisible(file.path(dest_dir, names(urls)))
}
