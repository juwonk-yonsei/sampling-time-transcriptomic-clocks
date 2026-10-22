#!/usr/bin/env bash
# Download the public input files into data/ and compare them with config/raw_checksums.tsv.
# Files that already exist are not downloaded again.
# GEO can regenerate series matrices and platform tables, so a later download may differ
# in bytes from the copy hashed on 2026-10-16; the comparison reports this as "mismatch".
# The data/locked/ files were opened only after the unlock commit (see README).
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

geo="https://ftp.ncbi.nlm.nih.gov/geo"
acc="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

files=(
  "data/raw/GSE39445_series_matrix.txt.gz|$geo/series/GSE39nnn/GSE39445/matrix/GSE39445_series_matrix.txt.gz"
  "data/raw/GSE56931_series_matrix.txt.gz|$geo/series/GSE56nnn/GSE56931/matrix/GSE56931_series_matrix.txt.gz"
  "data/raw/GSE113883/GSE113883_series_matrix.txt.gz|$geo/series/GSE113nnn/GSE113883/matrix/GSE113883_series_matrix.txt.gz"
  "data/raw/GSE113883/GSE113883_TPM_processed_data.txt.gz|$geo/series/GSE113nnn/GSE113883/suppl/GSE113883_TPM_processed_data.txt.gz"
  "data/raw/GSE196728/GSE196728_series_matrix.txt.gz|$geo/series/GSE196nnn/GSE196728/matrix/GSE196728_series_matrix.txt.gz"
  "data/raw/GSE196728/GSE196728_Peru_raw.xlsx|$geo/series/GSE196nnn/GSE196728/suppl/GSE196728_Peru_raw.xlsx"
  "data/locked/GSE48113_series_matrix.txt.gz|$geo/series/GSE48nnn/GSE48113/matrix/GSE48113_series_matrix.txt.gz"
  "data/locked/GSE122725_series_matrix.txt.gz|$geo/series/GSE122nnn/GSE122725/matrix/GSE122725_series_matrix.txt.gz"
  "data/locked/GSE295758/GSE295758_series_matrix.txt.gz|$geo/series/GSE295nnn/GSE295758/matrix/GSE295758_series_matrix.txt.gz"
  "data/locked/GSE295758/GSE295758_rawTranscriptCounts.csv.gz|$geo/series/GSE295nnn/GSE295758/suppl/GSE295758_rawTranscriptCounts.csv.gz"
  "data/locked/GSE295758/GSE295758_combatVstTranscriptCounts.csv.gz|$geo/series/GSE295nnn/GSE295758/suppl/GSE295758_combatVstTranscriptCounts.csv.gz"
  "data/raw/gpl/GPL10379.annot.gz|$geo/platforms/GPL10nnn/GPL10379/annot/GPL10379.annot.gz"
  "data/raw/gpl/GPL15331_platform.txt|$acc?acc=GPL15331&targ=self&view=full&form=text"
  "data/raw/gpl/GPL24324_platform.txt|$acc?acc=GPL24324&targ=self&view=full&form=text"
)

# PMC serves a browser check page to scripted requests, so this file is fetched by hand.
manual=(
  "data/raw/peters2015/ncomms9570-s2.xlsx|https://pmc.ncbi.nlm.nih.gov/articles/instance/4639797/bin/ncomms9570-s2.xlsx"
)

for entry in "${files[@]}"; do
  path="${entry%%|*}"
  url="${entry#*|}"
  if [[ -s "$path" ]]; then
    continue
  fi
  mkdir -p "$(dirname "$path")"
  echo "download $path"
  curl -L --fail --retry 3 -o "$path.part" "$url"
  mv "$path.part" "$path"
done

for entry in "${manual[@]}"; do
  path="${entry%%|*}"
  url="${entry#*|}"
  if [[ ! -s "$path" ]]; then
    echo "manual  $path  <- $url (Peters et al. 2015 Supplementary Data 1; open in a browser and save here)"
  fi
done

echo
printf "%s\t%s\t%s\n" status bytes path
tail -n +2 config/raw_checksums.tsv | while IFS=$'\t' read -r path bytes sha; do
  if [[ ! -f "$path" ]]; then
    printf "%s\t%s\t%s\n" missing "$bytes" "$path"
    continue
  fi
  got="$(sha256sum "$path" | cut -d' ' -f1)"
  if [[ "$got" == "$sha" ]]; then
    printf "%s\t%s\t%s\n" match "$bytes" "$path"
  else
    printf "%s\t%s\t%s\n" mismatch "$(stat -c %s "$path")" "$path"
  fi
done
