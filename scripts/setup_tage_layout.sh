#!/usr/bin/env bash
# Link the tAge resources into data/external/ in the layout read by third_party/clockspec
# and third_party/scmodage. Nothing is copied.
#   $1  folder with the tAge Zenodo v2 files (doi:10.5281/zenodo.22166800), module archives unpacked
#   $2  clone of https://github.com/Gladyshev-Lab/tAge at d048ba6016a7e5a114f45195e0b11071e2b8c424
# Afterwards export CLOCKSPEC_TAGE="$root/data/external" and PYTHONPATH="$root/third_party".
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
zenodo="$(cd "$1" && pwd)"
tage_repo="$(cd "$2" && pwd)"
ext="$root/data/external"

for need in "$zenodo/EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl" \
            "$zenodo/Rodent module clocks" "$zenodo/Multispecies module clocks" \
            "$tage_repo/inst/extdata/module_clocks_metadata.csv"; do
  [[ -e "$need" ]] || { echo "missing: $need" >&2; exit 1; }
done
commit="$(git -C "$tage_repo" rev-parse HEAD 2>/dev/null || echo unknown)"
[[ "$commit" == d048ba6016a7e5a114f45195e0b11071e2b8c424 ]] || echo "warning: tAge clone is at $commit" >&2

mkdir -p "$ext/tage_modules/rodent" "$ext/tage_modules/multispecies"
ln -sfn "$tage_repo" "$ext/tAge"
ln -sfn "$zenodo" "$ext/tage_small"
ln -sfn "$zenodo" "$ext/tage_zenodo"
ln -sfn "$zenodo/Rodent module clocks" "$ext/tage_modules/rodent/Rodent module clocks"
ln -sfn "$zenodo/Multispecies module clocks" "$ext/tage_modules/multispecies/Multispecies module clocks"
echo "linked tAge resources under $ext"
