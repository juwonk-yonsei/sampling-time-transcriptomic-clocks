#!/usr/bin/env bash
# Re-score the committed GSE48113 (array) and GSE295758 (counts) gene matrices with
# scripts/confirm/score_clocks.py and the packages in third_party/, and compare the outputs
# byte for byte with results/confirm/L1 and L3. Run scripts/setup_tage_layout.sh first.
#   PYTHON  interpreter of the clock environment (environment/python_clock_env_requirements.txt)
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
python="${PYTHON:-python}"
export CLOCKSPEC_TAGE="$root/data/external"
export PYTHONPATH="$root/third_party"
out="$(mktemp -d)"
trap 'rm -rf "$out"' EXIT

"$python" -c "import clockspec.clocks, scmodage.clocks, sys; [sys.exit('not using third_party: ' + m.__file__) for m in (clockspec.clocks, scmodage.clocks) if not m.__file__.startswith('$root/third_party/')]"

status=0
for spec in "L1|log2_genes.tsv.gz|array" "L3|gene_counts.tsv|counts"; do
  IFS='|' read -r set expr mode <<< "$spec"
  dir="$root/results/confirm/$set"
  "$python" "$root/scripts/confirm/score_clocks.py" --expr "$dir/$expr" --samples "$dir/samples_scoring.tsv" \
    --mode "$mode" --out "$out/scores.tsv" --components "$out/components.tsv" > /dev/null
  for kind in scores components; do
    if cmp -s "$out/$kind.tsv" "$dir/tage_$kind.tsv"; then
      echo "$set tage_$kind.tsv: identical"
    else
      echo "$set tage_$kind.tsv: DIFFERENT"; status=1
    fi
  done
done
exit $status
