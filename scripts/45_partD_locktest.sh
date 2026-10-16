#!/usr/bin/env bash
# Call load_locked_expression with no unlock file and record openat traces.
set -euo pipefail
root="/HDD8T1/WORK/GeroSci/TOD"
out="$root/results/stage6/rehearsal/lock"
mkdir -p "$out"
cd "$root"
rbin="/home/jkang/anaconda3/envs/r_env/bin/Rscript"
for which_data in L1 L2; do
  trace="$out/${which_data}.strace"
  set +e
  strace -f -e openat -o "$trace" "$rbin" -e "source('scripts/confirm/lib_lock.R'); load_locked_expression('${which_data}')" > "$out/${which_data}.stdout" 2> "$out/${which_data}.stderr"
  status=$?
  set -e
  echo "$status" > "$out/${which_data}.exit"
done
python3 - << 'PY'
from pathlib import Path
root = Path("/HDD8T1/WORK/GeroSci/TOD/results/stage6/rehearsal/lock")
targets = {
    "L1": "GSE48113_series_matrix.txt.gz",
    "L2": "GSE122725_series_matrix.txt.gz",
}
lines = ["dataset\texit_code\topenat_hits\tstderr_head"]
for which, name in targets.items():
    trace = (root / f"{which}.strace").read_text(errors="replace")
    hits = [line for line in trace.splitlines() if name in line and "openat" in line]
    err = (root / f"{which}.stderr").read_text(errors="replace").strip().splitlines()
    head = err[-1] if err else ""
    lines.append(f"{which}\t{(root / f'{which}.exit').read_text().strip()}\t{len(hits)}\t{head}")
(root / "lock_test.tsv").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
PY
