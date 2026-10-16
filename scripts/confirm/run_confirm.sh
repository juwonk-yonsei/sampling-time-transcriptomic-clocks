#!/usr/bin/env bash
# Confirmation entry point. Instruction 7 runs this after unlock.
# This script does not download L3 supplements.
set -euo pipefail
root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$root"
Rscript "$root/scripts/confirm/run_pipeline.R" --mode locked
