#!/usr/bin/env bash
# Step 06 driver: timing of the IS4 insertion in mutS against the LK PASS mutations of step 08.
set -euo pipefail
source config.sh
# config.sh must be the version shipped with this package (it defines the step 05-08 settings)
for v in DNA_SAMPLES VARIANT_ORGS KNOWN_VARIANTS RNA_EXCL_VARIANTS V_MINDEPTH MUTS_CONTIG LK_NEW; do
  if ! declare -p "$v" >/dev/null 2>&1; then
    echo "ERROR: config.sh does not define $v - this is an older config.sh."
    echo "       Replace it with the config.sh from LK_pipeline_scripts.zip (it also replaces align_config.sh),"
    echo "       then re-apply any paths you had changed. Nothing has been run."
    exit 1
  fi
done
[[ -s "$OUTDIR/08_variants.tsv" ]] || { echo "ERROR: run 08_run_variants.sh first"; exit 1; }
OUTDIR="$OUTDIR" NEW_CONTIG="$MUTS_CONTIG" REGION_START="$MUTS_START" REGION_END="$MUTS_END" VARIANTS="$OUTDIR/08_variants.tsv" \
  python3 06_muts_IS_timing.py refs/index/NEW.fna refs/index/OLD.fna "${DNA_SAMPLES[@]}" 2>&1 | tee "$OUTDIR/logs_06.txt"
