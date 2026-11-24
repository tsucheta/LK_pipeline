#!/usr/bin/env bash
# Runs the whole pipeline in order. Every step skips work already done (FORCE=1 bash run_all.sh redoes all).
# Run one step on its own with: bash <step script>
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
mkdir -p "$OUTDIR"
bash 00_align.sh        2>&1 | tee "$OUTDIR/logs_00.txt"
bash 01_library_qc.sh   2>&1 | tee "$OUTDIR/logs_01.txt"
bash 02_count.sh        2>&1 | tee "$OUTDIR/logs_02.txt"
bash 03_run_deseq2.sh   2>&1 | tee "$OUTDIR/logs_03.txt"
bash 08_run_variants.sh
bash 09_run_split_check.sh
# step 10 (polishing LK_1) changes the reference and is run separately: bash 10_run_polish.sh
bash 06_run_muts_timing.sh
bash 07_run_IS_timeline.sh
echo "Pipeline finished. Results in $OUTDIR/"
