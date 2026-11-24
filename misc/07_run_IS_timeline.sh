#!/usr/bin/env bash
# Step 07 driver: IS insertion timeline for every organism in VARIANT_ORGS.
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
GFF=()
for e in "${VARIANT_ORGS[@]}"; do IFS='|' read -r n fa gff rna <<< "$e"; [[ -s "$gff" ]] && GFF+=("$n=$gff"); done
python3 07_IS_timeline.py --new refs/index/NEW.fna --old refs/index/OLD.fna --gff "${GFF[@]}" \
  --samples "${DNA_SAMPLES[@]}" --out "$OUTDIR/07" 2>&1 | tee "$OUTDIR/logs_07.txt"
