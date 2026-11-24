#!/usr/bin/env bash
# Step 08 driver: genome-wide variant calling for every organism in VARIANT_ORGS (config.sh).
# Needs the DNA BAMs from step 00 (bams/DNA/<sample>_vs_LK1.bam). Output: $OUTDIR/08_*
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
if [[ -s "$OUTDIR/08_variants.tsv" && "$FORCE" != 1 ]]; then
  echo "08: $OUTDIR/08_variants.tsv exists, skipping (FORCE=1 to redo)"; exit 0
fi
[[ -s refs/index/NEW.fna ]] || { echo "ERROR: refs/index/NEW.fna missing - run 00_align.sh first"; exit 1; }
ORG=(); GFF=(); RNA=()
for e in "${VARIANT_ORGS[@]}"; do
  IFS='|' read -r n fa gff rna <<< "$e"
  [[ -s "$fa" ]] || { echo "ERROR: $fa not found"; exit 1; }
  ORG+=("$n=$fa"); [[ -s "$gff" ]] && GFF+=("$n=$gff"); [[ -d "$rna" ]] && RNA+=("$n=$rna")
done
python3 08_call_variants.py --ref refs/index/NEW.fna --bamsuffix LK1 \
  --org "${ORG[@]}" --gff "${GFF[@]}" --rna "${RNA[@]}" --rna-exclude "$RNA_EXCL_VARIANTS" \
  --known "$KNOWN_VARIANTS" --samples "${DNA_SAMPLES[@]}" --out "$OUTDIR/08" \
  --mindepth "$V_MINDEPTH" --maxearly "$V_MAXEARLY" --minlate "$V_MINLATE" --minalt "$V_MINALT" \
  --maxlowmq "$V_MAXLOWMQ" --maxdepthx "$V_MAXDEPTHX" --edge "$V_EDGE" --maxindel "$V_MAXINDEL" \
  2>&1 | tee "$OUTDIR/logs_08.txt"
