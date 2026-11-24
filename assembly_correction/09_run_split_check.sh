#!/usr/bin/env bash
# Step 09 driver: test every pseudogene/split gene of LK_1 (and the loci in SPLIT_LOCI) against the Illumina DNA
# reads. Needs the DNA BAMs from step 00 (bams/DNA/<sample>_vs_LK1.bam). If you only have SAM files, set
# DNA_SAM_DIR in config.sh: <sample>_vs_LK1.sam there is sorted and indexed into bams/DNA first.
# Output: $OUTDIR/09_loci.tsv, $OUTDIR/09_changes.tsv
set -euo pipefail
source config.sh
for v in DNA_SAMPLES LK_NEW LK_GFF SPLIT_CONTIG SPLIT_LOCI; do
  if ! declare -p "$v" >/dev/null 2>&1; then
    echo "ERROR: config.sh does not define $v - use the config.sh shipped with this package. Nothing has been run."; exit 1
  fi
done
mkdir -p "$OUTDIR" bams/DNA
if [[ -s "$OUTDIR/09_loci.tsv" && "$FORCE" != 1 ]]; then
  echo "09: $OUTDIR/09_loci.tsv exists, skipping (FORCE=1 to redo)"; exit 0
fi
REF=refs/index/NEW.fna; [[ -s $REF ]] || REF=$LK_NEW
[[ -s $REF.fai ]] || samtools faidx "$REF"
BAMS=()
for s in "${DNA_SAMPLES[@]}"; do
  b=bams/DNA/${s}_vs_LK1.bam
  if [[ ! -s $b && -n "${DNA_SAM_DIR:-}" && -s "$DNA_SAM_DIR/${s}_vs_LK1.sam" ]]; then
    echo "09: sorting $DNA_SAM_DIR/${s}_vs_LK1.sam -> $b"
    samtools sort -@ 4 -m 500M -o "$b.tmp.bam" "$DNA_SAM_DIR/${s}_vs_LK1.sam" && mv "$b.tmp.bam" "$b"
  fi
  [[ -s $b ]] || { echo "ERROR: $b not found (run 00_align.sh, or set DNA_SAM_DIR)"; exit 1; }
  [[ -s $b.bai ]] || samtools index "$b"
  BAMS+=("$s=$b")
done
python3 09_check_split_genes.py --ref "$REF" --gff "$LK_GFF" --contig "$SPLIT_CONTIG" --bam "${BAMS[@]}" \
  --locus "${SPLIT_LOCI[@]}" --minmapq "$MINMAPQ" --out "$OUTDIR/09" 2>&1 | tee "$OUTDIR/logs_09.txt"
