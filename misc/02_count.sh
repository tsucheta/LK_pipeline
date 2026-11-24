#!/usr/bin/env bash
# Step 02 - gene-level raw counts per organism (input for DESeq2). Never normalise these.
# Produces two tables per organism:
#   *_counts_unique.txt : MAPQ >= MINMAPQ (use for DESeq2; excludes reads from repeated IS copies,
#                         so near-identical transposase copies such as pgaptmp_000417 get ~0 counts
#                         and are dropped by the DESeq2 low-count filter)
#   *_counts_all.txt    : all primary alignments, multi-mappers placed at random by Bowtie2
#                         (use only for IS-family-level summaries)
# Organisms whose two count tables already exist are skipped (FORCE=1 to redo).
set -euo pipefail
source config.sh
source common.sh
[[ "$STRAND" == "AUTO" ]] && { echo "Set STRAND in config.sh from step 01 output first."; exit 1; }
mkdir -p "$OUTDIR/02_counts"
PE=""; [[ "$PAIRED" == "yes" ]] && PE="-p --countReadPairs -B -C"

for entry in "${ORGANISMS[@]}"; do
  IFS='|' read -r ORG GFF FTYPE ATTR BAMCOL EXCL <<< "$entry"
  echo "=== $ORG"
  U="$OUTDIR/02_counts/${ORG}_counts_unique.txt"; A="$OUTDIR/02_counts/${ORG}_counts_all.txt"
  if [[ -s "$U" && -s "$A" && "$FORCE" != 1 ]]; then
    echo "    already done, skipping ($U)"; continue
  fi
  ensure_bamcol "$BAMCOL"
  bam_list "$BAMCOL" > /dev/null          # stops here if a BAM is missing
  mapfile -t BAMS < <(bam_list "$BAMCOL")
  check_gff_contigs "$GFF" "${BAMS[0]}"
  featureCounts -T "$THREADS" $PE -s "$STRAND" -Q "$MINMAPQ" -t "$FTYPE" -g "$ATTR" \
    -a "$GFF" -o "$U" "${BAMS[@]}"
  featureCounts -T "$THREADS" $PE -s "$STRAND" -t "$FTYPE" -g "$ATTR" \
    -a "$GFF" -o "$A" "${BAMS[@]}"
  # a GFF that matches the contigs but not the feature type/attribute assigns nothing: catch it here
  assigned=$(awk -F'\t' '$1=="Assigned"{s=0; for(i=2;i<=NF;i++) s+=$i; print s}' "$U.summary")
  [[ "${assigned:-0}" -gt 0 ]] || { echo "ERROR: no reads assigned for $ORG - check feature type '$FTYPE' and attribute '$ATTR' in $GFF"; exit 1; }
  echo "    assigned read pairs (unique table, all libraries): $assigned"
done
echo "Counts written to $OUTDIR/02_counts/. Check the .summary files: Assigned should be similar across libraries."
