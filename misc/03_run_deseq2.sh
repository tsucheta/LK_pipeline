#!/usr/bin/env bash
# Step 03 driver - runs 03_deseq2.R for every organism in config.sh, with that organism's exclusions
# and a product table taken from its GFF.
# Organisms whose contrast_summary.csv already exists are skipped (FORCE=1 to redo all, or delete
# results/03_deseq2/<organism>/ to redo one).
# Usage: bash 03_run_deseq2.sh            (all organisms)
#        bash 03_run_deseq2.sh Inquilinaceae Nitratireductor     (only these)
set -euo pipefail
source config.sh
source common.sh
command -v Rscript >/dev/null || { echo "ERROR: Rscript not found (conda activate lk?)"; exit 1; }
OUT3="$OUTDIR/03_deseq2"; mkdir -p "$OUT3"
ONLY=("$@")

for entry in "${ORGANISMS[@]}"; do
  IFS='|' read -r ORG GFF FTYPE ATTR BAMCOL EXCL <<< "$entry"
  if [[ ${#ONLY[@]} -gt 0 ]] && [[ ! " ${ONLY[*]} " =~ " $ORG " ]]; then continue; fi
  echo "=== $ORG"
  if [[ -s "$OUT3/$ORG/contrast_summary.csv" && "$FORCE" != 1 ]]; then
    echo "    already done, skipping ($OUT3/$ORG/contrast_summary.csv)"; continue
  fi
  CNT="$OUTDIR/02_counts/${ORG}_counts_unique.txt"
  [[ -s "$CNT" ]] || { echo "ERROR: $CNT missing - run 02_count.sh first"; exit 1; }
  # too few counted genes -> DESeq2 cannot estimate dispersions; report and move on
  ngenes=$(python3 - "$CNT" "${EXCL:-}" <<'PY'
import sys
cnt, excl = sys.argv[1], set(x for x in sys.argv[2].split(",") if x)
rows = [l.rstrip("\n").split("\t") for l in open(cnt) if not l.startswith("#")]
hdr = rows[0]; keep = [i for i, h in enumerate(hdr) if i >= 6 and h.split("/")[-1].replace(".bam", "") not in excl]
print(sum(1 for r in rows[1:] if sum(int(r[i]) >= 10 for i in keep) >= 2))
PY
)
  if [[ "$ngenes" -lt 200 ]]; then
    echo "    only $ngenes genes have >=10 reads in >=2 libraries: too few for DESeq2, skipping $ORG"
    echo "    (see $OUTDIR/02_counts/${ORG}_counts_unique.txt.summary for where the reads went)"
    continue
  fi
  echo "    $ngenes genes pass the count filter"
  ensure_bamcol "$BAMCOL"
  write_products "$GFF" "$OUT3/${ORG}_products.tsv"
  echo "    excluding: ${EXCL:-none}"
  Rscript 03_deseq2.R --org "$ORG" --bamcol "$BAMCOL" --counts "$CNT" \
      --samples "$SAMPLE_SHEET" --contrasts "$CONTRASTS" --out "$OUT3" \
      --exclude "${EXCL:-}" --annot "$OUT3/${ORG}_products.tsv" 2>&1 | tee "$OUT3/${ORG}.log"
done
echo "Done. Per organism: $OUT3/<organism>/{library_qc,contrast_summary,DE_*,normalized_counts,size_factors}.csv and PCA.pdf"
