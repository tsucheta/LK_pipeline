#!/usr/bin/env bash
# Step 01 - library QC from BAMs: strandedness, antisense and intergenic read fractions.
# Flags genomic-DNA contamination or degraded RNA BEFORE any differential expression.
# Organisms whose *_library_qc.tsv already exists are skipped (FORCE=1 to redo).
# Requires: featureCounts (Subread >= 2.0.2), samtools
set -euo pipefail
source config.sh
source common.sh
mkdir -p "$OUTDIR/01_qc"
PE=""; [[ "$PAIRED" == "yes" ]] && PE="-p --countReadPairs -B -C"

for entry in "${ORGANISMS[@]}"; do
  IFS='|' read -r ORG GFF FTYPE ATTR BAMCOL EXCL <<< "$entry"
  echo "=== $ORG"
  if [[ -s "$OUTDIR/01_qc/${ORG}_library_qc.tsv" && "$FORCE" != 1 ]]; then
    echo "    already done, skipping ($OUTDIR/01_qc/${ORG}_library_qc.tsv)"; continue
  fi
  ensure_bamcol "$BAMCOL"
  bam_list "$BAMCOL" > /dev/null          # stops here if a BAM is missing
  mapfile -t BAMS < <(bam_list "$BAMCOL")
  check_gff_contigs "$GFF" "${BAMS[0]}"
  for s in 0 1 2; do
    featureCounts -T "$THREADS" $PE -s $s -t "$FTYPE" -g "$ATTR" -a "$GFF" \
      -o "$OUTDIR/01_qc/${ORG}_s${s}.txt" "${BAMS[@]}" > "$OUTDIR/01_qc/${ORG}_s${s}.log" 2>&1
  done
  # Summarise: assigned reads per strand setting + reads falling outside annotated features
  python3 - "$OUTDIR/01_qc" "$ORG" <<'PY'
import sys, csv
d, org = sys.argv[1], sys.argv[2]
def summ(s):
    rows = list(csv.reader(open(f"{d}/{org}_s{s}.txt.summary"), delimiter="\t"))
    names = [x.split("/")[-1].replace(".bam","") for x in rows[0][1:]]
    tab = {r[0]: [int(v) for v in r[1:]] for r in rows[1:]}
    return names, tab
names, s0 = summ(0); _, s1 = summ(1); _, s2 = summ(2)
out = open(f"{d}/{org}_library_qc.tsv", "w")
out.write("sample\ttotal\tassigned_unstranded\tpct_noFeature_(intergenic)\tassigned_fwd\tassigned_rev\tpct_antisense_if_stranded\tlikely_strand\n")
for i, n in enumerate(names):
    tot = sum(v[i] for v in s0.values())
    a0, a1, a2 = s0["Assigned"][i], s1["Assigned"][i], s2["Assigned"][i]
    nofeat = s0.get("Unassigned_NoFeatures", [0]*len(names))[i]
    hi, lo = max(a1, a2), min(a1, a2)
    strand = "1 (fwd)" if a1 > 3*a2 else ("2 (rev)" if a2 > 3*a1 else "0 (unstranded)")
    anti = 100*lo/(hi+lo) if (hi+lo) else 0
    out.write(f"{n}\t{tot}\t{a0}\t{100*nofeat/max(tot,1):.1f}\t{a1}\t{a2}\t{anti:.1f}\t{strand}\n")
out.close()
print(open(f"{d}/{org}_library_qc.tsv").read())
PY
done
cat <<'MSG'
How to read *_library_qc.tsv
  likely_strand        -> set STRAND in config.sh to this value (it should agree across samples)
  pct_antisense        -> stranded kits normally give <5%. Much higher in some libraries = genomic DNA.
  pct_noFeature        -> compare between libraries; an outlier (e.g. 2-3x the others) = gDNA or rRNA.
Libraries flagged here (expect LK-7R1/LK-7R2, and LK4-1/LK4-2 in the bacteria) are excluded in step 03 via config.sh.
MSG
