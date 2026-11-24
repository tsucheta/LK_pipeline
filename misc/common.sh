#!/usr/bin/env bash
# Helper functions shared by steps 01-03 (sourced, not run).

col_index () { head -1 "$SAMPLE_SHEET" | tr -d '\r' | tr '\t' '\n' | grep -nx "$1" | cut -d: -f1; }

# ensure_bamcol bam_Inquil : add the column to the sample sheet if missing, derived from bam_LK
# (bams/LK/LK_1.bam -> bams/Inquil/LK_1.bam). The original sheet is kept as <sheet>.bak.
ensure_bamcol () {
  local col=$1 dir=${1#bam_}
  sed -i 's/\r$//' "$SAMPLE_SHEET"
  [[ -n "$(col_index "$col")" ]] && return 0
  local lk; lk=$(col_index bam_LK)
  [[ -n "$lk" ]] || { echo "ERROR: $SAMPLE_SHEET has neither $col nor bam_LK; add a $col column by hand"; exit 1; }
  if ! tail -n +2 "$SAMPLE_SHEET" | cut -f"$lk" | grep -q '/LK/'; then
    echo "ERROR: cannot derive $col: bam_LK paths do not contain '/LK/'. Add a $col column to $SAMPLE_SHEET by hand."; exit 1
  fi
  [[ -f "$SAMPLE_SHEET.bak" ]] || cp "$SAMPLE_SHEET" "$SAMPLE_SHEET.bak"
  awk -F'\t' -v OFS='\t' -v c="$lk" -v name="$col" -v d="$dir" \
    'NR==1{print $0, name; next} {p=$c; sub("/LK/", "/" d "/", p); print $0, p}' "$SAMPLE_SHEET" > "$SAMPLE_SHEET.tmp"
  mv "$SAMPLE_SHEET.tmp" "$SAMPLE_SHEET"
  echo "    added column $col to $SAMPLE_SHEET (backup: $SAMPLE_SHEET.bak)"
}

# bam_list bam_Inquil : print the BAM paths of that column; stop if any BAM is missing
bam_list () {
  local c; c=$(col_index "$1")
  local missing=0
  while read -r b; do
    [[ -z "$b" ]] && continue
    [[ -s "$b" ]] || { echo "ERROR: BAM missing: $b (run 00_align.sh first)" >&2; missing=1; }
    echo "$b"
  done < <(tail -n +2 "$SAMPLE_SHEET" | cut -f"$c")
  [[ $missing -eq 0 ]] || exit 1
}

# check_gff_contigs GFF BAM : the GFF sequence names must match the BAM (= FASTA) contig names,
# otherwise featureCounts silently assigns nothing (common with Prokka GFFs and renamed contigs)
check_gff_contigs () {
  local gff=$1 bam=$2
  [[ -s "$gff" ]] || { echo "ERROR: GFF not found: $gff"; exit 1; }
  local n_gff n_hit
  n_gff=$(awk -F'\t' '!/^#/ && NF>=9{print $1}' "$gff" | sort -u | tee /tmp/gff_ctg.$$ | wc -l)
  n_hit=$(samtools view -H "$bam" | awk -F'\t' '$1=="@SQ"{sub("SN:","",$2); print $2}' | sort -u | comm -12 - /tmp/gff_ctg.$$ | wc -l)
  rm -f /tmp/gff_ctg.$$
  if [[ "$n_hit" -eq 0 ]]; then
    echo "ERROR: none of the $n_gff sequence names in $gff occur in $bam."
    echo "  GFF names : $(awk -F'\t' '!/^#/ && NF>=9{print $1}' "$gff" | sort -u | head -3 | tr '\n' ' ')"
    echo "  BAM names : $(samtools view -H "$bam" | awk -F'\t' '$1=="@SQ"{sub("SN:","",$2); print $2}' | grep -v -e QVFV -e LNDBAELA | head -3 | tr '\n' ' ')"
    echo "  Use a GFF made from the same FASTA that was put into refs/index/NEW.fna."
    exit 1
  fi
  echo "    GFF check: $n_hit of $n_gff GFF sequences found in the BAM header"
}

# write_products GFF OUT : two-column TSV (locus_tag, product) for the DESeq2 result tables
write_products () {
  python3 - "$1" "$2" <<'PY'
import sys
from urllib.parse import unquote
gff, out = sys.argv[1], sys.argv[2]
prod = {}
for l in open(gff):
    if l.startswith("##FASTA"): break
    if l.startswith("#") or l.count("\t") < 8: continue
    f = l.rstrip("\n").split("\t")
    a = {k: unquote(v) for k, v in (x.split("=", 1) for x in f[8].split(";") if "=" in x)}
    t = a.get("locus_tag")
    if t and (t not in prod or f[2] == "CDS") and a.get("product"):
        prod[t] = a["product"]
with open(out, "w") as o:
    for t, p in prod.items(): o.write(f"{t}\t{p}\n")
print(f"    {len(prod)} products written to {out}")
PY
}
