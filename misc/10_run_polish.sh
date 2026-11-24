#!/usr/bin/env bash
# Step 10: polish LK_1 with the Illumina DNA reads (Pilon), carry the annotation over, and re-check.
#   1. Segments the long-read assembly got wrong (~1 kb; SPLICE_LOCI) are replaced with the homologous sequence of
#      the older assembly (SPLICE_DONOR, LK_16) by 10_splice_segments.py; only differences >= 50 bp are imported.
#   2. The reads are re-mapped and Pilon corrects bases and short indels (--fix snps,indels; Pilon's local
#      reassembly is not used: in tests it introduced errors in repeats and at misassembled segments), up to
#      POLISH_ROUNDS rounds, re-mapping after each; stops early when a round changes nothing.
#   3. 10_lift_annotation.py: coordinate map old -> polished, lifted GFF, genes whose frame changed.
#   4. Re-check: step 09 on the polished assembly (should report no assembly errors left) and conversion of the
#      step-08 variant positions to polished coordinates.
# Needs: java + pilon jar (PILON_JAR), bowtie2, samtools, and the trimmed reads of step 00 (trimmed/<sample>_R1.fq.gz).
# Outputs in polish/: LK_1_polished.fna, LK_1_polished.gff, LK_1_polished_map.tsv, LK_1_polished_genes.tsv,
#   09_after_loci.tsv, 08_variants_polished.tsv, round<N>/pilon.changes
set -euo pipefail
source config.sh
for v in DNA_SAMPLES LK_NEW LK_GFF LK_OLD SPLIT_CONTIG SPLIT_LOCI; do
  declare -p "$v" >/dev/null 2>&1 || { echo "ERROR: config.sh does not define $v (add the step 09/10 block). Nothing has been run."; exit 1; }
done
PILON_JAR=${PILON_JAR:-pilon.jar}; PILON_MEM=${PILON_MEM:-24G}; POLISH_ROUNDS=${POLISH_ROUNDS:-2}
S=${POLISH_SAMPLE:-${DNA_SAMPLES[${#DNA_SAMPLES[@]}-1]}}      # newest DNA dataset by default
C=$SPLIT_CONTIG; P=polish
SPLICE_DONOR=${SPLICE_DONOR:-$LK_OLD}
declare -p SPLICE_LOCI >/dev/null 2>&1 || SPLICE_LOCI=()
[[ -s $PILON_JAR ]] || { echo "ERROR: Pilon jar not found: $PILON_JAR (set PILON_JAR in config.sh; conda install -c bioconda pilon also works, then PILON_JAR=\$(ls \$CONDA_PREFIX/share/pilon*/pilon*.jar))"; exit 1; }
command -v java >/dev/null || { echo "ERROR: java not found"; exit 1; }
[[ -s trimmed/${S}_R1.fq.gz && -s trimmed/${S}_R2.fq.gz ]] || { echo "ERROR: trimmed/${S}_R1.fq.gz / _R2.fq.gz not found (made by step 00)"; exit 1; }
[[ -s refs/index/NEW.fna ]] || { echo "ERROR: refs/index/NEW.fna not found (run step 00)"; exit 1; }
if [[ -s $P/LK_1_polished.gff && "$FORCE" != 1 ]]; then echo "10: $P/LK_1_polished.gff exists, skipping (FORCE=1 to redo)"; exit 0; fi
mkdir -p $P

# replace contig $C in a multi-FASTA: swap_contig in.fna contig.fna out.fna
swap_contig () {
  awk -v c="$C" -v f="$2" 'BEGIN{while((getline l < f)>0){if(l~/^>/) print ">" c; else print l}}
       /^>/{skip=(substr($1,2)==c)} !skip' "$1" > "$3"
}
map_reads () {   # ref.fna out.bam
  local idx=${1%.fna}
  [[ -s $idx.1.bt2 || -s $idx.1.bt2l ]] || bowtie2-build --threads "$THREADS" "$1" "$idx" > "$idx.build.log" 2>&1
  bowtie2 -p "$THREADS" --very-sensitive --no-unal -x "$idx" -1 trimmed/${S}_R1.fq.gz -2 trimmed/${S}_R2.fq.gz \
      --rg-id "$S" --rg "SM:$S" 2> "${2%.bam}.bowtie2.log" \
    | samtools fixmate -@ 4 -m - - | samtools sort -@ 4 - | samtools markdup -@ 4 -s - "$2" 2> "${2%.bam}.markdup.log"
  samtools index "$2"
}

REF=refs/index/NEW.fna; BAM=bams/DNA/${S}_vs_LK1.bam
[[ -s $BAM ]] || { echo "ERROR: $BAM not found (step 00)"; exit 1; }
[[ -s $BAM.bai ]] || samtools index "$BAM"
if [[ ${#SPLICE_LOCI[@]} -gt 0 ]]; then
  echo "=== replacing misassembled segments with $SPLICE_DONOR"
  mkdir -p $P/round0
  python3 10_splice_segments.py --asm "$LK_NEW" --contig "$C" --donor "$SPLICE_DONOR" --locus "${SPLICE_LOCI[@]}" \
      --out $P/round0/contig.fna 2>&1 | tee $P/logs_10_splice.txt
  swap_contig "$REF" $P/round0/contig.fna $P/round0/NEWP.fna; samtools faidx $P/round0/NEWP.fna
  cp $P/round0/contig.fna $P/LK_1_polished.fna
  REF=$P/round0/NEWP.fna; BAM=$P/round0/${S}_vs_LK1p.bam
  echo "    re-mapping $S to the spliced assembly"; map_reads "$REF" "$BAM"
fi
for r in $(seq 1 "$POLISH_ROUNDS"); do
  D=$P/round$r; mkdir -p $D
  echo "=== round $r: Pilon on $C with $BAM"
  java -Xmx"$PILON_MEM" -jar "$PILON_JAR" --genome "$REF" --frags "$BAM" --targets "$C" --fix snps,indels --minmq 10 \
       --changes --output pilon --outdir "$D" > "$D/pilon.log" 2>&1 || { tail -20 "$D/pilon.log"; exit 1; }
  # Pilon writes the polished target as ${C}_pilon
  awk -v c="$C" '/^>/{on=(substr($1,2)==c"_pilon" || substr($1,2)==c)} on' "$D/pilon.fasta" > "$D/contig.fna"
  [[ -s $D/contig.fna ]] || { echo "ERROR: contig $C not in $D/pilon.fasta"; exit 1; }
  n=$(grep -c . "$D/pilon.changes" || true)
  echo "    round $r: $n changes"
  [[ "$n" -eq 0 ]] && { echo "    no further changes, stopping (the previous assembly and BAM are final)"; break; }
  swap_contig "$REF" "$D/contig.fna" "$D/NEWP.fna"; samtools faidx "$D/NEWP.fna"
  awk -v c="$C" '/^>/{print ">" c; next} {print}' "$D/contig.fna" > $P/LK_1_polished.fna
  REF=$D/NEWP.fna
  echo "    re-mapping $S to the round-$r assembly"
  BAM=$D/${S}_vs_LK1p.bam; map_reads "$REF" "$BAM"
done

echo "=== lifting the annotation"
python3 10_lift_annotation.py --old "$LK_NEW" --new $P/LK_1_polished.fna --gff "$LK_GFF" --contig "$C" \
    --locus "${SPLIT_LOCI[@]}" --out $P/LK_1_polished 2>&1 | tee $P/logs_10_lift.txt

echo "=== re-check with step 09 on the polished assembly"
LOCI=()
for L in "${SPLIT_LOCI[@]}"; do
  nm=${L%%:*}; rng=${L##*:}; s=${rng%-*}; e=${rng#*-}
  printf '%s\t%s\n%s\t%s\n' "$C" "$s" "$C" "$e" > $P/.l.tsv
  mapfile -t q < <(python3 10_lift_annotation.py --map $P/LK_1_polished_map.tsv --convert $P/.l.tsv | cut -f3)
  LOCI+=("$nm:$C:${q[0]}-${q[1]}")
done
python3 09_check_split_genes.py --ref "$REF" --gff $P/LK_1_polished.gff --contig "$C" --bam "$S=$BAM" \
    --locus "${LOCI[@]}" --minmapq "$MINMAPQ" --out $P/09_after 2>&1 | tee $P/logs_10_check.txt
echo "    assembly errors left: $(grep -c 'assembly error' $P/09_after_loci.tsv || true) (see $P/09_after_loci.tsv)"

if [[ -s $OUTDIR/08_variants.tsv ]]; then
  echo "=== step-08 variant positions in polished coordinates"
  awk -F'\t' -v OFS='\t' 'NR==1{print "contig","position","organism","change"; next} $1=="LK"{print $2,$3,$1,$6}' \
      $OUTDIR/08_variants.tsv > $P/.v.tsv
  python3 10_lift_annotation.py --map $P/LK_1_polished_map.tsv --convert $P/.v.tsv > $P/08_variants_polished.tsv
  bad=$(awk -F'\t' 'NR>1 && $5=="in changed block"' $P/08_variants_polished.tsv | wc -l)
  echo "    $(($(wc -l < $P/08_variants_polished.tsv)-1)) LK variants converted; $bad fall inside a polished block (check them)"
fi
rm -f $P/.l.tsv $P/.v.tsv
echo "Done. Polished assembly: $P/LK_1_polished.fna; annotation: $P/LK_1_polished.gff; coordinate map: $P/LK_1_polished_map.tsv"
