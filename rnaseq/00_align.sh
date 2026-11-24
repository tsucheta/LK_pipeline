#!/usr/bin/env bash
# Step 00 - build all BAM files from FASTQ.
#   RNA : competitive mapping to LK_1 + all bacterial genomes, then split into one BAM per organism
#         (bams/LK, Phaeo, Ocean, Inquil, Nitrat)
#   DNA : mapped twice - to the OLD references (IS empty sites, steps 05-07)
#         and to the NEW references (variant calling, step 08; IS junctions); duplicates marked
# Use paired-end libraries only (never mate-pair: wrong orientation, chimeric junctions).
#
# Resuming: every step checks what already exists and only does what is missing.
#   - RNA sample whose per-organism BAMs all exist                     -> skipped
#   - RNA sample with <sample>.all.bam that holds every organism's contigs -> only the missing
#     per-organism BAMs are cut out (no re-mapping)
#   - RNA/DNA BAM mapped to an index that lacked Inquilinaceae or Nitratireductor contigs -> re-mapped
#     (its reads may sit on the wrong organism), trimmed reads are reused
# Requires: fastp, bowtie2, samtools (>= 1.7)
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
CUR_LOG=""
trap 'echo; echo "FAILED at line $LINENO: $BASH_COMMAND"; [[ -n "$CUR_LOG" && -f "$CUR_LOG" ]] && { echo "--- last lines of $CUR_LOG ---"; tail -n 15 "$CUR_LOG"; }' ERR
for tool in fastp bowtie2 bowtie2-build bowtie2-inspect samtools; do
  command -v $tool >/dev/null || { echo "ERROR: $tool not found in PATH (conda activate lk?)"; exit 1; }
done
# old fastp builds (< 0.19.5) lack --detect_adapter_for_pe; they still trim PE adapters by read-overlap
ADAPT_OPT=""
fastp --help 2>&1 | grep -q -- '--detect_adapter_for_pe' && ADAPT_OPT="--detect_adapter_for_pe" || \
  echo "WARNING: old fastp ($(fastp --version 2>&1 | head -1)); adapters trimmed by overlap only. Upgrade: conda install -c bioconda 'fastp>=0.23'"
# fastq.tsv saved from Excel/Windows has CRLF line endings: strip them
sed -i 's/\r$//' "$FASTQ_SHEET"
ORGS=(LK Phaeo Ocean Inquil Nitrat)
mkdir -p refs/index trimmed qc bams/{LK,Phaeo,Ocean,Inquil,Nitrat,DNA}
for f in "$LK_NEW" "$PHAEO_NEW" "$OCEAN_NEW" "$INQUIL_NEW" "$NITRAT" "$LK_OLD" "$OCEAN_OLD" "$PHAEO_OLD" ${INQUIL_OLD:+"$INQUIL_OLD"}; do
  [[ -s "$f" ]] || { echo "ERROR: reference not found: $f (edit config.sh)"; exit 1; }
done

# ---- 1. combined references (competitive mapping stops reads of one organism landing on another) ----
# The index is rebuilt only if it is missing or its contig list differs from the FASTA
# (e.g. an index built before Inquilinaceae/Nitratireductor were added).
make_ref () {   # name fasta...
  local name=$1; shift
  local fa="refs/index/${name}.fna"
  cat "$@" > "$fa.tmp"
  dup=$(grep '^>' "$fa.tmp" | cut -d' ' -f1 | sort | uniq -d | head -3)
  [[ -n "$dup" ]] && { echo "ERROR: duplicate contig names in $name: $dup"; exit 1; }
  if [[ -f "$fa" ]] && cmp -s "$fa" "$fa.tmp"; then rm "$fa.tmp"; else mv "$fa.tmp" "$fa"; rm -f "$fa.fai"; fi
  [[ -f "$fa.fai" ]] || samtools faidx "$fa"
  local idx_ok=no
  if [[ -f "refs/index/${name}.1.bt2" || -f "refs/index/${name}.1.bt2l" ]]; then
    if diff -q <(bowtie2-inspect -n "refs/index/${name}" | cut -d' ' -f1 | sort) <(cut -f1 "$fa.fai" | sort) >/dev/null; then
      idx_ok=yes
    else
      echo "    index ${name} does not match ${fa} (contigs added or removed): rebuilding"
      rm -f refs/index/${name}.*.bt2 refs/index/${name}.*.bt2l
    fi
  fi
  if [[ $idx_ok == no ]]; then
    echo "    building Bowtie2 index ${name} (takes a while)"
    bowtie2-build --threads "$THREADS" "$fa" "refs/index/${name}" > "refs/index/${name}.build.log" 2>&1
  fi
}
make_ref NEW "$LK_NEW" "$PHAEO_NEW" "$OCEAN_NEW" "$INQUIL_NEW" "$NITRAT"
# INQUIL_OLD (config.sh) is optional: if set, the 2014 Inquilinaceae bin replaces the 2023 one in OLD
# (this rebuilds the OLD index and re-maps every DNA library to it).
make_ref OLD "$LK_OLD" "$OCEAN_OLD" "$PHAEO_OLD" "${INQUIL_OLD:-$INQUIL_NEW}" "$NITRAT"

# BED of whole contigs per organism, used to split RNA BAMs
to_bed () { [[ -f "$1.fai" ]] || samtools faidx "$1"; awk 'BEGIN{OFS="\t"}{print $1,0,$2}' "$1.fai" > "$2"; }
to_bed "$LK_NEW" refs/index/LK.bed
to_bed "$PHAEO_NEW" refs/index/Phaeo.bed
to_bed "$OCEAN_NEW" refs/index/Ocean.bed
to_bed "$INQUIL_NEW" refs/index/Inquil.bed
to_bed "$NITRAT" refs/index/Nitrat.bed

# same_contigs BAM REF : true if the BAM header lists exactly the contigs of refs/index/REF.fna
# (false for BAMs mapped before Inquilinaceae/Nitratireductor were added to the index)
same_contigs () {
  [[ -s "$1" ]] || return 1
  diff -q <(samtools view -H "$1" 2>/dev/null | awk -F'\t' '$1=="@SQ"{sub("SN:","",$2); print $2}' | sort) \
          <(cut -f1 "refs/index/$2.fna.fai" | sort) >/dev/null
}

trim () {   # sample R1 R2 ; reuses trimmed reads if complete
  local S=$1 R1=$2 R2=$3
  CUR_LOG="qc/${S}.fastp.log"
  # fastp writes its JSON report only after it has finished both output files
  if [[ -s "qc/${S}.fastp.json" && -s "trimmed/${S}_R1.fq.gz" && -s "trimmed/${S}_R2.fq.gz" ]] && \
     gzip -t "trimmed/${S}_R1.fq.gz" "trimmed/${S}_R2.fq.gz" 2>/dev/null; then
    echo "    trimmed reads present, skipping fastp"
  else
    for f in "$R1" "$R2"; do
      [[ -s "$f" ]] || { echo "ERROR: $S: file not found or empty: '$f' (paths are relative to $(pwd))"; exit 1; }
    done
    fastp -i "$R1" -I "$R2" -o "trimmed/${S}_R1.fq.gz" -O "trimmed/${S}_R2.fq.gz" \
          $ADAPT_OPT -q 20 -l 36 -w "$THREADS" \
          -j "qc/${S}.fastp.json" -h "qc/${S}.fastp.html" 2> "qc/${S}.fastp.log"
  fi
}

split_rna () {   # sample ; cut per-organism BAMs out of <sample>.all.bam (only the missing ones)
  local S=$1
  for ORG in "${ORGS[@]}"; do
    [[ -s "bams/${ORG}/${S}.bam" && -f "bams/${ORG}/${S}.bam.bai" ]] && continue
    echo "    writing bams/${ORG}/${S}.bam"
    samtools view -@ 4 -b -L "refs/index/${ORG}.bed" -o "bams/${ORG}/${S}.bam" "bams/${S}.all.bam"
    samtools index "bams/${ORG}/${S}.bam"
  done
  samtools idxstats "bams/${S}.all.bam" > "qc/${S}.idxstats"
}

# ---- 2. per-sample trimming + mapping ----
tail -n +2 "$FASTQ_SHEET" | while IFS=$'\t' read -r S TYPE R1 R2; do
  [[ -z "$S" ]] && continue
  echo "=== $S ($TYPE)"
  [[ -n "${R2:-}" ]] || { echo "ERROR: $S has fewer than 4 columns - fastq.tsv must be TAB-separated"; exit 1; }

  if [[ "$TYPE" == "RNA" ]]; then
    missing=()
    for ORG in "${ORGS[@]}"; do [[ -s "bams/${ORG}/${S}.bam" ]] || missing+=("$ORG"); done
    if [[ ${#missing[@]} -eq 0 ]] && same_contigs "bams/LK/${S}.bam" NEW; then
      echo "    all organism BAMs present, skipping"; continue
    fi
    if same_contigs "bams/${S}.all.bam" NEW; then
      echo "    missing: ${missing[*]:-none} - cutting from existing bams/${S}.all.bam (no re-mapping)"
      [[ -f "bams/${S}.all.bam.bai" ]] || samtools index "bams/${S}.all.bam"
      split_rna "$S"; continue
    fi
    echo "    ${S}.all.bam absent or mapped without all organisms: mapping to refs/index/NEW"
    for ORG in "${ORGS[@]}"; do rm -f "bams/${ORG}/${S}.bam" "bams/${ORG}/${S}.bam.bai"; done
    trim "$S" "$R1" "$R2"
    CUR_LOG="qc/${S}.bowtie2.log"
    bowtie2 -p "$THREADS" --sensitive --no-unal -x refs/index/NEW \
        -1 "trimmed/${S}_R1.fq.gz" -2 "trimmed/${S}_R2.fq.gz" \
        --rg-id "$S" --rg "SM:$S" 2> "qc/${S}.bowtie2.log" \
      | samtools sort -@ 4 -o "bams/${S}.all.bam" -
    samtools index "bams/${S}.all.bam"
    split_rna "$S"

  elif [[ "$TYPE" == "DNA" ]]; then
    for REF in OLD NEW; do
      TAG=$([[ $REF == OLD ]] && echo LK16 || echo LK1)
      B="bams/DNA/${S}_vs_${TAG}.bam"
      if [[ -s "qc/${S}_vs_${TAG}.flagstat" ]] && same_contigs "$B" "$REF"; then
        echo "    ${S}_vs_${TAG} done, skipping"; continue
      fi
      [[ -s "$B" ]] && echo "    $B was mapped to an index without all organisms: re-mapping"
      trim "$S" "$R1" "$R2"
      CUR_LOG="qc/${S}_vs_${TAG}.bowtie2.log"
      bowtie2 -p "$THREADS" --very-sensitive --no-unal -x "refs/index/${REF}" \
          -1 "trimmed/${S}_R1.fq.gz" -2 "trimmed/${S}_R2.fq.gz" \
          --rg-id "$S" --rg "SM:$S" 2> "qc/${S}_vs_${TAG}.bowtie2.log" \
        | samtools fixmate -@ 4 -m - - \
        | samtools sort -@ 4 - \
        | samtools markdup -@ 4 -s - "$B" 2> "qc/${S}_vs_${TAG}.markdup.log"
      samtools index "$B"
      samtools flagstat "$B" > "qc/${S}_vs_${TAG}.flagstat"
    done
  else
    echo "Unknown type '$TYPE' for $S (use RNA or DNA)"; exit 1
  fi
done

# ---- 3. mapping summary ----
echo -e "sample\toverall_alignment_rate" > qc/alignment_summary.tsv
for f in qc/*.bowtie2.log; do
  echo -e "$(basename "$f" .bowtie2.log)\t$(grep 'overall alignment rate' "$f" | cut -d' ' -f1)" >> qc/alignment_summary.tsv
done
cat qc/alignment_summary.tsv
echo "Done. RNA BAMs: bams/<organism>/<sample>.bam ; DNA BAMs: bams/DNA/<sample>_vs_LK16.bam and _vs_LK1.bam"
