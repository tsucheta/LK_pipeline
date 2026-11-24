#!/usr/bin/env bash
# =====================================================================================================
#  LK pipeline - the ONLY file to edit. Every step (00-08) reads its input file names and settings here.
#  Paths are relative to the folder the scripts are run from (e.g. ~/LKPaper).
# =====================================================================================================
THREADS=16
OUTDIR=results                     # all tables and logs go here
FORCE=${FORCE:-0}                  # steps skip work that already exists; FORCE=1 bash <step> redoes it

# ---------------------------------------------------------------- input sheets
FASTQ_SHEET=fastq_new.tsv          # sample_id <tab> RNA|DNA <tab> R1.fastq.gz <tab> R2.fastq.gz   (step 00)
SAMPLE_SHEET=samples.tsv           # RNA libraries: sample, condition, bam_<organism> columns       (steps 01-03)
CONTRASTS=contrasts.tsv            # DESeq2 contrasts                                                (step 03)

# ---------------------------------------------------------------- reference genomes (FASTA)
# Current (2023, third-time) assemblies: RNA mapping, variant calling (combined into refs/index/NEW.fna)
LK_NEW=refs/LK_1.fna               # Leptolyngbya iicbica LK_1 (2023); matches the PGAP GFF below
PHAEO_NEW=refs/phaeo.fna           # Phaeodactylibacter sp. (bins.10)
OCEAN_NEW=refs/ocean.fna           # Oceanicaulis sp. (bins.9)
INQUIL_NEW=refs/inqui.fna          # Inquilinaceae sp. (bins.8)
NITRAT=refs/LK_bac_polished.fna    # Nitratireductor sp. isolate genome
# Old (first-time) assemblies: empty-site counts for dating IS insertions (combined into refs/index/OLD.fna)
LK_OLD=refs/LK_16.fna              # LK_16 (2014), GCA_004212215
OCEAN_OLD=refs/ocean_old.fna
PHAEO_OLD=refs/phaeo_old.fna
# Inquilinaceae and Nitratireductor have no old assembly: their 2023 genomes go into OLD, so step 07 reports
# their IS loci as "not testable". To date Inquilinaceae insertions, uncomment (re-maps all DNA to OLD, hours):
#INQUIL_OLD=refs/inqui_2014.fna

# ---------------------------------------------------------------- annotations (GFF, contig names as in the FASTA)
LK_GFF=refs/LK_1.gff
PHAEO_GFF=refs/phaeo.gff
OCEAN_GFF=refs/ocean.gff
INQUIL_GFF=refs/inqui.gff
NITRAT_GFF=refs/LK_bac_polished.gff

# ---------------------------------------------------------------- RNA-seq (steps 01-03)
# Libraries excluded in DESeq2. LK-7R1/7R2 are genomic DNA in every genome; LK4-1/LK4-2 (library LK_4) carry
# genomic DNA in the bacterial genomes, so they are excluded for the bacteria only (contrasts that then have one
# library per group are skipped automatically).
DNA_LIBS="LK-7R1_S78,LK-7R2_S79"
BAC_EXCL="$DNA_LIBS,LK4-1_S80,LK4-2_S81"
#   name | GFF | feature types | attribute | BAM column in the sample sheet | libraries excluded in DESeq2
ORGANISMS=(
  "LK|$LK_GFF|gene,pseudogene|locus_tag|bam_LK|$DNA_LIBS"
  "Phaeodactylibacter|$PHAEO_GFF|CDS|locus_tag|bam_Phaeo|$BAC_EXCL"
  "Oceanicaulis|$OCEAN_GFF|CDS|locus_tag|bam_Ocean|$BAC_EXCL"
  "Inquilinaceae|$INQUIL_GFF|CDS|locus_tag|bam_Inquil|$BAC_EXCL"
  "Nitratireductor|$NITRAT_GFF|CDS|locus_tag|bam_Nitrat|$BAC_EXCL"
)
STRAND=2                           # featureCounts -s (0 unstranded, 1 forward, 2 reverse); step 01 reports it
PAIRED=yes
MINMAPQ=10                         # MAPQ cut-off of the "unique" count table

# ---------------------------------------------------------------- DNA datasets and variant calling (steps 05-08)
# Paired-end DNA libraries, OLDEST FIRST (never mate-pair). BAMs: bams/DNA/<name>_vs_LK1.bam and _vs_LK16.bam
DNA_SAMPLES=(DNA2014_PE DNA2017_PE DNA2023_PE)
#   short name used in the variant tables | FASTA | GFF | folder of that organism's RNA BAMs (step 00)
VARIANT_ORGS=(
  "LK|$LK_NEW|$LK_GFF|bams/LK"
  "Inqui|$INQUIL_NEW|$INQUIL_GFF|bams/Inquil"
  "Phaeo|$PHAEO_NEW|$PHAEO_GFF|bams/Phaeo"
  "Ocean|$OCEAN_NEW|$OCEAN_GFF|bams/Ocean"
  "Nitrat|$NITRAT|$NITRAT_GFF|bams/Nitrat"
)
KNOWN_VARIANTS=known_variants_LK1.tsv     # labels only (column previously_reported); never used to call
RNA_EXCL_VARIANTS="$BAC_EXCL"             # RNA libraries left out of the RNA-support count
# Step 08 thresholds (defaults of 08_call_variants.py; change only with a reason, and report it)
V_MINDEPTH=15; V_MAXEARLY=0.05; V_MINLATE=0.20; V_MINALT=5; V_MAXLOWMQ=0.10; V_MAXDEPTHX=2.5; V_EDGE=150; V_MAXINDEL=3

# ---------------------------------------------------------------- mutS IS4 (step 06)
MUTS_CONTIG=LNDBAELA               # LK_1 contig carrying mutS
MUTS_START=562850; MUTS_END=564650 # region spanning pgaptmp_000416 - IS4 - pgaptmp_000418

# ---------------------------------------------------------------- split/pseudogene check (step 09)
SPLIT_CONTIG=LNDBAELA              # LK_1 contig whose PGAP pseudogenes are tested
# extra loci not annotated as pseudogenes (name:contig:start-end). valS (pgaptmp_000998-000999) is a PGAP
# pseudogene and is tested automatically; tmk is not annotated at all, so it is listed here.
SPLIT_LOCI=(
  "tmk_thymidylate_kinase:LNDBAELA:4131300-4132750"     # thymidylate kinase (34113at1117); LK_16 pgaptmp_000908
)
#DNA_SAM_DIR=sam                   # only if you have SAM files and no BAMs: folder with <sample>_vs_LK1.sam

# ---------------------------------------------------------------- polishing LK_1 (step 10)
PILON_JAR=/opt/pilon/pilon-1.24.jar		# path to the Pilon jar (v1.24; conda install -c bioconda pilon)
PILON_MEM=24G                      		# Java memory for Pilon (the combined NEW reference is ~20 Mb)
POLISH_ROUNDS=2					# Pilon rounds; stops earlier when a round changes nothing
#POLISH_SAMPLE=DNA2023_PE			# DNA dataset used for polishing (default: the newest in DNA_SAMPLES)
SPLICE_DONOR=$LK_OLD				# assembly whose sequence replaces misassembled segments (LK_16, 2014)

# LK_1 windows holding a misassembled segment (step 09: "N bp not covered by reads"); only differences >= 50 bp
# from the donor inside them are replaced. valS = pgaptmp_000998-000999; tmk = thymidylate kinase.
SPLICE_LOCI=(
  "LNDBAELA:1225361-1228552"       # valS
  "LNDBAELA:4131300-4132750"       # tmk
  "LNDBAELA:2643150-2643700"       # Clp protease N-terminal domain gene (pgaptmp_002201); gap in 2014/2017/2023 reads
  "LNDBAELA:4095050-4095220"       # pgaptmp_003455; gap in all three years
  "LNDBAELA:4320150-4320300"       # pgaptmp_003632; gap in all three years
)

