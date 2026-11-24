# LK pipeline

Scripts for: Tripathy et al., "Insertion-sequence disruption of mutS and a rising lineage with a
mismatch-repair-loss mutation spectrum in a decade-old laboratory-grown cyanobacterial culture",
Microbiological Research (submitted).

Edit `config/config.sh` (input paths) and `config/samples.tsv`, then run the steps in order.
Full description and run order: Supplementary Methods S2.14 and S3.

| Step | Script | Output |
|---|---|---|
| 00–03 | rnaseq/ | alignment, QC, featureCounts, DESeq2 (Tables S1.24–S1.27, S5) |
| 06 | 06_muts_IS_timing.py | IS4–mutS junction dating (Table S1.22) |
| 07 | 07_IS_timeline.py | IS timelines (Table S1.31) |
| 08 | 08_call_variants.py | variant calls (Tables S1.21, S1.23, S1.28) |
| 09–10 | assembly_correction/ | Illumina correction of LK_1 |
|  | tcmP_expression.py | IS5–tcmP read-through (Table S1.32) |

Requirements: Python 3, pysam, numpy, pandas, scipy; Bowtie2 2.3.4.1, samtools,
Subread 2.0.3, R with DESeq2 1.46.0, BLAST+, ISEScan 1.7.2.3.


Scripts used for the RNA-seq, variant-calling and IS-insertion analyses of the *Leptolyngbya iicbica* LK consortium
(DNA datasets 2014, 2017 and 2023; RNA-seq 2023). **All input file names and settings are in `config.sh`**;
no other file needs editing.

## Install

```bash
conda create -n lk -c bioconda -c conda-forge python=3.10 fastp bowtie2 samtools'>=1.7' subread \
      r-base bioconductor-deseq2 r-ggplot2 scipy pysam
conda activate lk
```
Steps 06-08 need only Python 3 (standard library) and samtools. Step 05 needs pysam.
Versions used for the paper: fastp 0.23.4, Bowtie2 2.3.4.1, samtools 1.7, featureCounts (Subread) 2.0.3, DESeq2 1.46.0.
Also tested end to end with samtools 1.19, Bowtie2 2.5.2, featureCounts 2.0.6 and DESeq2 1.42.

## Inputs (named in config.sh)

| What | config.sh variable | Default |
|---|---|---|
| FASTQ sheet: `sample_id  RNA/DNA  R1  R2` (tab-separated) | `FASTQ_SHEET` | fastq_new.tsv |
| RNA sample sheet (sample, condition, `bam_<organism>` columns) | `SAMPLE_SHEET` | samples.tsv |
| DESeq2 contrasts | `CONTRASTS` | contrasts.tsv |
| 2023 genomes (LK_1 + four bacteria) | `LK_NEW`, `PHAEO_NEW`, `OCEAN_NEW`, `INQUIL_NEW`, `NITRAT` | refs/*.fna |
| First-time assemblies | `LK_OLD`, `OCEAN_OLD`, `PHAEO_OLD` (optional `INQUIL_OLD`) | refs/*.fna |
| Annotations | `LK_GFF`, `PHAEO_GFF`, `OCEAN_GFF`, `INQUIL_GFF`, `NITRAT_GFF` | refs/*.gff |
| DNA datasets, oldest first | `DNA_SAMPLES` | DNA2014_PE DNA2017_PE DNA2023_PE |
| Earlier variant list (labels only) | `KNOWN_VARIANTS` | known_variants_LK1.tsv |

## Run

```bash
cd ~/LKPaper            # folder holding config.sh, the scripts, refs/ and the FASTQ files
conda activate lk
bash run_all.sh         # or one step at a time, in this order:
```

| Step | Script | Does | Main output (in `results/`) |
|---|---|---|---|
| 00 | `00_align.sh` | fastp trimming; Bowtie2 mapping of RNA and DNA to the combined 2023 genomes (NEW) and of DNA to the old assemblies (OLD); RNA BAMs split per organism; duplicates marked | `bams/` |
| 01 | `01_library_qc.sh` | strandedness and genomic-DNA check per library | `01_qc/` |
| 02 | `02_count.sh` | featureCounts (unique: MAPQ ≥ 10; all) | `02_counts/` |
| 03 | `03_run_deseq2.sh` → `03_deseq2.R` | DESeq2 per organism and contrast, with library exclusions | `03_deseq2/<organism>/` |
| 08 | `08_run_variants.sh` → `08_call_variants.py` | genome-wide calling of new variants, all organisms | `08_variants.tsv`, `08_<org>.vcf`, `08_spectrum.tsv`, `08_depth.tsv` |
| 09 | `09_run_split_check.sh` → `09_check_split_genes.py` | are the LK_1 pseudogenes/split genes real or assembly errors (Illumina base disagreements, coverage gaps, repeats) | `09_loci.tsv`, `09_changes.tsv` |
| 10 | `10_run_polish.sh` → `10_splice_segments.py`, Pilon, `10_lift_annotation.py` | replace the misassembled segments (valS, tmk) with the LK_16 sequence, polish LK_1 with the newest Illumina DNA reads, carry the annotation over (shift untouched genes, re-check the frame of touched ones, merge split genes that become intact), re-check with step 09, convert step-08 positions | `polish/LK_1_polished.fna`, `.gff`, `_map.tsv`, `_genes.tsv`, `09_after_loci.tsv`, `08_variants_polished.tsv` |
| 06 | `06_run_muts_timing.sh` → `06_muts_IS_timing.py` | IS4-in-mutS frequency per dataset, compared with the step-08 LK mutations | `06_mutS_IS_timing.tsv`, `06_timeline.tsv` |
| 07 | `07_run_IS_timeline.sh` → `07_IS_timeline.py` | new vs pre-existing IS insertions, all organisms | `07_IS_timeline.tsv`, `07_IS_summary.tsv` |
| (05) | `05_is_junction_support.py junctions.tsv out.tsv` | optional: junction/empty-site reads for insertions listed by hand | out.tsv |

Step 06 reads the step-08 table, so 08 runs before 06. Every step skips work already done;
`FORCE=1 bash <script>` redoes it. `fix_gff_names.py` is only needed if a GFF's contig names differ from its FASTA.

## Variant calling (step 08)

`08_call_variants.py` is the only variant-calling step: earlier GATK/SnpEff calls are not inputs and appear only as a
label (`previously_reported`) on calls that this step makes independently. The method and every threshold are described in
Supplementary Methods S2 (variant calling); the thresholds are the `V_*` variables in config.sh and are written
into the log. All candidates are kept in `08_variants.tsv`, with the reason in `FILTER` when a call is not PASS.

| FILTER | Rule |
|---|---|
| depth | ≥ 15 reads (MAPQ ≥ 20, base quality ≥ 20) in the first and the last dataset |
| new | derived allele ≤ 5% in the first dataset; ≥ 20% and ≥ 5 reads in the last |
| lowMQ | > 10% of reads at the site with MAPQ < 20 |
| highDP | depth > 2.5 × the organism's median in that dataset |
| strand | derived allele on one strand only, or Fisher strand-bias P < 0.001 (derived vs ancestral reads) |
| edge | within 150 bp of a contig end |
| cluster | ≥ 3 candidates within 100 bp |
| hpNoise | homopolymer indel below 3 × the slippage rate of the same slip at up to 100 runs of the same base and length |
| presentFirst / rareLast / lowDepthSpanning | indel that fails the depth or frequency rules after the spanning-read recount |

## Libraries excluded

LK-7R1/LK-7R2 are genomic DNA (excluded in every organism). LK4-1/LK4-2 contain genomic DNA in the bacterial genomes
(excluded for the bacteria, in DESeq2 and in the RNA support of step 08). Set in config.sh (`DNA_LIBS`, `BAC_EXCL`).

## Step 10: polishing LK_1 and converting coordinates

Step 10 is run on its own (not by `run_all.sh`), after steps 00-09, because it changes the reference:
```bash
bash 10_run_polish.sh
```
It needs Java and the Pilon jar (`PILON_JAR` in `config.sh`). The ~1-kb segments that the long-read assembly got
wrong (`SPLICE_LOCI`) are first replaced with the LK_16 sequence (`polish/logs_10_splice.txt` says what was replaced);
Pilon then corrects bases and short indels only. Read `polish/LK_1_polished_genes.tsv` (every gene a
correction touched and what happened to it) and `polish/09_after_loci.tsv` (should list no assembly errors).
Any table of LK_1 positions can then be converted to polished coordinates:
```bash
python3 10_lift_annotation.py --map polish/LK_1_polished_map.tsv --convert positions.tsv > positions_polished.tsv
```
(`positions.tsv`: tab-separated, contig in column 1 and position in column 2; other columns are kept.)
To re-run steps 00-09 on the polished genome, point `LK_NEW` and `LK_GFF` in `config.sh` at the polished files.
