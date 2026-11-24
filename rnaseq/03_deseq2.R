#!/usr/bin/env Rscript
# Step 03 - library QC + DESeq2 for ONE organism, all libraries in one model.
#
# Usage:
#   Rscript 03_deseq2.R --org LK --bamcol bam_LK \
#       --counts results/02_counts/LK_counts_unique.txt \
#       [--samples samples.tsv] [--contrasts contrasts.tsv] [--out results/03_deseq2] \
#       [--exclude LK-7R1_S78,LK-7R2_S79] [--annot LK_products.tsv]
#
# Input counts are RAW featureCounts output. DESeq2 normalises internally (median of ratios).
# --annot (optional): two-column TSV (gene id, product) to add gene names to result tables.

suppressPackageStartupMessages(library(DESeq2))

## ---------- arguments ----------
args <- commandArgs(trailingOnly = TRUE)
opt <- list(samples = "samples.tsv", contrasts = "contrasts.tsv", out = "results/03_deseq2",
            exclude = "", annot = "")
for (i in seq(1, length(args), by = 2)) opt[[sub("^--", "", args[i])]] <- args[i + 1]
for (k in c("org", "bamcol", "counts")) if (is.null(opt[[k]])) stop("missing --", k)
outdir <- file.path(opt$out, opt$org)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## ---------- read counts and sample sheet ----------
fc <- read.delim(opt$counts, comment.char = "#", check.names = FALSE)
ann_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
cts <- as.matrix(fc[, setdiff(colnames(fc), ann_cols)])
rownames(cts) <- fc$Geneid
gene_len <- setNames(fc$Length, fc$Geneid)

smp <- read.delim(opt$samples, check.names = FALSE, stringsAsFactors = FALSE)
bam <- smp[[opt$bamcol]]
idx <- match(colnames(cts), bam)
miss <- is.na(idx)
if (any(miss)) idx[miss] <- match(basename(colnames(cts))[miss], basename(bam))
if (any(is.na(idx))) stop("BAM columns not found in samples.tsv: ",
                          paste(colnames(cts)[is.na(idx)], collapse = ", "))
colnames(cts) <- smp$sample_id[idx]
smp <- smp[idx, ]
rownames(smp) <- smp$sample_id
storage.mode(cts) <- "integer"

## ---------- library-level QC (raw counts, before any exclusion) ----------
keepq <- rowSums(cts) > 0
lcpm <- log2(t(t(cts[keepq, ]) / colSums(cts[keepq, ])) * 1e6 + 1)
cormat <- cor(lcpm)
qc <- data.frame(
  sample = colnames(cts),
  group = smp$group, run = smp$run, library_id = smp$library_id,
  assigned_reads = colSums(cts),
  genes_ge10 = colSums(cts >= 10),
  top10_share_pct = apply(cts, 2, function(x) 100 * sum(sort(x, decreasing = TRUE)[1:10]) / sum(x)),
  rho_count_vs_length = apply(cts[keepq, ], 2, function(x)
    suppressWarnings(cor(log2(x + 1), log2(gene_len[keepq]), method = "spearman"))),
  cv_reads_per_kb = apply(cts[keepq, ], 2, function(x) {
    d <- (x / sum(x)) / gene_len[keepq]; sd(d) / mean(d) }),
  median_cor_to_others = sapply(seq_len(ncol(cormat)), function(j) median(cormat[j, -j]))
)
qc$flag <- ifelse(qc$rho_count_vs_length > 0.6 | qc$median_cor_to_others < 0.8,
                  "SUSPECT (gDNA/degradation?)", "ok")
write.csv(qc, file.path(outdir, "library_qc.csv"), row.names = FALSE)
pdf(file.path(outdir, "library_correlation_heatmap.pdf"), width = 9, height = 8)
heatmap(cormat, symm = TRUE, margins = c(10, 10), main = paste(opt$org, "log-CPM correlation"))
dev.off()
cat("\nLibrary QC (", opt$org, "):\n", sep = ""); print(qc[, c("sample", "group", "run",
    "rho_count_vs_length", "median_cor_to_others", "flag")], row.names = FALSE)

## ---------- exclusions and technical re-runs ----------
excl <- if (nzchar(opt$exclude)) strsplit(opt$exclude, ",")[[1]] else character(0)
if (length(excl)) {
  cat("\nExcluding:", paste(excl, collapse = ", "), "\n")
  keep_s <- !colnames(cts) %in% excl
  cts <- cts[, keep_s]; smp <- smp[keep_s, ]
}
smp$group <- factor(smp$group)
smp$run <- factor(smp$run)
dds <- DESeqDataSetFromMatrix(cts, smp, design = ~ group)
if (any(duplicated(smp$library_id))) {
  cat("Collapsing technical re-runs sharing a library_id:",
      paste(unique(smp$library_id[duplicated(smp$library_id)]), collapse = ", "), "\n")
  dds <- collapseReplicates(dds, groupby = dds$library_id, run = dds$sample_id)
}
reps <- table(colData(dds)$group)
cat("\nBiological replicates per group after collapsing:\n"); print(reps)

## ---------- fit ----------
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
dds <- DESeq(dds, quiet = TRUE)
norm <- counts(dds, normalized = TRUE)
write.csv(data.frame(gene = rownames(norm), norm, check.names = FALSE),
          file.path(outdir, "normalized_counts.csv"), row.names = FALSE)
write.csv(data.frame(library = colnames(dds), size_factor = sizeFactors(dds)),
          file.path(outdir, "size_factors.csv"), row.names = FALSE)

annot <- NULL
if (nzchar(opt$annot)) { annot <- read.delim(opt$annot, header = FALSE, quote = "")
                         annot <- setNames(annot[[2]], annot[[1]]) }
shrink_type <- if (requireNamespace("ashr", quietly = TRUE)) "ashr" else "normal"

## ---------- contrasts ----------
con <- read.delim(opt$contrasts, stringsAsFactors = FALSE)
cd <- as.data.frame(colData(dds))
summ <- list()
for (k in seq_len(nrow(con))) {
  num <- con$numerator[k]; den <- con$denominator[k]; nm <- con$name[k]
  if (!all(c(num, den) %in% levels(cd$group)) || min(reps[c(num, den)]) < 2) {
    cat("Skipping", nm, "(group missing or <2 replicates)\n"); next }
  res <- results(dds, contrast = c("group", num, den), alpha = 0.05)
  shr <- lfcShrink(dds, contrast = c("group", num, den), res = res, type = shrink_type, quiet = TRUE)
  cols <- rownames(cd)[cd$group %in% c(num, den)]
  tab <- data.frame(gene = rownames(res),
                    product = if (is.null(annot)) NA else unname(annot[rownames(res)]),
                    baseMean = res$baseMean, log2FC = res$log2FoldChange,
                    log2FC_shrunk = shr$log2FoldChange, lfcSE = res$lfcSE,
                    pvalue = res$pvalue, padj = res$padj,
                    norm[rownames(res), cols, drop = FALSE], check.names = FALSE)
  tab <- tab[order(tab$padj), ]
  write.csv(tab, file.path(outdir, paste0("DE_", nm, ".csv")), row.names = FALSE)
  # diagnostic: fold change should NOT be predicted by baseline expression
  base_den <- rowMeans(norm[rownames(res), rownames(cd)[cd$group == den], drop = FALSE])
  ok <- is.finite(res$log2FoldChange) & base_den > 0
  rho <- cor(res$log2FoldChange[ok], log2(base_den[ok]), method = "spearman")
  runs_num <- unique(as.character(cd$run[cd$group == num]))
  runs_den <- unique(as.character(cd$run[cd$group == den]))
  summ[[nm]] <- data.frame(
    contrast = nm, numerator = num, denominator = den,
    genes_tested = sum(!is.na(res$padj)),
    padj_lt_0.05 = sum(res$padj < 0.05, na.rm = TRUE),
    pct_DE = round(100 * mean(res$padj < 0.05, na.rm = TRUE), 1),
    up_lfc_gt1 = sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE),
    down_lfc_lt_minus1 = sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm = TRUE),
    rho_lfc_vs_baseline = round(rho, 2),
    runs = if (setequal(runs_num, runs_den)) "same run(s)" else
      paste0("CROSSES RUNS (", paste(runs_num, collapse = "+"), " vs ", paste(runs_den, collapse = "+"), ")"))
}
summ <- do.call(rbind, summ)
write.csv(summ, file.path(outdir, "contrast_summary.csv"), row.names = FALSE)
cat("\nContrast summary (", opt$org, "):\n", sep = ""); print(summ, row.names = FALSE)
cat("\nRed flags: pct_DE > 40%, |rho_lfc_vs_baseline| > 0.5, or 'CROSSES RUNS'.\n")

## ---------- PCA ----------
vsd <- if (ncol(dds) >= 3 && nrow(dds) > 1000) vst(dds, blind = FALSE) else varianceStabilizingTransformation(dds, blind = FALSE)
pc <- plotPCA(vsd, intgroup = c("group", "run"), returnData = TRUE)
pv <- round(100 * attr(pc, "percentVar"))
pdf(file.path(outdir, "PCA.pdf"), width = 8, height = 6)
grp <- factor(pc$group); rn <- factor(pc$run)
plot(pc$PC1, pc$PC2, col = rainbow(nlevels(grp))[grp], pch = c(16, 17, 15, 18)[as.integer(rn)],
     cex = 1.6, xlab = paste0("PC1 (", pv[1], "%)"), ylab = paste0("PC2 (", pv[2], "%)"),
     main = paste(opt$org, "- colour = condition, shape = sequencing run"))
text(pc$PC1, pc$PC2, labels = colnames(vsd), pos = 3, cex = 0.6)
legend("topright", legend = levels(grp), col = rainbow(nlevels(grp)), pch = 16, cex = 0.7, bty = "n")
legend("bottomright", legend = levels(rn), pch = c(16, 17, 15, 18)[seq_len(nlevels(rn))], cex = 0.7, bty = "n")
dev.off()
cat("\nDone. Outputs in", outdir, "\n")

