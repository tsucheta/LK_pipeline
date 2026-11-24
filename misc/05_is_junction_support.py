#!/usr/bin/env python3
"""Step 05 - read-level test of whether an IS insertion was present in older datasets.

For each insertion, the SAME reads (e.g. the 2014 or 2017 Illumina run) must be mapped to
  (a) the NEW assembly that contains the IS (LK_1_SCAFFOLD)      -> reads spanning the IS junctions
  (b) the OLD assembly without the IS   (LK_16 / LK_3 / bin)     -> reads spanning the empty target site
A read "spans" a position if it aligns continuously from pos-FLANK to pos+FLANK with little soft-clipping.

    frequency(with IS) = mean(left, right junction reads) / (mean(left, right) + empty-site reads)

Usage:
  python3 05_is_junction_support.py junctions.tsv out.tsv
junctions.tsv columns (tab-separated, header line required):
  name  dataset  bam_new  contig_new  left_junction  right_junction  bam_old  contig_old  target_site
  (left/right_junction = first/last base of the IS in the new assembly; target_site = insertion point
   in the old assembly; all 1-based)
Requires: pysam
"""
import sys, csv
import pysam

FLANK, MAXCLIP, MINMAPQ = 20, 5, 20

def spanning(bam, contig, pos1):
    """count reads aligned continuously across pos1 (1-based) +/- FLANK"""
    p = pos1 - 1
    n = 0
    with pysam.AlignmentFile(bam) as f:
        for r in f.fetch(contig, max(0, p - FLANK), p + FLANK + 1):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
                continue
            if r.mapping_quality < MINMAPQ:
                continue
            if r.reference_start > p - FLANK or r.reference_end < p + FLANK + 1:
                continue
            cig = r.cigartuples or []
            clip = sum(l for op, l in cig if op == 4)
            if clip > MAXCLIP:
                continue
            # no deletion/skip across the junction window
            blocks = r.get_blocks()
            if not any(s <= p - FLANK and e >= p + FLANK + 1 for s, e in blocks):
                continue
            n += 1
    return n

def main(inp, outp):
    rows = list(csv.DictReader(open(inp), delimiter="\t"))
    with open(outp, "w") as o:
        o.write("name\tdataset\tleft_junction_reads\tright_junction_reads\tempty_site_reads\tfreq_with_IS\tinterpretation\n")
        for r in rows:
            L = spanning(r["bam_new"], r["contig_new"], int(r["left_junction"]))
            R = spanning(r["bam_new"], r["contig_new"], int(r["right_junction"]))
            E = spanning(r["bam_old"], r["contig_old"], int(r["target_site"]))
            w = (L + R) / 2
            f = w / (w + E) if (w + E) else float("nan")
            if w + E < 10:
                msg = "too few reads to decide"
            elif f < 0.05:
                msg = "IS absent in this dataset"
            elif f > 0.95:
                msg = "IS present (fixed) in this dataset"
            else:
                msg = "IS polymorphic in this dataset"
            o.write(f"{r['name']}\t{r['dataset']}\t{L}\t{R}\t{E}\t{f:.3f}\t{msg}\n")
    print(open(outp).read())

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])

