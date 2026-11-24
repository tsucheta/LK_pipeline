#!/usr/bin/env python3
"""
fig1_IS_circos.py - Figure 1: insertion-sequence (IS) elements in the LK consortium.

Circos plots (pyCirclize) of each genome: inner ring = contigs with IS elements as ticks
coloured by family and labelled family-order (e.g. IS481-3); links = pairs of IS elements
whose BLASTn alignment covers at least --mincov of both elements (in the combined panel, grey
links show every BLASTn hit between the two genomes, --pair-mincov 0); outer tracks = number of
genes on the forward (green) and reverse (red) strand per window. A combined panel shows two
genomes together with the links between them in grey.

Inputs
  --genome LABEL|ISESCAN_TSV|ANNOTATION[|FASTA]   (repeat for each genome)
           ANNOTATION = GenBank (.gbk) or GFF3 (.gff). With a GFF3, add the genome FASTA as the
           4th field; contig lengths are taken from it (otherwise from a ##FASTA section in the
           GFF, or from its ##sequence-region lines).
  --blast  all-against-all BLASTn table, -outfmt "6 qseqid sseqid pident length qstart qend
           sstart send evalue bitscore qlen slen" (12 columns), with the
           ISEScan .is.fna sequence names (contig_begin_end_strand) as query/subject IDs
  --pair   LABEL1,LABEL2   genomes to draw together in the last panel (optional)

Contig names are matched after removing any "gnl|X|" prefix, so PGAP and Prokka names agree.

Example (as used for the paper):
  python3 fig1_IS_circos.py \
    --genome "LK|LK_1_polished.fna.tsv|LK_1_corrected.gff|LK_1_corrected.fna" \
    --genome "Oceanicaulis sp.|ocean.fna.tsv|Oceanicaulis.gbk" \
    --genome "Inquilinaceae sp.|inqui.fna.tsv|Inquilinaceae.gbk" \
    --genome "Phaeodactylibacter sp.|phaeo.fna.tsv|Phaeodactylibacter.gbk" \
    --pair "Inquilinaceae sp.,Oceanicaulis sp." \
    --blast all_IS_vs_all.tsv --out Figure1_IS_circos

Outputs: <out>.png, <out>.pdf, <out>.tif (300 dpi) and <out>_links.tsv (every link drawn).
Requires: pycirclize (tested 1.10), matplotlib, numpy, biopython (GenBank input only).
"""
import argparse
import re
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import numpy as np
from pycirclize import Circos

FAMILY_COLOURS = {
    "IS630": "#2a78d6", "IS4": "#eb6834", "IS5": "#1baf7a", "IS66": "#eda100",
    "IS30": "#e87ba4", "IS3": "#008b00", "IS481": "#6a4fd8", "IS110": "#8b5a00",
    "IS256": "#00a3b4", "IS21": "#b0305c", "ISL3": "#7f7f00", "ISNCY": "#555555",
    "IS1": "#c75b00",
}
OTHER = "#aaa9a3"
CONTIG = "#2a78d6"
CONTIG2 = "#eb6834"
PLUS, MINUS = "#5cb85c", "#e8603c"


def clean(name):
    return re.sub(r"^gnl\|[^|]*\|", "", name)


def fam_colour(f):
    return FAMILY_COLOURS.get(f, OTHER)


def read_isescan(path):
    els = []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(hdr):
                continue
            els.append(dict(contig=clean(f[ix["seqID"]]), family=f[ix["family"]],
                            begin=int(f[ix["isBegin"]]), end=int(f[ix["isEnd"]]),
                            strand=f[ix["strand"]] if "strand" in ix else ""))
    return els


def read_annotation(path, fasta=None):
    """Return {contig: length} and list of (contig, start, end, strand) for genes/CDS."""
    lengths, genes = {}, []
    if re.search(r"\.(gbk|gb|gbff|genbank)$", path, re.I):
        from Bio import SeqIO
        for rec in SeqIO.parse(path, "genbank"):
            c = clean(rec.name if rec.name else rec.id)
            lengths[c] = len(rec.seq)
            for ft in rec.features:
                if ft.type == "CDS":
                    s = "+" if ft.location.strand == 1 else "-"
                    genes.append((c, int(ft.location.start), int(ft.location.end), s))
    else:
        header_len, embedded, cds = {}, {}, {}
        in_fasta, name = False, None
        with open(path) as fh:
            for line in fh:
                if in_fasta:
                    if line.startswith(">"):
                        name = clean(line[1:].split()[0])
                        embedded[name] = 0
                    elif name:
                        embedded[name] += len(line.strip())
                    continue
                if line.startswith("##sequence-region"):
                    # normally "##sequence-region seqid start end"; tolerate missing/odd fields
                    f = line.split()
                    nums = [x for x in f[2:] if x.isdigit()]
                    if len(f) >= 3 and nums:
                        header_len[clean(f[1])] = int(nums[-1])
                elif line.startswith("##FASTA"):
                    in_fasta = True
                elif not line.startswith("#") and line.strip():
                    f = line.rstrip("\n").split("\t")
                    if len(f) > 8 and f[2] == "CDS":
                        m = re.search(r"(?:^|;)ID=([^;]+)", f[8])
                        key = m.group(1) if m else (f[0], f[3], f[4], f[6])
                        c, s, e, st = clean(f[0]), int(f[3]), int(f[4]), f[6]
                        if key in cds:  # split CDS: one gene spanning all its parts
                            c0, s0, e0, st0 = cds[key]
                            cds[key] = (c0, min(s0, s), max(e0, e), st0)
                        else:
                            cds[key] = (c, s, e, st)
        genes = list(cds.values())
        fasta_len = {}
        if fasta:
            name = None
            with open(fasta) as fh:
                for line in fh:
                    if line.startswith(">"):
                        name = clean(line[1:].split()[0])
                        fasta_len[name] = 0
                    elif name:
                        fasta_len[name] += len(line.strip())
        # priority: separate FASTA > FASTA inside the GFF > ##sequence-region header
        lengths = dict(header_len)
        lengths.update(embedded)
        lengths.update(fasta_len)
        for c, L in header_len.items():
            if c in lengths and lengths[c] != L:
                print(f"WARNING: {c}: ##sequence-region says {L:,} bp but the sequence is "
                      f"{lengths[c]:,} bp; using {lengths[c]:,} bp", file=sys.stderr)
    if not lengths:
        sys.exit(f"No contig lengths found in {path}: add the genome FASTA as the 4th field of "
                 f"--genome (LABEL|ISESCAN_TSV|GFF|FASTA).")
    return lengths, genes


def element_id(c, b, e):
    return (c, b, e)


def parse_blast_id(s):
    m = re.match(r"^(.*)_(\d+)_(\d+)_([+-]?)$", clean(s))
    if not m:
        return None
    return (m.group(1), int(m.group(2)), int(m.group(3)))


def read_blast(path, mincov):
    best = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            q, s = parse_blast_id(f[0]), parse_blast_id(f[1])
            if not q or not s or q == s:
                continue
            pid = float(f[2])
            qs, qe, ss, se = (int(x) for x in f[4:8])
            qlen, slen = int(f[10]), int(f[11])
            covq = (abs(qe - qs) + 1) / qlen
            covs = (abs(se - ss) + 1) / slen
            key = tuple(sorted([q, s]))
            bits = float(f[9])
            if key not in best or bits > best[key][3]:
                best[key] = (pid, covq, covs, bits)
    return {k: v for k, v in best.items() if min(v[1], v[2]) >= mincov}


def order_labels(els, contig_order):
    """Label IS elements family-N in genome order."""
    rank = {c: i for i, c in enumerate(contig_order)}
    counter = defaultdict(int)
    for e in sorted(els, key=lambda e: (rank.get(e["contig"], 1e9), e["begin"])):
        counter[e["family"]] += 1
        e["label"] = f"{e['family']}-{counter[e['family']]}"


def gene_counts(genes, contig, length, window):
    nb = max(1, int(np.ceil(length / window)))
    plus, minus = np.zeros(nb), np.zeros(nb)
    for c, s, e, st in genes:
        if c != contig:
            continue
        i = min(nb - 1, ((s + e) // 2) // window)
        (plus if st == "+" else minus)[i] += 1
    x = (np.arange(nb) + 0.5) * window
    x[-1] = min(x[-1], length - 1)
    return x, plus, minus, window


def draw(ax, genomes, links, window, title, prime_second=False, grey_links=False):
    sectors, colours, display = {}, {}, {}
    for gi, g in enumerate(genomes):
        for ci, c in enumerate(g["contigs"]):
            key = f"{g['label']}::{c}"
            sectors[key] = g["lengths"][c]
            colours[key] = CONTIG if gi == 0 else CONTIG2
            display[key] = f"{ci + 1}" + ("′" if (prime_second and gi == 1) else "")
    total = sum(sectors.values())
    nsec = len(sectors)
    space = 0 if nsec == 1 else min(2.0, 60 / nsec)
    circos = Circos(sectors, space=space)
    vmax = 1
    for gg in genomes:
        for cc in gg["contigs"]:
            _, p_, m_, _ = gene_counts(gg["genes"], cc, gg["lengths"][cc], window)
            vmax = max(vmax, p_.max(), m_.max())
    bymap = {}
    for g in genomes:
        for e in g["elements"]:
            bymap[(e["contig"], e["begin"], e["end"])] = (g, e)
    for sector in circos.sectors:
        gl, contig = sector.name.split("::", 1)
        g = next(x for x in genomes if x["label"] == gl)
        L = sector.size
        ring = sector.add_track((78, 82))
        ring.axis(fc=colours[sector.name], ec="none")
        ticks = sector.add_track((70, 78))
        for e in g["elements"]:
            if e["contig"] != contig:
                continue
            mid = (e["begin"] + e["end"]) / 2
            w = max(e["end"] - e["begin"], total * 0.0012)
            ticks.rect(max(0, mid - w / 2), min(L, mid + w / 2), fc=fam_colour(e["family"]), ec="none")
            ticks.text(e["label"], mid, r=63, size=2.6, orientation="vertical", color="#3d3d3a")
        x, p, m, wdw = gene_counts(g["genes"], contig, L, window)
        tp = sector.add_track((84, 91))
        tm = sector.add_track((92, 99))
        tp.bar(x, p, width=wdw * 0.95, color=PLUS, vmax=vmax)
        tm.bar(x, m, width=wdw * 0.95, color=MINUS, vmax=vmax)
        if nsec > 1:
            sector.text(display[sector.name], r=106, size=5, color="#6b6a64")
    n = 0
    for (a, b), (pid, cq, cs, bits) in links.items():
        if (a[0], a[1], a[2]) not in bymap or (b[0], b[1], b[2]) not in bymap:
            continue
        ga, ea = bymap[(a[0], a[1], a[2])]
        gb, eb = bymap[(b[0], b[1], b[2])]
        if ga["label"] not in [g["label"] for g in genomes] or gb["label"] not in [g["label"] for g in genomes]:
            continue
        if grey_links and ga is gb:
            continue
        if not grey_links and ga is not gb:
            continue
        col = "#9a9994" if grey_links else fam_colour(ea["family"])
        circos.link((f"{ga['label']}::{ea['contig']}", ea["begin"], ea["end"]),
                    (f"{gb['label']}::{eb['contig']}", eb["begin"], eb["end"]),
                    r1=69, r2=69, color=col, alpha=0.55 if not grey_links else 0.45, lw=0.3)
        n += 1
    circos.plotfig(ax=ax)
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", action="append", required=True, help="LABEL|ISESCAN_TSV|ANNOTATION[|FASTA]")
    ap.add_argument("--blast", required=True)
    ap.add_argument("--pair", help="LABEL1,LABEL2 for the combined panel")
    ap.add_argument("--mincov", type=float, default=0.5,
                    help="links within a genome: minimum alignment coverage of both elements (0.5)")
    ap.add_argument("--pair-mincov", type=float, default=0.0,
                    help="grey links between the paired genomes: minimum coverage (0 = every BLASTn hit)")
    ap.add_argument("--window", type=int, default=20000, help="window for gene-count tracks (bp)")
    ap.add_argument("--out", default="Figure1_IS_circos")
    a = ap.parse_args()

    genomes = []
    for spec in a.genome:
        parts = spec.split("|")
        if len(parts) not in (3, 4):
            sys.exit(f"--genome needs LABEL|ISESCAN_TSV|ANNOTATION[|FASTA], got: {spec}")
        label, tsv, ann = parts[:3]
        lengths, genes = read_annotation(ann, parts[3] if len(parts) == 4 else None)
        contigs = sorted(lengths, key=lambda c: -lengths[c])
        els = read_isescan(tsv)
        missing = {e["contig"] for e in els} - set(lengths)
        if missing:
            sys.exit(f"{label}: ISEScan contigs not in annotation: {sorted(missing)[:5]}")
        beyond = [e for e in els if e["end"] > lengths[e["contig"]]]
        if beyond:
            sys.exit(f"{label}: {len(beyond)} IS element(s) end beyond their contig "
                     f"(e.g. {beyond[0]['contig']}:{beyond[0]['end']:,} > {lengths[beyond[0]['contig']]:,} bp). "
                     f"The ISEScan table and the annotation/FASTA come from different assemblies.")
        order_labels(els, contigs)
        genomes.append(dict(label=label, lengths=lengths, contigs=contigs, genes=genes, elements=els))
        print(f"{label}: {len(contigs)} contigs, {sum(lengths.values()):,} bp, {len(els)} IS elements, "
              f"{len(genes)} CDS", file=sys.stderr)

    links = read_blast(a.blast, a.mincov)
    links_between = read_blast(a.blast, a.pair_mincov)
    letters = "ABCDEFGH"
    pair = a.pair.split(",") if a.pair else None
    npan = len(genomes) + (1 if pair else 0)
    rows = int(np.ceil(len(genomes) / 2)) + (1 if pair else 0)
    fig = plt.figure(figsize=(7.2, 3.55 * rows + 0.6), dpi=300)
    gsp = fig.add_gridspec(rows + 1, 2, height_ratios=[1] * rows + [0.17], hspace=0.12, wspace=0.05)
    linkrows = []
    for i, g in enumerate(genomes):
        ax = fig.add_subplot(gsp[i // 2, i % 2], projection="polar")
        n = draw(ax, [g], links, a.window, g["label"])
        linkrows.append((letters[i], g["label"], n))
        ax.set_title("")
        ax.text(-0.02, 1.04, letters[i], transform=ax.transAxes, fontsize=13, fontweight="bold",
                va="bottom", ha="left", bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="black", lw=1.2))
        ital = g["label"] if g["label"] == "LK" else g["label"].replace(" sp.", "")
        ax.text(0.09, 1.045, ital + (" sp." if g["label"].endswith(" sp.") else ""), transform=ax.transAxes,
                fontsize=10, va="bottom", ha="left", style="normal" if g["label"] == "LK" else "italic")
    if pair:
        gp = [next(g for g in genomes if g["label"] == p) for p in pair]
        ax = fig.add_subplot(gsp[rows - 1, :], projection="polar")
        n = draw(ax, gp, links_between, a.window, "", prime_second=True, grey_links=True)
        k = letters[len(genomes)]
        linkrows.append((k, " + ".join(pair), n))
        n1, n2 = len(gp[0]["contigs"]), len(gp[1]["contigs"])
        ax.text(-0.02, 1.04, k, transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom",
                bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="black", lw=1.2))
        ax.text(0.06, 1.045, f"{pair[0]} (1–{n1}) and {pair[1]} (1′–{n2}′)", transform=ax.transAxes,
                fontsize=10, style="italic", va="bottom")
    # legend
    lax = fig.add_subplot(gsp[rows, :]); lax.axis("off")
    used = sorted({e["family"] for g in genomes for e in g["elements"]},
                  key=lambda f: list(FAMILY_COLOURS).index(f) if f in FAMILY_COLOURS else 99)
    handles = [Patch(fc=fam_colour(f), label=f) for f in used if f in FAMILY_COLOURS]
    if any(f not in FAMILY_COLOURS for f in used):
        handles.append(Patch(fc=OTHER, label="other families"))
    handles += [Patch(fc=PLUS, label="genes, + strand"), Patch(fc=MINUS, label="genes, − strand")]
    lax.legend(handles=handles, ncol=8, loc="center", frameon=False, fontsize=6.5,
               title="IS family (ticks and links)", title_fontsize=7, handlelength=2.2)

    for ext in ("png", "pdf"):
        fig.savefig(f"{a.out}.{ext}", dpi=300, bbox_inches="tight")
    fig.savefig(f"{a.out}.tif", dpi=300, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})

    # links table
    with open(f"{a.out}_links.tsv", "w") as out:
        out.write("element_1\telement_2\tidentity_%\tcoverage_1\tcoverage_2\tbitscore\n")
        for (x, y), (pid, cq, cs, bits) in sorted(links.items()):
            out.write(f"{x[0]}:{x[1]}-{x[2]}\t{y[0]}:{y[1]}-{y[2]}\t{pid:.1f}\t{cq:.2f}\t{cs:.2f}\t{bits}\n")
    for k, lab, n in linkrows:
        print(f"Panel {k} ({lab}): {n} links", file=sys.stderr)


if __name__ == "__main__":
    main()
