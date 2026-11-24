#!/usr/bin/env python3
"""Rename the contigs in a Prokka GFF so they match the FASTA used for mapping.

Prokka renames contigs (e.g. k141_123 -> gnl|X|BOIBFNEL_1) and drops contigs < 200 bp, so names cannot be
matched directly. This pairs the GFF ##sequence-region entries with the FASTA (.fai) contigs IN ORDER,
skipping FASTA contigs whose length does not match, and refuses to write anything unless every GFF contig
finds a partner of identical length.

Usage: python3 fix_gff_names.py refs/inqui.gff refs/inqui.fna.fai refs/inqui.fixed.gff
Then set refs/inqui.fixed.gff in config.sh (or move it over refs/inqui.gff after keeping a copy).
"""
import sys

def main(gff, fai, out):
    regions = []          # (name, length) in GFF order
    for l in open(gff):
        if l.startswith("##sequence-region"):
            _, name, s, e = l.split()[:4]
            regions.append((name, int(e) - int(s) + 1))
        elif l.startswith("##FASTA"):
            break
    if not regions:
        sys.exit("ERROR: no ##sequence-region lines in the GFF; cannot pair contigs by length")
    fa = [(l.split("\t")[0], int(l.split("\t")[1])) for l in open(fai) if l.strip()]
    fa_names = {n for n, _ in fa}
    if all(n in fa_names for n, _ in regions):
        sys.exit("Nothing to do: every GFF contig name is already in the FASTA. The problem is elsewhere.")
    ren, j = {}, 0
    for name, ln in regions:
        while j < len(fa) and fa[j][1] != ln:
            j += 1
        if j == len(fa):
            sys.exit(f"ERROR: no FASTA contig of length {ln} left for GFF contig {name}; the GFF was not made from this FASTA. Nothing written.")
        ren[name] = fa[j][0]; j += 1
    n = 0
    with open(out, "w") as o:
        in_fasta = False
        for l in open(gff):
            if l.startswith("##FASTA"):
                in_fasta = True
            if in_fasta:
                if l.startswith(">") and l[1:].split()[0] in ren:
                    l = ">" + ren[l[1:].split()[0]] + "\n"
                o.write(l); continue
            if l.startswith("##sequence-region"):
                p = l.split(); p[1] = ren.get(p[1], p[1]); l = " ".join(p) + "\n"
            elif not l.startswith("#") and "\t" in l:
                f = l.split("\t", 1)
                if f[0] in ren: l = ren[f[0]] + "\t" + f[1]; n += 1
            o.write(l)
    print(f"{len(ren)} contigs renamed ({n} feature lines), e.g. {next(iter(ren))} -> {next(iter(ren.values()))}. Written: {out}")

if __name__ == "__main__":
    if len(sys.argv) != 4: sys.exit(__doc__)
    main(*sys.argv[1:])
