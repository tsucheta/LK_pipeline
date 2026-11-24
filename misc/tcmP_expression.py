#!/usr/bin/env python3
"""
tcmP_expression.py - is tcmP (pgaptmp_004127) expressed out of line with its neighbours after the IS5 insertion?

Reads a featureCounts table (the standard output: a '#' command line, then Geneid Chr Start End Strand Length
<sample columns>) and, for every library:
  - TPM of every gene (counts / length, scaled to 1e6)
  - TPM of tcmP, the IS5 transposase in front of it (pgaptmp_004128) and the neighbouring genes
  - the percentile of tcmP among all genes with length >= 300 bp
  - tcmP / pgaptmp_004126 (radical SAM, the next gene downstream on the same strand)
Libraries known to contain genomic DNA (LK-7R1, LK-7R2) are left out by default (--exclude).

Optional --bam: for each RNA BAM, counts the reads that cross the IS5 -> tcmP junction (start in the IS5, continue
into tcmP) against the reads inside the first 300 bp of tcmP. Reads crossing the junction show transcription that
runs from the IS into tcmP, which is what a promoter inside the IS would produce. Needs samtools.

Usage:
  python3 tcmP_expression.py counts.txt
  python3 tcmP_expression.py counts.txt --bam rna_bams/*.bam --contig LNDBAELA --tcmp 4859836-4861230
(--tcmp: tcmP CDS in the coordinates of the assembly the BAMs were mapped to; default = uncorrected LK_1,
 use 4859963-4861357 for BAMs mapped to the corrected assembly)
"""
import argparse, os, re, subprocess, sys

FOCUS = ['pgaptmp_004124', 'pgaptmp_004125', 'pgaptmp_004126', 'pgaptmp_004127', 'pgaptmp_004128',
         'pgaptmp_004129', 'pgaptmp_004130', 'pgaptmp_004131']
NAMES = {'pgaptmp_004126': 'radical SAM (next gene, same strand)', 'pgaptmp_004127': 'tcmP',
         'pgaptmp_004128': 'IS5 transposase (inserted)', 'pgaptmp_004129': 'Uma2 endonuclease'}

def tag(gid):
    m = re.search(r'pgaptmp_\d+', gid)
    return m.group(0) if m else gid

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('counts')
    ap.add_argument('--exclude', default='LK-7R1,LK-7R2', help='comma-separated library name fragments to leave out')
    ap.add_argument('--bam', nargs='*', default=[])
    ap.add_argument('--contig', default='LNDBAELA')
    ap.add_argument('--tcmp', default='4859836-4861230')
    a = ap.parse_args()
    excl = [e for e in a.exclude.split(',') if e]

    rows, header = [], None
    for line in open(a.counts):
        if line.startswith('#'): continue
        c = line.rstrip('\n').split('\t')
        if header is None: header = c; continue
        rows.append(c)
    si = [i for i in range(6, len(header)) if not any(e in header[i] for e in excl)]
    left = [header[i] for i in range(6, len(header)) if i not in si]
    names = [os.path.basename(header[i]).replace('.bam', '').replace('.sorted', '') for i in si]
    genes = [(tag(r[0]), int(r[5]), [float(r[i]) for i in si]) for r in rows]
    if left: print('left out (genomic DNA):', ', '.join(os.path.basename(x) for x in left), file=sys.stderr)

    # TPM
    tpm = {}
    for j in range(len(si)):
        rpk = [(g, c[j] / (L / 1000)) for g, L, c in genes if L > 0]
        tot = sum(v for _, v in rpk) / 1e6
        for g, v in rpk: tpm.setdefault(g, []).append(v / tot if tot else 0)
    raw = {g: c for g, L, c in genes}
    length = {g: L for g, L, c in genes}

    print('gene\tlength\tdescription\t' + '\t'.join(f'{n} counts\t{n} TPM' for n in names))
    for g in FOCUS:
        if g not in tpm: continue
        print(f'{g}\t{length[g]}\t{NAMES.get(g, "")}\t' + '\t'.join(f'{raw[g][j]:.0f}\t{tpm[g][j]:.1f}' for j in range(len(si))))

    t, n = 'pgaptmp_004127', 'pgaptmp_004126'
    if t in tpm:
        print('\nlibrary\ttcmP TPM\ttcmP percentile (genes >= 300 bp)\ttcmP / pgaptmp_004126')
        long_genes = [g for g in tpm if length[g] >= 300]
        for j, nm in enumerate(names):
            v = tpm[t][j]
            pct = 100 * sum(1 for g in long_genes if tpm[g][j] < v) / len(long_genes)
            ratio = (v / tpm[n][j]) if n in tpm and tpm[n][j] > 0 else float('nan')
            print(f'{nm}\t{v:.1f}\t{pct:.0f}\t{ratio:.2f}')

    if a.bam:
        s, e = map(int, a.tcmp.split('-'))
        start = e                      # tcmP is on the minus strand: its 5' end (start codon) is the high coordinate
        print(f'\nIS5 -> tcmP read-through (tcmP start codon at {a.contig}:{start}; minus strand)')
        # reverse-stranded paired-end libraries (featureCounts -s 2): for a minus-strand gene, a sense fragment has
        # read 1 on the plus strand or read 2 on the minus strand. DNA contamination gives ~50% each way.
        print('BAM\tcrossing junction: sense\tcrossing junction: antisense\tfirst 300 bp of tcmP: sense\tfirst 300 bp: antisense\tsense share crossing')
        for b in a.bam:
            if any(x in os.path.basename(b) for x in excl): continue
            n = {('x', 's'): 0, ('x', 'a'): 0, ('i', 's'): 0, ('i', 'a'): 0}
            out = subprocess.run(['samtools', 'view', '-F', '0xD04', '-q', '10', b, f'{a.contig}:{start - 300}-{start + 300}'],
                                 capture_output=True, text=True, check=True).stdout
            for line in out.splitlines():
                c = line.split('\t'); flag = int(c[1]); pos = int(c[3]); cig = c[5]
                ref_len = sum(int(k) for k, op in re.findall(r'(\d+)([MDN=X])', cig))
                end = pos + ref_len - 1
                rev = bool(flag & 16); r1 = bool(flag & 64) or not (flag & 1)
                sense = (r1 and not rev) or (not r1 and rev)
                where = 'x' if (pos <= start - 20 and end >= start + 20) else ('i' if (pos >= start - 300 and end <= start) else None)
                if where: n[(where, 's' if sense else 'a')] += 1
            xs, xa = n[('x', 's')], n[('x', 'a')]
            share = f'{xs / (xs + xa):.2f}' if xs + xa else 'NA'
            print(f"{os.path.basename(b)}\t{xs}\t{xa}\t{n[('i', 's')]}\t{n[('i', 'a')]}\t{share}")

if __name__ == '__main__':
    main()
