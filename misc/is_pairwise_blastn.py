#!/usr/bin/env python3
"""
is_pairwise_blastn.py - pairwise nucleotide identity and coverage between IS elements (replaces pyani ANIb).

ANIb cuts the query into 1,020-bp fragments, which is poorly defined for 1-2 kb IS elements. This script aligns
every IS element against every other with BLASTn and reports, for each pair, the identity of the best alignment and
the fraction of each element it covers. Input: the IS nucleotide sequences written by ISEScan (<genome>.is.fna).

Needs python3 and BLAST+ (makeblastdb, blastn).
Usage:
  python3 is_pairwise_blastn.py LK=isescan/LK_1.fna.is.fna Inquilinaceae=isescan/inqui.fna.is.fna \
      Oceanicaulis=isescan/ocean.fna.is.fna Phaeodactylibacter=isescan/phaeo.fna.is.fna \
      Nitratireductor=isescan/LK_bac_polished.fna.is.fna --out is_pairwise --threads 8
Outputs:
  <out>.tsv          one row per pair with an alignment (identity, alignment length, coverage of each element)
  <out>_summary.tsv  per family: pairs within and between organisms, how many are near-identical
  <out>_circos_links.txt  pairs with >= --link-id identity over >= --link-cov of both elements (Circos link format)
Options: --evalue 1e-10, --link-id 99, --link-cov 0.9
"""
import argparse, os, re, subprocess, sys, collections

def read_fasta(p):
    seqs, n, buf = [], None, []
    for l in open(p):
        if l.startswith('>'):
            if n: seqs.append((n, ''.join(buf)))
            n, buf = l[1:].strip(), []
        else: buf.append(l.strip())
    if n: seqs.append((n, ''.join(buf)))
    return seqs

def family(name):
    m = re.search(r'(IS\d+|ISL\d|ISAs1|ISKra4|IS200_IS605|IS200/IS605|IS1595|IS1380|IS1182|IS1634|IS481|IS66|IS630|IS256|Tn3|new)', name)
    return m.group(1) if m else ''

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('inputs', nargs='+', help='Organism=path/to/<genome>.is.fna')
    ap.add_argument('--out', default='is_pairwise'); ap.add_argument('--threads', default='4')
    ap.add_argument('--evalue', default='1e-10')
    ap.add_argument('--link-id', type=float, default=99.0); ap.add_argument('--link-cov', type=float, default=0.9)
    a = ap.parse_args()

    allseq, meta = [], {}
    for i, arg in enumerate(a.inputs):
        org, path = arg.split('=', 1)
        for j, (hdr, s) in enumerate(read_fasta(path)):
            sid = f'is{i}_{j}'
            name = re.sub(r'\s+', '_', hdr)
            meta[sid] = (name, org, family(hdr), len(s), hdr.split()[0])
            allseq.append((sid, s))
    fa = a.out + '_all.fna'
    with open(fa, 'w') as F:
        for sid, s in allseq: F.write(f'>{sid}\n{s}\n')
    print(f'{len(allseq)} IS elements from {len(a.inputs)} genomes', file=sys.stderr)
    subprocess.run(['makeblastdb', '-in', fa, '-dbtype', 'nucl', '-out', a.out + '_db'], check=True, stdout=subprocess.DEVNULL)
    res = subprocess.run(['blastn', '-task', 'blastn', '-query', fa, '-db', a.out + '_db', '-evalue', a.evalue,
                          '-num_threads', a.threads, '-max_target_seqs', '100000', '-dust', 'no',
                          '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore'],
                         check=True, capture_output=True, text=True).stdout
    best = {}
    for l in res.splitlines():
        q, s, pid, ln, qs, qe, ss, se, ev, bit = l.split('\t')
        if q == s: continue
        key = tuple(sorted((q, s)))
        bit = float(bit)
        if key not in best or bit > best[key][0]:
            # coverage of each element by this alignment
            if key[0] == q: c1, c2 = (abs(int(qe) - int(qs)) + 1), (abs(int(se) - int(ss)) + 1)
            else: c1, c2 = (abs(int(se) - int(ss)) + 1), (abs(int(qe) - int(qs)) + 1)
            best[key] = (bit, float(pid), int(ln), c1 / meta[key[0]][3], c2 / meta[key[1]][3], ev)
    rows = []
    for (x, y), (bit, pid, ln, cx, cy, ev) in best.items():
        rows.append((meta[x], meta[y], pid, ln, min(cx, 1.0), min(cy, 1.0), ev))
    rows.sort(key=lambda r: (-r[2] * min(r[4], r[5]), -r[2]))
    with open(a.out + '.tsv', 'w') as O:
        O.write('IS element 1\tOrganism 1\tFamily 1\tLength 1 (bp)\tIS element 2\tOrganism 2\tFamily 2\tLength 2 (bp)\t'
                'BLASTn identity (%)\tAlignment length (bp)\tCoverage of element 1\tCoverage of element 2\tE-value\n')
        for m1, m2, pid, ln, c1, c2, ev in rows:
            O.write(f'{m1[0]}\t{m1[1]}\t{m1[2]}\t{m1[3]}\t{m2[0]}\t{m2[1]}\t{m2[2]}\t{m2[3]}\t{pid:.2f}\t{ln}\t{c1:.3f}\t{c2:.3f}\t{ev}\n')
    # per family summary
    S = collections.defaultdict(lambda: [0, 0, 0, 0])
    for m1, m2, pid, ln, c1, c2, ev in rows:
        fam = m1[2] if m1[2] == m2[2] else f'{m1[2]}/{m2[2]}'
        same = m1[1] == m2[1]
        ident = pid >= a.link_id and min(c1, c2) >= a.link_cov
        S[fam][0 if same else 2] += 1
        if ident: S[fam][1 if same else 3] += 1
    with open(a.out + '_summary.tsv', 'w') as O:
        O.write(f'Family\tPairs within an organism\tof which near-identical (>= {a.link_id}% over >= {a.link_cov:.0%})\t'
                'Pairs between organisms\tof which near-identical\n')
        for fam in sorted(S): O.write(f'{fam}\t' + '\t'.join(map(str, S[fam])) + '\n')
    # circos links (contig and coordinates parsed from the ISEScan header: <seqID>_<begin>_<end>_<strand> ...)
    n = 0
    with open(a.out + '_circos_links.txt', 'w') as O:
        for m1, m2, pid, ln, c1, c2, ev in rows:
            if pid < a.link_id or min(c1, c2) < a.link_cov: continue
            p1 = re.match(r'(.+)_(\d+)_(\d+)_[+-]$', m1[4]); p2 = re.match(r'(.+)_(\d+)_(\d+)_[+-]$', m2[4])
            if p1 and p2:
                O.write(f'{p1.group(1)} {p1.group(2)} {p1.group(3)} {p2.group(1)} {p2.group(2)} {p2.group(3)} id={pid:.1f}\n'); n += 1
    print(f'{len(rows)} pairs with an alignment; {n} near-identical pairs written as Circos links', file=sys.stderr)

if __name__ == '__main__':
    main()
