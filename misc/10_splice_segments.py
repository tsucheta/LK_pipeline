#!/usr/bin/env python3
"""
10_splice_segments.py - replace misassembled segments of LK_1 with the homologous sequence of an older assembly.

Pilon corrects bases and short indels but cannot rebuild a segment of ~1 kb that the long-read assembly got wrong
(step 09 reports it as "N bp not covered by reads"). For each window given with --locus, the LK_1 window (plus
--pad bp each side) is aligned to the donor assembly (e.g. LK_16, 2014): the window is located in the donor by
unique flanking anchors (either strand), the two sequences are walked in parallel, and every difference block of
>= --minblock bp is replaced by the donor sequence. Shorter differences (single bases, small indels) are NOT
imported - they may be real changes since the donor year - and are left to Pilon and the Illumina reads.
The result must be checked with the 2023 reads (10_run_polish.sh does this with step 09).

Needs python3 (standard library) and 10_lift_annotation.py in the same folder.
Usage:
  python3 10_splice_segments.py --asm refs/LK_1.fna --contig LNDBAELA --donor refs/LK_16.fna \
      --locus LNDBAELA:1225361-1228552 LNDBAELA:4131300-4132750 --out polish/LK_1_spliced.fna
Outputs: the spliced FASTA (one contig) and <out>.splices.tsv (what was replaced, with what).
"""
import argparse, importlib.util, os, sys
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('lift', os.path.join(here, '10_lift_annotation.py'))
lift = importlib.util.module_from_spec(spec); spec.loader.exec_module(lift)

def read_all(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name: seqs[name] = ''.join(buf).upper()
                name, buf = line[1:].split()[0], []
            else: buf.append(line.strip())
    if name: seqs[name] = ''.join(buf).upper()
    return seqs

def locate(donor, left, right):
    """find left and right anchors (40 bp) on the same donor contig and strand; returns (contig, seq, a, b)
    with seq oriented like LK_1 and [a, b) spanning anchor to anchor"""
    hits = []
    for c, s in donor.items():
        for strand, t in (('+', s), ('-', lift.rc(s))):
            if t.count(left) == 1 and t.count(right) == 1:
                a, b = t.index(left), t.index(right) + len(right)
                if 0 <= b - a < 50000: hits.append((c, strand, t, a, b))
    return hits

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--asm', required=True); ap.add_argument('--contig', required=True); ap.add_argument('--donor', required=True)
    ap.add_argument('--locus', nargs='+', required=True); ap.add_argument('--out', required=True)
    ap.add_argument('--pad', type=int, default=500); ap.add_argument('--minblock', type=int, default=50)
    a = ap.parse_args()
    seq = read_all(a.asm)[a.contig] if a.contig in read_all(a.asm) else sys.exit(f'ERROR: {a.contig} not in {a.asm}')
    donor = read_all(a.donor)
    edits, rows = [], []
    for L in a.locus:
        c, rng = L.rsplit(':', 1)
        if c != a.contig: continue
        s, e = map(int, rng.replace(',', '').split('-'))
        # anchors just outside the padded window, shifted outwards until unique in LK_1 and found in the donor
        found = None
        for shift in range(0, 2000, 50):
            w0, w1 = max(0, s - 1 - a.pad - shift), min(len(seq), e + a.pad + shift)
            left, right = seq[w0:w0 + 40], seq[w1 - 40:w1]
            if seq.count(left) != 1 or seq.count(right) != 1: continue
            h = locate(donor, left, right)
            if len(h) == 1: found = (w0, w1, h[0]); break
        if not found:
            rows.append([L, '-', '-', '-', '-', 'not spliced: window not found once in the donor']); continue
        w0, w1, (dc, strand, t, da, db) = found
        mine, theirs = seq[w0:w1], t[da:db]
        for o0, o1, n0, n1 in lift.align(mine, theirs):
            if max(o1 - o0, n1 - n0) < a.minblock: continue
            edits.append((w0 + o0, w0 + o1, theirs[n0:n1]))
            rows.append([L, f'{w0 + o0 + 1}-{w0 + o1}', str(o1 - o0), f'{dc}:{strand}', str(n1 - n0),
                         f'replaced {o1 - o0} bp of {a.contig} with {n1 - n0} bp of the donor'])
        if not any(r[0] == L and r[5].startswith('replaced') for r in rows):
            rows.append([L, '-', '-', f'{dc}:{strand}', '-', f'no difference >= {a.minblock} bp from the donor; nothing replaced'])
    for o0, o1, new in sorted(edits, reverse=True): seq = seq[:o0] + new + seq[o1:]
    with open(a.out, 'w') as O:
        O.write(f'>{a.contig}\n' + '\n'.join(seq[i:i + 80] for i in range(0, len(seq), 80)) + '\n')
    with open(a.out + '.splices.tsv', 'w') as T:
        T.write('locus\treplaced_LK_1_segment\tLK_1_bp\tdonor_contig:strand\tdonor_bp\tresult\n')
        for r in rows: T.write('\t'.join(r) + '\n')
    for r in rows: print('\t'.join(r), file=sys.stderr)

if __name__ == '__main__':
    main()
