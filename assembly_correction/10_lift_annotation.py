#!/usr/bin/env python3
"""
10_lift_annotation.py - carry the LK_1 annotation over to the polished assembly and convert coordinates.

  1. Aligns the old and the polished contig (they differ only at the corrected bases and segments) and writes the
     coordinate map: <out>_map.tsv, one row per changed block (old start-end -> new start-end, old/new sequence).
  2. Lifts every GFF feature. A gene that no change touches is only shifted. A gene that a change touches has its
     reading frame re-checked on the polished sequence:
       - intact gene, frame still runs from its start codon to a stop at its (lifted) end  -> shifted, "frame ok"
       - pseudogene / split gene (fragments of one gene merged as in step 09) whose start codon now reads through
         in one open frame over >= 90% of the fragments' span                          -> one intact gene
         (locus tag of the first fragment, pseudo removed, Note lists the merged tags)
       - anything else                                                                  -> shifted, flagged "check"
     Loci given with --locus (genes missing from the annotation, e.g. tmk) are searched for their longest ORF on
     the polished sequence and added as new genes when it is >= --minorf bp.
  3. --convert FILE converts any table of positions (columns: contig, position; header allowed, extra columns kept)
     to polished coordinates; a position inside a replaced segment is reported as "in changed block".

Needs python3 (standard library only).
Usage:
  python3 10_lift_annotation.py --old refs/LK_1.fna --new polish/LK_1_polished.fna --gff refs/LK_1.gff \
      --contig LNDBAELA --locus tmk_dTMP_kinase:LNDBAELA:4131300-4132750 --out polish/LK_1_polished
  python3 10_lift_annotation.py --map polish/LK_1_polished_map.tsv --convert positions.tsv > positions_polished.tsv
Outputs: <out>.gff, <out>_map.tsv, <out>_genes.tsv (every gene a change touched, and what was done)
"""
import argparse, sys
from urllib.parse import unquote, quote

def read_contig(path, contig):
    seq, on, names = [], False, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if on: break
                name = line[1:].split()[0]; names.append(name)
                on = name == contig or (contig + '_pilon') == name
            elif on: seq.append(line.strip())
    if not seq:
        if len(names) == 1:  # a one-contig FASTA under another name (e.g. LNDBAELA_pilon) is accepted
            return read_contig(path, names[0])
        sys.exit(f'ERROR: contig {contig} not found in {path}')
    return ''.join(seq).upper()

# ------------------------------------------------------------------ 1. old -> new alignment (blocks)
def align(old, new, k=40, maxgap=6000):
    """walk both sequences; at a mismatch find the nearest point where they agree again over k bp.
    Returns change blocks (o0, o1, n0, n1) as 0-based half-open intervals."""
    i = j = 0; blocks = []
    no, nn = len(old), len(new)
    while i < no and j < nn:
        # fast skip over identical stretches
        step = 1 << 16
        while step >= 1:
            if old[i:i + step] == new[j:j + step] and i + step <= no and j + step <= nn: i += step; j += step
            else: step >>= 1
        if i >= no or j >= nn: break
        best = None
        for d in range(0, maxgap):
            if best and d >= best[0]: break
            probe = old[i + d:i + d + k]
            if len(probe) < k: break
            lim = maxgap if best is None else best[0] - d
            p = new.find(probe, j, j + lim + k)
            if p >= 0:
                cost = d + (p - j)
                if best is None or cost < best[0]: best = (cost, d, p - j)
        if best is None:
            sys.exit(f'ERROR: assemblies do not re-synchronise after old position {i + 1}; is --new a polished version of --old?')
        _, d, e = best
        blocks.append((i, i + d, j, j + e)); i += d; j += e
    if i < no or j < nn: blocks.append((i, no, j, nn))
    return [b for b in blocks if b[1] > b[0] or b[3] > b[2]]

def make_lift(blocks):
    """returns lift(p, side) for a 1-based old position; side 'L'/'R' decides where a position inside a
    changed block goes (start -> left edge, end -> right edge); None if inside and side is None"""
    import bisect
    starts = [b[0] for b in blocks]
    def lift(p, side=None):
        x = p - 1
        k = bisect.bisect_right(starts, x) - 1
        off = 0
        if k >= 0:
            o0, o1, n0, n1 = blocks[k]
            if o0 <= x < o1:
                if o1 - o0 == n1 - n0: return n0 + (x - o0) + 1         # substitution: same position
                if side is None: return None
                return n0 + 1 if side == 'L' else max(n0 + 1, n1)
            off = n1 - o1
        return x + off + 1
    return lift

# ------------------------------------------------------------------ 2. GFF
COMP = str.maketrans('ACGTN', 'TGCAN')
STARTS = {'ATG', 'GTG', 'TTG'}; STOPS = {'TAA', 'TAG', 'TGA'}
def rc(s): return s.translate(COMP)[::-1]
def orf_from(seq, start, strand, limit):
    """1-based start of a start codon (+: first base, -: last base); returns the 1-based stop end coordinate of
    the open frame (including the stop codon) or None if no stop within limit bp"""
    if strand == '+':
        for i in range(start - 1, min(len(seq) - 2, start - 1 + limit), 3):
            if seq[i:i + 3] in STOPS: return i + 3
    else:
        for i in range(start, max(2, start - limit), -3):
            if rc(seq[i - 3:i]) in STOPS: return i - 2
    return None

def longest_orf(seq, a, b, minlen):
    best = None
    for strand in '+-':
        s = seq[a - 1:b] if strand == '+' else rc(seq[a - 1:b])
        for fr in range(3):
            st = None
            for i in range(fr, len(s) - 2, 3):
                c = s[i:i + 3]
                if st is None and c in STARTS: st = i
                elif st is not None and c in STOPS:
                    L = i + 3 - st
                    if L >= minlen and (best is None or L > best[0]):
                        g0, g1 = (a + st, a + i + 2) if strand == '+' else (b - i - 2, b - st)
                        best = (L, g0, g1, strand)
                    st = None
    return best

def parse_attr(s): return dict(kv.split('=', 1) for kv in s.split(';') if '=' in kv)
def fmt_attr(d): return ';'.join(f'{k}={v}' for k, v in d.items())

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--old'); ap.add_argument('--new'); ap.add_argument('--gff'); ap.add_argument('--contig', default='LNDBAELA')
    ap.add_argument('--locus', nargs='*', default=[]); ap.add_argument('--out'); ap.add_argument('--minorf', type=int, default=300)
    ap.add_argument('--map', help='existing <out>_map.tsv (for --convert only)'); ap.add_argument('--convert')
    a = ap.parse_args()

    if a.map:
        blocks = []
        for line in open(a.map):
            if line.startswith('#') or line.startswith('old_start'): continue
            c = line.rstrip('\n').split('\t'); blocks.append((int(c[0]) - 1, int(c[1]), int(c[2]) - 1, int(c[3])))
        contig = open(a.map).readline().lstrip('#').split()[0]
    else:
        old, new = read_contig(a.old, a.contig), read_contig(a.new, a.contig)
        blocks = align(old, new); contig = a.contig
        with open(a.out + '_map.tsv', 'w') as M:
            M.write(f'#{contig} old length {len(old)} -> polished length {len(new)}; {len(blocks)} changed blocks; coordinates 1-based inclusive\n')
            M.write('old_start\told_end\tnew_start\tnew_end\told_seq\tnew_seq\tlength_change\tcumulative_shift_after\n')
            cum = 0
            for o0, o1, n0, n1 in blocks:
                cum += (n1 - n0) - (o1 - o0)
                M.write(f'{o0 + 1}\t{o1}\t{n0 + 1}\t{n1}\t{old[o0:o1][:60] or "-"}{"..." if o1 - o0 > 60 else ""}\t'
                        f'{new[n0:n1][:60] or "-"}{"..." if n1 - n0 > 60 else ""}\t{(n1 - n0) - (o1 - o0)}\t{cum}\n')
        print(f'{len(blocks)} changed blocks; length {len(old)} -> {len(new)}', file=sys.stderr)
    lift = make_lift(blocks)

    if a.convert:
        out = sys.stdout
        for line in open(a.convert):
            c = line.rstrip('\n').split('\t')
            if len(c) < 2 or not c[1].replace(',', '').isdigit():
                out.write(line.rstrip('\n') + '\tpolished_position\n'); continue
            if c[0] != contig: out.write(line.rstrip('\n') + '\tother contig (unchanged)\n'); continue
            p = lift(int(c[1].replace(',', '')))
            out.write(line.rstrip('\n') + '\t' + (str(p) if p else 'in changed block') + '\n')
        return
    if not a.gff: return

    # ---- read the GFF, group fragments of split genes (as in step 09)
    lines = [l.rstrip('\n') for l in open(a.gff)]
    feats = []
    for idx, l in enumerate(lines):
        if l.startswith('#') or not l.strip(): continue
        c = l.split('\t')
        if len(c) < 9 or c[0] != contig: continue
        feats.append((idx, c))
    def touched(s, e): return [b for b in blocks if b[0] < e and b[1] >= s - 1]
    genes = {}
    for idx, c in feats:
        A = parse_attr(c[8]); lt = A.get('locus_tag')
        if not lt or c[2] not in ('gene', 'pseudogene', 'CDS'): continue
        g = genes.setdefault(lt, {'s': int(c[3]), 'e': int(c[4]), 'strand': c[6], 'pseudo': False, 'product': '', 'rows': []})
        g['s'] = min(g['s'], int(c[3])); g['e'] = max(g['e'], int(c[4])); g['rows'].append(idx)
        if A.get('pseudo', '').lower() == 'true' or c[2] == 'pseudogene': g['pseudo'] = True
        if A.get('product'): g['product'] = unquote(A['product'])
    order = sorted(genes, key=lambda t: genes[t]['s'])
    groups = []
    for lt in order:
        g = genes[lt]
        if g['pseudo'] and groups and groups[-1]['pseudo'] and groups[-1]['strand'] == g['strand'] \
           and g['s'] - groups[-1]['e'] < 2000 and groups[-1]['product'] == g['product']:
            groups[-1]['tags'].append(lt); groups[-1]['e'] = max(groups[-1]['e'], g['e'])
        else:
            groups.append({'tags': [lt], 's': g['s'], 'e': g['e'], 'strand': g['strand'], 'pseudo': g['pseudo'], 'product': g['product']})

    new = read_contig(a.new, contig)
    newrow = {}            # line index -> replacement line ('' = drop)
    report = []
    for G in groups:
        T = touched(G['s'], G['e'])
        if not T: continue
        s1, e1 = lift(G['s'], 'L'), lift(G['e'], 'R')
        tags = G['tags']; name = tags[0] if len(tags) == 1 else tags[0] + '..' + tags[-1]
        st = G['strand']; span = e1 - s1 + 1
        start = s1 if st == '+' else e1
        codon = new[start - 1:start + 2] if st == '+' else rc(new[start - 3:start])
        stop = orf_from(new, start, st, span + 3000) if codon in STARTS else None
        L = (stop - start + 1) if (stop and st == '+') else ((start - stop + 1) if stop else 0)
        what = ''
        if not G['pseudo']:
            endpos = e1 if st == '+' else s1
            what = 'frame ok (shifted)' if stop == endpos else f'CHECK: frame from the start codon now ends at {stop} (annotation end {endpos})'
        elif stop and L >= 0.9 * span:
            g0, g1 = (start, stop) if st == '+' else (stop, start)
            what = f'now one intact gene {g0}-{g1} ({L} bp, {L // 3 - 1} aa)'
            # rewrite: keep the first fragment's gene and CDS lines with new coordinates, drop the others
            first = genes[tags[0]]['rows']
            cds_seen = False
            for lt in tags:
                for r in genes[lt]['rows']:
                    if lt != tags[0]: newrow[r] = ''; continue
                    c = lines[r].split('\t'); A = parse_attr(c[8])
                    if c[2] == 'CDS':                 # a CDS split over several lines becomes one line
                        if cds_seen: newrow[r] = ''; continue
                        cds_seen = True
                    for k in ('pseudo', 'pseudogene', 'partial', 'start_range', 'end_range', 'exception'): A.pop(k, None)
                    if A.get('gene_biotype') == 'pseudogene': A['gene_biotype'] = 'protein_coding'
                    note = f'merged from {",".join(tags)} after Illumina polishing' if len(tags) > 1 else 'frameshift corrected by Illumina polishing'
                    A['Note'] = quote(note, safe=' ,')
                    if c[2] == 'pseudogene': c[2] = 'gene'
                    c[3], c[4] = str(g0), str(g1); c[7] = '0' if c[2] == 'CDS' else c[7]
                    c[8] = fmt_attr(A); newrow[r] = '\t'.join(c)
            # a pseudogene annotated without a CDS line gets one
            if not any(lines[r].split('\t')[2] == 'CDS' for r in first):
                gl = lines[first[0]].split('\t'); A = parse_attr(gl[8])
                cds = gl[:]; cds[2] = 'CDS'; cds[3], cds[4] = str(g0), str(g1); cds[7] = '0'
                B = {k: v for k, v in A.items() if k in ('locus_tag', 'product', 'Name', 'gene')}
                B['ID'] = 'cds-' + tags[0]; B['Parent'] = A.get('ID', 'gene-' + tags[0]); B['Note'] = A.get('Note', '')
                if A.get('product') is None and G['product']: B['product'] = quote(G['product'], safe=' ,')
                newrow[first[-1]] = (newrow.get(first[-1]) or lines[first[-1]]) + '\n' + '\t'.join(cds)
        else:
            what = 'still disrupted after polishing (shifted)' + (f'; open frame from the start covers {L} of {span} bp' if L else '')
        report.append([name, ','.join(tags), G['product'], str(G['s']), str(G['e']), st, 'pseudo' if G['pseudo'] else 'gene',
                       str(len(T)), f'{s1}-{e1}', what])

    added = []
    for Lc in a.locus:
        nm, c0, rng = Lc.split(':'); s, e = map(int, rng.replace(',', '').split('-'))
        if c0 != contig: continue
        s1, e1 = lift(s, 'L'), lift(e, 'R')
        o = longest_orf(new, s1, e1, a.minorf)
        if o:
            Lb, g0, g1, strand = o
            # the longest ORF runs from the most upstream start codon; the true start may lie further in
            # (NCBI PGAP re-calls start codons at deposition)
            attr = f'ID=gene-{nm};locus_tag={nm};Note=longest ORF in the corrected locus after Illumina polishing (not annotated in LK_1; verify the product)'
            added.append(f'{contig}\tpolish\tgene\t{g0}\t{g1}\t.\t{strand}\t.\t{attr}')
            added.append(f'{contig}\tpolish\tCDS\t{g0}\t{g1}\t.\t{strand}\t0\tID=cds-{nm};Parent=gene-{nm};locus_tag={nm}')
            what = f'added: longest ORF {g0}-{g1} ({Lb} bp, {Lb // 3 - 1} aa, {strand} strand)'
        else: what = f'no ORF >= {a.minorf} bp in {s1}-{e1}'
        report.append([nm, nm, 'user locus', str(s), str(e), '.', 'locus', str(len(touched(s, e))), f'{s1}-{e1}', what])

    with open(a.out + '.gff', 'w') as O:
        for idx, l in enumerate(lines):
            if l.startswith('##sequence-region') and contig in l:
                p = l.split(); O.write(f'{p[0]} {p[1]} 1 {len(new)}\n'); continue
            if idx in newrow:
                if newrow[idx]: O.write(newrow[idx] + '\n')
                continue
            c = l.split('\t')
            if len(c) >= 9 and c[0] == contig and not l.startswith('#'):
                c[3], c[4] = str(lift(int(c[3]), 'L')), str(lift(int(c[4]), 'R')); l = '\t'.join(c)
            O.write(l + '\n')
        for l in added: O.write(l + '\n')
    with open(a.out + '_genes.tsv', 'w') as R:
        R.write('locus\tlocus_tags\tproduct\told_start\told_end\tstrand\ttype\tchanges_in_gene\tpolished_span\tresult\n')
        for r in report: R.write('\t'.join(r) + '\n')
    for r in report: print(f'{r[0]}\t{r[2]}\t{r[9]}', file=sys.stderr)

if __name__ == '__main__':
    main()
