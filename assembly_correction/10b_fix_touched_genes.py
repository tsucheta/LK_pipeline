#!/usr/bin/env python3
"""
10b_fix_touched_genes.py - re-predict the genes that a polishing correction touched, and put them into the lifted GFF.

10_lift_annotation.py shifts every gene and re-checks the reading frame of the genes a correction falls in. Genes whose
old model no longer fits the corrected sequence are marked "CHECK" (and user loci such as tmk are added as the longest
ORF). This script runs Prodigal (pyrodigal; trained on the whole corrected contig) and, for each such gene, takes the
predicted gene(s) overlapping it:
  - one predicted gene      -> the old locus tag gets the new coordinates
  - one predicted gene covering several old genes (e.g. fragments of one gene) -> the first tag keeps it, the others
    are removed and listed as merged into it
  - several predicted genes -> the first keeps the tag, the others get <tag>_b, <tag>_c …
Genes that merged or were repaired by 10_lift_annotation.py ("now one intact gene") are checked against Prodigal too.
Every new model is translated and must have no internal stop codon.

Needs python3 and pyrodigal (pip install pyrodigal).
Usage:
  python3 10b_fix_touched_genes.py --fna polish/LK_1_polished.fna --gff polish/LK_1_polished.gff \
      --genes polish/LK_1_polished_genes.tsv --out polish/LK_1_polished_fixed.gff
Outputs: <out> (GFF) and <out>.report.tsv (old model -> new model for every gene handled).
"""
import argparse, sys, re
import pyrodigal

def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith('>'):
            if name: seqs[name] = ''.join(buf).upper()
            name, buf = line[1:].split()[0], []
        else: buf.append(line.strip())
    if name: seqs[name] = ''.join(buf).upper()
    return seqs

def attrs(s): return dict(kv.split('=', 1) for kv in s.split(';') if '=' in kv)
def fmt(d): return ';'.join(f'{k}={v}' for k, v in d.items())

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fna', required=True); ap.add_argument('--gff', required=True); ap.add_argument('--genes', required=True)
    ap.add_argument('--out', required=True); ap.add_argument('--minov', type=float, default=0.3)
    a = ap.parse_args()
    G = read_fasta(a.fna)
    # targets from the lift report
    targets = []
    for line in open(a.genes):
        c = line.rstrip('\n').split('\t')
        if c[0] == 'locus' or len(c) < 10: continue
        res = c[9]
        if res.startswith('CHECK') or res.startswith('added') or res.startswith('now one intact gene'):
            s, e = map(int, c[8].split('-'))
            targets.append({'tags': c[1].split(','), 'start': s, 'end': e, 'strand': c[5], 'why': res.split(':')[0].split(' (')[0]})
    if not targets: sys.exit('no genes to fix in ' + a.genes)
    lines = [l.rstrip('\n') for l in open(a.gff)]
    contig = next(l.split('\t')[0] for l in lines if not l.startswith('#') and l.strip())
    seq = G[contig] if contig in G else next(iter(G.values()))
    print(f'training Prodigal on {contig} ({len(seq):,} bp)', file=sys.stderr)
    gf = pyrodigal.GeneFinder(meta=False, closed=False)
    gf.train(seq.encode())
    pred = [(g.begin, g.end, '+' if g.strand == 1 else '-', g.translate(translation_table=11).rstrip('*'), g.partial_begin or g.partial_end)
            for g in gf.find_genes(seq.encode())]
    # assign predicted genes to targets: best overlap first; named loci (e.g. tmk) before repaired genes before CHECK genes
    report, replace, remove, add, rename = [], {}, set(), [], {}
    prio = {'added': 0, 'now one intact gene': 1, 'CHECK': 2}
    cand = []
    for ti, T in enumerate(targets):
        tl = T['end'] - T['start'] + 1
        for k, (b_, e_, st, prot, partial) in enumerate(pred):
            ov = min(e_, T['end']) - max(b_, T['start']) + 1
            if ov <= 0 or (T['strand'] in '+-' and st != T['strand']): continue
            if ov / (e_ - b_ + 1) >= a.minov or ov / tl >= a.minov:
                cand.append((prio.get(T['why'], 3), -ov / max(tl, e_ - b_ + 1), ti, k))
    owner, got = {}, {}
    for _, _, ti, k in sorted(cand):
        if k in owner or ti in got: continue
        owner[k] = ti; got[ti] = k
    extra = {}
    for _, _, ti, k in sorted(cand):          # predicted genes left over: extra genes of the target they overlap most
        if k not in owner: owner[k] = ti; extra.setdefault(ti, []).append(k)
    for ti, T in enumerate(targets):
        tag0 = T['tags'][0]; span = f"{T['start']}-{T['end']}"
        if ti in got:
            k = got[ti]; b_, e_, st, prot, partial = pred[k]
            replace[tag0] = (b_, e_, st, len(prot)); remove.update(T['tags'][1:])
            report.append([','.join(T['tags']), T['why'], span, f'{tag0}: {b_}-{e_} ({st})', f'{len(prot)} aa' + (' (partial)' if partial else ''),
                           'internal stop!' if '*' in prot else 'ok'])
            for j, k2 in enumerate(sorted(extra.get(ti, []), key=lambda k: pred[k][0])):
                b2, e2, st2, prot2, partial2 = pred[k2]; tag = f'{tag0}_{"bcdefgh"[j]}'
                add.append((tag, b2, e2, st2, len(prot2)))
                report.append([','.join(T['tags']), T['why'], span, f'{tag}: {b2}-{e2} ({st2}) (additional gene)', f'{len(prot2)} aa',
                               'internal stop!' if '*' in prot2 else 'ok'])
        else:
            ks = [k for _, _, t2, k in sorted(cand, key=lambda x: x[1]) if t2 == ti]   # best overlap first
            if ks:
                k = ks[0]; O = targets[owner[k]]; into = O['tags'][0]
                remove.update(T['tags'])
                if O['why'] == 'added' and into not in rename:   # a named locus (e.g. tmk) takes over the PGAP tag it replaces
                    rename[into] = T['tags'][0]
                    report.append([','.join(T['tags']), T['why'], span, f'replaced by {into}, which keeps the locus tag {T["tags"][0]} ({pred[k][0]}-{pred[k][1]})', '', ''])
                else:
                    report.append([','.join(T['tags']), T['why'], span, f'merged into {into} ({pred[k][0]}-{pred[k][1]})', '', ''])
            else:
                report.append([','.join(T['tags']), T['why'], span, 'no Prodigal gene here; old model kept', '', 'check'])
    # rewrite the GFF
    was_pseudo = set()
    for l in lines:
        c = l.split('\t')
        if len(c) >= 9 and not l.startswith('#'):
            A = attrs(c[8])
            if c[2] == 'pseudogene' or A.get('pseudo', '').lower() == 'true' or A.get('gene_biotype') == 'pseudogene': was_pseudo.add(A.get('locus_tag'))
    cds_done = set()
    out = []
    for l in lines:
        if l.startswith('#') or not l.strip(): out.append(l); continue
        c = l.split('\t'); A = attrs(c[8]) if len(c) >= 9 else {}
        lt = A.get('locus_tag')
        if lt in remove: continue
        if lt in replace and c[2] in ('gene', 'CDS', 'pseudogene'):
            b, e, st, L = replace[lt]
            if c[2] == 'CDS':                       # a split pseudogene CDS becomes one CDS line
                if lt in cds_done: continue
                cds_done.add(lt); c[7] = '0'
            c[3], c[4], c[6] = str(b), str(e), st
            if c[2] == 'pseudogene': c[2] = 'gene'
            for k in ('pseudo', 'pseudogene', 'partial', 'start_range', 'end_range', 'exception'): A.pop(k, None)
            if A.get('gene_biotype') == 'pseudogene': A['gene_biotype'] = 'protein_coding'
            A['Note'] = 're-predicted with Prodigal on the corrected assembly (intact; previously a pseudogene in the long-read assembly)' \
                if lt in was_pseudo else 're-predicted with Prodigal on the corrected assembly'
            if lt in rename:                        # e.g. tmk_dTMP_kinase -> locus_tag pgaptmp_003484, gene=tmk, product=dTMP kinase
                nt = rename[lt]; gname, _, prod = lt.partition('_')
                A['ID'] = ('cds-' if c[2] == 'CDS' else 'gene-') + nt; A['locus_tag'] = nt; A['gene'] = gname
                if 'Parent' in A: A['Parent'] = 'gene-' + nt
                if c[2] == 'CDS': A['product'] = prod.replace('_', ' ') or 'hypothetical protein'
                else: A.pop('Name', None); A['gene_biotype'] = 'protein_coding'
                A['Note'] += f', {gname} not annotated by PGAP in the long-read assembly, where {nt} was a fragment model of this locus'
            c[8] = fmt(A); l = '\t'.join(c)
        out.append(l)
    for tag, b, e, st, L in add:
        out.append(f'{contig}\tProdigal\tgene\t{b}\t{e}\t.\t{st}\t.\tID=gene-{tag};locus_tag={tag};Note=new gene predicted with Prodigal on the corrected assembly')
        out.append(f'{contig}\tProdigal\tCDS\t{b}\t{e}\t.\t{st}\t0\tID=cds-{tag};Parent=gene-{tag};locus_tag={tag};product=hypothetical protein')
    open(a.out, 'w').write('\n'.join(out) + '\n')
    with open(a.out + '.report.tsv', 'w') as R:
        R.write('old locus tag(s)\treason\tlifted span\tnew model\tprotein\tcheck\n')
        for r in report: R.write('\t'.join(r) + '\n')
    for r in report: print('\t'.join(r), file=sys.stderr)

if __name__ == '__main__':
    main()
