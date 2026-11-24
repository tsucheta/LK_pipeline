#!/usr/bin/env python3
"""
09_check_split_genes.py - are the split/pseudogenes of the long-read LK_1 assembly real, or assembly errors?

For every PGAP pseudogene (pseudo=true) in the GFF, grouped into loci (fragments of one gene merged), plus any
loci given with --locus, the Illumina DNA reads of each dataset are piled up on the assembly and three tests run:

  1. Base disagreement: a position where most Illumina reads (>= --major, default 0.5) carry a different base or
     an indel than the assembly. If such a change lies in the locus the assembly is wrong there; applying it is
     tested for whether it restores the open reading frame.
  2. Coverage gap: a run of >= 50 bp with depth < 0.2x the genome median (a stretch of assembly that the reads
     do not support, e.g. a misassembled segment).
  3. Clipping pile-up: >= max(5, 0.2 x depth) reads soft-clipped (>= 20 bp) at the same position (a breakpoint).

Verdict per locus
  assembly error       most 2023 Illumina reads disagree with the assembly in the locus (base or indel), or
                       the locus has a coverage gap or a majority clipping pile-up. An error of the 2023 long-read
                       assembly is normally also seen in the 2014 and 2017 reads.
  split supported      the 2023 reads agree with the assembly, depth normal. If the older reads disagree, the
                       split is a recent mutation (it should then also be in the step-08 variant table).
  minority alleles     reads agree in the majority but 20-50% carry another base, indel or clip (mixed culture)
  unresolved           depth below --mindepth; inspect in IGV

Needs python3 (standard library) + samtools (>= 1.7). BAMs must be sorted and indexed and aligned to a FASTA that
contains the assembly contigs (refs/index/NEW.fna from step 00 is fine).

Usage:
  python3 09_check_split_genes.py --ref refs/index/NEW.fna --gff refs/LK_1.gff --contig LNDBAELA \
      --bam DNA2014_PE=bams/DNA/DNA2014_PE_vs_LK1.bam DNA2017_PE=... DNA2023_PE=... \
      --locus tmk:LNDBAELA:4131300-4132750 valS:LNDBAELA:1225361-1228552 --out results/09
Outputs: <out>_loci.tsv (one row per locus) and <out>_changes.tsv (every majority disagreement).
"""
import argparse, re, subprocess, sys
from urllib.parse import unquote

def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True).stdout

def read_fasta_contig(path, contig):
    seq, on = [], False
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if on: break
                on = line[1:].split()[0] == contig
            elif on: seq.append(line.strip())
    if not seq: sys.exit(f'ERROR: contig {contig} not in {path}')
    return ''.join(seq).upper()

def read_pseudo(gff, contig):
    out = []
    with open(gff) as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 9 or c[0] != contig: continue
            a = dict(kv.split('=', 1) for kv in c[8].split(';') if '=' in kv)
            pseudo = a.get('pseudo', '').lower() == 'true' or c[2] == 'pseudogene'
            if not pseudo or c[2] not in ('gene', 'pseudogene', 'CDS'): continue
            lt = a.get('locus_tag') or a.get('ID', '')
            out.append((int(c[3]), int(c[4]), c[6], lt, unquote(a.get('product', a.get('Note', '')))))
    # one entry per locus tag (gene and CDS lines repeat it)
    prod = {}
    for s, e, st, lt, p in out:
        if p and not prod.get(lt): prod[lt] = p
    seen, uniq = set(), []
    for s, e, st, lt, p in sorted(out):
        if lt in seen: continue
        seen.add(lt); uniq.append([s, e, st, lt, prod.get(lt, '')])
    # merge fragments of one gene: same strand, same product, gap < 2 kb
    loci = []
    for s, e, st, lt, p in uniq:
        if loci and loci[-1]['strand'] == st and s - loci[-1]['end'] < 2000 and loci[-1]['product'] == p:
            loci[-1]['end'] = max(loci[-1]['end'], e); loci[-1]['tags'].append(lt)
        else:
            loci.append({'start': s, 'end': e, 'strand': st, 'tags': [lt], 'product': p})
    return loci

PILE = re.compile(r'[+-](\d+)')
def pileup(bam, ref, region, minmapq, minbq):
    """per position: depth, ref-agreeing count, base counts, indel counts {'+ACG':n,'-2':n}"""
    cmd = f'samtools mpileup -aa -B -q {minmapq} -Q {minbq} -d 1000000 -f {ref} -r {region} {bam} 2>/dev/null'
    res = {}
    for line in sh(cmd).splitlines():
        c = line.split('\t'); pos = int(c[1]); s = c[4]; i = 0
        base = {'A': 0, 'C': 0, 'G': 0, 'T': 0}; agree = 0; indel = {}
        while i < len(s):
            ch = s[i]
            if ch == '^': i += 2; continue
            if ch == '$': i += 1; continue
            if ch in '+-':
                m = PILE.match(s, i); n = int(m.group(1)); j = m.end()
                key = ch + s[j:j + n].upper(); indel[key] = indel.get(key, 0) + 1; i = j + n; continue
            if ch in '.,': agree += 1
            elif ch.upper() in base: base[ch.upper()] += 1
            i += 1
        depth = agree + sum(base.values())
        res[pos] = (depth, agree, base, indel)
    return res

def clips(bam, region, minmapq, minclip=20):
    cnt = {}
    for line in sh(f'samtools view -q {minmapq} -F 0x904 {bam} {region}').splitlines():
        c = line.split('\t'); pos = int(c[3]); cig = re.findall(r'(\d+)([MIDNSHP=X])', c[5])
        if cig and cig[0][1] == 'S' and int(cig[0][0]) >= minclip: cnt[pos] = cnt.get(pos, 0) + 1
        end = pos + sum(int(n) for n, op in cig if op in 'MDN=X') - 1
        if cig and cig[-1][1] == 'S' and int(cig[-1][0]) >= minclip: cnt[end + 1] = cnt.get(end + 1, 0) + 1
    return cnt

def depth_all(bam, region):
    d = {}
    for line in sh(f'samtools depth -a -r {region} {bam}').splitlines():
        c = line.split('\t'); d[int(c[1])] = int(c[2])
    return d

def median_depth(bam, contig, minmapq):
    hist = {}
    out = subprocess.Popen(f'samtools depth -a -Q {minmapq} -r {contig} {bam}', shell=True, stdout=subprocess.PIPE,
                           universal_newlines=True)
    n = 0
    for line in out.stdout:
        d = int(line.rsplit('\t', 1)[1]); hist[d] = hist.get(d, 0) + 1; n += 1
    out.wait(); k = 0
    for d in sorted(hist):
        k += hist[d]
        if k >= n / 2: return d
    return 0

COMP = str.maketrans('ACGTN', 'TGCAN')
STOP = {'TAA', 'TAG', 'TGA'}
def longest_orf(seq, strand):
    s = seq if strand == '+' else seq.translate(COMP)[::-1]
    best = 0
    for fr in range(3):
        run = 0
        for i in range(fr, len(s) - 2, 3):
            if s[i:i + 3] in STOP: best = max(best, run); run = 0
            else: run += 3
        best = max(best, run)
    return best

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ref', required=True); ap.add_argument('--gff', required=True); ap.add_argument('--contig', required=True)
    ap.add_argument('--bam', nargs='+', required=True, help='name=path, oldest first')
    ap.add_argument('--locus', nargs='*', default=[], help='name:contig:start-end (extra loci, e.g. genes missing from the annotation)')
    ap.add_argument('--out', required=True); ap.add_argument('--pad', type=int, default=300)
    ap.add_argument('--minmapq', type=int, default=10); ap.add_argument('--minbq', type=int, default=13)
    ap.add_argument('--major', type=float, default=0.5); ap.add_argument('--mindepth', type=int, default=10)
    a = ap.parse_args()
    bams = [b.split('=', 1) for b in a.bam]
    seq = read_fasta_contig(a.ref, a.contig)
    loci = read_pseudo(a.gff, a.contig)
    for L in a.locus:
        name, c, rng = L.split(':'); s, e = map(int, rng.replace(',', '').split('-'))
        loci.append({'start': s, 'end': e, 'strand': '.', 'tags': [name], 'product': 'user locus'})
    print(f'{len(loci)} loci ({sum(1 for l in loci if l["product"] != "user locus")} from PGAP pseudogenes)', file=sys.stderr)
    med = {}
    for n, b in bams:
        med[n] = median_depth(b, a.contig, a.minmapq); print(f'{n}: median depth on {a.contig} = {med[n]}', file=sys.stderr)
    LT = open(a.out + '_loci.tsv', 'w'); CH = open(a.out + '_changes.tsv', 'w')
    LT.write('\t'.join(['locus', 'locus_tags', 'product', 'start', 'end', 'strand'] +
                       [f'{n}_{k}' for n, _ in bams for k in ('mean_depth_x_median', 'majority_changes', 'minority_changes', 'gap_bp', 'clip_pileups')] +
                       ['longest_ORF_assembly_bp', 'longest_ORF_after_newest_changes_bp', 'verdict']) + '\n')
    CH.write('locus\tdataset\tcontig\tpos\tassembly\tillumina_majority\tfraction\tdepth\n')
    for L in loci:
        s0, e0 = max(1, L['start'] - a.pad), min(len(seq), L['end'] + a.pad)
        region = f'{a.contig}:{s0}-{e0}'; name = L['tags'][0] if len(L['tags']) == 1 else L['tags'][0] + '..' + L['tags'][-1]
        row = [name, ','.join(L['tags']), L['product'], str(L['start']), str(L['end']), L['strand']]
        per = {}
        for n, b in bams:
            P = pileup(b, a.ref, region, a.minmapq, a.minbq); C = clips(b, region, a.minmapq)
            inl = [p for p in P if L['start'] <= p <= L['end']]
            md = sum(P[p][0] for p in inl) / max(1, len(inl))
            maj, mino, changes = 0, 0, []
            for p in inl:
                d, ag, base, ind = P[p]
                if d < a.mindepth: continue
                cand = [(v, k) for k, v in base.items()] + [(v, k) for k, v in ind.items()]
                if not cand: continue
                v, k = max(cand)
                f = v / d
                if f >= a.major: maj += 1; changes.append((p, k, f, d))
                elif f >= 0.2: mino += 1
            # gap: >= 50 bp with MAPQ-filtered depth < 0.2x median. If reads of any MAPQ cover it (>= 0.5x median)
            # it is a repeat (reads map equally well elsewhere), not unsupported sequence.
            DA = depth_all(b, region)
            gap, rep, run, runrep = 0, 0, 0, 0
            for p in list(range(L['start'], L['end'] + 1)) + [None]:
                low = p is not None and P.get(p, (0,))[0] < 0.2 * med[n]
                if low:
                    run += 1; runrep += DA.get(p, 0) >= 0.5 * med[n]
                else:
                    if run >= 50:
                        if runrep > run / 2: rep += run
                        else: gap += run
                    run, runrep = 0, 0
            # clipping pile-ups: >= half the reads clipped = the assembly is wrong there; 0.2-0.5 = a minority structural variant
            cp = sum(1 for p, v in C.items() if s0 <= p <= e0 and v >= max(5, 0.5 * (P.get(p, (0,))[0] + v)))
            cpm = sum(1 for p, v in C.items() if s0 <= p <= e0 and max(5, 0.2 * (P.get(p, (0,))[0] + v)) <= v < 0.5 * (P.get(p, (0,))[0] + v))
            per[n] = (maj, mino, gap, cp, md, changes, cpm, rep)
            row += [f'{md / max(1, med[n]):.2f}', str(maj), str(mino), f'{gap}' + (f' (+{rep} repeat)' if rep else ''), f'{cp} (+{cpm} minority)']
            for p, k, f, d in changes:
                asm = seq[p - 1] if k[0] not in '+-' else seq[p - 1] + ('' if k[0] == '+' else seq[p:p + len(k) - 1])
                CH.write(f'{name}\t{n}\t{a.contig}\t{p}\t{asm}\t{k if k[0] not in "+-" else ("ins " + k[1:] if k[0] == "+" else "del " + k[1:])}\t{f:.2f}\t{d}\n')
        # apply the newest dataset's majority changes and re-measure the longest ORF
        new = bams[-1][0]; loc = list(seq[L['start'] - 1:L['end']])
        for p, k, f, d in sorted(per[new][5], reverse=True):
            i = p - L['start']
            if k[0] == '+': loc[i] = loc[i] + k[1:]
            elif k[0] == '-': del loc[i + 1:i + 1 + len(k) - 1]
            else: loc[i] = k
        st = L['strand'] if L['strand'] in '+-' else '+'
        orf0 = max(longest_orf(seq[L['start'] - 1:L['end']], s) for s in ((st,) if L['strand'] in '+-' else '+-'))
        orf1 = max(longest_orf(''.join(loc), s) for s in ((st,) if L['strand'] in '+-' else '+-'))
        old = [per[n] for n, _ in bams[:-1]]; nw = per[new]
        lowdepth = nw[4] < a.mindepth
        if nw[2] > 0 or nw[3] > 0:
            v = f'assembly error: {nw[2]} bp not covered by reads' + (f', {nw[0]} base/indel disagreements' if nw[0] else '') + ' (misassembled or error-dense segment)'
        elif nw[0] > 0: v = f'assembly error: {nw[0]} base/indel disagreements'
        elif nw[7] > 0: v = 'unresolved: repeat, reads map equally well elsewhere (inspect in IGV)'
        elif lowdepth: v = 'unresolved: low depth (inspect in IGV)'
        elif old and any(x[0] > 0 for x in old): v = 'split supported; recent mutation (absent in older reads, check step 08)'
        elif nw[1] > 0 or nw[6] > 0: v = 'split supported in the majority; minority alleles present (inspect in IGV)'
        else: v = 'split supported'
        row += [str(orf0), str(orf1), v]
        LT.write('\t'.join(row) + '\n'); print(f'{name}\t{v}', file=sys.stderr)
    LT.close(); CH.close()

if __name__ == '__main__':
    main()
