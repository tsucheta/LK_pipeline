#!/usr/bin/env python3
"""
08_call_variants.py - genome-wide calling of variants that arose between DNA datasets, for every organism.
This is the only variant-calling step of the pipeline; it replaces the earlier GATK/SnpEff lists and
08_scan_new_mutations.py. Earlier calls are NOT used as input; a --known file can be given only to label
calls that were also reported before (a known site that is not re-detected here is not reported).

What is called
  A site where one allele is absent from the FIRST dataset and present in the LAST one (SNVs and 1-3 bp
  indels), in reads mapped competitively to all 2023 genomes (bams/DNA/<sample>_vs_LK1.bam, refs/index/NEW.fna).
  The 2023 assemblies may already carry the new allele, so both directions are tested and every call is
  written as ancestral>derived with the frequency of the DERIVED allele in each dataset.

High-confidence rules (FILTER = PASS only if all hold; failing calls are kept with the reason)
  depth          >= --mindepth (15) reads with MAPQ >= 20 and base quality >= 20 in the first AND last dataset
  first dataset  derived allele <= --maxearly (5%)
  last dataset   derived allele >= --minlate (20%) and >= --minalt (5) reads
  lowMQ          <= 10% of the reads at the site have MAPQ < 20 (repeats, IS copies, rRNA operons)
  highDP         depth <= 2.5 x the organism's median depth in that dataset (collapsed repeats)
  strand         derived reads on both strands, and Fisher strand-bias P >= 0.001 (derived vs ancestral reads;
                 for indels from the spanning reads)
  edge           >= 150 bp from a contig end
  cluster        fewer than 3 candidate calls within 100 bp (misassembly, divergent strain segments)
  presentFirst / rareLast / lowDepthSpanning
                 indels that met the rules in the pileup but not after the spanning-read recount (kept, not PASS)
  hpNoise        homopolymer indels: last-dataset frequency >= 3x the slippage rate of the same +/-1 slip at up to
                 100 other runs of the same base and length in that genome (reported per dataset)
Duplicates, secondary and QC-fail reads are ignored. Indels are found in the pileup (left-normalised, so a
homopolymer slip is counted once whatever position the aligner gave it), then RECOUNTED from reads that span
the whole run plus 8 bp on each side, one read per fragment: each read's net indel length is taken from its
CIGAR, so indel frequencies are exact and need no correction in either direction.

Needs python3 (standard library) + samtools (>= 1.7).

Usage
  python3 08_call_variants.py --ref refs/index/NEW.fna --bamsuffix LK1 \\
      --org LK=refs/LK_1.fna Inqui=refs/inqui.fna Phaeo=refs/phaeo.fna Ocean=refs/ocean.fna Nitrat=refs/LK_bac_polished.fna \\
      --gff LK=refs/LK_1.gff Inqui=refs/inqui.gff Phaeo=refs/phaeo.gff Ocean=refs/ocean.gff Nitrat=refs/LK_bac_polished.gff \\
      --rna LK=bams/LK Inqui=bams/Inquil Phaeo=bams/Phaeo Ocean=bams/Ocean Nitrat=bams/Nitrat \\
      --rna-exclude LK-7R1_S78,LK-7R2_S79,LK4-1_S80,LK4-2_S81 \\
      --known known_variants_LK1.tsv \\
      --samples DNA2014_PE DNA2017_PE DNA2023_PE
  (samples from oldest to newest; each --org FASTA only supplies the contig names of that organism;
   --gff, --rna and --known are optional)
Outputs (prefix --out, default results/08)
  _variants.tsv   every candidate with FILTER, gene/effect, per-dataset derived/depth and frequency, RNA support
  _<org>.vcf      the same per organism (VCF 4.2, REF/ALT relative to the 2023 assembly)
  _spectrum.tsv   PASS calls per organism: transitions, transversions, homopolymer and other indels
  _depth.tsv      median depth per organism and dataset (organisms below --mindepth cannot be assessed)
"""
import argparse, math, os, re, subprocess, sys
from collections import Counter, defaultdict, deque
from urllib.parse import unquote

MINMQ, MINBQ = 20, 20
COMP = str.maketrans("ACGTN", "TGCAN")
CODON = {a + b + c: aa for (a, b, c), aa in zip(
    [(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}

def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True).stdout

def read_fasta(path, want=None):
    seqs, name, buf = {}, None, []
    for l in open(path):
        if l.startswith(">"):
            if name is not None and (want is None or name in want): seqs[name] = "".join(buf).upper()
            name, buf = l[1:].split()[0], []
        else: buf.append(l.strip())
    if name is not None and (want is None or name in want): seqs[name] = "".join(buf).upper()
    return seqs

def headers(path):
    return [l[1:].split()[0] for l in open(path) if l.startswith(">")]

# ------------------------------------------------------------------ pileup parsing
def parse_reads(bases, quals, mqs, refb):
    """yield (allele, is_reverse, indel or None, mapq) per read; allele is a base, or '*' for a deletion placeholder"""
    i, k, n = 0, 0, len(bases)
    out = []
    while i < n:
        c = bases[i]
        if c == '^': i += 2; continue
        if c == '$': i += 1; continue
        if c in '+-':
            m = re.match(r'\d+', bases[i + 1:]); L = int(m.group()); s = bases[i + 1 + len(m.group()): i + 1 + len(m.group()) + L]
            if out: out[-1][2] = (c, s.upper())
            i += 1 + len(m.group()) + L; continue
        mq = ord(mqs[k]) - 33 if k < len(mqs) else 60
        if c in '.,': out.append([refb, c == ',', None, mq])
        elif c in '*#': out.append(['*', c == '#', None, mq])
        elif c in '<>': pass
        else: out.append([c.upper(), c.islower(), None, mq])
        k += 1; i += 1
    return out

def left_norm(seq, a, typ, s):
    """left-normalise an indel. a = 1-based anchor position; returns (anchor, inserted/deleted sequence)"""
    L = len(s)
    if typ == '+':
        while a >= 2 and seq[a - 1] == s[-1]:
            s = s[-1] + s[:-1]; a -= 1
    else:
        s = seq[a:a + L]
        while a >= 2 and seq[a - 1] == seq[a + L - 1]:
            a -= 1; s = seq[a:a + L]
    return a, s

def homopolymer_len(seq, a, s):
    """length of the reference run of base s[0] starting right after the anchor (only for mono-base indels)"""
    if len(set(s)) != 1: return 0
    b, k, n = s[0], a, 0
    while k < len(seq) and seq[k] == b: n += 1; k += 1
    return n

def fisher_p(a, b, c, d):
    """two-sided Fisher exact test for [[a, b], [c, d]]"""
    n = a + b + c + d; r1, c1 = a + b, a + c
    def lp(x): return (math.lgamma(r1 + 1) + math.lgamma(n - r1 + 1) + math.lgamma(c1 + 1) + math.lgamma(n - c1 + 1)
                       - math.lgamma(n + 1) - math.lgamma(x + 1) - math.lgamma(r1 - x + 1) - math.lgamma(c1 - x + 1) - math.lgamma(n - r1 - c1 + x + 1))
    lo, hi = max(0, r1 + c1 - n), min(r1, c1)
    p0 = lp(a)
    return min(1.0, sum(math.exp(lp(x)) for x in range(lo, hi + 1) if lp(x) <= p0 + 1e-9))

def span_counts(bam, ctg, rs, re_, flank=8, minmq=MINMQ):
    """reads (one per fragment) that span rs-flank .. re_+flank; returns Counter(net indel length inside rs..re_)
    and per-delta strand counts. Exact for homopolymer slips: no read ends inside the run."""
    cnt, strand, seen = Counter(), defaultdict(lambda: [0, 0]), set()
    out = sh(f"samtools view -F 0xF04 -q {minmq} {bam} '{ctg}:{max(1, rs - flank)}-{re_ + flank}'")
    for line in out.splitlines():
        f = line.split("\t")
        if f[0] in seen: continue
        pos, cig, flag = int(f[3]), f[5], int(f[1])
        ops = re.findall(r"(\d+)([MIDNSHP=X])", cig)
        end = pos - 1 + sum(int(l) for l, o in ops if o in "MDN=X")
        if pos > rs - flank or end < re_ + flank: continue
        p, d = pos, 0
        for l, o in ops:
            l = int(l)
            if o in "M=X": p += l
            elif o in "DN":
                if p + l - 1 >= rs and p <= re_: d -= l
                p += l
            elif o == "I":
                if rs - 1 <= p - 1 <= re_: d += l
        seen.add(f[0]); cnt[d] += 1; strand[d][1 if flag & 16 else 0] += 1
    return cnt, strand

# ------------------------------------------------------------------ annotation
def parse_gff(path):
    feats = defaultdict(list)
    for l in open(path):
        if l.startswith("##FASTA"): break
        if l.startswith("#") or l.count("\t") < 8: continue
        f = l.rstrip("\n").split("\t")
        if f[2] not in ("CDS", "gene", "pseudogene", "tRNA", "rRNA"): continue
        a = {k: unquote(v) for k, v in (x.split("=", 1) for x in f[8].split(";") if "=" in x)}
        feats[f[0]].append(dict(type=f[2], s=int(f[3]), e=int(f[4]), strand=f[6], tag=a.get("locus_tag", a.get("ID", "")),
                                product=a.get("product", ""), pseudo="pseudo" in a or f[2] == "pseudogene"))
    for c in feats: feats[c].sort(key=lambda x: x["s"])
    return feats

def annotate(feats, seq, ctg, pos, ref, alt, anc, der):
    """effect of ancestral -> derived. ref/alt are VCF alleles at pos (1-based); anc/der the same pair in ancestral order"""
    F = feats.get(ctg, [])
    cds = [f for f in F if f["type"] == "CDS" and f["s"] <= pos + max(len(ref), len(alt)) - 1 and f["e"] >= pos]
    hit = [f for f in F if f["s"] <= pos and f["e"] >= pos]
    if not cds:
        if hit:
            f = hit[0]; return f["tag"], f["product"], f"non-coding feature ({f['type']})"
        near = min(F, key=lambda f: min(abs(f["s"] - pos), abs(f["e"] - pos)), default=None)
        if not near: return "", "", "intergenic"
        d = min(abs(near["s"] - pos), abs(near["e"] - pos))
        return "", "", f"intergenic; {d} bp from {near['tag']} ({near['product'][:40]})"
    f = cds[0]; tag, prod = f["tag"], f["product"]
    pseudo = " (pseudogene)" if f["pseudo"] else ""
    if len(ref) != len(alt):
        L = abs(len(ref) - len(alt))
        off = (pos + 1 - f["s"]) if f["strand"] == "+" else (f["e"] - pos)
        cod = max(1, off // 3 + 1)
        return tag, prod, ("frameshift" if L % 3 else "in-frame indel") + f" near codon {cod}{pseudo}"
    gs = seq[f["s"] - 1:f["e"]]
    if f["strand"] == "+": cdsq, off, a_, d_ = gs, pos - f["s"], anc, der
    else: cdsq, off, a_, d_ = gs.translate(COMP)[::-1], f["e"] - pos, anc.translate(COMP), der.translate(COMP)
    ci = off // 3; c0 = list(cdsq[ci * 3: ci * 3 + 3])
    if len(c0) < 3: return tag, prod, "coding (partial codon)" + pseudo
    c_anc = c0[:]; c_anc[off % 3] = a_; c_der = c0[:]; c_der[off % 3] = d_
    aa0, aa1 = CODON.get("".join(c_anc), "X"), CODON.get("".join(c_der), "X")
    kind = "synonymous" if aa0 == aa1 else ("stop gained" if aa1 == "*" else ("stop lost" if aa0 == "*" else "missense"))
    return tag, prod, f"{kind} p.{aa0}{ci + 1}{aa1} ({''.join(c_anc)}>{''.join(c_der)}){pseudo}"

# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True); ap.add_argument("--org", nargs="+", required=True)
    ap.add_argument("--samples", nargs="+", required=True, help="DNA datasets, oldest first")
    ap.add_argument("--bamsuffix", default="LK1")
    ap.add_argument("--gff", nargs="*", default=[]); ap.add_argument("--rna", nargs="*", default=[])
    ap.add_argument("--rna-exclude", default=""); ap.add_argument("--known", default="")
    ap.add_argument("--mindepth", type=int, default=15); ap.add_argument("--maxearly", type=float, default=0.05)
    ap.add_argument("--minlate", type=float, default=0.20); ap.add_argument("--minalt", type=int, default=5)
    ap.add_argument("--maxlowmq", type=float, default=0.10); ap.add_argument("--maxdepthx", type=float, default=2.5)
    ap.add_argument("--edge", type=int, default=150); ap.add_argument("--maxindel", type=int, default=3)
    ap.add_argument("--out", default="results/08")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    S = a.samples; NS = len(S)
    bams = [f"bams/DNA/{s}_vs_{a.bamsuffix}.bam" for s in S]
    for b in bams:
        if not os.path.exists(b): sys.exit(f"ERROR: missing {b}")
    bamref = set(l.split("\t")[1][3:] for l in sh(f"samtools view -H {bams[0]}").splitlines() if l.startswith("@SQ"))
    org_of = {}
    for spec in a.org:
        org, path = spec.split("=", 1)
        n = 0
        for h in headers(path):
            for c in (h, "NZ_" + h, h[3:] if h.startswith("NZ_") else ""):
                if c and c in bamref: org_of[c] = org; n += 1; break
        print(f"{org}: {n} contigs found in the BAM header")
        if n == 0: print(f"  WARNING: none of the contigs of {path} are in {bams[0]}; {org} is skipped")
    if not org_of: sys.exit("ERROR: none of the --org contig names are in the BAM header")
    seqs = read_fasta(a.ref, set(org_of))
    gffs = {o: parse_gff(p) for o, p in (x.split("=", 1) for x in a.gff)}

    bed = a.out + "_contigs.bed"
    with open(bed, "w") as o:
        for c in org_of:
            if c in seqs: o.write(f"{c}\t0\t{len(seqs[c])}\n")
    cmd = f"samtools mpileup -B -q 0 -Q {MINBQ} -s -d 100000 -l {bed} -f {a.ref} " + " ".join(bams)
    print(f"Scanning {sum(len(seqs[c]) for c in org_of if c in seqs):,} bp in {NS} datasets (this takes a while) ...")
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True)

    hist = defaultdict(Counter)                       # (org, k) -> Counter(depth) for medians
    snv = {}                                          # (ctg, pos, allele) -> per-sample stats
    ind = defaultdict(lambda: [[0, 0] for _ in range(NS)])   # (ctg, anchor, typ, seq) -> per-sample [fwd, rev]
    ind_lowmq = defaultdict(float)
    recent = {}                                       # ctg -> deque of (pos, [hq depth per sample], [lowmq frac])
    depth_at = {}                                     # (ctg, pos) -> ([depth], [lowmq]) kept for indel anchors only
    nline = 0
    for line in p.stdout:
        nline += 1
        if nline % 2_000_000 == 0: print(f"  {nline:,} positions", flush=True)
        f = line.rstrip("\n").split("\t")
        ctg, pos, refb = f[0], int(f[1]), f[2].upper()
        org = org_of.get(ctg)
        if org is None: continue
        cols = [f[3 + 4 * k: 7 + 4 * k] for k in range(NS)]
        for k in range(NS): hist[(org, k)][min(int(cols[k][0]), 100000)] += 1
        raw = [c[1] if len(c) == 4 else "" for c in cols]
        has_ind = any(('+' in r or '-' in r) for r in raw)
        nalt = [sum(r.count(x) for x in "ACGTacgt") for r in raw]
        dq = recent.setdefault(ctg, deque(maxlen=64))
        if not has_ind and max(nalt[0], nalt[-1]) < a.minalt:
            dq.append((pos, None)); continue
        per, dep, lmq = [], [], []
        for k, c in enumerate(cols):
            if len(c) < 4 or c[0] == '0': per.append([]); dep.append(0); lmq.append(0.0); continue
            rr = parse_reads(c[1], c[2], c[3], refb)
            hq = [r for r in rr if r[3] >= MINMQ]
            per.append(hq); dep.append(sum(1 for r in hq if r[0] != '*'))
            lmq.append(1 - len(hq) / len(rr) if rr else 0.0)
        dq.append((pos, (dep, lmq)))
        # SNV candidates
        if max(nalt[0], nalt[-1]) >= a.minalt:
            alle = set()
            for k in (0, NS - 1):
                alle |= {r[0] for r in per[k] if r[0] in "ACGT" and r[0] != refb}
            for al in alle:
                st = []
                for k in range(NS):
                    af = sum(1 for r in per[k] if r[0] == al and not r[1]); ar = sum(1 for r in per[k] if r[0] == al and r[1])
                    rf = sum(1 for r in per[k] if r[0] == refb and not r[1]); rr_ = sum(1 for r in per[k] if r[0] == refb and r[1])
                    st.append((af, ar, rf, rr_, dep[k], lmq[k]))
                snv[(ctg, pos, al)] = st
        # indel evidence (left-normalised)
        if has_ind:
            seq = seqs[ctg]
            for k in range(NS):
                for r in per[k]:
                    if r[2] is None: continue
                    typ, s = r[2]
                    if len(s) > a.maxindel or 'N' in s: continue
                    na, ns_ = left_norm(seq, pos, typ, s)
                    key = (ctg, na, typ, ns_)
                    ind[key][k][1 if r[1] else 0] += 1
                    ind_lowmq[key] = max(ind_lowmq[key], lmq[k])
                    if na == pos: depth_at[(ctg, na)] = (dep, lmq)
                    elif (ctg, na) not in depth_at:
                        # depth at the normalised anchor if that line was parsed, else the depth here (a few bp away)
                        got = next((x[1] for x in dq if x[0] == na and x[1] is not None), None)
                        depth_at[(ctg, na)] = got or (dep, lmq)
    p.wait()
    if p.returncode not in (0, None): sys.exit("ERROR: samtools mpileup failed")

    med = {}
    for (org, k), h in hist.items():
        tot, acc = sum(h.values()), 0
        for d in sorted(h):
            acc += h[d]
            if acc >= tot / 2: med[(org, k)] = d; break
    orgs = sorted(set(org_of.values()))
    with open(a.out + "_depth.tsv", "w") as o:
        o.write("organism\t" + "\t".join(f"{s}_median_depth" for s in S) + "\tassessable\n")
        for org in orgs:
            ms = [med.get((org, k), 0) for k in range(NS)]
            o.write(f"{org}\t" + "\t".join(map(str, ms)) + "\t" + ("yes" if ms[0] >= a.mindepth and ms[-1] >= a.mindepth else "NO (median depth below --mindepth)") + "\n")
    print(open(a.out + "_depth.tsv").read())

    # ------------------------------------------------------------------ evaluate candidates
    calls = []
    def basic(fr):
        d0, dl = fr[0][1], fr[-1][1]
        return (d0 >= a.mindepth and dl >= a.mindepth and fr[0][2] is not None and fr[0][2] <= a.maxearly
                and fr[-1][0] >= a.minalt and fr[-1][2] >= a.minlate)

    for (ctg, pos, al), st in snv.items():
        refb = seqs[ctg][pos - 1]
        alt_n = [x[0] + x[1] for x in st]; ref_n = [x[2] + x[3] for x in st]; dep = [x[4] for x in st]; lmq = [x[5] for x in st]
        for direction in ("ref=ancestral", "ref=derived"):
            der_n = alt_n if direction == "ref=ancestral" else ref_n
            # reference = derived: this allele must BE the ancestral state, i.e. the majority allele in the first dataset
            # (otherwise a minor error allele at the same site would be reported as a second, spurious call)
            if direction == "ref=derived" and (not dep[0] or alt_n[0] / dep[0] < 0.5):
                continue
            fr = [(der_n[k], dep[k], der_n[k] / dep[k] if dep[k] else None) for k in range(NS)]
            if not basic(fr): continue
            L = st[-1]
            if direction == "ref=ancestral": df, dr, af_, ar_ = L[0], L[1], L[2], L[3]
            else: df, dr, af_, ar_ = L[2], L[3], L[0], L[1]
            anc, der = (refb, al) if direction == "ref=ancestral" else (al, refb)
            cls = "transition" if {refb, al} in ({"A", "G"}, {"C", "T"}) else "transversion"
            calls.append(dict(ctg=ctg, pos=pos, ref=refb, alt=al, anc=anc, der=der, cls=cls, direction=direction, fr=fr,
                              strand=(df, dr, af_, ar_), lowmq=max(lmq[0], lmq[-1]), indel=False))
    for (ctg, na, typ, s), cnt in ind.items():
        if (ctg, na) not in depth_at: continue
        dep, lmq = depth_at[(ctg, na)]
        seq = seqs[ctg]; refb = seq[na - 1]
        alt_n = [c[0] + c[1] for c in cnt]
        if typ == '+': ref_a, alt_a = refb, refb + s
        else: ref_a, alt_a = refb + s, refb
        # the 2023 assembly may carry the derived state: then the indel allele (vs reference) is ANCESTRAL.
        # Indels are under-counted near read ends, so a fixed ancestral indel shows at ~85-95%; scale by the
        # first-dataset level before testing that direction.
        for direction in ("ref=ancestral", "ref=derived"):
            if direction == "ref=ancestral":
                der_n = alt_n; fr = [(der_n[k], dep[k], der_n[k] / dep[k] if dep[k] else None) for k in range(NS)]
            else:
                f0 = alt_n[0] / dep[0] if dep[0] else 0
                if f0 < 0.70: continue
                fr = []
                for k in range(NS):
                    if not dep[k]: fr.append((0, 0, None)); continue
                    fk = max(0.0, 1 - (alt_n[k] / dep[k]) / f0)
                    fr.append((round(fk * dep[k]), dep[k], fk))
            if not basic(fr): continue
            df, dr = (cnt[-1][0], cnt[-1][1]) if direction == "ref=ancestral" else (None, None)
            anc, der = (ref_a, alt_a) if direction == "ref=ancestral" else (alt_a, ref_a)
            hp = homopolymer_len(seq, na, s)
            gained = (len(der) > len(anc))
            cls = ("ins" if gained else "del") + str(len(s)) + (f"_homopolymer{hp + (0 if typ == '+' else 0)}" if hp >= 5 else "")
            calls.append(dict(ctg=ctg, pos=na, ref=ref_a, alt=alt_a, anc=anc, der=der, cls=cls, direction=direction, fr=fr,
                              strand=(df, dr, None, None), lowmq=max(lmq[0], lmq[-1], ind_lowmq[(ctg, na, typ, s)]), indel=True))

    # ------------------------------------------------------------------ indels: recount with spanning reads
    # Pileup-based indel counts miss reads that end inside a run, and the 'reference = derived' direction needs a
    # correction that depends on read length. Every indel candidate is therefore recounted from reads that span the
    # whole run (+8 bp each side): each read's net indel length in the run is read directly from its CIGAR.
    kept = []
    for c in calls:
        if not c["indel"]:
            kept.append(c); continue
        seq = seqs[c["ctg"]]; L = abs(len(c["alt"]) - len(c["ref"])); s_ = (c["alt"] if len(c["alt"]) > len(c["ref"]) else c["ref"])[1:]
        hp = homopolymer_len(seq, c["pos"], s_)
        rs, re_ = c["pos"], c["pos"] + max(hp, L)          # anchor base .. end of the run / deleted bases
        d_alt = L if len(c["alt"]) > len(c["ref"]) else -L # ALT relative to REF
        d_der, d_anc = (d_alt, 0) if c["direction"] == "ref=ancestral" else (0, d_alt)
        fr, st_last, st_anc = [], (0, 0), (0, 0)
        for k, b in enumerate(bams):
            cnt, strand = span_counts(b, c["ctg"], rs, re_)
            der, anc = cnt[d_der], cnt[d_anc]
            n = der + anc
            fr.append((der, n, der / n if n else None))
            if k == NS - 1: st_last = tuple(strand[d_der]); st_anc = tuple(strand[d_anc])
        # strand counts of derived AND ancestral spanning reads: the same Fisher strand-bias test as for SNVs
        c["fr"] = fr; c["strand"] = (st_last[0], st_last[1], st_anc[0], st_anc[1]); c["span"] = True
        c["pileup_fr"] = None
        if not basic(fr):
            # keep it, with the reason, so a site that only looked new in the pileup counts is visible
            why = []
            if fr[0][1] < a.mindepth or fr[-1][1] < a.mindepth: why.append("lowDepthSpanning")
            if fr[0][2] is not None and fr[0][2] > a.maxearly: why.append("presentFirst")
            if fr[-1][2] is not None and (fr[-1][2] < a.minlate or fr[-1][0] < a.minalt): why.append("rareLast")
            c["pre"] = why or ["notNew"]
        kept.append(c)
    print(f"Indels recounted from spanning reads: {sum(1 for c in calls if c['indel'])} candidates, "
          f"{sum(1 for c in kept if c['indel'] and not c.get('pre'))} still meet the depth/frequency rules")
    calls = kept

    # ------------------------------------------------------------------ filters
    bypos = defaultdict(list)
    for c in calls: bypos[c["ctg"]].append(c["pos"])
    for c in bypos.values(): c.sort()
    import bisect
    for c in calls:
        org = org_of[c["ctg"]]; F = list(c.get("pre", []))
        if c["lowmq"] > a.maxlowmq: F.append("lowMQ")
        for k in (0, NS - 1):
            m = med.get((org, k), 0)
            if m and c["fr"][k][1] > a.maxdepthx * m: F.append("highDP"); break
        df, dr, af_, ar_ = c["strand"]
        if df is not None:
            if min(df, dr) < 1: F.append("strand")
            elif af_ is not None and fisher_p(af_, ar_, df, dr) < 0.001: F.append("strand")
        L = len(seqs[c["ctg"]])
        if c["pos"] <= a.edge or c["pos"] > L - a.edge: F.append("edge")
        P = bypos[c["ctg"]]; lo = bisect.bisect_left(P, c["pos"] - 100); hi = bisect.bisect_right(P, c["pos"] + 100)
        if hi - lo >= 3: F.append("cluster")
        c["filter"] = ";".join(F) or "PASS"
        c["org"] = org

    # ------------------------------------------------------------------ homopolymer slippage background
    # For each homopolymer indel call, measure the same slip (same base, run length, +1/-1) at up to 100 other runs of
    # that length in the same genome. Polymerase slippage in long runs can reach several percent; a call is flagged
    # hpNoise unless its last-dataset frequency is >= 3x that background and its first-dataset background is low.
    import random
    rng = random.Random(1)
    hp_calls = [c for c in calls if "homopolymer" in c["cls"]]
    groups = defaultdict(list)
    for c in hp_calls:
        L = int(re.search(r"homopolymer(\d+)", c["cls"]).group(1))
        b = (c["alt"][1:] if len(c["alt"]) > len(c["ref"]) else c["ref"][1:])[0]
        typ = '+' if len(c["alt"]) > len(c["ref"]) else '-'
        groups[(c["org"], b, L, typ)].append(c)
    for (org, b, L, typ), cs in groups.items():
        runs = []
        for ctg in [x for x in org_of if org_of[x] == org and x in seqs]:
            q = seqs[ctg]
            for base in {b, b.translate(COMP)}:
                for m in re.finditer(f"(?<!{base}){base}{{{L}}}(?!{base})", q):
                    st = m.start() + 1                          # 1-based first base of the run
                    if st - 1 < a.edge or m.end() > len(q) - a.edge: continue
                    if any(cc["ctg"] == ctg and abs(cc["pos"] - st) < 50 for cc in cs): continue
                    runs.append((ctg, st, base))
        rng.shuffle(runs); runs = runs[:100]
        if not runs:
            for c in cs: c["hpbg"] = None
            continue
        delta = 1 if typ == '+' else -1
        bg = []
        for k, bm in enumerate(bams):
            tot = sl = 0
            for ctg, st, base in runs:
                cnt, _ = span_counts(bm, ctg, st - 1, st + L - 1)
                tot += sum(cnt.values()); sl += cnt[delta]
            bg.append(sl / tot if tot else None)
        for c in cs: c["hpbg"] = (bg, len(runs))
    for c in calls:
        hb = c.get("hpbg")
        if not hb: continue
        bg, n = hb; fl = c["fr"][-1][2]
        if bg[-1] is not None and fl is not None and fl < 3 * bg[-1]:
            c["filter"] = "hpNoise" if c["filter"] == "PASS" else c["filter"] + ";hpNoise"

    # ------------------------------------------------------------------ known sites (labels only)
    known = set()
    if a.known and os.path.exists(a.known):
        for l in open(a.known):
            if l.startswith("#") or not l.strip(): continue
            x = l.rstrip("\n").split("\t")
            known.add((x[0], int(x[1])))
    for c in calls:
        L = max(len(c["ref"]), len(c["alt"]))
        c["known"] = "yes" if any((c["ctg"], q) in known for q in range(c["pos"] - 10, c["pos"] + L + 10)) else ""

    # ------------------------------------------------------------------ annotation
    for c in calls:
        G = gffs.get(c["org"])
        c["tag"], c["product"], c["effect"] = annotate(G, seqs[c["ctg"]], c["ctg"], c["pos"], c["ref"], c["alt"], c["anc"], c["der"]) if G else ("", "", "")

    # ------------------------------------------------------------------ RNA support (pooled RNA reads per organism)
    rna_dirs = dict(x.split("=", 1) for x in a.rna)
    rexcl = set(x for x in a.rna_exclude.split(",") if x)
    for c in calls: c["rna"] = "not run"
    for org, d in rna_dirs.items():
        cs = [c for c in calls if c["org"] == org]
        rb = sorted(os.path.join(d, x) for x in os.listdir(d) if x.endswith(".bam") and x[:-4] not in rexcl) if os.path.isdir(d) else []
        if not cs or not rb: continue
        rbed = a.out + f"_{org}_rna.bed"
        with open(rbed, "w") as o:
            for c in sorted(cs, key=lambda c: (c["ctg"], c["pos"])):
                o.write(f"{c['ctg']}\t{max(0, c['pos'] - 15)}\t{c['pos'] + 15}\n")
        out = sh(f"samtools mpileup -B -q 0 -Q {MINBQ} -s -d 1000000 -l {rbed} -f {a.ref} " + " ".join(rb))
        pile = {}
        for line in out.splitlines():
            f = line.split("\t"); pos = int(f[1]); refb = f[2].upper()
            reads = []
            for k in range(len(rb)):
                col = f[3 + 4 * k: 7 + 4 * k]
                if len(col) == 4 and col[0] != '0': reads += [r for r in parse_reads(col[1], col[2], col[3], refb) if r[3] >= MINMQ]
            pile[(f[0], pos)] = reads
        for c in cs:
            reads = pile.get((c["ctg"], c["pos"]), [])
            dep = sum(1 for r in reads if r[0] != '*')
            if not c["indel"]:
                der = sum(1 for r in reads if r[0] == c["der"])
            else:
                typ = '+' if len(c["alt"]) > len(c["ref"]) else '-'
                s = (c["alt"][1:] if typ == '+' else c["ref"][1:])
                n_ind = 0
                for q in range(c["pos"] - 15, c["pos"] + 16):
                    for r in pile.get((c["ctg"], q), []):
                        if r[2] and r[2][0] == typ and len(r[2][1]) == len(s) and left_norm(seqs[c["ctg"]], q, typ, r[2][1]) == (c["pos"], s):
                            n_ind += 1
                der = n_ind if c["direction"] == "ref=ancestral" else max(0, dep - n_ind)
            c["rna"] = f"{der}/{dep}" + (f"={der / dep:.2f}" if dep else "")

    # ------------------------------------------------------------------ write
    calls.sort(key=lambda c: (c["org"], c["ctg"], c["pos"]))
    cols = (["organism", "contig", "pos", "REF", "ALT", "change(ancestral>derived)", "class", "reference_allele", "FILTER",
             "first_seen", "locus_tag", "product", "effect", "previously_reported"]
            + [f"{s}_derived/depth" for s in S] + [f"{s}_derived_freq" for s in S] + ["RNA_derived/depth_2023", "lowMQ_fraction", "derived_strand_fwd/rev", "homopolymer_slip_background_per_dataset"])
    with open(a.out + "_variants.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for c in calls:
            first = next((S[k] for k, x in enumerate(c["fr"]) if x[2] is not None and x[2] >= 0.05), "?")
            df, dr = c["strand"][0], c["strand"][1]
            o.write("\t".join(map(str, [c["org"], c["ctg"], c["pos"], c["ref"], c["alt"], f"{c['anc']}>{c['der']}", c["cls"], c["direction"],
                                        c["filter"], first, c["tag"], c["product"], c["effect"], c["known"]]
                                  + [f"{x[0]}/{x[1]}" for x in c["fr"]] + ["NA" if x[2] is None else f"{x[2]:.3f}" for x in c["fr"]]
                                  + [c["rna"], f"{c['lowmq']:.2f}", "NA" if df is None else f"{df}/{dr}",
                                     "" if not c.get("hpbg") else ",".join("NA" if x is None else f"{x:.4f}" for x in c["hpbg"][0]) + f" ({c['hpbg'][1]} runs)"])) + "\n")
    for org in orgs:
        with open(f"{a.out}_{org}.vcf", "w") as o:
            o.write("##fileformat=VCFv4.2\n##source=08_call_variants.py\n##reference=" + a.ref + "\n")
            for c in sorted(org_of):
                if org_of[c] == org and c in seqs: o.write(f"##contig=<ID={c},length={len(seqs[c])}>\n")
            for fid, desc in [("lowMQ", f">{a.maxlowmq:.0%} of reads MAPQ<20"), ("highDP", f"depth >{a.maxdepthx}x median"),
                              ("strand", "derived allele on one strand or strand bias"), ("edge", f"within {a.edge} bp of a contig end"),
                              ("cluster", ">=3 candidates within 100 bp"),
                              ("hpNoise", "homopolymer indel below 3x the slippage rate at other runs of the same length"),
                              ("presentFirst", "indel: derived allele >5% of spanning reads in the first dataset (not new)"),
                              ("rareLast", "indel: derived allele <20% or <5 spanning reads in the last dataset"),
                              ("lowDepthSpanning", "indel: <15 reads span the run in the first or last dataset")]:
                o.write(f'##FILTER=<ID={fid},Description="{desc}">\n')
            o.write('##INFO=<ID=DER,Number=1,Type=String,Description="derived allele: REF or ALT">\n'
                    '##INFO=<ID=CLASS,Number=1,Type=String,Description="transition, transversion, insN/delN[_homopolymerL]">\n'
                    '##INFO=<ID=DF,Number=.,Type=Float,Description="derived-allele frequency per dataset (oldest first)">\n'
                    '##INFO=<ID=EFF,Number=1,Type=String,Description="locus|effect">\n')
            o.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for c in calls:
                if c["org"] != org: continue
                dfs = ",".join("." if x[2] is None else f"{x[2]:.3f}" for x in c["fr"])
                eff = (c["tag"] + "|" + c["effect"]).replace(";", ",").replace(" ", "_")
                o.write(f"{c['ctg']}\t{c['pos']}\t.\t{c['ref']}\t{c['alt']}\t.\t{c['filter']}\t"
                        f"DER={'ALT' if c['direction'] == 'ref=ancestral' else 'REF'};CLASS={c['cls']};DF={dfs};EFF={eff}\n")
    with open(a.out + "_spectrum.tsv", "w") as o:
        o.write("organism\tPASS_calls\ttransitions\ttransversions\thomopolymer_indels\tother_indels\tTs_fraction_of_SNVs\tfailed_filters\tassessable\n")
        for org in orgs:
            H = [c for c in calls if c["org"] == org and c["filter"] == "PASS"]
            ts = sum(c["cls"] == "transition" for c in H); tv = sum(c["cls"] == "transversion" for c in H)
            hp = sum("homopolymer" in c["cls"] for c in H); oi = sum(c["indel"] and "homopolymer" not in c["cls"] for c in H)
            nf = sum(1 for c in calls if c["org"] == org and c["filter"] != "PASS")
            ok = med.get((org, 0), 0) >= a.mindepth and med.get((org, NS - 1), 0) >= a.mindepth
            o.write(f"{org}\t{len(H)}\t{ts}\t{tv}\t{hp}\t{oi}\t" + (f"{ts / (ts + tv):.2f}" if ts + tv else "NA") + f"\t{nf}\t{'yes' if ok else 'NO (low depth)'}\n")
    print(open(a.out + "_spectrum.tsv").read())
    print(f"All calls with filters: {a.out}_variants.tsv ; VCF per organism: {a.out}_<org>.vcf")

if __name__ == "__main__":
    main()
