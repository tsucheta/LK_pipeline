#!/usr/bin/env python3
"""
07_IS_timeline.py - genome-wide IS insertion timeline for every organism.

For every transposase gene annotated in the NEW assemblies it asks:
  1. Is this IS copy NEW relative to the OLD assembly (the flanking host DNA is contiguous in the old genome)?
  2. In which DNA dataset does it first appear, and at what frequency (reads across the IS junctions in the
     NEW-reference BAM vs reads across the empty site in the OLD-reference BAM)?
  3. Which genes does it disrupt or sit next to, and is any of them a DNA-repair gene (a mutator candidate)?

Needs python3 (standard library) + samtools. BAMs from 00_align.sh: bams/DNA/<sample>_vs_LK1.bam / _vs_LK16.bam

Usage:
  python3 07_IS_timeline.py --new refs/index/NEW.fna --old refs/index/OLD.fna \
      --gff LK=refs/LK_1.gff Phaeo=refs/phaeo.gff Ocean=refs/ocean.gff \
            Inquil=refs/inqui.gff Nitrat=refs/LK_bac_polished.gff \
      --samples DNA2014_PE DNA2017_PE DNA2023_PE
An organism whose contigs are identical in --new and --old (Inquilinaceae and Nitratireductor, unless INQUIL_OLD
was set in config.sh) cannot be dated: its IS loci are reported as "not testable" instead of "present in old".
Outputs: results/07_IS_timeline.tsv (one row per IS copy) and results/07_IS_summary.tsv (per organism).
Options: --window 600 (host DNA taken beyond each transposase CDS), --flank 20, --minmapq 10
"""
import argparse, os, re, subprocess, sys
from urllib.parse import unquote

REPAIR = re.compile(r"\b(mut[SLHYMT]|uvr[ABCD]|dnaQ|rec[AFGNOQR]|ruv[ABC]|mfd|ung|nth|nei|alkA|ada|xth|nfo|radA|"
                    r"umu[CD]|dinB|polA|lig[AD]|sbc[BCD]|addA|addB|recBCD)\b|mismatch repair|DNA repair|"
                    r"exonuclease|endonuclease (III|IV|V|VIII)|glycosylase|photolyase|recombinase RecA", re.I)
IS_RX = re.compile(r"transposase|insertion (element|sequence)|\bIS\d+|\bIS[A-Z]?\d* family", re.I)

def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True).stdout

def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name: seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else: buf.append(line.strip())
    if name: seqs[name] = "".join(buf).upper()
    return seqs

def rc(s): return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]

def parse_gff(path):
    feats = []
    for l in open(path):
        if l.startswith("#") or l.count("\t") < 8: continue
        f = l.rstrip("\n").split("\t")
        if f[2] not in ("CDS", "gene", "pseudogene", "tRNA", "rRNA"): continue
        a = {k: unquote(v) for k, v in (x.split("=", 1) for x in f[8].split(";") if "=" in x)}
        feats.append(dict(ctg=f[0], s=int(f[3]), e=int(f[4]), type=f[2], tag=a.get("locus_tag", a.get("ID", "")),
                          gene=a.get("gene", ""), product=a.get("product", a.get("Note", ""))))
    return feats

def walk(a, b):
    """length of the shared start of a and b, tolerating isolated mismatches (stop when a 20-bp window has >=6)"""
    n, i = min(len(a), len(b)), 0
    while i < n:
        w = sum(1 for k in range(i, min(i + 20, n)) if a[k] != b[k])
        if w >= 6: break
        i += 1
    # the window test stops up to ~15 bp early; extend through the exact matches up to the real junction
    while i < n and a[i] == b[i]: i += 1
    return i

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--new", required=True); ap.add_argument("--old", required=True)
    ap.add_argument("--gff", nargs="+", required=True, help="ORG=path.gff (contig names as in --new)")
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--window", type=int, default=600); ap.add_argument("--flank", type=int, default=20)
    ap.add_argument("--minmapq", type=int, default=10); ap.add_argument("--out", default="results/07")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)

    print("Loading references ..."); new, old = read_fasta(a.new), read_fasta(a.old)
    oldidx = []
    for n, s in old.items(): oldidx += [(n, "+", s), (n, "-", rc(s))]

    bam_refs = {}
    def resolve(bam, name):
        if bam not in bam_refs:
            bam_refs[bam] = set(l.split("\t")[1][3:] for l in sh(f"samtools view -H {bam}").splitlines() if l.startswith("@SQ"))
        for c in (name, "NZ_" + name, name[3:] if name.startswith("NZ_") else ""):
            if c and c in bam_refs[bam]: return c
        return None

    def spanning(bam, ctg, x, y):
        n = 0
        for line in sh(f"samtools view -F 0xF04 -q {a.minmapq} {bam} '{ctg}:{x}-{y}'").splitlines():
            f = line.split("\t"); pos = int(f[3]); ops = re.findall(r"(\d+)([MIDNSHP=X])", f[5])
            span = sum(int(l) for l, o in ops if o in "MDN=X"); indel = sum(int(l) for l, o in ops if o in "ID")
            if pos <= x - a.flank and pos + span - 1 >= y + a.flank and indel <= 3: n += 1
        return n

    for smp in a.samples:
        for suf in ("LK1", "LK16"):
            if not os.path.exists(f"bams/DNA/{smp}_vs_{suf}.bam"): sys.exit(f"ERROR: bams/DNA/{smp}_vs_{suf}.bam missing")

    W, A = a.window, 60
    rows, summary = [], {}
    for spec in a.gff:
        org, path = spec.split("=", 1)
        feats = parse_gff(path)
        loci = sorted([f for f in feats if f["type"] == "CDS" and IS_RX.search(f["product"])], key=lambda f: (f["ctg"], f["s"]))
        merged = []
        for f in loci:   # merge transposase CDSs within 200 bp (multi-ORF IS)
            if merged and merged[-1]["ctg"] == f["ctg"] and f["s"] - merged[-1]["e"] < 200:
                merged[-1]["e"] = max(merged[-1]["e"], f["e"]); merged[-1]["tag"] += "," + f["tag"]
            else: merged.append(dict(f))
        print(f"{org}: {len(merged)} transposase loci")
        S = summary.setdefault(org, {"loci": len(merged), "new": 0, "in_old": 0, "unresolved": 0, "not_testable": 0, "first_seen": {}, "repair_hits": []})
        ctgs = {L["ctg"] for L in merged}
        missing = [c for c in ctgs if c not in new]
        if missing:
            print(f"  WARNING: {len(missing)} {org} GFF contigs are not in --new (e.g. {missing[0]}); GFF and FASTA must use the same names")
        same = bool(ctgs) and all(c in old and c in new and old[c] == new[c] for c in ctgs - set(missing))
        if same:
            print(f"  {org}: the same assembly is in --new and --old, so insertions cannot be dated (not testable)")
        for L in merged:
            ctg = L["ctg"]
            if ctg not in new: S["unresolved"] += 1; continue
            if same:
                rows.append(dict(org=org, locus=L["tag"], ctg=ctg, cds=f"{L['s']}-{L['e']}", product=L["product"],
                                 status="not testable (no older assembly of this organism in --old)"))
                S["not_testable"] += 1; continue
            seq = new[ctg]; rs, re_ = L["s"] - W, L["e"] + W
            if rs - A < 1 or re_ + A > len(seq): S["unresolved"] += 1; continue
            reg = seq[rs - 1 - A: re_ + A]
            la, ra = reg[:A], reg[-A:]
            hits = []
            for n, st, s in oldidx:
                i = s.find(la)
                while i >= 0:
                    j = s.find(ra, i + A)
                    if j >= 0 and j - i < len(reg) + 5000: hits.append((n, st, s, i, j))
                    i = s.find(la, i + 1)
            base = dict(org=org, locus=L["tag"], ctg=ctg, cds=f"{L['s']}-{L['e']}", product=L["product"])
            if len(hits) != 1:
                rows.append({**base, "status": "unresolved (flanks not unique/absent in old)"}); S["unresolved"] += 1; continue
            n, st, s, i, j = hits[0]
            oreg = s[i: j + A]
            ins = len(reg) - len(oreg)
            if ins < 300:
                rows.append({**base, "status": "present in old assembly", "ins_len": ins}); S["in_old"] += 1; continue
            p = walk(reg, oreg); q = walk(reg[::-1], oreg[::-1])
            J1 = rs - A + p                     # first inserted base (NEW coords)
            J2 = rs - A + len(reg) - q - 1      # last inserted base
            T = (i + p) if st == "+" else (len(s) - (i + p))
            # genes hit
            dis = [f for f in feats if f["ctg"] == ctg and f["type"] in ("CDS", "pseudogene") and f["e"] >= J1 - 30 and f["s"] <= J2 + 30
                   and not IS_RX.search(f["product"])]
            adj = [f for f in feats if f["ctg"] == ctg and f["type"] == "CDS" and not IS_RX.search(f["product"])
                   and (0 < J1 - f["e"] <= 300 or 0 < f["s"] - J2 <= 300)]
            lab = lambda fs: "; ".join(f"{f['tag']}:{f['gene'] or ''}:{f['product']}" for f in fs)[:300]
            rep = list({f["tag"]: f for f in dis + adj if REPAIR.search(f["product"] + " " + f["gene"])}.values())
            res = {**base, "status": "NEW vs old assembly", "ins_len": J2 - J1 + 1, "new_coords": f"{ctg}:{J1}-{J2}",
                   "old_site": f"{n}:{T}", "disrupted": lab(dis), "adjacent": lab(adj), "repair_gene": lab(rep) or ""}
            first = None
            for smp in a.samples:
                bn, bo = f"bams/DNA/{smp}_vs_LK1.bam", f"bams/DNA/{smp}_vs_LK16.bam"
                cn, co = resolve(bn, ctg), resolve(bo, n)
                if not cn or not co: res[smp] = "NA (contig not in BAM)"; continue
                l_ = spanning(bn, cn, J1 - 1, J1); r_ = spanning(bn, cn, J2, J2 + 1); e_ = spanning(bo, co, T, T + 1)
                isr = (l_ + r_) / 2; tot = isr + e_
                fq = isr / tot if tot >= 10 else None
                res[smp] = f"{l_}/{r_}/{e_}=" + ("NA" if fq is None else f"{fq:.2f}")
                if fq is not None and fq >= 0.05 and first is None: first = smp
            res["first_seen"] = first or "not detected"
            rows.append(res); S["new"] += 1
            S["first_seen"][res["first_seen"]] = S["first_seen"].get(res["first_seen"], 0) + 1
            if rep: S["repair_hits"].append(f"{L['tag']}->{lab(rep)} (first {res['first_seen']})")
            print(f"  {org} {L['tag']}: NEW insertion {ins} bp, first seen {res['first_seen']}" + (f"  ** repair gene: {lab(rep)}" if rep else ""))

    cols = ["org", "locus", "product", "status", "ins_len", "new_coords", "old_site", "disrupted", "adjacent", "repair_gene"] + a.samples + ["first_seen", "cds", "ctg"]
    with open(a.out + "_IS_timeline.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows: o.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    with open(a.out + "_IS_summary.tsv", "w") as o:
        o.write("organism\ttransposase_loci\tnew_insertions\tpresent_in_old\tunresolved\tnot_testable\t" + "\t".join(f"first_seen_{s}" for s in a.samples) + "\tnot_detected\trepair_gene_insertions\n")
        for org, S in summary.items():
            o.write(f"{org}\t{S['loci']}\t{S['new']}\t{S['in_old']}\t{S['unresolved']}\t{S['not_testable']}\t" +
                    "\t".join(str(S['first_seen'].get(s, 0)) for s in a.samples) + f"\t{S['first_seen'].get('not detected', 0)}\t" +
                    ("; ".join(S["repair_hits"]) or "none") + "\n")
    print("\n" + open(a.out + "_IS_summary.tsv").read())
    print(f"Per-insertion table: {a.out}_IS_timeline.tsv  (cells = IS_left/IS_right/empty=frequency)")

if __name__ == "__main__":
    main()
