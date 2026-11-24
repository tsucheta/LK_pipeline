#!/usr/bin/env python3
"""
06_mutS_IS_timing.py - when did the IS4 insertion into mutS happen?

For each DNA dataset it counts
  IS reads    : reads mapped to the NEW assembly (LK_1) that run straight across an IS4/mutS junction
  empty reads : reads mapped to the OLD assembly (LK_16) that run straight across the insertion site
and reports  IS frequency = mean(left, right IS reads) / (that + empty reads).

The insertion point is found automatically: the mutS sequence flanking the IS in LK_1
(end of pgaptmp_000416 and start of pgaptmp_000418) is searched in the old reference,
and the two sequences are compared base by base to locate the exact junctions.

Needs only python3 (standard library) and samtools.

Usage:
  python3 06_mutS_IS_timing.py refs/index/NEW.fna refs/index/OLD.fna \
      DNA2014_PE DNA2017_PE DNA2023_PE            # sample names; BAMs are bams/DNA/<name>_vs_LK1.bam / _vs_LK16.bam
Options (env vars): NEW_CONTIG (default LNDBAELA), REGION_START=562850, REGION_END=564650, FLANK=20, MINMAPQ=10,
                    VARIANTS=results/08_variants.tsv (PASS calls of ORG=LK from step 08 form the timeline)
"""
import os, sys, subprocess, re, math

def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True).stdout

def read_fasta(path, want=None):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None and (want is None or name in want): seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None and (want is None or name in want): seqs[name] = "".join(buf).upper()
    return seqs

def rc(s): return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]

def main():
    if len(sys.argv) < 4: sys.exit(__doc__)
    new_fa, old_fa, samples = sys.argv[1], sys.argv[2], sys.argv[3:]
    ctg = os.environ.get("NEW_CONTIG", "LNDBAELA")
    rs, re_ = int(os.environ.get("REGION_START", 562850)), int(os.environ.get("REGION_END", 564650))
    FL, MQ = int(os.environ.get("FLANK", 20)), int(os.environ.get("MINMAPQ", 10))

    new = read_fasta(new_fa, {ctg})
    if ctg not in new: sys.exit(f"ERROR: contig {ctg} not in {new_fa}. Set NEW_CONTIG=<name> (see: grep '>' {new_fa} | head)")
    reg = new[ctg][rs - 1:re_]
    left, right = reg[:60], reg[-60:]           # 60 bp of mutS on each side, well outside the IS

    old = read_fasta(old_fa)
    hit = None
    for name, s in old.items():
        for strand, q in (("+", s), ("-", None)):
            if strand == "-": q = rc(s)
            i, j = q.find(left), q.find(right)
            if i >= 0 and j > i:
                hit = (name, strand, q, i, j); break
        if hit: break
    if not hit: sys.exit("ERROR: could not find the mutS flanks in the old reference (sequence differs?)")
    oname, strand, q, i, j = hit
    oreg = q[i:j + 60]
    print(f"Old reference: mutS region found on {oname} ({strand} strand), length {len(oreg)} bp vs {len(reg)} bp in LK_1")
    ins = len(reg) - len(oreg)
    if ins < 500:
        print(f"!! Old reference region is only {ins} bp shorter: the IS is ALREADY present in the old assembly (or absent from both).")
        print("   -> the insertion predates the old (LK_16) assembly; reads below still show its frequency per year.")
    # exact junctions: longest common prefix / suffix
    p = 0
    while p < len(oreg) and reg[p] == oreg[p]: p += 1
    s_ = 0
    while s_ < len(oreg) - p and reg[-1 - s_] == oreg[-1 - s_]: s_ += 1
    J1 = rs + p                    # first inserted base in LK_1 (left junction is between J1-1 and J1)
    J2 = rs + len(reg) - s_ - 1    # last inserted base in LK_1 (right junction between J2 and J2+1)
    tsd = p + s_ - len(oreg)       # >0 = target-site duplication length
    # insertion point in old reference coordinates (forward strand of the old contig)
    T_local = p                    # junction between oreg[p-1] and oreg[p]
    if strand == "+": T = i + T_local          # 0-based index of first base after the junction -> 1-based = T+1
    else: T = len(q) - (i + T_local)           # map back to forward coordinate
    print(f"LK_1 insertion: {ctg}:{J1}-{J2} ({J2-J1+1} bp inserted; target-site duplication ~{max(tsd,0)} bp)")
    print(f"Old insertion point: {oname} between {T} and {T+1}\n")

    def spanning(bam, contig, a, b):
        """reads aligned continuously from <= a-FL to >= b+FL (no clipping at the junction, indels <= 3 bp)"""
        n = 0
        out = sh(f"samtools view -F 0xF04 -q {MQ} {bam} {contig}:{a}-{b}")
        for line in out.splitlines():
            f = line.split("\t"); pos, cig = int(f[3]), f[5]
            ops = re.findall(r"(\d+)([MIDNSHP=X])", cig)
            span = sum(int(l) for l, o in ops if o in "MDN=X")
            indel = sum(int(l) for l, o in ops if o in "ID")
            if pos <= a - FL and pos + span - 1 >= b + FL and indel <= 3: n += 1
        return n

    os.makedirs("results", exist_ok=True)
    rows, isf, isc = [], {}, {}
    def resolve(bam, name):
        """contig name as written in the BAM header (handles NZ_ prefix differences)"""
        refs = [l.split("\t")[1][3:] for l in sh(f"samtools view -H {bam}").splitlines() if l.startswith("@SQ")]
        for cand in (os.environ.get("OLD_CONTIG", ""), name, "NZ_" + name, name[3:] if name.startswith("NZ_") else ""):
            if cand and cand in refs: return cand
        sys.exit(f"ERROR: {name} not in the header of {bam}. Names there: {', '.join(r for r in refs if 'QVFV' in r or name[-6:] in r)[:300]}")

    print("sample\tIS_left\tIS_right\tempty_site\tIS_frequency\tverdict")
    for smp in samples:
        bn, bo = f"bams/DNA/{smp}_vs_LK1.bam", f"bams/DNA/{smp}_vs_LK16.bam"
        for b in (bn, bo):
            if not os.path.exists(b): sys.exit(f"ERROR: missing {b}")
        oname = resolve(bo, oname)
        L = spanning(bn, ctg, J1 - 1, J1)
        R = spanning(bn, ctg, J2, J2 + 1)
        E = spanning(bo, oname, T, T + 1)
        isr = (L + R) / 2
        tot = isr + E
        if tot < 10: v, fq = "too few reads", "NA"
        else:
            f = isr / tot; fq = f"{f:.3f}"
            v = "absent" if f < 0.05 else ("fixed" if f > 0.95 else "polymorphic")
        print(f"{smp}\t{L}\t{R}\t{E}\t{fq}\t{v}")
        rows.append(f"{smp}\t{L}\t{R}\t{E}\t{fq}\t{v}"); isf[smp] = None if fq == "NA" else float(fq); isc[smp] = (L, R, E)
    with open(os.environ.get("OUTDIR", "results") + "/06_mutS_IS_timing.tsv", "w") as o:
        o.write("sample\tIS_left\tIS_right\tempty_site\tIS_frequency\tverdict\n" + "\n".join(rows) + "\n")

    # ---- combined timeline with the PASS point mutations from step 08 ----
    af = os.environ.get("VARIANTS", "results/08_variants.tsv")
    if not os.path.exists(af):
        print(f"\n(no {af}; run step 08 to get the combined timeline)"); return
    lines = [l.rstrip("\n").split("\t") for l in open(af) if l.strip()]
    hdr, body = lines[0], lines[1:]
    H = {h: i for i, h in enumerate(hdr)}
    body = [r for r in body if r[H["organism"]] == os.environ.get("ORG", "LK") and r[H["FILTER"]] == "PASS"]

    def wilson(k, n, z=1.96):
        """95% Wilson interval for k of n"""
        if n == 0: return (0.0, 1.0)
        p = k / n; d = 1 + z * z / n; c = p + z * z / (2 * n); h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
        return (max(0.0, (c - h) / d), min(1.0, (c + h) / d))
    def detected(k, n):
        """present, not an error: >= 3 reads and the lower 95% bound above 1%"""
        return n > 0 and k >= 3 and wilson(k, n)[0] > 0.01

    # IS4 counts per dataset: insertion 'reads' = mean of the two junctions
    is_kn = {s_: (round((isc[s_][0] + isc[s_][1]) / 2), round((isc[s_][0] + isc[s_][1]) / 2) + isc[s_][2]) for s_ in samples}
    is_first = next((s_ for s_ in samples if detected(*is_kn[s_])), None)
    muts = []
    for r in body:
        kn = {}
        for s_ in samples:
            c = r[H.get(f"{s_}_derived/depth", -1)] if f"{s_}_derived/depth" in H else "0/0"
            k, n = (int(x) for x in c.split("/"))
            kn[s_] = (k, n)
        name = f"{r[H['contig']]}:{r[H['pos']]} {r[H['change(ancestral>derived)']]} ({r[H['effect']] or r[H['class']]})"
        muts.append(dict(name=name, kn=kn, first=next((s_ for s_ in samples if detected(*kn[s_])), None)))

    fmt = lambda k, n: "NA" if n == 0 else f"{k}/{n}={k / n:.2f}"
    out = ["event\t" + "\t".join(samples) + "\tfirst_detected\tlineage_vs_IS4_in_" + samples[-1]]
    out.append("mutS IS4 insertion\t" + "\t".join(fmt(*is_kn[s_]) for s_ in samples) + f"\t{is_first or 'not detected'}\t-")
    last = samples[-1]
    kI, nI = is_kn[last]
    for m in muts:
        k, n = m["kn"][last]
        # same lineage as the IS4 if the last-dataset frequencies do not differ (chi-square, 1 df, P > 0.05)
        tot = k + kI; N = n + nI
        if n and nI and 0 < tot < N:
            e = [n * tot / N, n * (N - tot) / N, nI * tot / N, nI * (N - tot) / N]
            chi = sum((o - x) ** 2 / x for o, x in zip([k, n - k, kI, nI - kI], e))
            m["same"] = chi < 3.841
        else:
            m["same"] = None
        m["lin"] = {True: "same frequency as IS4", False: ("lower" if n and k / n < kI / nI else "higher") + " than IS4", None: "NA"}[m["same"]]
        out.append(m["name"] + "\t" + "\t".join(fmt(*m["kn"][s_]) for s_ in samples) + f"\t{m['first'] or 'not detected'}\t{m['lin']}")
    with open(os.environ.get("OUTDIR", "results") + "/06_timeline.tsv", "w") as o: o.write("\n".join(out) + "\n")
    print("\nTimeline (derived reads/depth = frequency; detected = >=3 reads and lower 95% bound > 1%):")
    for l in out: print(l)

    print("\nInterpretation:")
    if is_first is None:
        print("- IS4 not detected in any dataset: check the junction coordinates."); return
    i_is = samples.index(is_first)
    same = [m for m in muts if m["same"]]
    other = [m for m in muts if m["same"] is False]
    print(f"- IS4 first detected in {is_first}. {len(muts)} PASS mutations: {len(same)} at the same {last} frequency as the IS4, "
          f"{len(other)} at a different frequency.")
    early = [m for m in muts if m["first"] and samples.index(m["first"]) < i_is]
    for grp, label in ((same, "same-frequency (candidate IS4-lineage) mutations"), (other, "other mutations")):
        if not grp: continue
        e = [m for m in grp if m["first"] and samples.index(m["first"]) < i_is]
        if not e:
            print(f"- {label}: none detected before {is_first} -> consistent with arising after, or together with, the IS4; "
                  "the order within that interval is not resolved by population frequencies.")
        else:
            print(f"- {label}: {len(e)} detected before the IS4:")
            for m in e:
                print(f"    {m['name']}: " + ", ".join(f"{s_} {fmt(*m['kn'][s_])}" for s_ in samples))
                # the IS4 may be present below detection: compare with its upper 95% bound in that dataset
                s0 = m["first"]; kk, nn = m["kn"][s0]; up = wilson(*is_kn[s0])[1]
                if kk / nn <= up:
                    print(f"      {s0}: {kk / nn:.3f} is within the IS4's upper 95% bound ({up:.3f}; {is_kn[s0][0]}/{is_kn[s0][1]} reads),"
                          " so the IS4 could already have been present at that level: no conflict.")
                else:
                    print(f"      {s0}: {kk / nn:.3f} exceeds the IS4's upper 95% bound ({up:.3f}): this mutation was more common than the IS4"
                          " could have been. If real (check IGV; homopolymer slippage can mimic it), it arose outside or before the IS4 lineage.")
            if not label.startswith("same"):
                print("  These do not share the IS4 frequency (a nested subclone, a different lineage or standing variation); "
                      "only a mutation above the IS4's upper bound in the earlier dataset argues against MutS loss coming first.")
    if same:
        hi = max(m["kn"][last][0] / m["kn"][last][1] for m in same if m["kn"][last][1])
        print(f"- In {last}: IS4 {kI / nI:.2f}; same-frequency mutations up to {hi:.2f} -> consistent with these mutations sitting in the IS4 lineage.")

if __name__ == "__main__":
    main()
