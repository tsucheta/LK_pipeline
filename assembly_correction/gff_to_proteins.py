#!/usr/bin/env python3
"""
gff_to_proteins.py - translate the CDS features of a GFF3 annotation into a protein FASTA.

Used to make the LK_1 protein set for BUSCO protein mode (Supplementary Table S3.3).
Pseudogenes are skipped. Split CDS (several CDS lines with the same ID / locus tag) are
joined in genome order before translation. No external libraries are needed.

Usage:
  python3 gff_to_proteins.py --gff LK_1_corrected.gff --fasta LK_1_corrected.fna \
      --out LK_1_corrected.faa [--table 11] [--keep-pseudo]
"""
import argparse
import sys
from collections import OrderedDict
from urllib.parse import unquote

BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
AA11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
TABLE = dict(zip(CODONS, AA11))
STARTS11 = {"TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"}
COMP = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs


def parse_attrs(field):
    attrs = {}
    for item in field.strip().strip(";").split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k.strip()] = unquote(v.strip())
    return attrs


def translate(nt, table_starts=True):
    prot = []
    for i in range(0, len(nt) - 2, 3):
        prot.append(TABLE.get(nt[i:i + 3], "X"))
    p = "".join(prot)
    # alternative start codons are read as Met (translation table 11)
    if table_starts and nt[:3] in STARTS11 and p:
        p = "M" + p[1:]
    if p.endswith("*"):
        p = p[:-1]
    return p


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--fasta", required=True, help="genome FASTA; sequence names must match GFF column 1")
    ap.add_argument("--out", required=True)
    ap.add_argument("--table", type=int, default=11, help="only table 11 is implemented")
    ap.add_argument("--keep-pseudo", action="store_true", help="also translate pseudogene CDS")
    a = ap.parse_args()
    if a.table != 11:
        sys.exit("Only translation table 11 is implemented.")

    genome = read_fasta(a.fasta)
    pseudo_genes = set()
    cds = OrderedDict()   # key -> dict(seqid, strand, parts[], phase, tag, product)

    with open(a.gff) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            seqid, _, ftype, start, end, _, strand, phase, attr = f
            at = parse_attrs(attr)
            is_pseudo = (at.get("pseudo", "").lower() == "true"
                         or at.get("gene_biotype", "") == "pseudogene"
                         or ftype == "pseudogene")
            if ftype in ("gene", "pseudogene"):
                if is_pseudo and "ID" in at:
                    pseudo_genes.add(at["ID"])
                continue
            if ftype != "CDS":
                continue
            if is_pseudo or at.get("Parent") in pseudo_genes:
                if not a.keep_pseudo:
                    continue
            tag = at.get("locus_tag") or at.get("ID") or at.get("Parent")
            key = at.get("ID") or tag
            rec = cds.setdefault(key, {"seqid": seqid, "strand": strand, "parts": [],
                                       "tag": tag, "product": at.get("product", ""),
                                       "phase": None})
            rec["parts"].append((int(start), int(end), phase))

    n_out = n_skip = n_internal_stop = 0
    with open(a.out, "w") as out:
        for key, r in cds.items():
            if r["seqid"] not in genome:
                n_skip += 1
                continue
            parts = sorted(r["parts"])
            nt = "".join(genome[r["seqid"]][s - 1:e] for s, e, _ in parts)
            if r["strand"] == "-":
                nt = nt.translate(COMP)[::-1]
                ph = parts[-1][2]
            else:
                ph = parts[0][2]
            if ph in ("1", "2"):
                nt = nt[int(ph):]
            prot = translate(nt)
            if "*" in prot:
                n_internal_stop += 1
            out.write(f">{r['tag']} {r['product']}\n")
            for i in range(0, len(prot), 60):
                out.write(prot[i:i + 60] + "\n")
            n_out += 1

    print(f"{n_out} proteins written to {a.out}", file=sys.stderr)
    if n_skip:
        print(f"{n_skip} CDS skipped: sequence name not in FASTA", file=sys.stderr)
    if n_internal_stop:
        print(f"WARNING: {n_internal_stop} proteins contain internal stop codons (*)", file=sys.stderr)


if __name__ == "__main__":
    main()
