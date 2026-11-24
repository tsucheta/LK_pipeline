#!/usr/bin/env python3
"""
compare_busco.py - compare BUSCO runs side by side (Supplementary Table S3.3).

Reads the full_table.tsv of each BUSCO run (BUSCO v5), prints one summary line per run
(C, S, D, F, M, n) and writes a table of every BUSCO whose status differs between runs,
e.g. the two genes missing from the uncorrected LK_1 assembly.

Usage:
  python3 compare_busco.py --runs LK_16=busco_LK16 LK_3=busco_LK3 \
      LK_1_raw=busco_LK1_raw LK_1_corrected=busco_LK1_corr --out busco_differences.tsv

Each run is NAME=PATH, where PATH is a BUSCO output folder or a full_table.tsv file.
"""
import argparse
import glob
import os
import sys

ORDER = {"Complete": 0, "Duplicated": 1, "Fragmented": 2, "Missing": 3}
SHORT = {"Complete": "S", "Duplicated": "D", "Fragmented": "F", "Missing": "M"}


def find_table(path):
    if os.path.isfile(path):
        return path
    hits = glob.glob(os.path.join(path, "**", "full_table.tsv"), recursive=True)
    if not hits:
        sys.exit(f"No full_table.tsv found under {path}")
    if len(hits) > 1:
        print(f"Note: several full_table.tsv under {path}; using {hits[0]}", file=sys.stderr)
    return hits[0]


def read_table(path):
    status, desc, genes = {}, {}, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            bid, st = f[0], f[1]
            # Duplicated BUSCOs appear on several lines; keep the strongest status
            if bid not in status or ORDER[st] < ORDER[status[bid]]:
                status[bid] = st
            if st != "Missing" and len(f) > 2:
                genes.setdefault(bid, []).append(f[2])
            if len(f) > 6 and f[-1] and not f[-1].replace(".", "").isdigit():
                desc[bid] = f[-1]
    return status, desc, genes


def summary(status):
    n = len(status)
    c = {k: sum(1 for v in status.values() if v == k) for k in ORDER}
    s, d, f, m = c["Complete"], c["Duplicated"], c["Fragmented"], c["Missing"]
    pct = lambda x: 100 * x / n if n else 0
    return (f"C:{s + d} [S:{s}, D:{d}], F:{f}, M:{m}, n:{n}  "
            f"(C {pct(s + d):.1f}% [S {pct(s):.1f}%, D {pct(d):.1f}%], F {pct(f):.1f}%, M {pct(m):.1f}%)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", nargs="+", required=True, help="NAME=PATH pairs")
    ap.add_argument("--out", default="busco_differences.tsv")
    a = ap.parse_args()

    runs = []
    for item in a.runs:
        if "=" not in item:
            sys.exit(f"Use NAME=PATH, got: {item}")
        name, path = item.split("=", 1)
        runs.append((name, *read_table(find_table(path))))

    print("Run\tSummary")
    for name, status, _, _ in runs:
        print(f"{name}\t{summary(status)}")

    all_ids = sorted(set().union(*[r[1].keys() for r in runs]))
    desc = {}
    for _, _, d, _ in runs:
        desc.update(d)
    n_diff = 0
    with open(a.out, "w") as out:
        out.write("BUSCO_id\t" + "\t".join(r[0] for r in runs) + "\t"
                  + "\t".join(f"{r[0]}_gene" for r in runs) + "\tDescription\n")
        for bid in all_ids:
            sts = [r[1].get(bid, "Missing") for r in runs]
            if len(set(sts)) == 1:
                continue
            n_diff += 1
            gns = [",".join(sorted(set(r[3].get(bid, [])))) or "-" for r in runs]
            out.write(f"{bid}\t" + "\t".join(SHORT[s] for s in sts) + "\t"
                      + "\t".join(gns) + f"\t{desc.get(bid, '')}\n")
    print(f"\n{n_diff} BUSCOs differ between runs; written to {a.out}")


if __name__ == "__main__":
    main()
