#!/usr/bin/env python3
"""
fig4_IS4_mutation_frequencies.py - Figure 4: rise of a lineage carrying the IS4 insertion in mutS.

(A) Frequency of the IS4 insertion and of the 12 new mutations in the 2014, 2017 and 2023 DNA
    reads, coloured by 2023 frequency group; dashed line = the pre-existing variant (2,398,457).
(B) 2023 frequency of each mutation with its 95% Clopper-Pearson interval; vertical line and
    grey band = insertion frequency and its interval.
(C) Mutation spectrum.

All numbers are read from Supplementary Table S1 (sheets S1.21_LK_mutations and
S1.22_mutS_IS4_junctions), so the figure always matches the table.

Usage:
  python3 fig4_IS4_mutation_frequencies.py --table Supplementary_Table_S1.xlsx --out Figure4

Outputs <out>.png, <out>.pdf and <out>.tif (300 dpi). Requires matplotlib, scipy, openpyxl.
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import openpyxl
from scipy.stats import beta

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#888780", "axes.linewidth": 0.8,
                     "xtick.color": "#3d3d3a", "ytick.color": "#3d3d3a"})
INK, MUTED, GRID, BAND = "#1a1a19", "#6b6a64", "#e6e5e0", "#e3e3dd"
COL = {"A": "#2a78d6", "B": "#eb6834", "C": "#1baf7a"}
MRK = {"A": "o", "B": "s", "C": "^"}
YRS = [2014, 2017, 2023]


def cp(k, n):
    lo = 0.0 if k == 0 else beta.ppf(0.025, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(0.975, k + 1, n - k)
    return lo, hi


def kn(cell):
    k, n = str(cell).split("/")
    return int(k), int(n)


def find_header(rows, first):
    for i, r in enumerate(rows):
        if r and r[0] == first:
            return i
    raise SystemExit(f"Header starting with '{first}' not found")


def load(path):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    rows = list(wb["S1.21_LK_mutations"].iter_rows(values_only=True))
    h = find_header(rows, "#")
    H = rows[h]
    ix = {name: H.index(name) for name in H if name}
    muts = []
    for r in rows[h + 1:]:
        if r[0] is None or not str(r[0]).isdigit():
            break
        muts.append(dict(pos=int(r[ix["LK_1 position (LNDBAELA)"]]),
                         g=str(r[ix["Group (2023 frequency)"]])[0],
                         cls=str(r[ix["Class"]]),
                         kn=[kn(r[ix[f"{y} reads (derived/total)"]]) for y in YRS]))
    # pre-existing variant: first row of the "Earlier calls" block whose status says pre-existing
    pre = None
    j = find_header(rows, "Site")
    J = rows[j]
    for r in rows[j + 1:]:
        if r[0] and "pre-existing" in str(r[1]):
            pre = dict(pos=int(r[2]), kn=[kn(r[J.index(f"{y} reads")]) for y in YRS])
            break
    # IS4 insertion: mean junction reads / (mean junction reads + intact-site reads)
    rows = list(wb["S1.22_mutS_IS4_junctions"].iter_rows(values_only=True))
    h = find_header(rows, "Dataset")
    is4 = {}
    for r in rows[h + 1:]:
        if r[0] is None:
            break
        if str(r[0]).strip() in map(str, YRS):
            j = (float(r[1]) + float(r[2])) / 2
            is4[int(r[0])] = (j, j + float(r[3]))
    return muts, pre, [is4[y] for y in YRS]


def spectrum(muts):
    cats = {"Transitions": 0, "±1-nt\nhomopolymer slips": 0, "Other 1-nt\ndeletion": 0, "Transversions": 0}
    for m in muts:
        c = m["cls"].lower()
        if "transition" in c:
            cats["Transitions"] += 1
        elif "transversion" in c:
            cats["Transversions"] += 1
        elif "homopolymer" in c:
            cats["±1-nt\nhomopolymer slips"] += 1
        else:
            cats["Other 1-nt\ndeletion"] += 1
    return cats


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True, help="Supplementary_Table_S1.xlsx")
    ap.add_argument("--out", default="Figure4_IS4_mutation_frequencies")
    a = ap.parse_args()
    muts, pre, is4 = load(a.table)
    is4f = [j / t for j, t in is4]
    j23, t23 = is4[-1]
    is4_lo, is4_hi = cp(round(j23), round(t23))

    fig = plt.figure(figsize=(7.2, 3.6), dpi=300)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.25, 0.8], wspace=0.62)

    # A
    ax = fig.add_subplot(gs[0])
    for y in (0.25, 0.5, 0.75):
        ax.axhline(y, color=GRID, lw=0.6, zorder=0)
    for m in muts:
        f = [k / n for k, n in m["kn"]]
        ax.plot(YRS, f, color=COL[m["g"]], lw=1.0, alpha=0.75, marker=MRK[m["g"]], ms=3.2,
                mec="white", mew=0.5, zorder=2)
    if pre:
        f = [k / n for k, n in pre["kn"]]
        ax.plot(YRS, f, color="#9a9994", lw=1.1, ls="--", marker="D", ms=3, mec="white", mew=0.5, zorder=2)
        ax.annotate(f"{pre['pos']:,}\n(pre-existing)", xy=(2017, f[1]), xytext=(2015.0, 0.43),
                    fontsize=7, color=MUTED, arrowprops=dict(arrowstyle="-", color="#b4b2a9", lw=0.6))
    ax.plot(YRS, is4f, color=INK, lw=2.2, marker="o", ms=4.5, mec="white", mew=0.7, zorder=3)
    ax.annotate("IS4 in mutS", xy=(2023, is4f[-1]), xytext=(2020.6, 0.76), fontsize=7.5,
                arrowprops=dict(arrowstyle="-", color=INK, lw=0.6))
    for g in "ABC":
        fs = [m["kn"][2][0] / m["kn"][2][1] for m in muts if m["g"] == g]
        ax.text(2023.7, sum(fs) / len(fs), g, fontsize=9, fontweight="bold", va="center", color=INK)
    ax.set_xticks(YRS)
    ax.set_xlim(2013.3, 2024.6)
    ax.set_ylim(-0.02, 0.8)
    ax.set_ylabel("Frequency in DNA reads")
    ax.set_title("A  Rise of the IS4 lineage", loc="left", fontweight="bold", fontsize=9.5)

    # B
    ax = fig.add_subplot(gs[1])
    ax.axvspan(is4_lo, is4_hi, color=BAND, zorder=0)
    ax.axvline(is4f[-1], color=INK, lw=1.0, zorder=1)
    ax.text(is4f[-1] + 0.02, len(muts) - 0.6, "IS4", fontsize=7, color=INK)
    order = sorted(muts, key=lambda m: ("ABC".index(m["g"]), -m["kn"][2][0] / m["kn"][2][1]))
    for i, m in enumerate(order):
        k, n = m["kn"][2]
        f = k / n
        lo, hi = cp(k, n)
        ax.plot([lo, hi], [i, i], color=COL[m["g"]], lw=1.4, solid_capstyle="round", zorder=2)
        ax.plot(f, i, marker=MRK[m["g"]], color=COL[m["g"]], ms=4, mec="white", mew=0.5, zorder=3)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{m['pos']:,}" for m in order], fontsize=7, color=MUTED)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.invert_yaxis()
    ax.set_xlim(0, 0.85)
    ax.set_xticks([0, 0.25, 0.5, 0.75])
    ax.set_xlabel("2023 frequency (95% CI)")
    ax.set_title("B  Three frequency groups", loc="left", fontweight="bold", fontsize=9.5)

    # C
    ax = fig.add_subplot(gs[2])
    sp = spectrum(muts)
    labs = list(sp)
    ax.barh(range(len(labs)), [sp[l] for l in labs], color="#5f5e5a", height=0.55)
    for i, l in enumerate(labs):
        ax.text(sp[l] + 0.2, i, str(sp[l]), va="center", fontsize=8)
    ax.set_yticks(range(len(labs)))
    ax.set_yticklabels(labs, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlim(0, max(sp.values()) + 1.5)
    ax.set_xticks(range(0, max(sp.values()) + 1, 2))
    ax.set_xlabel("New mutations")
    ax.set_title("C  Spectrum", loc="left", fontweight="bold", fontsize=9.5)

    def rng(g):
        fs = [m["kn"][2][0] / m["kn"][2][1] for m in muts if m["g"] == g]
        return f"{min(fs):.2f}–{max(fs):.2f}"
    handles = [Line2D([], [], color=INK, lw=2.2, marker="o", ms=4.5, mec="white", label="IS4 insertion in mutS"),
               Line2D([], [], color=COL["A"], marker="o", ms=4, label=f"Group A (same as IS4, {rng('A')})"),
               Line2D([], [], color=COL["B"], marker="s", ms=4, label=f"Group B ({rng('B')})"),
               Line2D([], [], color=COL["C"], marker="^", ms=4, label=f"Group C ({rng('C')})"),
               Line2D([], [], color="#9a9994", ls="--", marker="D", ms=3.5, label="Pre-existing variant")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=7.5,
               bbox_to_anchor=(0.5, -0.06))
    fig.subplots_adjust(bottom=0.27)
    for ext in ("png", "pdf"):
        fig.savefig(f"{a.out}.{ext}", dpi=300, bbox_inches="tight")
    fig.savefig(f"{a.out}.tif", dpi=300, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    print(f"IS4 frequency 2014/2017/2023: {', '.join(f'{x:.3f}' for x in is4f)} "
          f"(2023 95% CI {is4_lo:.2f}–{is4_hi:.2f}); spectrum: {dict((k.replace(chr(10), ' '), v) for k, v in sp.items())}")


if __name__ == "__main__":
    main()
