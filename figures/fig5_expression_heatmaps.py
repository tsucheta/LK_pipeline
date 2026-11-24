#!/usr/bin/env python3
"""
fig5_expression_heatmaps.py - Figure 5: expression of the genes hit by the new mutations.

(A) LK genes affected by the new mutations and the IS4 insertion in mutS (Table S1.25).
(B) Oceanicaulis sp. split genes (Table S1.27).
Cells show the DESeq2 log2 fold change (unshrunken) for each contrast; * padj < 0.05;
grey "n.d." = no estimate (gene absent from the table or log2FC empty). Contrasts that compare
libraries from different sequencing runs are marked with a dagger.

All values are read from Supplementary Table S1, so the figure always matches the table.

Usage:
  python3 fig5_expression_heatmaps.py --table Supplementary_Table_S1.xlsx --out Figure5

Edit GENES_A / GENES_B below to change the rows. Outputs <out>.png, .pdf, .tif (300 dpi).
Requires matplotlib, numpy, openpyxl.
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import openpyxl

GENES_A = [
    ("pgaptmp_002725", "PG-binding fragment (p.Gln93*)"),
    ("pgaptmp_002724", "HEAT-repeat fragment (p.Gln93*)"),
    ("pgaptmp_000955", "adenylate/guanylate cyclase (frameshift)"),
    ("pgaptmp_003298", "AAA-ATPase fragment (c.2818dupC)"),
    ("pgaptmp_003297", "WD40 fragment (c.2818dupC)"),
    ("pgaptmp_002003", "glycosyltransferase (c.123dupC)"),
    ("pgaptmp_000416", "MutS upstream fragment (pseudogene)"),
    ("pgaptmp_000418", "MutS downstream fragment"),
    ("pgaptmp_004128", "IS5 transposase"),
]
GENES_B = [
    ("JIFKDEIN_00421", "transmembrane protein, small fragment"),
    ("JIFKDEIN_00422", "transmembrane protein, large fragment"),
    ("JIFKDEIN_03191", "vitamin B12 transporter BtuB, fragment"),
    ("JIFKDEIN_03190", "vitamin B12 transporter BtuB, fragment"),
]
# contrasts whose libraries come from different sequencing runs (Supplementary Table S3.1)
CROSS_RUN = {"Dark_24h", "Time_in_light", "Nitrogen_2h_light"}
NICE = {"Dark_2h": "Dark, 2 h", "Dark_24h": "Dark, 24 h", "Time_in_light": "Time in light",
        "Time_in_dark": "Time in dark", "Nitrogen_2h_light": "N removal, 2 h light",
        "Nitrogen_2h_dark": "N removal, 2 h dark", "Nitrogen_24h_dark": "N removal, 24 h dark"}
# vertical separators after these contrasts (darkness | time | nitrogen)
BREAK_AFTER = {"Dark_24h", "Time_in_dark"}


def read_deseq(wb, sheet):
    rows = list(wb[sheet].iter_rows(values_only=True))
    h = next(i for i, r in enumerate(rows) if r and r[0] and str(r[0]).startswith(("Gene", "Locus")))
    top, sub = rows[h - 1], rows[h]
    contrasts, cur = [], None
    for j, (t, s) in enumerate(zip(top, sub)):
        if t:
            name, rest = str(t).split(" ", 1) if " " in str(t) else (str(t), "")
            cur = dict(name=name, desc=rest.strip("()"), cols={})
            contrasts.append(cur)
        if cur is not None and s in ("log2FC", "padj") and s not in cur["cols"]:
            cur["cols"][s] = j
    data = {}
    for r in rows[h + 1:]:
        if r and r[0]:
            data[str(r[0])] = {c["name"]: (r[c["cols"]["log2FC"]], r[c["cols"]["padj"]]) for c in contrasts}
    return contrasts, data


def heat(ax, cax, genes, contrasts, data, title, vmax, note):
    n, m = len(genes), len(contrasts)
    M = np.full((n, m), np.nan)
    P = np.full((n, m), np.nan)
    for i, (g, _) in enumerate(genes):
        for j, c in enumerate(contrasts):
            v = data.get(g, {}).get(c["name"])
            if v and v[0] is not None:
                M[i, j] = float(v[0])
                P[i, j] = float(v[1]) if v[1] is not None else np.nan
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#e5e5e3")
    im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")
    for i in range(n):
        for j in range(m):
            if np.isnan(M[i, j]):
                ax.text(j, i, "n.d.", ha="center", va="center", fontsize=7.5, color="#8a8984")
                continue
            star = "*" if (not np.isnan(P[i, j]) and P[i, j] < 0.05) else ""
            col = "white" if abs(M[i, j]) > 0.6 * vmax else "#1a1a19"
            ax.text(j, i, f"{M[i, j]:.2f}{star}", ha="center", va="center", fontsize=7.5, color=col)
    ax.set_xticks(np.arange(m) + 0.5, minor=True)
    ax.set_yticks(np.arange(n) + 0.5, minor=True)
    ax.grid(which="minor", color="white", lw=1.5)
    ax.tick_params(which="minor", length=0)
    for j, c in enumerate(contrasts):
        if c["name"] in BREAK_AFTER and j < m - 1:
            ax.axvline(j + 0.5, color="#3d3d3a", lw=0.9)
    labels = []
    for c in contrasts:
        d = c["desc"].replace(" vs ", " vs\n")
        labels.append(f"{NICE.get(c['name'], c['name'])}\n{d}{' †' if c['name'] in CROSS_RUN else ''}")
    ax.xaxis.tick_top()
    ax.set_xticks(range(m))
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_yticks(range(n))
    ax.set_yticklabels([f"{g}  {l}" for g, l in genes], fontsize=7.5)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = plt.colorbar(im, cax=cax)
    cb.set_label("log2 fold change (DESeq2)", fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    ax.text(0, 1, title, transform=ax.figure.transFigure, fontsize=9)  # placeholder, repositioned below
    return note


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True, help="Supplementary_Table_S1.xlsx")
    ap.add_argument("--out", default="Figure5_expression_heatmaps")
    a = ap.parse_args()
    wb = openpyxl.load_workbook(a.table, read_only=True, data_only=True)
    cA, dA = read_deseq(wb, "S1.25_LK_DESeq2")
    cB, dB = read_deseq(wb, "S1.27_Oceanicaulis_DESeq2")

    def vmax_of(genes, cons, data):
        vals = [abs(float(data[g][c["name"]][0])) for g, _ in genes for c in cons
                if g in data and data[g][c["name"]][0] is not None]
        return max(4.0, np.ceil(max(vals) * 10) / 10) if vals else 4.0

    fig = plt.figure(figsize=(12, 12), dpi=300)
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 0.025], height_ratios=[len(GENES_A) + 1.2, len(GENES_B) + 1.6],
                          hspace=0.62, wspace=0.03, left=0.32, right=0.92, top=0.86, bottom=0.06)
    axA, caxA = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    axB, caxB = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    heat(axA, caxA, GENES_A, cA, dA, "", vmax_of(GENES_A, cA, dA), "")
    heat(axB, caxB, GENES_B, cB, dB, "", 4.0, "")
    for t in list(fig.texts):
        t.remove()
    for ax, letter, title, note in [
        (axA, "A", "LK genes affected by the new mutations and the IS4 insertion in mutS",
         "* padj < 0.05.  † Contrast compares libraries from different sequencing runs.  Grey, no estimate."),
        (axB, "B", "Oceanicaulis sp.: split genes",
         "* padj < 0.05.  n.d., no estimate (too few reads).  † Different sequencing runs.")]:
        bb = ax.get_position()
        fig.text(0.015, bb.y1 + 0.075, letter, fontsize=20, fontweight="bold", va="bottom",
                 bbox=dict(boxstyle="square,pad=0.12", fc="white", ec="black", lw=1.5))
        fig.text(0.08, bb.y1 + 0.085, title, fontsize=9, va="bottom",
                 style="italic" if letter == "B" else "normal")
        fig.text(0.08, bb.y0 - 0.012, note, fontsize=6.5, va="top")
    for ext in ("png", "pdf"):
        fig.savefig(f"{a.out}.{ext}", dpi=300, bbox_inches="tight")
    fig.savefig(f"{a.out}.tif", dpi=300, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    missing = [g for g, _ in GENES_A if g not in dA] + [g for g, _ in GENES_B if g not in dB]
    if missing:
        print("No DESeq2 estimate (shown as n.d.):", ", ".join(missing))


if __name__ == "__main__":
    main()
