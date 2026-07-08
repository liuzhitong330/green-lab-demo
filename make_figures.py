"""
Two SVG figures for the Green lab demo:

1. muscarinic_heatmap.svg
   All 13 genes × 8 regions, grouped by function.
   Muscarinic M1 (Chrm1) vs OPC lineage markers highlighted.

2. m1_opc_scatter.svg
   Per-region scatter: X = Chrm1 (M1 brake on OPC differentiation),
   Y = Olig2 (OL lineage cell density proxy).
   Dot size = Plp1 (mature myelin content).
   Interpretation: upper-right quadrant = high OL density + high M1 brake
   = regions where M1 antagonism (clemastine) has the most therapeutic leverage.
"""
import csv, math
from pathlib import Path

OUT = Path(__file__).parent

STRUCTURES = ["Isocortex","Hippocampus","Cerebellum","Striatum",
              "Hypothalamus","Midbrain","Thalamus","Olfact_bulb"]
STRUCT_LABELS = ["Isocortex","Hippocampus","Cerebellum","Striatum",
                 "Hypothalamus","Midbrain","Thalamus","Olf. bulb"]

GROUPS = [
    ("muscarinic ACh receptors",  ["Chrm1","Chrm2","Chrm3","Chrm4"],      "#c0392b"),
    ("OPC / OL lineage",          ["Olig2","Pdgfra","Cspg4","Myrf","Sox10"], "#2980b9"),
    ("mature myelin",             ["Mbp","Plp1","Cnp"],                    "#27ae60"),
    ("neuroaxonal health",        ["Nefl"],                                "#8e44ad"),
]
ALL_GENES  = [g for _, gs, _ in GROUPS for g in gs]
GENE_COLOR = {g: c for _, gs, c in GROUPS for g in gs}

rows_by_gene = {}
with open(OUT / "expression.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        rows_by_gene[r["gene"]] = [float(r[s]) for s in STRUCTURES]


def row_norm(vals):
    mn, mx = min(vals), max(vals)
    return [(v - mn) / (mx - mn) if mx > mn else 0.0 for v in vals]


def red_blue(v):
    """0 → pale blue, 1 → deep red (hot = high expression)"""
    if v <= 0.5:
        t = v / 0.5
        r = int(220 + (255 - 220) * t)
        g = int(235 + (255 - 235) * t)
        b = 255
    else:
        t = (v - 0.5) / 0.5
        r = 255
        g = int(255 + (50 - 255) * t)
        b = int(255 + (50 - 255) * t)
    return f"rgb({r},{g},{b})"


# ── Figure 1: heatmap ─────────────────────────────────────────────────────────
CELL_W = 70
CELL_H = 26
LEFT   = 100
TOP    = 115
n_genes   = len(ALL_GENES)
n_structs = len(STRUCTURES)
W = LEFT + n_structs * CELL_W + 70
H = TOP + n_genes * CELL_H + 90

# Highlight Isocortex and Hippocampus (highest M1 + highest OL density)
ISO_I  = STRUCTURES.index("Isocortex")
HIPP_I = STRUCTURES.index("Hippocampus")
highlights = [
    (ISO_I,  "#fff0f0", "#a93226", "M1 peak"),
    (HIPP_I, "#f0f8ff", "#1a5c8a", "OPC-rich"),
]

col_bg = ""
for si, bg, _, _ in highlights:
    x = LEFT + si * CELL_W
    col_bg += (f'<rect x="{x}" y="{TOP-60}" width="{CELL_W}" '
               f'height="{n_genes*CELL_H+64}" fill="{bg}" opacity="0.65"/>')

cells = ""
dividers = ""
gi_abs = 0
for _, gs, _ in GROUPS:
    for gene in gs:
        if gene not in rows_by_gene:
            gi_abs += 1; continue
        nv = row_norm(rows_by_gene[gene])
        for si, v in enumerate(nv):
            x = LEFT + si * CELL_W
            y = TOP + gi_abs * CELL_H
            cells += (f'<rect x="{x}" y="{y}" width="{CELL_W}" height="{CELL_H}" '
                      f'fill="{red_blue(v)}" stroke="white" stroke-width="1.5"/>')
            raw = rows_by_gene[gene][si]
            if v > 0.78 and raw > 0.3:
                fc = "white" if v > 0.91 else "#333"
                cells += (f'<text x="{x+CELL_W/2:.0f}" y="{y+CELL_H/2+4:.0f}" '
                          f'text-anchor="middle" font-size="8" fill="{fc}">{raw:.1f}</text>')
        gi_abs += 1
    dividers += (f'<line x1="{LEFT}" y1="{TOP+gi_abs*CELL_H}" '
                 f'x2="{LEFT+n_structs*CELL_W}" y2="{TOP+gi_abs*CELL_H}" '
                 f'stroke="#bbb" stroke-width="1.5"/>')

# Column headers
col_hdrs = ""
for si, lbl in enumerate(STRUCT_LABELS):
    x = LEFT + si * CELL_W + CELL_W // 2
    hi = next(((fc, sub) for idx, _, fc, sub in highlights if idx == si), None)
    fc = hi[0] if hi else "#444"
    fw = "bold" if hi else "normal"
    col_hdrs += (f'<text x="{x}" y="{TOP-16}" text-anchor="middle" font-size="10" '
                 f'fill="{fc}" font-weight="{fw}">{lbl}</text>')
    if hi:
        col_hdrs += (f'<text x="{x}" y="{TOP-4}" text-anchor="middle" '
                     f'font-size="8.5" fill="{hi[0]}" font-style="italic">↑ {hi[1]}</text>')

# Row labels
row_lbls = ""
gi_abs = 0
for _, gs, _ in GROUPS:
    for gene in gs:
        y = TOP + gi_abs * CELL_H + CELL_H // 2 + 4
        row_lbls += (f'<text x="{LEFT-5}" y="{y}" text-anchor="end" font-size="10.5" '
                     f'fill="{GENE_COLOR[gene]}" font-weight="600" '
                     f'font-family="ui-monospace,Menlo,monospace">{gene}</text>')
        gi_abs += 1

# Group labels on the left
gi_abs = 0
grp_lbls = ""
for grp_name, gs, grp_col in GROUPS:
    mid_y = TOP + (gi_abs + len(gs)/2) * CELL_H
    grp_lbls += (f'<text transform="rotate(-90,18,{mid_y:.0f})" '
                 f'x="18" y="{mid_y:.0f}" text-anchor="middle" font-size="8.5" '
                 f'fill="{grp_col}" font-weight="600">{grp_name}</text>')
    gi_abs += len(gs)

# Colorbar
cb_x = LEFT + n_structs * CELL_W + 10
cb_h = n_genes * CELL_H
colorbar = (f'<defs><linearGradient id="cb1" x1="0" y1="1" x2="0" y2="0">'
            f'<stop offset="0%" stop-color="{red_blue(0)}"/>'
            f'<stop offset="50%" stop-color="{red_blue(0.5)}"/>'
            f'<stop offset="100%" stop-color="{red_blue(1)}"/>'
            f'</linearGradient></defs>'
            f'<rect x="{cb_x}" y="{TOP}" width="13" height="{cb_h}" '
            f'fill="url(#cb1)" stroke="#ccc" stroke-width="0.5"/>'
            f'<text x="{cb_x+6}" y="{TOP-4}" text-anchor="middle" font-size="8" fill="#555">high</text>'
            f'<text x="{cb_x+6}" y="{TOP+cb_h+10}" text-anchor="middle" font-size="8" fill="#555">low</text>')

svg1 = f"""<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif;background:white;">
  <text x="{W//2}" y="22" text-anchor="middle" font-size="13" font-weight="600" fill="#222">
    Muscarinic Receptors, OPC Lineage, and Myelin Genes Across Brain Regions
  </text>
  <text x="{W//2}" y="38" text-anchor="middle" font-size="10" fill="#666">
    Allen Mouse Brain Atlas · ISH expression energy (row-normalized per gene)
  </text>
  <text x="{W//2}" y="53" text-anchor="middle" font-size="10" fill="#a93226">
    Chrm1 (M1 receptor) and Olig2 (OL lineage) co-peak in isocortex and hippocampus
  </text>
  {col_bg}{colorbar}{col_hdrs}{cells}{dividers}{row_lbls}{grp_lbls}
</svg>"""

with open(OUT / "muscarinic_heatmap.svg", "w") as f:
    f.write(svg1)
print("Wrote muscarinic_heatmap.svg")


# ── Figure 2: M1 brake vs OL density scatter ─────────────────────────────────
# X = Chrm1 (M1 brake on OPC differentiation)
# Y = Olig2 (OL lineage density)
# Dot size = Plp1 (mature myelin — how differentiated OLs already are)
# Color = Plp1 value (warm = already myelinated, cool = less)

FW2, FH2 = 620, 480
PAD_L = 75
PAD_B = 70
PAD_T = 90
PAD_R = 180
AW = FW2 - PAD_L - PAD_R
AH = FH2 - PAD_T - PAD_B

chrm1 = rows_by_gene["Chrm1"]
olig2 = rows_by_gene["Olig2"]
plp1  = rows_by_gene["Plp1"]

max_x = max(chrm1) * 1.12
max_y = max(olig2) * 1.12
min_x = 0
min_y = 0

REGION_COLORS = {
    "Isocortex":   "#a93226",
    "Hippocampus": "#1a5c8a",
    "Cerebellum":  "#27ae60",
    "Striatum":    "#e67e22",
    "Hypothalamus":"#8e44ad",
    "Midbrain":    "#2c3e50",
    "Thalamus":    "#16a085",
    "Olfact_bulb": "#7f8c8d",
}

LABEL_OFFSETS = {
    # (lx_off, ly_off, text-anchor)  — tuned to avoid cluster overlap
    "Isocortex":    ( 28, -14, "start"),
    "Hippocampus":  (-28,   5, "end"),
    "Cerebellum":   (  0,  38, "middle"),
    "Striatum":     (-31, -14, "end"),
    "Hypothalamus": ( 31,  20, "start"),
    "Midbrain":     (-34,   5, "end"),
    "Thalamus":     ( 31, -15, "start"),
    "Olfact_bulb":  (  0,  38, "middle"),
}

dots = ""
labels = ""
for i, region in enumerate(STRUCTURES):
    cx = PAD_L + (chrm1[i] - min_x) / (max_x - min_x) * AW
    cy = PAD_T + AH - (olig2[i] - min_y) / (max_y - min_y) * AH
    r_dot = 10 + plp1[i] / max(plp1) * 18
    col = REGION_COLORS.get(region, "#888")

    dots += (f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r_dot:.1f}" '
             f'fill="{col}" opacity="0.82" stroke="white" stroke-width="1.5"/>')
    dots += (f'<text x="{cx:.1f}" y="{cy+4:.1f}" text-anchor="middle" '
             f'font-size="8" fill="white" font-weight="600">{plp1[i]:.0f}</text>')

    lx, ly, ta = LABEL_OFFSETS.get(region, (28, 5, "start"))
    lbl = region.replace("Olfact_bulb", "Olf. bulb")
    labels += (f'<text x="{cx+lx:.1f}" y="{cy+ly:.1f}" '
               f'text-anchor="{ta}" font-size="10" fill="{col}" font-weight="600">{lbl}</text>')

# Axes
ax_lines = (
    f'<line x1="{PAD_L}" y1="{PAD_T}" x2="{PAD_L}" y2="{PAD_T+AH}" stroke="#ccc" stroke-width="1"/>'
    f'<line x1="{PAD_L}" y1="{PAD_T+AH}" x2="{PAD_L+AW}" y2="{PAD_T+AH}" stroke="#ccc" stroke-width="1"/>'
)

# Axis ticks and labels
x_ticks = ""
for val in [0, 4, 8, 12, 16]:
    if val > max_x: continue
    tx = PAD_L + val / max_x * AW
    x_ticks += (f'<line x1="{tx:.1f}" y1="{PAD_T+AH}" x2="{tx:.1f}" y2="{PAD_T+AH+5}" stroke="#999" stroke-width="1"/>'
                f'<text x="{tx:.1f}" y="{PAD_T+AH+16}" text-anchor="middle" font-size="9" fill="#666">{val}</text>')

y_ticks = ""
for val in [0, 4, 8, 12, 16]:
    if val > max_y: continue
    ty = PAD_T + AH - val / max_y * AH
    y_ticks += (f'<line x1="{PAD_L-5}" y1="{ty:.1f}" x2="{PAD_L}" y2="{ty:.1f}" stroke="#999" stroke-width="1"/>'
                f'<text x="{PAD_L-8}" y="{ty+4:.1f}" text-anchor="end" font-size="9" fill="#666">{val}</text>')

# Legend + annotation in the right-side column (keeps chart area clean)
leg_x = FW2 - PAD_R + 15
leg_y = PAD_T + 10
legend = (
    f'<text x="{leg_x}" y="{leg_y}" font-size="9.5" fill="#444" font-weight="600">Dot size = Plp1</text>'
    f'<text x="{leg_x}" y="{leg_y+13}" font-size="9" fill="#666">(mature myelin level)</text>'
    f'<text x="{leg_x}" y="{leg_y+26}" font-size="9" fill="#666">value shown inside dot</text>'
)

# Annotation box — placed in legend column below dot-size legend
q_x = leg_x
q_y = leg_y + 60
quadrant_note = (
    f'<rect x="{q_x-4}" y="{q_y-14}" width="155" height="54" rx="4" '
    f'fill="#fff5f5" stroke="#c0392b" stroke-width="1"/>'
    f'<text x="{q_x+2}" y="{q_y}" font-size="9" fill="#c0392b" font-weight="600">'
    f'Upper-right quadrant:</text>'
    f'<text x="{q_x+2}" y="{q_y+13}" font-size="9" fill="#c0392b">'
    f'High OL density +</text>'
    f'<text x="{q_x+2}" y="{q_y+26}" font-size="9" fill="#c0392b">'
    f'High M1 brake</text>'
    f'<text x="{q_x+2}" y="{q_y+39}" font-size="9" fill="#c0392b">'
    f'→ most M1 antagonist leverage</text>'
)

svg2 = f"""<svg viewBox="0 0 {FW2} {FH2}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif;background:white;">
  <text x="{FW2//2}" y="22" text-anchor="middle" font-size="13" font-weight="600" fill="#222">
    M1 Receptor Brake vs. OL Lineage Density — Therapeutic Opportunity Map
  </text>
  <text x="{FW2//2}" y="38" text-anchor="middle" font-size="10" fill="#666">
    X = Chrm1 (M1 receptor expression — higher = more brake on OPC differentiation)
  </text>
  <text x="{FW2//2}" y="52" text-anchor="middle" font-size="10" fill="#666">
    Y = Olig2 (OL lineage density) · Dot size = Plp1 (mature myelin level, value in dot)
  </text>
  <text x="{PAD_L+AW/2:.0f}" y="{PAD_T+AH+38}" text-anchor="middle" font-size="10" fill="#444">
    Chrm1 expression (ISH energy) — M1 receptor brake intensity
  </text>
  <text transform="rotate(-90,22,{PAD_T+AH/2:.0f})" x="22" y="{PAD_T+AH/2:.0f}"
    text-anchor="middle" font-size="10" fill="#444">
    Olig2 expression — OL lineage density
  </text>
  {ax_lines}{x_ticks}{y_ticks}{dots}{labels}{quadrant_note}{legend}
</svg>"""

with open(OUT / "m1_opc_scatter.svg", "w") as f:
    f.write(svg2)
print("Wrote m1_opc_scatter.svg")
