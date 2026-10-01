#!/usr/bin/env python3
"""Generate publication-quality figures for the IGH locus manuscript."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import numpy as np

OUT = '/mnt/user-data/outputs'

# ── Colour scheme ──
COL = {
    'VH': '#4C72B0',
    'VHH': '#DD5555',
    'D': '#55A868',
    'J': '#C4A000',
    'C': '#8172B2',
    'C_HCAb': '#C44E52',
    'C_conv': '#4C72B0',
    'C_other': '#8172B2',
}
SP_COL = {
    'L_guanicoe': '#2563eb',
    'L_glama': '#1d4ed8',
    'V_vicugna': '#dc2626',
    'V_pacos': '#ea580c',
}

# ── Gene coordinates (from gene_coords_aln.tsv, mapped to ~233 kb reference) ──
# Using approximate kb positions on AM773729.1 for illustration
genes_ref = [
    ('vhh3-1', 0.5, 1.1, 'VHH'),
    ('vh3-1', 3.0, 3.6, 'VH'),
    ('vh3-2', 6.0, 6.6, 'VH'),
    ('vh1-1', 9.5, 10.1, 'VH'),
    ('ighd-1', 18, 18.1, 'D'),
    ('ighd-2', 19, 19.1, 'D'),
    ('ighd-3', 20, 20.1, 'D'),
    ('ighd-4', 21, 21.1, 'D'),
    ('ighd-5', 22, 22.1, 'D'),
    ('ighd-6', 23, 23.1, 'D'),
    ('ighd-7', 24, 24.1, 'D'),
    ('ighJ-1', 28, 28.1, 'J'),
    ('ighJ-2', 29, 29.1, 'J'),
    ('ighJ-3', 30, 30.1, 'J'),
    ('ighJ-4', 31, 31.1, 'J'),
    ('ighJ-5', 32, 32.1, 'J'),
    ('ighJ-6', 33, 33.1, 'J'),
    ('ighJ-7', 34, 34.1, 'J'),
    ('ighmu', 40, 50.6, 'C'),
    ('ighdelta', 55, 65.1, 'C'),
    ('ighgamma2b', 75, 85.1, 'C_HCAb'),
    ('ighgamma1a', 95, 106.5, 'C_conv'),
    ('ighgamma1b', 115, 131.1, 'C_conv'),
    ('ighgamma2c', 140, 150.2, 'C_HCAb'),
    ('ighepsilon', 160, 170.3, 'C'),
    ('ighalpha', 180, 191.0, 'C'),
]

# ═══════════════════════════════════════════════════════════════════════
# FIGURE 1 — Schematic map of the IGH locus
# ═══════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(14, 3.5))

# Draw backbone
ax.plot([0, 222.8], [0, 0], color='#333', lw=2, zorder=1)

type_colors = {
    'VHH': '#C44E52', 'VH': '#4C72B0', 'D': '#55A868',
    'J': '#C4A000', 'C': '#8172B2', 'C_HCAb': '#C44E52', 'C_conv': '#4C72B0',
}

for gene, start, end, gtype in genes_ref:
    w = max(end - start, 1.5)
    col = type_colors.get(gtype, '#999')
    h = 3 if gtype in ('C', 'C_HCAb', 'C_conv') else 2
    rect = mpatches.FancyBboxPatch(
        (start, -h/2), w, h,
        boxstyle="round,pad=0.15", facecolor=col, edgecolor='white', lw=0.8, zorder=2
    )
    ax.add_patch(rect)
    # Label — short name
    short = gene.replace('igh', '').replace('gamma', 'γ').replace('delta', 'δ').replace('mu', 'μ').replace('epsilon', 'ε').replace('alpha', 'α')
    fs = 6.5 if gtype in ('C', 'C_HCAb', 'C_conv') else 5.5
    rot = 0 if gtype in ('C', 'C_HCAb', 'C_conv') else 45
    yoff = -4.5 if gtype in ('D', 'J') else 4 if gtype in ('VHH', 'VH') else 0
    if gtype in ('C', 'C_HCAb', 'C_conv'):
        ax.text(start + w/2, 0, short, ha='center', va='center', fontsize=fs, color='white', fontweight='bold', zorder=3)
    else:
        ax.text(start + w/2, yoff, short, ha='center', va='center', fontsize=fs, rotation=rot, color='#333', zorder=3)

# Region labels
ax.annotate('V region', xy=(5, 5.5), fontsize=9, ha='center', color='#555', style='italic')
ax.annotate('D region', xy=(21, -7), fontsize=9, ha='center', color='#555', style='italic')
ax.annotate('J region', xy=(31, 5.5), fontsize=9, ha='center', color='#555', style='italic')
ax.annotate('C region', xy=(130, 5.5), fontsize=9, ha='center', color='#555', style='italic')

# HCAb markers
for gene, start, end, gtype in genes_ref:
    if gtype == 'C_HCAb':
        w = max(end - start, 1.5)
        ax.annotate('HCAb', xy=(start + w/2, -3.5), fontsize=6, ha='center', color='#C44E52', fontweight='bold')

# Scale bar
ax.plot([190, 210], [-8, -8], color='black', lw=1.5)
ax.text(200, -9.5, '20 kb', ha='center', va='top', fontsize=8)

# Legend
legend_elements = [
    mpatches.Patch(facecolor='#C44E52', label='VHH / HCAb constant'),
    mpatches.Patch(facecolor='#4C72B0', label='VH / conventional constant'),
    mpatches.Patch(facecolor='#55A868', label='D segments'),
    mpatches.Patch(facecolor='#C4A000', label='J segments'),
    mpatches.Patch(facecolor='#8172B2', label='Other constant (μ, δ, ε, α)'),
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=7, framealpha=0.9, ncol=2)

ax.set_xlim(-5, 228)
ax.set_ylim(-12, 9)
ax.set_xlabel('Position (kb)', fontsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_yticks([])

plt.tight_layout()
plt.savefig(f'{OUT}/Fig1_IGH_locus_map.png', dpi=300, bbox_inches='tight')
plt.close()
print("Fig 1 done")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE — Tajima's D by gene region (all SAC)
# ═══════════════════════════════════════════════════════════════════════
genes_d = [
    ('vhh3-1', 'VHH', -1.6231),
    ('vh3-1', 'VH', -2.0455),
    ('vh3-2', 'VH', -1.5329),
    ('vh1-1', 'VH', -1.3982),
    ('ighd-1', 'D', 0.5543),
    ('ighd-2', 'D', -1.4263),
    ('ighd-3', 'D', 0.3271),
    ('ighd-4', 'D', 0.5230),
    ('ighd-5', 'D', -1.5320),
    ('ighd-6', 'D', 0.2919),
    ('ighd-7', 'D', -1.0257),
    ('ighJ-1', 'J', -0.5666),
    ('ighJ-2', 'J', 0.2776),
    ('ighJ-3', 'J', -1.3755),
    ('ighJ-4', 'J', -1.4284),
    ('ighJ-5', 'J', -1.5834),
    ('ighJ-6', 'J', -1.4879),
    ('ighJ-7', 'J', -1.5834),
    ('ighmu', 'C', -1.8822),
    ('ighdelta', 'C', -1.1018),
    ('ighgamma2b', 'C_HCAb', -1.1356),
    ('ighgamma1a', 'C_conv', -1.3457),
    ('ighgamma1b', 'C_conv', -0.9064),
    ('ighgamma2c', 'C_HCAb', -0.5854),
    ('ighepsilon', 'C', -0.5535),
    ('ighalpha', 'C', -0.8588),
]

names = [g[0].replace('igh', '').replace('gamma', 'γ').replace('delta', 'δ').replace('mu', 'μ').replace('epsilon', 'ε').replace('alpha', 'α') for g in genes_d]
d_vals = [g[2] for g in genes_d]
colors = [type_colors.get(g[1], '#999') for g in genes_d]

fig, ax = plt.subplots(figsize=(14, 5))
bars = ax.bar(range(len(names)), d_vals, color=colors, edgecolor='white', lw=0.5, width=0.75)

# Significance thresholds
ax.axhline(y=-1.5, color='#888', ls='--', lw=0.8, alpha=0.7)
ax.axhline(y=-2.0, color='#444', ls='--', lw=0.8, alpha=0.7)
ax.axhline(y=0, color='black', lw=0.5)
ax.text(len(names)-0.3, -1.45, 'p < 0.05', fontsize=7, ha='right', color='#888')
ax.text(len(names)-0.3, -1.95, 'p < 0.01', fontsize=7, ha='right', color='#444')

ax.set_xticks(range(len(names)))
ax.set_xticklabels(names, rotation=45, ha='right', fontsize=8)
ax.set_ylabel("Tajima's D", fontsize=11)
# Title removed — goes in caption
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

legend_elements2 = [
    mpatches.Patch(facecolor='#C44E52', label='VHH / HCAb constant'),
    mpatches.Patch(facecolor='#4C72B0', label='VH / conventional constant'),
    mpatches.Patch(facecolor='#55A868', label='D segments'),
    mpatches.Patch(facecolor='#C4A000', label='J segments'),
    mpatches.Patch(facecolor='#8172B2', label='Other constant'),
]
ax.legend(handles=legend_elements2, fontsize=7, loc='upper right', framealpha=0.9)

plt.tight_layout()
plt.savefig(f'{OUT}/Fig_Tajima_D_per_gene.png', dpi=300, bbox_inches='tight')
plt.close()
print("Fig Tajima D per gene done")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE — Tajima's D by functional category per species
# ═══════════════════════════════════════════════════════════════════════
categories = ['VHH', 'VH_conv', 'D_seg', 'J_seg', 'C_HCAb', 'C_conv', 'C_other', 'Full']
cat_labels = ['VHH', 'VH\nconv.', 'D\nseg.', 'J\nseg.', 'C\nHCAb', 'C\nconv.', 'C\nother', 'Full\nlocus']

# Data: species × category
data = {
    'L. guanicoe': [-1.0132, -1.3356, 0.0447, -0.8539, -0.1452, -0.5659, -0.5888, -0.5136],
    'L. glama':    [-1.0575, -1.0508, -0.2606, -0.9870, -0.5731, -0.6319, -0.8896, -0.6478],
    'V. vicugna':  [-1.2180, -1.5341, -0.8579, -1.6519, -0.8626, -0.7392, -0.5567, -0.6308],
    'V. pacos':    [-1.0992, -1.3032, 0.5148, -1.2974, -0.8809, -0.6252, -0.5240, -0.3090],
}
sp_colors = ['#2563eb', '#1d4ed8', '#dc2626', '#ea580c']

x = np.arange(len(categories))
w = 0.18

fig, ax = plt.subplots(figsize=(12, 5.5))
for i, (sp, vals) in enumerate(data.items()):
    ax.bar(x + i*w - 1.5*w, vals, w, label=sp, color=sp_colors[i], edgecolor='white', lw=0.5)

ax.axhline(y=-1.5, color='#888', ls='--', lw=0.8, alpha=0.7)
ax.axhline(y=0, color='black', lw=0.5)
ax.text(len(categories)-0.8, -1.45, 'p < 0.05', fontsize=7, ha='right', color='#888')

ax.set_xticks(x)
ax.set_xticklabels(cat_labels, fontsize=9)
ax.set_ylabel("Tajima's D", fontsize=11)
# Title removed — goes in caption
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
leg6 = ax.legend(fontsize=8, loc='upper right', framealpha=0.9)
for t in leg6.get_texts():
    t.set_fontstyle('italic')

plt.tight_layout()
plt.savefig(f'{OUT}/Fig_Tajima_D_by_species.png', dpi=300, bbox_inches='tight')
plt.close()
print("Fig Tajima D by species done")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE — F_ST heatmap
# ═══════════════════════════════════════════════════════════════════════
species_labels = ['L. guanicoe', 'L. glama', 'V. vicugna', 'V. pacos']
fst_matrix = np.array([
    [0.0000, 0.1094, 0.3221, 0.1785],
    [0.1094, 0.0000, 0.2328, 0.0645],
    [0.3221, 0.2328, 0.0000, 0.1935],
    [0.1785, 0.0645, 0.1935, 0.0000],
])

fig, ax = plt.subplots(figsize=(6, 5))
im = ax.imshow(fst_matrix, cmap='YlOrRd', vmin=0, vmax=0.35, aspect='equal')

for i in range(4):
    for j in range(4):
        val = fst_matrix[i, j]
        color = 'white' if val > 0.2 else 'black'
        ax.text(j, i, f'{val:.4f}', ha='center', va='center', fontsize=10, color=color, fontweight='bold')

ax.set_xticks(range(4))
ax.set_yticks(range(4))
ax.set_xticklabels(species_labels, fontsize=9, rotation=30, ha='right', style='italic')
ax.set_yticklabels(species_labels, fontsize=9, style='italic')
# Title removed — goes in caption

cbar = plt.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label('$F_{ST}$', fontsize=10)

plt.tight_layout()
plt.savefig(f'{OUT}/Fig_FST_heatmap.png', dpi=300, bbox_inches='tight')
plt.close()
print("Fig FST heatmap done")

print("\nAll figures generated.")
