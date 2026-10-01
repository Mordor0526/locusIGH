#!/usr/bin/env python3
"""Plot IGH phylogeny — rectangular, compact, Word-friendly."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import Phylo
import io, re

TREE_FILE = '/home/claude/igh_tree_103.treefile'
OUT = '/mnt/user-data/outputs'

# Species colors
SP_COL = {
    'g': '#2563eb',   # L. guanicoe
    'l': '#1d4ed8',   # L. glama
    'v': '#dc2626',   # V. vicugna
    'a': '#ea580c',   # V. pacos
}
SP_LABEL = {
    'g': 'L. guanicoe',
    'l': 'L. glama',
    'v': 'V. vicugna',
    'a': 'V. pacos',
}

def get_sp_key(name):
    if name and ('cbac' in name.lower()):
        return 'out'
    if name:
        last = name.strip()[-1].lower()
        if last in SP_COL:
            return last
    return 'out'

# Clean the Newick: IQ-TREE uses "SH-aLRT/UFBoot" as internal labels
# Bio.Phylo can't parse "/" in names well, strip internal labels for clean reading
with open(TREE_FILE) as f:
    nwk = f.read().strip()

# Remove internal node labels (numbers like 55.1/75 before the colon)
nwk_clean = re.sub(r'\)[\d./]+:', '):', nwk)

tree = Phylo.read(io.StringIO(nwk_clean), 'newick')

# Root on cbactrianus
outgroup = None
for tip in tree.get_terminals():
    if 'cbac' in tip.name.lower():
        outgroup = tip
        break
if outgroup:
    tree.root_with_outgroup(outgroup)

# Rename tips for display
name_map = {}
for tip in tree.get_terminals():
    sp = get_sp_key(tip.name)
    if sp == 'out':
        name_map[tip.name] = 'C. bactrianus'
    else:
        name_map[tip.name] = f'{SP_LABEL[sp]} {tip.name}'

# Color tips
def get_color(name):
    sp = get_sp_key(name)
    return SP_COL.get(sp, '#6b7280')

# Plot
fig, ax = plt.subplots(figsize=(12, 28))

# Draw tree
Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False,
           label_func=lambda c: name_map.get(c.name, '') if c.is_terminal() else '',
           label_colors=lambda name: get_color(name) if name else '#333')

# Increase tip label font size and set italic
for text in ax.texts:
    text.set_fontsize(9)
    text.set_fontstyle('italic')

# Style
ax.set_ylabel('')
ax.set_xlabel('Substitutions per site', fontsize=11)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_yticks([])

# Legend
import matplotlib.patches as mpatches
legend_elements = [
    mpatches.Patch(facecolor='#2563eb', label='L. guanicoe (n = 38)'),
    mpatches.Patch(facecolor='#1d4ed8', label='L. glama (n = 14)'),
    mpatches.Patch(facecolor='#dc2626', label='V. vicugna (n = 33)'),
    mpatches.Patch(facecolor='#ea580c', label='V. pacos (n = 17)'),
    mpatches.Patch(facecolor='#6b7280', label='C. bactrianus (n = 1)'),
]
leg = ax.legend(handles=legend_elements, loc='upper right', fontsize=10, framealpha=0.9)
for t in leg.get_texts():
    # Extract species name part and set italic
    t.set_fontstyle('italic')

plt.tight_layout()
plt.savefig(f'{OUT}/Fig2_IGH_phylogeny_rect.png', dpi=300, bbox_inches='tight')
plt.close()
print("Fig 2 rectangular done")
