#!/usr/bin/env python3
"""PCA scatter plot of IgH locus - 102 CSA + outgroup."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Parse data
data = []
with open('/dev/stdin') as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 8:
            continue
        data.append({
            'name': parts[0], 'species': parts[1], 'genus': parts[2],
            'pc1': float(parts[3]), 'pc2': float(parts[4]), 'pc3': float(parts[5])
        })

# Species config: color, marker, label
config = {
    'L_guanicoe':   {'color': '#2563eb', 'marker': 'o', 'label': 'L. guanicoe (n=38)'},
    'L_glama':      {'color': '#1d4ed8', 'marker': 's', 'label': 'L. glama (n=14)'},
    'V_vicugna':    {'color': '#dc2626', 'marker': 'o', 'label': 'V. vicugna (n=33)'},
    'V_pacos':      {'color': '#ea580c', 'marker': 's', 'label': 'V. pacos (n=17)'},
    'C_bactrianus': {'color': '#6b7280', 'marker': 'D', 'label': 'C. bactrianus (n=1)'},
}

fig, ax = plt.subplots(figsize=(10, 8))

# Plot each species
for sp, cfg in config.items():
    pts = [d for d in data if d['species'] == sp]
    if not pts:
        continue
    x = [p['pc1'] for p in pts]
    y = [p['pc2'] for p in pts]
    ax.scatter(x, y, c=cfg['color'], marker=cfg['marker'], s=60,
               edgecolors='white', linewidths=0.5, label=cfg['label'], zorder=3)

# Draw convex hulls for each genus
from matplotlib.patches import Polygon
from scipy.spatial import ConvexHull

for genus, color, alpha in [('Lama', '#2563eb', 0.08), ('Vicugna', '#dc2626', 0.08)]:
    pts = [(d['pc1'], d['pc2']) for d in data if d['genus'] == genus]
    if len(pts) < 3:
        continue
    arr = np.array(pts)
    hull = ConvexHull(arr)
    hull_pts = arr[hull.vertices]
    poly = Polygon(hull_pts, alpha=alpha, color=color, zorder=1)
    ax.add_patch(poly)

ax.set_xlabel('PC1 (19.0%)', fontsize=12)
ax.set_ylabel('PC2 (8.5%)', fontsize=12)
ax.set_title('PCA of IgH Locus — South American Camelids', fontsize=14, fontweight='bold')

# Legend outside plot, to the right
legend = ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=10,
                   framealpha=0.9, borderaxespad=0)
for text in legend.get_texts():
    text.set_fontstyle('italic')

ax.axhline(0, color='gray', linewidth=0.5, linestyle='--', zorder=0)
ax.axvline(0, color='gray', linewidth=0.5, linestyle='--', zorder=0)
ax.grid(True, alpha=0.2)
ax.tick_params(labelsize=10)

plt.savefig('/mnt/user-data/outputs/IgH_PCA_103.png', dpi=300, bbox_inches='tight')
print("Saved PCA figure")
