#!/usr/bin/env python3
"""Rebuild Supplementary Figure 1 with non-overlapping PRISMA-style boxes."""
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/研究升级/Communications_Biology_v6'
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    'font.family': 'Arial', 'font.size': 8,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

fig, ax = plt.subplots(figsize=(7.5, 6.1))
ax.set_xlim(0, 12)
ax.set_ylim(0, 12)
ax.axis('off')
edge = '#244A66'

def box(cx, cy, w, h, text, *, dashed=False, fontsize=7.5):
    patch = FancyBboxPatch(
        (cx-w/2, cy-h/2), w, h,
        boxstyle='round,pad=0.06,rounding_size=0.10',
        facecolor='white', edgecolor=edge, linewidth=1.2,
        linestyle='--' if dashed else '-')
    ax.add_patch(patch)
    ax.text(cx, cy, text, ha='center', va='center', fontsize=fontsize,
            linespacing=1.18)
    return {'cx':cx, 'cy':cy, 'w':w, 'h':h}

def arrow(x1, y1, x2, y2, *, dashed=False):
    ax.annotate('', xy=(x2,y2), xytext=(x1,y1),
                arrowprops=dict(arrowstyle='->', color='#555555', lw=1.0,
                                linestyle='--' if dashed else '-',
                                shrinkA=0, shrinkB=0))

ax.text(.25, 11.68,
        'Supplementary Fig. 1 | Public-data identification and eligibility flow',
        fontsize=10.5, fontweight='bold', va='top')

# Identification
db = box(3.0, 10.25, 5.2, 1.05,
         'Records identified through database searches\nGEO n = 42; PubMed n = 3')
add = box(9.0, 10.25, 4.4, 1.05,
          'Additional record from citation and\naccession checking; n = 1')

# Screening and eligibility. Main flow is left; exclusions are right.
screen = box(3.8, 8.05, 4.2, .95, 'Records screened\nn = 46')
ex_title = box(9.0, 8.05, 4.7, .95, 'Excluded at title/summary\nn = 34')
full = box(3.8, 5.72, 4.2, .95, 'Full records assessed for eligibility\nn = 12')
ex_full = box(9.0, 5.72, 5.1, 1.45,
              'Excluded after full assessment; n = 5\n'
              'One donor (1); split series (2)\n'
              'Glass only (1); strain/vibration (1)', fontsize=7.1)
included = box(3.8, 3.35, 4.8, 1.0,
               'Studies included in quantitative synthesis\nn = 7')

# This dataset was not part of the seven-study synthesis.
validation = box(9.0, 3.35, 5.1, 1.35,
                 'Separate independent validation source\n'
                 'GSE147390; one donor\n'
                 'n = 5,329 osteoblasts', dashed=True, fontsize=7.1)

# Keep a visible gap at both ends of every connector. Arrow shafts and heads do
# not touch the box outlines, which prevents apparent overlap after reduction.
gap = .14
arrow(db['cx'], db['cy']-db['h']/2-gap,
      screen['cx']-.55, screen['cy']+screen['h']/2+gap)
arrow(add['cx'], add['cy']-add['h']/2-gap,
      screen['cx']+.55, screen['cy']+screen['h']/2+gap)
arrow(screen['cx']+screen['w']/2+gap, screen['cy'],
      ex_title['cx']-ex_title['w']/2-gap, ex_title['cy'])
arrow(screen['cx'], screen['cy']-screen['h']/2-gap,
      full['cx'], full['cy']+full['h']/2+gap)
arrow(full['cx']+full['w']/2+gap, full['cy'],
      ex_full['cx']-ex_full['w']/2-gap, ex_full['cy'])
arrow(full['cx'], full['cy']-full['h']/2-gap,
      included['cx'], included['cy']+included['h']/2+gap)

fig.subplots_adjust(left=.025, right=.985, top=.98, bottom=.03)
fig.savefig(OUT/'Supplementary_Figure_1_PRISMA.png', dpi=600,
            bbox_inches='tight', pad_inches=.06, facecolor='white')
fig.savefig(OUT/'Supplementary_Figure_1_PRISMA.tif', dpi=600,
            bbox_inches='tight', pad_inches=.06, facecolor='white',
            pil_kwargs={'compression':'tiff_lzw'})
fig.savefig(OUT/'Supplementary_Figure_1_PRISMA.pdf',
            bbox_inches='tight', pad_inches=.06, facecolor='white')
plt.close(fig)
print(OUT)
