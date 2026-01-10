#!/usr/bin/env python3
"""Plot VMS extended applications."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Arc
import matplotlib.gridspec as gridspec
from pathlib import Path

FIGURES_DIR = Path('data/figures')

plt.style.use('dark_background')

fig = plt.figure(figsize=(16, 12))
fig.suptitle('VMS TOPOLOGICAL ANALYSIS - EXTENDED APPLICATIONS',
             fontsize=16, fontweight='bold', color='#00ff88')

gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.25)

# 1. Application suitability chart
ax1 = fig.add_subplot(gs[0, 0])

applications = ['Luna\n(Apollo)', 'Ionosfera\n(Flares)', 'Marte\n(InSight)',
                'Hielo\n(Glaciar)', 'Tierra\n(Sismica)', 'Oceano\n(Acustica)', 'Audio\n(Aire)']
coherence = [83, 70, 65, 55, 24, 40, 10]
colors = ['#ffe66d', '#ff6b6b', '#e74c3c', '#3498db', '#4ecdc4', '#9b59b6', '#95a5a6']

bars = ax1.barh(applications, coherence, color=colors, edgecolor='white', linewidth=2)
ax1.set_xlabel('Coherencia Topologica Estimada (%)', fontsize=11)
ax1.set_title('Idoneidad para VMS por Aplicacion', fontweight='bold')
ax1.set_xlim(0, 100)

for bar, val in zip(bars, coherence):
    ax1.text(val + 2, bar.get_y() + bar.get_height()/2, f'{val}%',
             va='center', fontsize=10, fontweight='bold')

ax1.axvline(50, color='white', linestyle='--', alpha=0.5)
ax1.text(52, 6.5, 'Umbral efectivo', fontsize=9, color='white', alpha=0.7)

# 2. Ionospheric monitoring concept
ax2 = fig.add_subplot(gs[0, 1])
ax2.set_xlim(0, 10)
ax2.set_ylim(0, 10)
ax2.axis('off')
ax2.set_title('MONITOREO IONOSFERICO - Llamaradas Solares', fontweight='bold', color='#ff6b6b')

# Sun
sun = Circle((1, 8), 0.8, color='#FFD700', ec='#FFA500', linewidth=3)
ax2.add_patch(sun)
ax2.text(1, 8, 'SOL', fontsize=12, ha='center', va='center', fontweight='bold')

# Earth
earth = Circle((8, 5), 0.6, color='#4169E1', ec='white', linewidth=2)
ax2.add_patch(earth)
ax2.text(8, 5, 'TIERRA', fontsize=8, ha='center', va='center', fontweight='bold')

# Ionosphere arc
arc = Arc((8, 5), 2.5, 2.5, angle=0, theta1=40, theta2=140, color='#00ff88', linewidth=4)
ax2.add_patch(arc)
ax2.text(8, 7, 'Ionosfera', fontsize=10, ha='center', color='#00ff88', fontweight='bold')

# Solar flare arrow
ax2.annotate('', xy=(6.5, 5.5), xytext=(2, 7.5),
            arrowprops=dict(arrowstyle='->', color='#ff6b6b', lw=4))
ax2.text(3.5, 7.2, 'Llamarada\nSolar', fontsize=11, ha='center', color='#ff6b6b', fontweight='bold')

# Effects list
ax2.text(0.5, 4, 'Efectos detectables:', fontsize=11, fontweight='bold', color='white')
effects = ['TEC (contenido electronico)', 'VLF/ELF propagacion',
           'Campo geomagnetico', 'Scintilacion GPS']
for i, eff in enumerate(effects):
    ax2.text(0.5, 3.3 - i*0.6, f'> {eff}', fontsize=10, color='#aaaaaa')

# VMS result box
vms_box = FancyBboxPatch((0.3, 0.3), 4, 1.2, boxstyle='round,pad=0.1',
                          facecolor='#1a1a1a', edgecolor='#00ff88', linewidth=2)
ax2.add_patch(vms_box)
ax2.text(2.3, 0.9, 'VMS detecta "montanas"\nde perturbacion ionosferica',
         fontsize=10, ha='center', color='#00ff88')

# 3. Why topology works
ax3 = fig.add_subplot(gs[1, 0])

t = np.linspace(0, 100, 1000)
np.random.seed(42)

# High coherence (lunar/ionospheric-like)
lunar = np.exp(-((t-30)/20)**2) * np.sin(2*np.pi*0.1*t)
for i in range(20):
    lunar += np.exp(-((t-30-i*2)/15)**2) * np.sin(2*np.pi*(0.1+i*0.01)*t) / (i+1)
lunar = lunar / np.max(np.abs(lunar))

# Low coherence (audio-like)
audio = np.exp(-((t-30)/5)**2) * np.sin(2*np.pi*0.5*t)
audio += np.random.randn(len(t)) * 0.5
audio = audio / np.max(np.abs(audio))

ax3.plot(t, lunar + 1.5, color='#ffe66d', linewidth=1.5, label='Alta coherencia (Luna/Ionosfera)')
ax3.fill_between(t, lunar + 1.5, 1.5, alpha=0.3, color='#ffe66d')

ax3.plot(t, audio - 1.5, color='#95a5a6', linewidth=0.5, label='Baja coherencia (Audio/Aire)')
ax3.fill_between(t, audio - 1.5, -1.5, alpha=0.3, color='#95a5a6')

ax3.set_xlabel('Tiempo')
ax3.set_ylabel('Amplitud')
ax3.set_title('Por Que la Topologia Importa', fontweight='bold')
ax3.legend(loc='upper right')
ax3.axhline(0, color='gray', linestyle='--', alpha=0.3)

# Annotations
ax3.annotate('Estructura\nCONECTADA', xy=(60, 2.2), fontsize=11, ha='center', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='#ffe66d', alpha=0.5, edgecolor='white'))
ax3.annotate('Estructura\nDISPERSA', xy=(60, -0.8), fontsize=11, ha='center', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='#95a5a6', alpha=0.5, edgecolor='white'))

# 4. Priority applications table
ax4 = fig.add_subplot(gs[1, 1])
ax4.axis('off')

# Create table manually
table_data = [
    ['Prioridad', 'Aplicacion', 'Coherencia', 'Datos'],
    ['1', 'IONOSFERA / FLARES', '~70%', 'GOES/GNSS'],
    ['2', 'MARTE / InSight', '~65%', 'NASA PDS'],
    ['3', 'HIELO / Glaciares', '~55%', 'Seismic arrays'],
    ['4', 'VOLCANES / Tremor', '~40%', 'Observatorios'],
    ['5', 'OCEANO / SOFAR', '~40%', 'Hydrophones'],
]

table = ax4.table(cellText=table_data, loc='upper center', cellLoc='center',
                 colWidths=[0.15, 0.35, 0.2, 0.25])
table.auto_set_font_size(False)
table.set_fontsize(11)
table.scale(1.2, 2.5)

# Style header
for i in range(4):
    table[(0, i)].set_facecolor('#333333')
    table[(0, i)].set_text_props(fontweight='bold', color='#00ff88')

# Style priority column
for i in range(1, 6):
    table[(i, 0)].set_facecolor('#222222')
    table[(i, 0)].set_text_props(fontweight='bold', color='#ffe66d')

ax4.set_title('Aplicaciones Prioritarias', fontweight='bold', y=0.95)

# Add key insight text
insight = '"Cualquier senal que viaja por un medio de baja dispersion\ncon estructura preservada es candidata para VMS topologico"'
ax4.text(0.5, 0.15, insight, transform=ax4.transAxes, fontsize=10,
        ha='center', va='center', style='italic',
        bbox=dict(boxstyle='round', facecolor='#1a1a1a', edgecolor='#00ff88', linewidth=2))

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'vms_extended_applications.png', dpi=150, bbox_inches='tight',
            facecolor='#0a0a0a', edgecolor='none')
print(f'Saved: {FIGURES_DIR}/vms_extended_applications.png')
