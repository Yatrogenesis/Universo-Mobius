#!/usr/bin/env python3
"""
VMS Topological Analysis - Visualizations
==========================================

Create comprehensive plots of all analysis results.

Author: Francisco Molina Burgos
Date: January 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle, Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.gridspec as gridspec
from pathlib import Path
import json
import sys

# Paths
DATA_DIR = Path(__file__).parent.parent / 'data'
RESULTS_DIR = DATA_DIR / 'results'
FIGURES_DIR = DATA_DIR / 'figures'
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Style
plt.style.use('dark_background')
plt.rcParams['figure.facecolor'] = '#0a0a0a'
plt.rcParams['axes.facecolor'] = '#0a0a0a'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.labelcolor'] = '#ffffff'
plt.rcParams['text.color'] = '#ffffff'
plt.rcParams['xtick.color'] = '#ffffff'
plt.rcParams['ytick.color'] = '#ffffff'
plt.rcParams['grid.color'] = '#333333'
plt.rcParams['font.size'] = 10


def plot_earth_vs_moon_comparison():
    """Compare Earth and Moon seismic topology."""
    print("Creating Earth vs Moon comparison...")

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('VMS TOPOLOGICAL COMPARISON: EARTH vs MOON', fontsize=16, fontweight='bold', color='#00ff88')

    # Data
    earth_data = {
        'coherent': 214,
        'scattered': 1642,
        'coda_decay_s': 346,
        'coherence_ratio': 0.245,
        'dominant_freq': 0.5
    }

    moon_data = {
        'coherent': 3,
        'scattered': 7,
        'coda_decay_s': 14400,  # 4 hours
        'coherence_ratio': 0.833,
        'dominant_freq': 0.5
    }

    # 1. Component counts
    ax = axes[0, 0]
    x = np.arange(2)
    width = 0.35
    ax.bar(x - width/2, [earth_data['coherent'], moon_data['coherent']], width,
           label='Coherent', color='#00ff88', alpha=0.8)
    ax.bar(x + width/2, [earth_data['scattered'], moon_data['scattered']], width,
           label='Scattered', color='#ff6b6b', alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(['Earth', 'Moon'])
    ax.set_ylabel('Number of Components')
    ax.set_title('Topological Components', fontweight='bold')
    ax.legend()
    ax.set_yscale('log')

    # 2. Coherence ratio pie charts
    ax = axes[0, 1]
    ax.pie([earth_data['coherence_ratio'], 1 - earth_data['coherence_ratio']],
           labels=['Coherent', 'Scattered'], colors=['#00ff88', '#ff6b6b'],
           autopct='%1.1f%%', startangle=90)
    ax.set_title('Earth Coherence: 24.5%', fontweight='bold')

    ax = axes[0, 2]
    ax.pie([moon_data['coherence_ratio'], 1 - moon_data['coherence_ratio']],
           labels=['Coherent', 'Scattered'], colors=['#00ff88', '#ff6b6b'],
           autopct='%1.1f%%', startangle=90)
    ax.set_title('Moon Coherence: 83.3%', fontweight='bold')

    # 3. Coda decay comparison
    ax = axes[1, 0]
    categories = ['Earth\n(Japan M7.5)', 'Moon\n(S-IVB Impact)']
    values = [earth_data['coda_decay_s'] / 60, moon_data['coda_decay_s'] / 60]  # Minutes
    colors = ['#4ecdc4', '#ffe66d']
    bars = ax.bar(categories, values, color=colors, edgecolor='white', linewidth=2)
    ax.set_ylabel('Coda Decay Time (minutes)')
    ax.set_title('Coda Duration', fontweight='bold')
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                f'{val:.0f} min', ha='center', fontsize=12, fontweight='bold')

    # 4. Synthetic seismograms
    ax = axes[1, 1]
    t_earth = np.linspace(0, 600, 1000)  # 10 minutes
    t_moon = np.linspace(0, 600, 1000)

    # Earth: sharp arrivals, quick decay
    earth_signal = (np.exp(-((t_earth - 30) / 3)**2) * np.sin(2*np.pi*2*t_earth) +
                   np.exp(-((t_earth - 60) / 8)**2) * np.sin(2*np.pi*0.8*t_earth) * 2 +
                   np.exp(-(t_earth - 80) / 100) * (t_earth > 80) * np.random.randn(len(t_earth)) * 0.3)
    earth_signal = earth_signal / np.max(np.abs(earth_signal))

    ax.plot(t_earth, earth_signal + 1.5, color='#4ecdc4', linewidth=0.5, label='Earth')

    # Moon: emergent, long ringing
    moon_signal = np.exp(-((t_moon - 200) / 100)**2) * np.exp(-np.maximum(t_moon - 200, 0) / 300)
    for i in range(50):
        freq = np.random.uniform(0.3, 1.0)
        moon_signal += moon_signal * np.sin(2*np.pi*freq*t_moon + np.random.uniform(0, 2*np.pi)) / np.sqrt(i+1)
    moon_signal = moon_signal / np.max(np.abs(moon_signal))

    ax.plot(t_moon, moon_signal - 1.5, color='#ffe66d', linewidth=0.5, label='Moon')
    ax.set_xlabel('Time (seconds)')
    ax.set_ylabel('Amplitude')
    ax.set_title('Waveform Comparison (10 min window)', fontweight='bold')
    ax.legend(loc='upper right')
    ax.set_xlim(0, 600)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.3)

    # 5. Key metrics table
    ax = axes[1, 2]
    ax.axis('off')

    table_data = [
        ['Metric', 'Earth', 'Moon', 'Ratio'],
        ['Coherent', f'{earth_data["coherent"]}', f'{moon_data["coherent"]}', f'{moon_data["coherent"]/earth_data["coherent"]:.2f}x'],
        ['Scattered', f'{earth_data["scattered"]}', f'{moon_data["scattered"]}', f'{moon_data["scattered"]/earth_data["scattered"]:.2f}x'],
        ['Coda (min)', f'{earth_data["coda_decay_s"]/60:.0f}', f'{moon_data["coda_decay_s"]/60:.0f}', f'{moon_data["coda_decay_s"]/earth_data["coda_decay_s"]:.0f}x'],
        ['Coherence', f'{earth_data["coherence_ratio"]*100:.1f}%', f'{moon_data["coherence_ratio"]*100:.1f}%', '-'],
    ]

    table = ax.table(cellText=table_data, loc='center', cellLoc='center',
                     colWidths=[0.3, 0.2, 0.2, 0.2])
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.2, 2)

    for i in range(len(table_data[0])):
        table[(0, i)].set_facecolor('#333333')
        table[(0, i)].set_text_props(fontweight='bold', color='#00ff88')

    ax.set_title('Summary Metrics', fontweight='bold', y=0.95)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'earth_vs_moon_comparison.png', dpi=150, bbox_inches='tight',
                facecolor='#0a0a0a', edgecolor='none')
    print(f"  Saved: {FIGURES_DIR / 'earth_vs_moon_comparison.png'}")
    plt.close()


def plot_lunar_interior():
    """Visualize lunar interior model."""
    print("Creating lunar interior model...")

    # Load model
    model_file = RESULTS_DIR / 'lunar_interior_model.json'
    if model_file.exists():
        with open(model_file) as f:
            model = json.load(f)
    else:
        # Use default
        model = {
            'layers': [
                {'name': 'Megaregolith', 'depth_range_km': [0, 25], 'vp_km_s': 2.0, 'Q': 100},
                {'name': 'Upper Crust', 'depth_range_km': [25, 45], 'vp_km_s': 5.5, 'Q': 4000},
                {'name': 'Lower Crust', 'depth_range_km': [45, 60], 'vp_km_s': 6.8, 'Q': 5000},
                {'name': 'Upper Mantle', 'depth_range_km': [60, 500], 'vp_km_s': 7.7, 'Q': 6000},
                {'name': 'Middle Mantle', 'depth_range_km': [500, 1000], 'vp_km_s': 8.0, 'Q': 1500},
                {'name': 'Lower Mantle', 'depth_range_km': [1000, 1400], 'vp_km_s': 8.2, 'Q': 500},
                {'name': 'Core', 'depth_range_km': [1400, 1737], 'vp_km_s': 4.5, 'Q': 1000},
            ],
            'total_radius_km': 1737
        }

    fig = plt.figure(figsize=(16, 8))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.2, 1])

    # Left: Cross-section
    ax1 = fig.add_subplot(gs[0])

    colors = ['#8B4513', '#A0522D', '#CD853F', '#2E8B57', '#3CB371', '#90EE90', '#FFD700']
    R = model['total_radius_km']

    # Draw layers as concentric circles
    for i, layer in enumerate(model['layers']):
        r_outer = R - layer['depth_range_km'][0]
        r_inner = R - layer['depth_range_km'][1]

        circle = plt.Circle((0, 0), r_outer, color=colors[i], alpha=0.8)
        ax1.add_patch(circle)

        # Label
        r_mid = (r_outer + r_inner) / 2
        angle = np.pi / 4 + i * 0.15
        x = r_mid * np.cos(angle) * 0.7
        y = r_mid * np.sin(angle) * 0.7
        ax1.annotate(layer['name'], (x, y), fontsize=9, ha='center', va='center',
                    fontweight='bold', color='white',
                    bbox=dict(boxstyle='round', facecolor='black', alpha=0.5))

    # Core (special)
    core_r = R - 1400
    ax1.add_patch(plt.Circle((0, 0), core_r, color='#FFD700', alpha=1))
    ax1.text(0, 0, 'CORE\n350 km', ha='center', va='center', fontsize=10, fontweight='bold')

    # Deep moonquake zone
    theta = np.linspace(0, 2*np.pi, 100)
    r_inner = R - 1100
    r_outer = R - 700
    ax1.plot(r_inner * np.cos(theta), r_inner * np.sin(theta), 'r--', linewidth=2, alpha=0.7)
    ax1.plot(r_outer * np.cos(theta), r_outer * np.sin(theta), 'r--', linewidth=2, alpha=0.7)
    ax1.text(R - 900, -200, 'Deep\nMoonquake\nZone', color='red', fontsize=9, ha='center')

    ax1.set_xlim(-R*1.1, R*1.1)
    ax1.set_ylim(-R*1.1, R*1.1)
    ax1.set_aspect('equal')
    ax1.axis('off')
    ax1.set_title('LUNAR INTERIOR MODEL\n(Based on Apollo + VMS Topology)', fontsize=14, fontweight='bold', color='#00ff88')

    # Right: Velocity and Q profiles
    ax2 = fig.add_subplot(gs[1])

    depths = []
    vps = []
    Qs = []

    for layer in model['layers']:
        d1, d2 = layer['depth_range_km']
        depths.extend([d1, d2])
        vps.extend([layer['vp_km_s'], layer['vp_km_s']])
        Qs.extend([layer['Q'], layer['Q']])

    ax2.plot(vps, depths, 'c-', linewidth=3, label='Vp (km/s)')
    ax2.set_xlabel('P-wave Velocity (km/s)', color='cyan')
    ax2.set_ylabel('Depth (km)')
    ax2.invert_yaxis()
    ax2.set_xlim(0, 10)
    ax2.tick_params(axis='x', colors='cyan')

    ax2_twin = ax2.twiny()
    ax2_twin.plot(Qs, depths, 'y-', linewidth=3, label='Q factor')
    ax2_twin.set_xlabel('Q Factor (attenuation)', color='yellow')
    ax2_twin.set_xlim(0, 7000)
    ax2_twin.tick_params(axis='x', colors='yellow')

    # Annotations
    ax2.axhspan(700, 1100, alpha=0.2, color='red', label='Moonquake zone')
    ax2.axhline(45, color='white', linestyle=':', alpha=0.5)
    ax2.text(9.5, 45, 'Moho', fontsize=8, va='center')

    ax2.set_title('Velocity & Attenuation Profile', fontsize=12, fontweight='bold')
    ax2.legend(loc='lower left')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'lunar_interior_model.png', dpi=150, bbox_inches='tight',
                facecolor='#0a0a0a', edgecolor='none')
    print(f"  Saved: {FIGURES_DIR / 'lunar_interior_model.png'}")
    plt.close()


def plot_apollo_waveforms():
    """Plot Apollo seismogram waveforms."""
    print("Creating Apollo waveforms plot...")

    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    fig.suptitle('APOLLO LUNAR SEISMOGRAMS (VMS Analysis)', fontsize=14, fontweight='bold', color='#00ff88')

    apollo_files = [
        ('sivb_impact_realistic', 'S-IVB Impact (Artificial)', '#ff6b6b'),
        ('deep_moonquake_realistic', 'Deep Moonquake (700-900 km)', '#4ecdc4'),
        ('meteorite_impact_realistic', 'Meteorite Impact', '#ffe66d'),
    ]

    for ax, (name, title, color) in zip(axes, apollo_files):
        filepath = DATA_DIR / 'apollo' / f'{name}.npy'
        sr_file = DATA_DIR / 'apollo' / f'{name}_sr.npy'

        if filepath.exists():
            data = np.load(filepath)
            sr = np.load(sr_file)[0]
            t = np.arange(len(data)) / sr / 3600  # Hours

            ax.plot(t, data, color=color, linewidth=0.3, alpha=0.8)
            ax.fill_between(t, data, 0, alpha=0.3, color=color)
            ax.set_ylabel('Amplitude')
            ax.set_title(title, fontweight='bold', color=color)
            ax.set_xlim(0, max(t))
            ax.grid(True, alpha=0.3)

            # Envelope
            from scipy.signal import hilbert
            envelope = np.abs(hilbert(data))
            # Smooth
            window = int(sr * 60)  # 1 minute window
            if window > 1:
                envelope = np.convolve(envelope, np.ones(window)/window, mode='same')
            ax.plot(t, envelope, 'w-', linewidth=1.5, alpha=0.7, label='Envelope')
            ax.plot(t, -envelope, 'w-', linewidth=1.5, alpha=0.7)
        else:
            ax.text(0.5, 0.5, 'Data not found', transform=ax.transAxes, ha='center')

    axes[-1].set_xlabel('Time (hours)')

    # Add annotation
    axes[0].annotate('Moon "rings like a bell"\nfor HOURS after impact',
                    xy=(2, 0.5), fontsize=11, color='white',
                    bbox=dict(boxstyle='round', facecolor='black', alpha=0.7))

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'apollo_waveforms.png', dpi=150, bbox_inches='tight',
                facecolor='#0a0a0a', edgecolor='none')
    print(f"  Saved: {FIGURES_DIR / 'apollo_waveforms.png'}")
    plt.close()


def plot_iris_earthquakes():
    """Plot IRIS earthquake catalog."""
    print("Creating IRIS earthquake map...")

    catalog_file = RESULTS_DIR / 'iris_2024_m6plus.json'
    if not catalog_file.exists():
        print("  Catalog not found, skipping...")
        return

    with open(catalog_file) as f:
        events = json.load(f)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle('IRIS EARTHQUAKE CATALOG 2024 (M6.0+)', fontsize=14, fontweight='bold', color='#00ff88')

    # Left: Map
    ax = axes[0]

    lats = [e['lat'] for e in events]
    lons = [e['lon'] for e in events]
    mags = [e['mag'] for e in events]
    depths = [e['depth_km'] for e in events]

    # Size by magnitude
    sizes = [(m - 5.5) ** 3 * 50 for m in mags]

    # Color by depth
    scatter = ax.scatter(lons, lats, s=sizes, c=depths, cmap='plasma',
                        alpha=0.7, edgecolors='white', linewidth=0.5)

    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Earthquake Locations', fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Depth (km)')

    # Annotate largest
    largest = max(events, key=lambda e: e['mag'])
    ax.annotate(f"M{largest['mag']:.1f}\n{largest['region'][:20]}",
               (largest['lon'], largest['lat']),
               xytext=(10, 10), textcoords='offset points',
               fontsize=8, color='red',
               arrowprops=dict(arrowstyle='->', color='red'))

    # Right: Statistics
    ax = axes[1]

    # Magnitude distribution
    ax.hist(mags, bins=np.arange(6, 8, 0.2), color='#4ecdc4', edgecolor='white', alpha=0.8)
    ax.set_xlabel('Magnitude')
    ax.set_ylabel('Count')
    ax.set_title('Magnitude Distribution', fontweight='bold')
    ax.axvline(np.mean(mags), color='red', linestyle='--', label=f'Mean: {np.mean(mags):.1f}')
    ax.legend()

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'iris_earthquakes_2024.png', dpi=150, bbox_inches='tight',
                facecolor='#0a0a0a', edgecolor='none')
    print(f"  Saved: {FIGURES_DIR / 'iris_earthquakes_2024.png'}")
    plt.close()


def plot_vms_summary():
    """Create summary dashboard."""
    print("Creating VMS summary dashboard...")

    fig = plt.figure(figsize=(16, 12))
    fig.suptitle('VMS TOPOLOGICAL ANALYSIS - COMPLETE RESULTS', fontsize=18, fontweight='bold', color='#00ff88', y=0.98)

    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.3)

    # 1. Coherence comparison bar chart
    ax1 = fig.add_subplot(gs[0, 0])
    categories = ['Audio\n(air)', 'LIGO\n(GW)', 'Earth\n(seismic)', 'Moon\n(Apollo)']
    coherence = [10, 6, 24.5, 83.3]
    colors = ['#ff6b6b', '#9b59b6', '#4ecdc4', '#ffe66d']
    bars = ax1.bar(categories, coherence, color=colors, edgecolor='white', linewidth=2)
    ax1.set_ylabel('Coherence Ratio (%)')
    ax1.set_title('Topological Coherence by Signal Type', fontweight='bold')
    ax1.set_ylim(0, 100)
    for bar, val in zip(bars, coherence):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                f'{val:.0f}%', ha='center', fontsize=10, fontweight='bold')

    # 2. Coda duration
    ax2 = fig.add_subplot(gs[0, 1])
    categories = ['Audio', 'Earth', 'Moon']
    coda_min = [0.001, 6, 240]  # minutes
    colors = ['#ff6b6b', '#4ecdc4', '#ffe66d']
    bars = ax2.bar(categories, coda_min, color=colors, edgecolor='white', linewidth=2)
    ax2.set_ylabel('Coda Duration (minutes)')
    ax2.set_title('Signal Persistence', fontweight='bold')
    ax2.set_yscale('log')
    for bar, val in zip(bars, coda_min):
        label = f'{val:.0f}' if val >= 1 else f'{val*60:.0f}ms'
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.5,
                label, ha='center', fontsize=10, fontweight='bold')

    # 3. VMS suitability radar
    ax3 = fig.add_subplot(gs[0, 2], projection='polar')
    categories = ['Coherence', 'Duration', 'Structure', 'Separability', 'SNR']

    # Values for each signal type (0-1 scale)
    audio_vals = [0.1, 0.01, 0.3, 0.5, 0.3]
    earth_vals = [0.25, 0.3, 0.7, 0.6, 0.7]
    moon_vals = [0.85, 1.0, 0.9, 0.8, 0.5]

    angles = np.linspace(0, 2*np.pi, len(categories), endpoint=False).tolist()
    angles += angles[:1]

    audio_vals += audio_vals[:1]
    earth_vals += earth_vals[:1]
    moon_vals += moon_vals[:1]

    ax3.plot(angles, audio_vals, 'o-', color='#ff6b6b', linewidth=2, label='Audio')
    ax3.fill(angles, audio_vals, color='#ff6b6b', alpha=0.2)
    ax3.plot(angles, earth_vals, 'o-', color='#4ecdc4', linewidth=2, label='Earth')
    ax3.fill(angles, earth_vals, color='#4ecdc4', alpha=0.2)
    ax3.plot(angles, moon_vals, 'o-', color='#ffe66d', linewidth=2, label='Moon')
    ax3.fill(angles, moon_vals, color='#ffe66d', alpha=0.2)

    ax3.set_xticks(angles[:-1])
    ax3.set_xticklabels(categories, size=8)
    ax3.set_title('VMS Suitability', fontweight='bold', y=1.1)
    ax3.legend(loc='upper right', bbox_to_anchor=(1.3, 1))

    # 4. Apollo waveform example
    ax4 = fig.add_subplot(gs[1, :2])

    filepath = DATA_DIR / 'apollo' / 'sivb_impact_realistic.npy'
    if filepath.exists():
        data = np.load(filepath)
        sr = 6.625
        t = np.arange(len(data)) / sr / 3600
        ax4.plot(t, data, color='#ffe66d', linewidth=0.3, alpha=0.8)
        ax4.fill_between(t, data, 0, alpha=0.2, color='#ffe66d')
        ax4.set_xlabel('Time (hours)')
        ax4.set_ylabel('Amplitude')
        ax4.set_title('Apollo S-IVB Impact: Moon Rings for 4 HOURS', fontweight='bold', color='#ffe66d')
        ax4.set_xlim(0, 4)
        ax4.grid(True, alpha=0.3)

    # 5. Lunar interior mini
    ax5 = fig.add_subplot(gs[1, 2])

    R = 1737
    layers = [
        (0, 60, '#8B4513', 'Crust'),
        (60, 500, '#2E8B57', 'Upper Mantle'),
        (500, 1400, '#3CB371', 'Lower Mantle'),
        (1400, 1737, '#FFD700', 'Core'),
    ]

    for d1, d2, color, name in layers:
        r = R - d1
        circle = plt.Circle((0, 0), r, color=color, alpha=0.8)
        ax5.add_patch(circle)

    ax5.set_xlim(-R*1.2, R*1.2)
    ax5.set_ylim(-R*1.2, R*1.2)
    ax5.set_aspect('equal')
    ax5.axis('off')
    ax5.set_title('Lunar Interior', fontweight='bold')

    # 6. Key findings text
    ax6 = fig.add_subplot(gs[2, :])
    ax6.axis('off')

    findings = """
╔═══════════════════════════════════════════════════════════════════════════════════════════════════╗
║                                    KEY FINDINGS - VMS TOPOLOGICAL ANALYSIS                          ║
╠═══════════════════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                                     ║
║  1. LUNAR SEISMOLOGY: VMS ideal → 83% coherence, 4-hour coda, connected topological structures     ║
║                                                                                                     ║
║  2. EARTH SEISMOLOGY: VMS effective → 24% coherence, clear P/S/surface wave "mountains"           ║
║                                                                                                     ║
║  3. GRAVITATIONAL WAVES: VMS limited → 6% coherence, detector noise dominates                      ║
║                                                                                                     ║
║  4. AUDIO CLEANING: VMS moderate → Best for solid-medium transmission (laser microphones)          ║
║                                                                                                     ║
║  RECOMMENDATION: Apply VMS topological methods to full Apollo PSE dataset (1969-1977)              ║
║                  for refined lunar interior model and undiscovered event detection                  ║
║                                                                                                     ║
╚═══════════════════════════════════════════════════════════════════════════════════════════════════╝
"""
    ax6.text(0.5, 0.5, findings, transform=ax6.transAxes, fontsize=9,
            fontfamily='monospace', ha='center', va='center',
            bbox=dict(boxstyle='round', facecolor='#1a1a1a', edgecolor='#00ff88', linewidth=2))

    plt.savefig(FIGURES_DIR / 'vms_summary_dashboard.png', dpi=150, bbox_inches='tight',
                facecolor='#0a0a0a', edgecolor='none')
    print(f"  Saved: {FIGURES_DIR / 'vms_summary_dashboard.png'}")
    plt.close()


def main():
    print("=" * 70)
    print("GENERATING VMS ANALYSIS VISUALIZATIONS")
    print("=" * 70)

    plot_earth_vs_moon_comparison()
    plot_lunar_interior()
    plot_apollo_waveforms()
    plot_iris_earthquakes()
    plot_vms_summary()

    print("\n" + "=" * 70)
    print(f"All figures saved to: {FIGURES_DIR}")
    print("=" * 70)

    # List files
    for f in sorted(FIGURES_DIR.glob('*.png')):
        size_kb = f.stat().st_size / 1024
        print(f"  {f.name}: {size_kb:.0f} KB")


if __name__ == "__main__":
    main()
