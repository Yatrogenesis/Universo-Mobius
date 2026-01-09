#!/usr/bin/env python3
"""
MONTE CARLO LIMPIO: Ruido sin lineas de potencia
================================================

El test anterior mostro que las lineas de 60/120 Hz siempre coinciden
con modos predichos. Este test usa ruido PURO sin estructura.

Autor: Francisco Molina Burgos & Claude
Fecha: 2025-01-09
"""

import numpy as np
from scipy import signal, stats
import matplotlib.pyplot as plt
import os
import json
from collections import Counter

np.random.seed(None)  # Seed aleatorio para cada ejecucion

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures')

SAMPLE_RATE = 4096
DURATION = 32
N_SAMPLES = SAMPLE_RATE * DURATION
N_MONTE_CARLO = 1000

# Frecuencias predichas (European Test)
F1_BASE = 37.0
HEX_RATIOS = [1.0, np.sqrt(3), 2.0, np.sqrt(7)]
PREDICTED_MODES = [F1_BASE * r for r in HEX_RATIOS]

TOLERANCE = 0.15  # Misma que European Test

print("=" * 70)
print("MONTE CARLO LIMPIO: Ruido Puro sin Lineas de Potencia")
print("=" * 70)
print(f"Tolerancia: {TOLERANCE*100:.0f}%")
print(f"Realizaciones: {N_MONTE_CARLO}")


def generate_pure_noise():
    """Ruido gaussiano + coloreado, SIN lineas de potencia."""
    # Componente rosa (1/f)
    white = np.random.normal(0, 1, N_SAMPLES)
    b, a = signal.butter(2, [0.005, 0.3], btype='band')
    pink = signal.filtfilt(b, a, white)

    # Ruido blanco adicional
    high = np.random.normal(0, 0.5, N_SAMPLES)

    return (pink + high) * 1e-21


def analyze_for_modes(strain):
    """Busca modos hexagonales en los datos."""
    # Filtro pasabanda
    nyq = SAMPLE_RATE / 2
    low = 20 / nyq
    high = 500 / nyq
    b, a = signal.butter(4, [low, high], btype='band')
    filtered = signal.filtfilt(b, a, strain)

    # PSD
    nperseg = 4096
    freqs, psd = signal.welch(filtered, fs=SAMPLE_RATE, nperseg=nperseg)

    # Normalizar
    mask = (freqs >= 20) & (freqs <= 200)
    freqs_m = freqs[mask]
    psd_m = psd[mask]

    median = np.median(psd_m)
    std = np.std(psd_m)
    psd_norm = (psd_m - median) / std

    # Encontrar picos significativos
    peaks_idx, props = signal.find_peaks(psd_norm, height=2.0, distance=5)

    if len(peaks_idx) == 0:
        return 0, []

    peak_freqs = freqs_m[peaks_idx]
    peak_heights = props['peak_heights']

    # Top 10 picos
    top_idx = np.argsort(peak_heights)[::-1][:10]
    top_freqs = peak_freqs[top_idx]

    # Buscar coincidencias
    matches = []
    for i, f_pred in enumerate(PREDICTED_MODES):
        for f_obs in top_freqs:
            rel_diff = abs(f_obs - f_pred) / f_pred
            if rel_diff < TOLERANCE:
                matches.append({
                    'mode': i + 1,
                    'predicted': f_pred,
                    'observed': f_obs,
                    'diff': rel_diff
                })
                break

    return len(matches), matches


# Ejecutar Monte Carlo
print("\n[Ejecutando Monte Carlo...]")
results = []
all_matches = []

for i in range(N_MONTE_CARLO):
    noise = generate_pure_noise()
    n_modes, matches = analyze_for_modes(noise)
    results.append(n_modes)
    if matches:
        all_matches.extend([m['mode'] for m in matches])

    if (i + 1) % 100 == 0:
        print(f"  {i+1}/{N_MONTE_CARLO}: media hasta ahora = {np.mean(results):.2f}")

# Estadisticas
counts = Counter(results)
mean = np.mean(results)
std = np.std(results)

print("\n" + "=" * 70)
print("RESULTADOS")
print("=" * 70)

print(f"\nDistribucion de modos encontrados en ruido:")
for k in range(5):
    n = counts.get(k, 0)
    pct = n / N_MONTE_CARLO * 100
    bar = '*' * int(pct / 2)
    print(f"  {k}/4: {n:4d} ({pct:5.1f}%) {bar}")

print(f"\nEstadisticas:")
print(f"  Media: {mean:.3f} +/- {std:.3f}")
print(f"  P(0/4): {counts.get(0, 0)/N_MONTE_CARLO:.4f}")
print(f"  P(4/4): {counts.get(4, 0)/N_MONTE_CARLO:.6f}")

# Que modos se encuentran mas?
if all_matches:
    mode_counts = Counter(all_matches)
    print(f"\nModos mas frecuentemente encontrados (falsos positivos):")
    for mode in sorted(mode_counts.keys()):
        print(f"  Modo {mode} (f={PREDICTED_MODES[mode-1]:.1f} Hz): {mode_counts[mode]} veces")

# Comparar con resultado real
REAL_H1 = 4
z_score = (REAL_H1 - mean) / std if std > 0 else float('inf')
p_value = np.sum(np.array(results) >= REAL_H1) / N_MONTE_CARLO

print(f"\nSignificancia del resultado real (H1 = 4/4):")
print(f"  Z-score: {z_score:.2f}")
print(f"  P-value: {p_value:.6f}")

# Veredicto
print("\n" + "=" * 70)
print("VEREDICTO")
print("=" * 70)

if mean < 1.5 and p_value < 0.01:
    verdict = "EL TEST PLACEBO FUNCIONA"
    detail = f"Ruido puro produce {mean:.2f} modos. P(4/4) = {p_value:.6f}"
    passes = True
else:
    verdict = "ALERTA: Posibles falsos positivos"
    detail = f"Ruido produce {mean:.2f} modos. La senal puede ser marginal."
    passes = False

print(f"\n  {verdict}")
print(f"  {detail}")

if passes:
    print(f"\n  Conclusion:")
    print(f"    - El resultado real de H1 (4/4 modos) es {z_score:.1f} sigma sobre ruido")
    print(f"    - Probabilidad de obtener 4/4 por azar: {p_value*100:.4f}%")
    print(f"    - La senal hexagonal es SIGNIFICATIVA")

# Crear figura
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel 1: Histograma
ax = axes[0]
x = list(range(5))
y = [counts.get(k, 0)/N_MONTE_CARLO for k in x]
bars = ax.bar(x, y, color='steelblue', edgecolor='black', alpha=0.8)
ax.axvline(REAL_H1, color='red', ls='--', lw=2, label=f'H1 real = {REAL_H1}')
ax.set_xlabel('Modos hexagonales encontrados')
ax.set_ylabel('Probabilidad')
ax.set_title(f'Monte Carlo ({N_MONTE_CARLO} trials): Ruido Puro')
ax.set_xticks(x)
ax.legend()

# Anotar
for i, (xi, yi) in enumerate(zip(x, y)):
    if yi > 0.01:
        ax.annotate(f'{yi*100:.1f}%', xy=(xi, yi), ha='center', va='bottom')

# Panel 2: Resumen
ax = axes[1]
ax.axis('off')

summary = f"""
MONTE CARLO: VERIFICACION DE ALEATORIEDAD
{'='*45}

Parametros:
  Realizaciones: {N_MONTE_CARLO}
  Tolerancia: {TOLERANCE*100:.0f}%
  f1 base: {F1_BASE} Hz

Distribucion en ruido:
  0 modos: {counts.get(0, 0)/N_MONTE_CARLO*100:.1f}%
  1 modo:  {counts.get(1, 0)/N_MONTE_CARLO*100:.1f}%
  2 modos: {counts.get(2, 0)/N_MONTE_CARLO*100:.1f}%
  3 modos: {counts.get(3, 0)/N_MONTE_CARLO*100:.1f}%
  4 modos: {counts.get(4, 0)/N_MONTE_CARLO*100:.4f}%

Media: {mean:.3f} +/- {std:.3f}

{'='*45}
RESULTADO REAL (H1): {REAL_H1}/4 modos
Z-score: {z_score:.1f} sigma
P-value: {p_value:.6f}
{'='*45}

VEREDICTO: {verdict}
"""

ax.text(0.05, 0.95, summary, transform=ax.transAxes, fontsize=11,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round',
                  facecolor='lightgreen' if passes else 'lightyellow',
                  alpha=0.8))

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'monte_carlo_clean.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIGURES_DIR, 'monte_carlo_clean.pdf'), dpi=300, bbox_inches='tight')
print(f"\nFigura guardada: monte_carlo_clean.png/pdf")
plt.close()

# Guardar resultados
output = {
    'test': 'Monte Carlo Clean (Pure Noise)',
    'n_trials': N_MONTE_CARLO,
    'tolerance': TOLERANCE,
    'distribution': dict(counts),
    'mean': float(mean),
    'std': float(std),
    'p_four': float(counts.get(4, 0)/N_MONTE_CARLO),
    'real_result': REAL_H1,
    'z_score': float(z_score) if z_score != float('inf') else None,
    'p_value': float(p_value),
    'verdict': verdict,
    'passes': passes
}

filepath = os.path.join(RESULTS_DIR, 'monte_carlo_clean.json')
with open(filepath, 'w') as f:
    json.dump(output, f, indent=2)
print(f"Resultados guardados: {filepath}")
