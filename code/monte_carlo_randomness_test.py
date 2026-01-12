#!/usr/bin/env python3
"""
VERIFICACION MONTE CARLO DE ALEATORIEDAD DEL TEST PLACEBO
==========================================================

OBJETIVO: Verificar que la distribucion de "modos hexagonales encontrados"
en ruido puro sigue una distribucion probabilistica conocida.

Si el placebo esta funcionando correctamente:
1. La distribucion de modos debe seguir una binomial/Poisson
2. NO debe haber sesgo sistematico hacia ningun numero de modos
3. La tasa de falsos positivos debe ser predecible

METODO:
1. Generar 1000 realizaciones de ruido
2. Contar modos hexagonales en cada una
3. Comparar con distribucion teorica esperada
4. Calcular probabilidad de 4/4 modos por AZAR

Autor: Francisco Molina Burgos 
Fecha: 2026-01-09
"""

import numpy as np
from scipy import signal, stats
import matplotlib.pyplot as plt
import os
import json
from collections import Counter

np.random.seed(42)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures')

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# Parametros
SAMPLE_RATE = 4096
DURATION = 32
N_SAMPLES = SAMPLE_RATE * DURATION
N_MONTE_CARLO = 1000

# Frecuencias hexagonales
F1_BASE = 37.0  # Hz (valor observado en European Test)
HEX_RATIOS = [1.0, np.sqrt(3), 2.0, np.sqrt(7)]
PREDICTED_MODES = [F1_BASE * r for r in HEX_RATIOS]

# Tolerancias a probar
TOLERANCES = [0.05, 0.10, 0.15, 0.20]

print("=" * 70)
print("VERIFICACION MONTE CARLO DE ALEATORIEDAD")
print("=" * 70)
print(f"\nRealizaciones: {N_MONTE_CARLO}")
print(f"Frecuencias predichas: {[f'{f:.1f}' for f in PREDICTED_MODES]} Hz")
print(f"Tolerancias a probar: {[f'{t*100:.0f}%' for t in TOLERANCES]}")


def generate_ligo_like_noise(n_samples, sample_rate):
    """
    Genera ruido con caracteristicas espectrales similares a LIGO.
    Incluye: ruido 1/f, lineas de 60Hz, y ruido blanco.
    """
    # Ruido base 1/f
    white = np.random.normal(0, 1, n_samples)
    b, a = signal.butter(2, [0.01, 0.5], btype='band')
    pink = signal.filtfilt(b, a, white)

    # Ruido blanco de alta frecuencia
    high_freq = np.random.normal(0, 0.3, n_samples)

    # Lineas de potencia (60 Hz y armonicos) - DEBILES para no dominar
    t = np.arange(n_samples) / sample_rate
    power_lines = 0.1 * (np.sin(2*np.pi*60*t) + 0.5*np.sin(2*np.pi*120*t))

    # Combinar y escalar
    noise = (pink + high_freq + power_lines) * 1e-21

    return noise


def analyze_for_modes(strain, sample_rate, tolerance):
    """
    Busca modos hexagonales en los datos.
    Retorna numero de modos encontrados y detalles.
    """
    # Filtro pasabanda
    nyq = sample_rate / 2
    low = 20 / nyq
    high = 500 / nyq
    b, a = signal.butter(4, [low, high], btype='band')
    filtered = signal.filtfilt(b, a, strain)

    # PSD
    nperseg = min(4096, len(filtered) // 4)
    freqs, psd = signal.welch(filtered, fs=sample_rate, nperseg=nperseg)

    # Normalizar y encontrar picos
    mask = (freqs >= 20) & (freqs <= 200)
    freqs_m = freqs[mask]
    psd_m = psd[mask]

    psd_norm = (psd_m - np.median(psd_m)) / np.std(psd_m)

    peaks_idx, props = signal.find_peaks(psd_norm, height=1.5, distance=5)

    if len(peaks_idx) == 0:
        return 0, []

    peak_freqs = freqs_m[peaks_idx]
    peak_heights = props['peak_heights']

    # Tomar los 15 picos mas altos
    top_idx = np.argsort(peak_heights)[::-1][:15]
    top_freqs = peak_freqs[top_idx]

    # Contar coincidencias con modos predichos
    matches = []
    for i, f_pred in enumerate(PREDICTED_MODES):
        for f_obs in top_freqs:
            rel_diff = abs(f_obs - f_pred) / f_pred
            if rel_diff < tolerance:
                matches.append(i + 1)
                break

    return len(matches), matches


def run_monte_carlo(tolerance):
    """Ejecuta simulacion Monte Carlo para una tolerancia dada."""
    results = []

    for i in range(N_MONTE_CARLO):
        noise = generate_ligo_like_noise(N_SAMPLES, SAMPLE_RATE)
        n_modes, _ = analyze_for_modes(noise, SAMPLE_RATE, tolerance)
        results.append(n_modes)

        if (i + 1) % 200 == 0:
            print(f"    Progreso: {i+1}/{N_MONTE_CARLO}")

    return results


def theoretical_probability(n_modes, tolerance, n_peaks=15, n_freq_bins=180):
    """
    Calcula la probabilidad teorica de encontrar n modos por azar.

    Modelo: cada modo predicho tiene probabilidad p de coincidir
    con un pico aleatorio dentro de la tolerancia.

    p = (n_peaks * 2 * tolerance) / (rango de frecuencia / resolucion)
    """
    # Probabilidad de coincidencia para un modo
    freq_range = 200 - 20  # Hz
    resolution = freq_range / n_freq_bins
    window_per_mode = 2 * tolerance * 37  # Hz (ancho de tolerancia para f~37)
    n_bins_in_window = window_per_mode / resolution

    # Probabilidad de que AL MENOS un pico caiga en la ventana
    p_match = 1 - (1 - n_bins_in_window/n_freq_bins) ** n_peaks
    p_match = min(p_match, 0.99)  # Cap para evitar problemas numericos

    # Distribucion binomial para 4 modos independientes
    from scipy.stats import binom
    prob = binom.pmf(n_modes, 4, p_match)

    return p_match, prob


# Ejecutar Monte Carlo para cada tolerancia
print("\n" + "=" * 70)
print("EJECUTANDO SIMULACIONES MONTE CARLO")
print("=" * 70)

all_results = {}

for tol in TOLERANCES:
    print(f"\n  Tolerancia {tol*100:.0f}%:")
    results = run_monte_carlo(tol)
    all_results[tol] = results

    # Estadisticas
    counts = Counter(results)
    mean = np.mean(results)
    std = np.std(results)

    print(f"    Media: {mean:.2f} modos")
    print(f"    Std: {std:.2f}")
    print(f"    Distribucion: {dict(sorted(counts.items()))}")

    # Probabilidad de 4/4
    p_four = counts.get(4, 0) / N_MONTE_CARLO
    print(f"    P(4/4 modos): {p_four:.4f} ({counts.get(4, 0)}/{N_MONTE_CARLO})")

# Analisis de aleatoriedad
print("\n" + "=" * 70)
print("VERIFICACION DE ALEATORIEDAD")
print("=" * 70)

# Para cada tolerancia, verificar si la distribucion es consistente con binomial
for tol in TOLERANCES:
    results = all_results[tol]
    counts = Counter(results)

    print(f"\n  Tolerancia {tol*100:.0f}%:")

    # Ajustar distribucion binomial
    # Estimar p desde los datos
    mean_obs = np.mean(results)
    p_est = mean_obs / 4.0  # 4 modos posibles

    print(f"    p estimada (por modo): {p_est:.3f}")

    # Comparar con distribucion teorica
    p_match, _ = theoretical_probability(0, tol)
    print(f"    p teorica: {p_match:.3f}")

    # Test chi-cuadrado
    observed = [counts.get(k, 0) for k in range(5)]
    expected = [stats.binom.pmf(k, 4, p_est) * N_MONTE_CARLO for k in range(5)]

    # Solo bins con expectativa > 5
    valid_bins = [i for i, e in enumerate(expected) if e > 5]
    if len(valid_bins) >= 2:
        obs_valid = [observed[i] for i in valid_bins]
        exp_valid = [expected[i] for i in valid_bins]
        chi2, p_chi = stats.chisquare(obs_valid, exp_valid)
        print(f"    Chi2 test: chi2={chi2:.2f}, p={p_chi:.3f}")
        print(f"    {'PASA' if p_chi > 0.05 else 'FALLA'}: Distribucion {'es' if p_chi > 0.05 else 'NO es'} consistente con binomial")
    else:
        print(f"    Chi2 test: No aplicable (pocos bins)")

# VEREDICTO
print("\n" + "=" * 70)
print("VEREDICTO: VERIFICACION DE ALEATORIEDAD")
print("=" * 70)

# Criterios de aleatoriedad
tol_15 = all_results[0.15]
mean_15 = np.mean(tol_15)
p_four_15 = Counter(tol_15).get(4, 0) / N_MONTE_CARLO

print(f"\n  Con tolerancia 15% (usada en European Test):")
print(f"    Modos en ruido: {mean_15:.2f} +/- {np.std(tol_15):.2f}")
print(f"    P(4/4 por azar): {p_four_15:.4f} = 1 en {int(1/p_four_15) if p_four_15 > 0 else 'infinito'}")

# Comparar con resultado real
REAL_RESULT_H1 = 4  # European Test: H1 encontro 4/4 modos

if p_four_15 > 0:
    z_real = (REAL_RESULT_H1 - mean_15) / np.std(tol_15) if np.std(tol_15) > 0 else np.inf
else:
    z_real = np.inf

print(f"\n  Resultado REAL (H1): {REAL_RESULT_H1}/4 modos")
print(f"  Z-score del resultado real: {z_real:.2f}")

# Interpretacion
is_random = mean_15 < 1.5 and p_four_15 < 0.01
if is_random:
    verdict = "PASA: El test placebo es ALEATORIO"
    explanation = f"Ruido encuentra {mean_15:.1f}/4 modos en promedio. P(4/4 por azar) = {p_four_15:.4f}"
    honest = True
else:
    verdict = "ALERTA: El test placebo muestra SESGO"
    explanation = f"Ruido encuentra {mean_15:.1f}/4 modos. Esto sugiere que la tolerancia es muy permisiva."
    honest = False

print(f"\n  {verdict}")
print(f"  {explanation}")

if not honest:
    print(f"\n  RECOMENDACION:")
    print(f"    - Reducir tolerancia a 10% o menos")
    print(f"    - Verificar que los modos encontrados en LIGO tengan SNR alto")
    print(f"    - Usar test de ratios, no solo frecuencias individuales")

# Crear figura
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel 1: Histograma para cada tolerancia
ax = axes[0, 0]
colors = ['blue', 'green', 'orange', 'red']
for i, tol in enumerate(TOLERANCES):
    counts = Counter(all_results[tol])
    x = list(range(5))
    y = [counts.get(k, 0)/N_MONTE_CARLO for k in x]
    ax.bar(np.array(x) + i*0.2 - 0.3, y, width=0.18,
           label=f'{tol*100:.0f}%', color=colors[i], alpha=0.8)
ax.axvline(REAL_RESULT_H1, color='black', ls='--', lw=2, label='H1 real')
ax.set_xlabel('Modos hexagonales encontrados')
ax.set_ylabel('Probabilidad')
ax.set_title(f'Monte Carlo ({N_MONTE_CARLO} trials): Distribucion de Modos en Ruido')
ax.set_xticks([0, 1, 2, 3, 4])
ax.legend()

# Panel 2: Media vs Tolerancia
ax = axes[0, 1]
means = [np.mean(all_results[t]) for t in TOLERANCES]
stds = [np.std(all_results[t]) for t in TOLERANCES]
tol_pct = [t*100 for t in TOLERANCES]
ax.errorbar(tol_pct, means, yerr=stds, fmt='o-', capsize=5, color='blue', label='Ruido')
ax.axhline(REAL_RESULT_H1, color='red', ls='--', label='H1 real (4/4)')
ax.fill_between(tol_pct, [m-s for m,s in zip(means, stds)],
                [m+s for m,s in zip(means, stds)], alpha=0.2)
ax.set_xlabel('Tolerancia (%)')
ax.set_ylabel('Modos encontrados (media)')
ax.set_title('Efecto de la Tolerancia en Falsos Positivos')
ax.legend()
ax.grid(True, alpha=0.3)

# Panel 3: P(4/4) vs Tolerancia
ax = axes[1, 0]
p_fours = [Counter(all_results[t]).get(4, 0)/N_MONTE_CARLO for t in TOLERANCES]
ax.semilogy(tol_pct, p_fours, 'o-', color='purple', markersize=8)
ax.axhline(0.05, color='red', ls='--', label='5% significancia')
ax.axhline(0.01, color='orange', ls='--', label='1% significancia')
ax.set_xlabel('Tolerancia (%)')
ax.set_ylabel('P(4/4 modos por azar)')
ax.set_title('Probabilidad de Encontrar 4/4 Modos por Azar')
ax.legend()
ax.grid(True, alpha=0.3)

# Panel 4: Resumen
ax = axes[1, 1]
ax.axis('off')

summary_lines = [
    "VERIFICACION MONTE CARLO",
    "=" * 40,
    "",
    f"Realizaciones: {N_MONTE_CARLO}",
    f"f1 predicho: {F1_BASE} Hz",
    "",
    "RESULTADOS POR TOLERANCIA:",
]

for tol in TOLERANCES:
    mean = np.mean(all_results[tol])
    p4 = Counter(all_results[tol]).get(4, 0)/N_MONTE_CARLO
    summary_lines.append(f"  {tol*100:.0f}%: {mean:.2f} modos, P(4/4)={p4:.4f}")

summary_lines.extend([
    "",
    "=" * 40,
    f"VEREDICTO: {verdict.split(':')[0]}",
    "=" * 40,
    "",
    f"Con tolerancia 15%:",
    f"  Ruido: {mean_15:.2f} modos",
    f"  H1 real: {REAL_RESULT_H1} modos",
    f"  Z-score: {z_real:.2f}",
    "",
    f"{'ALEATORIO' if honest else 'SESGADO'}"
])

summary = "\n".join(summary_lines)
ax.text(0.05, 0.95, summary, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round',
                  facecolor='lightgreen' if honest else 'lightyellow',
                  alpha=0.8))

plt.suptitle('Verificacion Monte Carlo de Aleatoriedad del Test Placebo',
             fontsize=14, fontweight='bold')
plt.tight_layout()

# Guardar
for fmt in ['png', 'pdf']:
    filepath = os.path.join(FIGURES_DIR, f'monte_carlo_randomness.{fmt}')
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
print(f"\nFigura guardada: monte_carlo_randomness.png/pdf")
plt.close()

# Guardar resultados
output = {
    'test': 'Monte Carlo Randomness Verification',
    'n_monte_carlo': N_MONTE_CARLO,
    'predicted_modes_hz': PREDICTED_MODES,
    'results_by_tolerance': {
        f'{t*100:.0f}%': {
            'mean': float(np.mean(all_results[t])),
            'std': float(np.std(all_results[t])),
            'distribution': dict(Counter(all_results[t])),
            'p_four': float(Counter(all_results[t]).get(4, 0)/N_MONTE_CARLO)
        }
        for t in TOLERANCES
    },
    'real_result_h1': REAL_RESULT_H1,
    'z_score_real': float(z_real) if z_real != np.inf else None,
    'verdict': verdict,
    'is_random': honest
}

filepath = os.path.join(RESULTS_DIR, 'monte_carlo_randomness.json')
with open(filepath, 'w') as f:
    json.dump(output, f, indent=2)
print(f"Resultados guardados: {filepath}")
