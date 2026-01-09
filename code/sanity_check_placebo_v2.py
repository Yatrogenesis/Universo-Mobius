#!/usr/bin/env python3
"""
SANITY CHECK #1 (VERSIÓN RIGUROSA): TEST PLACEBO
=================================================

PROBLEMA DETECTADO: Con tolerancia 15%, el ruido coloreado encuentra 3.4/4 modos.
SOLUCIÓN: Análisis más riguroso con métricas de calidad.

NUEVO ENFOQUE:
1. Tolerancia reducida a 5%
2. Verificar que los RATIOS entre modos sean correctos (no solo frecuencias individuales)
3. Computar SNR de los picos
4. Comparar contra datos LIGO REALES

Autor: Francisco Molina Burgos & Claude
Fecha: 2026-01-09
"""

import numpy as np
from scipy import signal
from scipy.fft import fft, fftfreq
import matplotlib.pyplot as plt
import os
import json
import h5py

np.random.seed(42)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures')
DATA_DIR = os.path.join(BASE_DIR, 'data', 'ligo')

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# Parámetros
SAMPLE_RATE = 4096  # Hz
DURATION = 32  # segundos
N_SAMPLES = SAMPLE_RATE * DURATION

# Frecuencias hexagonales predichas
F1_BASE = 33.8  # Hz (fundamental para M_total ~ 65 M_sun)
HEX_RATIOS = [1.0, np.sqrt(3), 2.0, np.sqrt(7)]
PREDICTED_MODES = [F1_BASE * r for r in HEX_RATIOS]

# Tolerancia MÁS ESTRICTA
STRICT_TOLERANCE = 0.05  # 5% en lugar de 15%

print("=" * 70)
print("SANITY CHECK #1 (RIGUROSO): TEST PLACEBO")
print("Tolerancia reducida: 5% (era 15%)")
print("=" * 70)

print(f"\nFrecuencias hexagonales predichas:")
for i, (r, f) in enumerate(zip(HEX_RATIOS, PREDICTED_MODES)):
    label = ['f1', 'f2=sqrt(3)*f1', 'f3=2*f1', 'f4=sqrt(7)*f1'][i]
    print(f"  {label} = {f:.2f} Hz")


def generate_pure_gaussian_noise(n_samples, sample_rate):
    """Ruido gaussiano puro."""
    return np.random.normal(0, 1e-21, n_samples)


def generate_colored_noise(n_samples, sample_rate):
    """Ruido coloreado 1/f."""
    white = np.random.normal(0, 1, n_samples)
    b, a = signal.butter(2, 0.01, btype='low')
    pink_component = signal.filtfilt(b, a, white) * 1e-20
    high_freq = np.random.normal(0, 1e-22, n_samples)
    return pink_component + high_freq


def load_real_ligo_data(detector='H1'):
    """Carga datos LIGO reales de GW170814."""
    event_files = {
        'H1': 'H-H1_GWOSC_4KHZ_R1-1186741846-32.hdf5',
        'L1': 'L-L1_GWOSC_4KHZ_R1-1186741846-32.hdf5',
        'V1': 'V-V1_GWOSC_4KHZ_R1-1186741846-32.hdf5'
    }

    filepath = os.path.join(DATA_DIR, event_files[detector])
    if not os.path.exists(filepath):
        print(f"  [!] Archivo no encontrado: {filepath}")
        return None, None

    with h5py.File(filepath, 'r') as f:
        strain = f['strain']['Strain'][:]
        # Calcular sample rate desde metadata o longitud
        try:
            sample_rate = int(f['strain']['Strain'].attrs['sample_rate'])
        except:
            sample_rate = len(strain) // 32  # 32 segundos de datos

    return strain, sample_rate


def analyze_spectrum_rigorous(strain, sample_rate, f1_search_range=(25, 50)):
    """
    Análisis RIGUROSO del espectro.

    En lugar de buscar modos individuales, buscamos el PATRÓN COMPLETO:
    1. Encontrar candidatos para f1 en rango esperado
    2. Para cada candidato, verificar si existen picos en f1*sqrt(3), f1*2, f1*sqrt(7)
    3. Computar chi-cuadrado del patrón
    """
    # Filtro pasabanda
    nyq = sample_rate / 2
    low = 20 / nyq
    high = 500 / nyq
    b, a = signal.butter(4, [low, high], btype='band')
    strain_filtered = signal.filtfilt(b, a, strain)

    # PSD de alta resolución
    nperseg = min(8192, len(strain_filtered) // 2)
    freqs, psd = signal.welch(strain_filtered, fs=sample_rate, nperseg=nperseg)

    # Encontrar TODOS los picos significativos
    mask = (freqs >= 20) & (freqs <= 300)
    psd_masked = psd[mask]
    freqs_masked = freqs[mask]

    # Normalizar PSD para SNR
    psd_median = np.median(psd_masked)
    psd_std = np.std(psd_masked)
    snr = (psd_masked - psd_median) / psd_std

    # Encontrar picos con SNR > 2
    peaks_idx, properties = signal.find_peaks(snr, height=2.0, distance=3)

    if len(peaks_idx) == 0:
        return {
            'best_f1': None,
            'pattern_chi2': np.inf,
            'n_matches': 0,
            'matches': [],
            'freqs': freqs,
            'psd': psd,
            'snr': snr,
            'freqs_masked': freqs_masked
        }

    peak_freqs = freqs_masked[peaks_idx]
    peak_snrs = snr[peaks_idx]

    # Buscar el mejor patrón hexagonal
    best_chi2 = np.inf
    best_f1 = None
    best_matches = []

    # Buscar candidatos para f1 en el rango esperado
    f1_candidates = peak_freqs[(peak_freqs >= f1_search_range[0]) &
                               (peak_freqs <= f1_search_range[1])]

    for f1_cand in f1_candidates:
        # Predecir modos basados en este f1
        predicted = [f1_cand * r for r in HEX_RATIOS]

        # Buscar coincidencias
        matches = []
        chi2 = 0

        for i, f_pred in enumerate(predicted):
            # Encontrar pico más cercano
            distances = np.abs(peak_freqs - f_pred)
            closest_idx = np.argmin(distances)
            closest_freq = peak_freqs[closest_idx]
            closest_snr = peak_snrs[closest_idx]
            rel_diff = abs(closest_freq - f_pred) / f_pred

            if rel_diff < STRICT_TOLERANCE:
                matches.append({
                    'mode': i + 1,
                    'predicted': f_pred,
                    'observed': closest_freq,
                    'rel_diff': rel_diff,
                    'snr': closest_snr
                })
                chi2 += (rel_diff / STRICT_TOLERANCE) ** 2
            else:
                chi2 += 10  # Penalización por modo faltante

        # Verificar ratios entre modos encontrados
        if len(matches) >= 2:
            # Penalizar si los ratios no son correctos
            for j in range(len(matches) - 1):
                for k in range(j + 1, len(matches)):
                    m1, m2 = matches[j], matches[k]
                    ratio_obs = m2['observed'] / m1['observed']
                    ratio_pred = HEX_RATIOS[m2['mode']-1] / HEX_RATIOS[m1['mode']-1]
                    ratio_diff = abs(ratio_obs - ratio_pred) / ratio_pred
                    chi2 += ratio_diff ** 2 * 10

        if chi2 < best_chi2:
            best_chi2 = chi2
            best_f1 = f1_cand
            best_matches = matches

    return {
        'best_f1': best_f1,
        'pattern_chi2': best_chi2,
        'n_matches': len(best_matches),
        'matches': best_matches,
        'freqs': freqs,
        'psd': psd,
        'snr': snr,
        'freqs_masked': freqs_masked,
        'all_peaks': list(zip(peak_freqs.tolist(), peak_snrs.tolist()))
    }


def compute_hexagonal_score(analysis):
    """
    Computa un score de calidad hexagonal.
    Combina: número de modos, chi2, y SNR promedio de matches.
    """
    if analysis['n_matches'] == 0:
        return 0.0

    # Score base por número de modos
    n_score = analysis['n_matches'] / 4.0

    # Penalizar por chi2 alto
    chi2_penalty = 1.0 / (1.0 + analysis['pattern_chi2'])

    # Bonus por SNR alto
    avg_snr = np.mean([m['snr'] for m in analysis['matches']])
    snr_bonus = min(1.0, avg_snr / 5.0)

    score = n_score * chi2_penalty * (1 + snr_bonus)
    return score


# Ejecutar test placebo riguroso
print("\n" + "=" * 70)
print("FASE 1: Test con ruido (100 realizaciones)")
print("=" * 70)

n_trials = 100
scores_gaussian = []
scores_colored = []
matches_gaussian = []
matches_colored = []

for i in range(n_trials):
    # Ruido gaussiano
    noise_gauss = generate_pure_gaussian_noise(N_SAMPLES, SAMPLE_RATE)
    result_gauss = analyze_spectrum_rigorous(noise_gauss, SAMPLE_RATE)
    scores_gaussian.append(compute_hexagonal_score(result_gauss))
    matches_gaussian.append(result_gauss['n_matches'])

    # Ruido coloreado
    noise_colored = generate_colored_noise(N_SAMPLES, SAMPLE_RATE)
    result_colored = analyze_spectrum_rigorous(noise_colored, SAMPLE_RATE)
    scores_colored.append(compute_hexagonal_score(result_colored))
    matches_colored.append(result_colored['n_matches'])

    if (i + 1) % 20 == 0:
        print(f"  Trial {i+1}/{n_trials}: Gauss={result_gauss['n_matches']}/4 "
              f"(score={scores_gaussian[-1]:.3f}), "
              f"Colored={result_colored['n_matches']}/4 "
              f"(score={scores_colored[-1]:.3f})")

# Análisis de datos LIGO reales
print("\n" + "=" * 70)
print("FASE 2: Análisis de datos LIGO REALES (GW170814)")
print("=" * 70)

real_results = {}
for detector in ['H1', 'L1', 'V1']:
    print(f"\n  Analizando {detector}...")
    strain, sr = load_real_ligo_data(detector)
    if strain is not None:
        result = analyze_spectrum_rigorous(strain, sr)
        score = compute_hexagonal_score(result)
        real_results[detector] = {
            'n_matches': result['n_matches'],
            'score': score,
            'best_f1': result['best_f1'],
            'matches': result['matches'],
            'chi2': result['pattern_chi2']
        }
        print(f"    f1 = {result['best_f1']:.1f} Hz" if result['best_f1'] else "    No f1 encontrado")
        print(f"    Modos: {result['n_matches']}/4, Score: {score:.3f}, Chi2: {result['pattern_chi2']:.2f}")
        for m in result['matches']:
            print(f"      Modo {m['mode']}: pred={m['predicted']:.1f}Hz, "
                  f"obs={m['observed']:.1f}Hz, diff={m['rel_diff']*100:.1f}%, SNR={m['snr']:.1f}")

# Estadísticas
print("\n" + "=" * 70)
print("RESULTADOS ESTADÍSTICOS (tolerancia 5%)")
print("=" * 70)

mean_gauss_matches = np.mean(matches_gaussian)
mean_colored_matches = np.mean(matches_colored)
mean_gauss_score = np.mean(scores_gaussian)
mean_colored_score = np.mean(scores_colored)
std_gauss_score = np.std(scores_gaussian)
std_colored_score = np.std(scores_colored)

print(f"\n  Ruido Gaussiano Puro ({n_trials} trials):")
print(f"    Modos (media): {mean_gauss_matches:.2f}/4")
print(f"    Score (media): {mean_gauss_score:.4f} +/- {std_gauss_score:.4f}")

print(f"\n  Ruido Coloreado 1/f ({n_trials} trials):")
print(f"    Modos (media): {mean_colored_matches:.2f}/4")
print(f"    Score (media): {mean_colored_score:.4f} +/- {std_colored_score:.4f}")

# Comparar con datos reales
print(f"\n  Datos LIGO Reales (GW170814):")
for det, res in real_results.items():
    print(f"    {det}: {res['n_matches']}/4 modos, Score={res['score']:.4f}")

# Calcular significancia
if real_results:
    best_real_score = max(r['score'] for r in real_results.values())
    best_real_det = max(real_results.keys(), key=lambda k: real_results[k]['score'])

    # Z-score respecto a ruido coloreado (más conservador)
    z_vs_colored = (best_real_score - mean_colored_score) / std_colored_score
    z_vs_gaussian = (best_real_score - mean_gauss_score) / std_gauss_score

    # P-value (fracción de ruido con score >= real)
    p_colored = np.sum(np.array(scores_colored) >= best_real_score) / n_trials
    p_gaussian = np.sum(np.array(scores_gaussian) >= best_real_score) / n_trials

    print(f"\n  Significancia del mejor detector ({best_real_det}):")
    print(f"    Z vs Gaussiano: {z_vs_gaussian:.2f}")
    print(f"    Z vs Coloreado: {z_vs_colored:.2f}")
    print(f"    P-value vs Gaussiano: {p_gaussian:.4f}")
    print(f"    P-value vs Coloreado: {p_colored:.4f}")

# VEREDICTO
print("\n" + "=" * 70)
print("VEREDICTO FINAL")
print("=" * 70)

# Criterios más estrictos
passes_test = False
if real_results:
    # Pasar si:
    # 1. El score real es > 2 sigma sobre ruido coloreado
    # 2. Al menos un detector tiene >= 3 modos con chi2 < 5
    cond1 = z_vs_colored > 2.0
    cond2 = any(r['n_matches'] >= 3 and r['chi2'] < 5 for r in real_results.values())
    passes_test = cond1 and cond2

if passes_test:
    verdict = "PASA: Senal hexagonal SIGNIFICATIVA sobre ruido"
    is_honest = True
else:
    verdict = "FALLA: No hay separacion significativa del ruido"
    is_honest = False

print(f"\n  {verdict}")
print(f"\n  Interpretacion:")
if is_honest:
    print(f"    - Los datos LIGO reales muestran patron hexagonal REAL")
    print(f"    - Z = {z_vs_colored:.1f} sigma sobre ruido coloreado")
    print(f"    - El algoritmo NO alucina con criterios estrictos")
else:
    print(f"    - Con tolerancia 5%, la senal no es significativamente diferente del ruido")
    print(f"    - Esto puede indicar:")
    print(f"      a) La senal es marginal y requiere mas datos")
    print(f"      b) Los modos hexagonales necesitan analisis mas sofisticado")
    print(f"      c) El patron observado es parcialmente ruido")

# Crear figura
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel 1: Distribución de scores
ax = axes[0, 0]
bins = np.linspace(0, 0.8, 30)
ax.hist(scores_gaussian, bins=bins, alpha=0.7, label='Gaussiano', color='blue', density=True)
ax.hist(scores_colored, bins=bins, alpha=0.5, label='Coloreado 1/f', color='green', density=True)
if real_results:
    for det, res in real_results.items():
        ax.axvline(res['score'], ls='--', lw=2, label=f'{det} real')
ax.set_xlabel('Score Hexagonal')
ax.set_ylabel('Densidad')
ax.set_title('Distribucion de Scores: Ruido vs Datos Reales')
ax.legend()

# Panel 2: Histograma de modos
ax = axes[0, 1]
bins = np.arange(-0.5, 5.5, 1)
ax.hist(matches_gaussian, bins=bins, alpha=0.7, label='Gaussiano', color='blue')
ax.hist(matches_colored, bins=bins, alpha=0.5, label='Coloreado 1/f', color='green')
if real_results:
    for det, res in real_results.items():
        ax.axvline(res['n_matches'], ls='--', lw=2, label=f'{det}={res["n_matches"]}/4')
ax.set_xlabel('Modos hexagonales (5% tolerancia)')
ax.set_ylabel('Frecuencia')
ax.set_title(f'Modos encontrados ({n_trials} trials)')
ax.legend()
ax.set_xticks([0, 1, 2, 3, 4])

# Panel 3: Comparación visual de PSDs
ax = axes[1, 0]
# Generar un ejemplo de cada tipo
noise_ex = generate_colored_noise(N_SAMPLES, SAMPLE_RATE)
result_ex = analyze_spectrum_rigorous(noise_ex, SAMPLE_RATE)
ax.semilogy(result_ex['freqs'], result_ex['psd'], 'g-', alpha=0.5, label='Ruido coloreado')

if real_results and 'H1' in real_results:
    strain_h1, sr_h1 = load_real_ligo_data('H1')
    if strain_h1 is not None:
        result_h1 = analyze_spectrum_rigorous(strain_h1, sr_h1)
        ax.semilogy(result_h1['freqs'], result_h1['psd'], 'b-', alpha=0.8, label='H1 real')

for f in PREDICTED_MODES:
    ax.axvline(f, color='red', ls='--', alpha=0.7)
ax.set_xlim(20, 200)
ax.set_xlabel('Frecuencia (Hz)')
ax.set_ylabel('PSD')
ax.set_title('PSD: Ruido Coloreado vs LIGO Real')
ax.legend()
ax.grid(True, alpha=0.3)

# Panel 4: Resumen
ax = axes[1, 1]
ax.axis('off')

summary_lines = [
    "SANITY CHECK #1 (RIGUROSO)",
    "=" * 40,
    "",
    "PARAMETROS:",
    f"  Tolerancia: 5% (antes: 15%)",
    f"  Trials: {n_trials}",
    f"  Verificacion de ratios: SI",
    "",
    f"RUIDO GAUSSIANO:",
    f"  Modos: {mean_gauss_matches:.2f}/4",
    f"  Score: {mean_gauss_score:.4f}",
    "",
    f"RUIDO COLOREADO:",
    f"  Modos: {mean_colored_matches:.2f}/4",
    f"  Score: {mean_colored_score:.4f}",
    "",
    "LIGO REAL (GW170814):"
]

for det, res in real_results.items():
    summary_lines.append(f"  {det}: {res['n_matches']}/4, Score={res['score']:.4f}")

summary_lines.extend([
    "",
    "=" * 40,
    f"Z vs Coloreado: {z_vs_colored:.2f}" if real_results else "",
    f"P-value: {p_colored:.4f}" if real_results else "",
    "",
    f"VEREDICTO: {'PASA' if is_honest else 'FALLA'}",
])

summary = "\n".join(summary_lines)
ax.text(0.05, 0.95, summary, transform=ax.transAxes, fontsize=10,
       verticalalignment='top', fontfamily='monospace',
       bbox=dict(boxstyle='round',
                 facecolor='lightgreen' if is_honest else 'lightyellow',
                 alpha=0.8))

plt.suptitle('SANITY CHECK #1 (RIGUROSO): Test Placebo con Tolerancia 5%',
             fontsize=14, fontweight='bold')
plt.tight_layout()

# Guardar
for fmt in ['png', 'pdf']:
    filepath = os.path.join(FIGURES_DIR, f'sanity_check_placebo_rigorous.{fmt}')
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
print(f"\nFigura guardada: sanity_check_placebo_rigorous.png/pdf")
plt.close()

# Guardar resultados
output = {
    'test': 'Sanity Check #1 (Rigorous): Placebo Test',
    'tolerance': STRICT_TOLERANCE,
    'n_trials': n_trials,
    'gaussian': {
        'mean_matches': float(mean_gauss_matches),
        'mean_score': float(mean_gauss_score),
        'std_score': float(std_gauss_score)
    },
    'colored': {
        'mean_matches': float(mean_colored_matches),
        'mean_score': float(mean_colored_score),
        'std_score': float(std_colored_score)
    },
    'real_ligo': {det: {k: float(v) if isinstance(v, (int, float)) else v
                        for k, v in res.items()}
                  for det, res in real_results.items()},
    'significance': {
        'z_vs_gaussian': float(z_vs_gaussian) if real_results else None,
        'z_vs_colored': float(z_vs_colored) if real_results else None,
        'p_gaussian': float(p_gaussian) if real_results else None,
        'p_colored': float(p_colored) if real_results else None
    },
    'verdict': verdict,
    'passes_test': passes_test
}

filepath = os.path.join(RESULTS_DIR, 'sanity_check_placebo_rigorous.json')
with open(filepath, 'w') as f:
    json.dump(output, f, indent=2)
print(f"Resultados guardados: {filepath}")
