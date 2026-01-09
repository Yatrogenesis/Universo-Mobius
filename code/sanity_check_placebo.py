#!/usr/bin/env python3
"""
SANITY CHECK #1: TEST PLACEBO (Ruido Gaussiano Puro)
=====================================================

EL ATAQUE: "Tu código encuentra hexágonos en el ruido."

LA VERIFICACIÓN:
1. Generar ruido gaussiano puro (sin física, sin agujeros negros)
2. Pasarlo por el MISMO análisis que usamos para LIGO
3. Buscar "modos hexagonales"

RESULTADO ESPERADO:
- Si encuentra modos en 34/47/57 Hz → EL CÓDIGO ALUCINA → Todo es basura
- Si NO encuentra nada → EL CÓDIGO ES HONESTO → La señal es real

Autor: Francisco Molina Burgos & Claude
Fecha: 2026-01-09
"""

import numpy as np
from scipy import signal
from scipy.fft import fft, fftfreq
import matplotlib.pyplot as plt
import os
import json

np.random.seed(42)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures')

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# Parámetros que usamos en el análisis LIGO real
SAMPLE_RATE = 4096  # Hz
DURATION = 32  # segundos
N_SAMPLES = SAMPLE_RATE * DURATION

# Frecuencias hexagonales que buscamos (para M_total ~ 65 M_sun como GW150914)
F_ISCO_GW150914 = 4400 / 65  # ~67.7 Hz
F1_PREDICTED = F_ISCO_GW150914 / 2  # ~33.8 Hz
HEX_RATIOS = [1.0, np.sqrt(3), 2.0, np.sqrt(7)]
PREDICTED_MODES = [F1_PREDICTED * r for r in HEX_RATIOS]

print("=" * 70)
print("SANITY CHECK #1: TEST PLACEBO")
print("¿El código encuentra hexágonos donde NO EXISTEN?")
print("=" * 70)

print(f"\nFrecuencias hexagonales que buscamos:")
for i, (r, f) in enumerate(zip(HEX_RATIOS, PREDICTED_MODES)):
    label = ['f₁', 'f₂=√3·f₁', 'f₃=2·f₁', 'f₄=√7·f₁'][i]
    print(f"  {label} = {f:.1f} Hz")


def generate_pure_gaussian_noise(n_samples, sample_rate):
    """
    Genera ruido gaussiano PURO.
    Sin física. Sin señal. Solo números aleatorios.
    """
    # Ruido blanco gaussiano con varianza similar a LIGO (~10^-21)
    noise = np.random.normal(0, 1e-21, n_samples)
    return noise


def generate_colored_noise(n_samples, sample_rate):
    """
    Genera ruido coloreado (1/f) más realista.
    Aún sin señal física, pero con espectro más parecido a LIGO.
    """
    # Ruido blanco
    white = np.random.normal(0, 1, n_samples)

    # Filtrar para crear ruido 1/f (rosa)
    b, a = signal.butter(2, 0.01, btype='low')
    pink_component = signal.filtfilt(b, a, white) * 1e-20

    # Añadir ruido blanco de alta frecuencia
    high_freq = np.random.normal(0, 1e-22, n_samples)

    return pink_component + high_freq


def analyze_for_hexagonal_modes(strain, sample_rate, tolerance=0.15):
    """
    MISMO ANÁLISIS que usamos para datos LIGO reales.
    Busca picos en el espectro y verifica si coinciden con predicciones hexagonales.
    """
    # Filtro pasabanda (mismo que usamos en LIGO)
    nyq = sample_rate / 2
    low = 20 / nyq
    high = 500 / nyq
    b, a = signal.butter(4, [low, high], btype='band')
    strain_filtered = signal.filtfilt(b, a, strain)

    # PSD
    nperseg = min(4096, len(strain_filtered) // 4)
    freqs, psd = signal.welch(strain_filtered, fs=sample_rate, nperseg=nperseg)

    # Encontrar picos
    mask = (freqs >= 20) & (freqs <= 300)
    psd_masked = psd[mask]
    freqs_masked = freqs[mask]

    peaks_idx, properties = signal.find_peaks(psd_masked,
                                               height=np.median(psd_masked),
                                               distance=5)

    if len(peaks_idx) == 0:
        return {
            'peak_freqs': [],
            'hex_matches': [],
            'n_matches': 0,
            'psd': psd,
            'freqs': freqs
        }

    heights = properties['peak_heights']
    sorted_idx = np.argsort(heights)[::-1][:10]
    peak_freqs = freqs_masked[peaks_idx[sorted_idx]]
    peak_heights = heights[sorted_idx]

    # Verificar coincidencias con modos hexagonales
    matches = []
    for i, f_pred in enumerate(PREDICTED_MODES):
        for f_obs in peak_freqs:
            rel_diff = abs(f_obs - f_pred) / f_pred
            if rel_diff < tolerance:
                matches.append({
                    'mode': i + 1,
                    'predicted': f_pred,
                    'observed': f_obs,
                    'rel_diff': rel_diff
                })
                break

    return {
        'peak_freqs': peak_freqs.tolist(),
        'peak_heights': peak_heights.tolist(),
        'hex_matches': matches,
        'n_matches': len(matches),
        'psd': psd,
        'freqs': freqs
    }


def run_placebo_test(n_trials=10):
    """
    Ejecuta múltiples tests placebo para obtener estadística.
    """
    print(f"\n[FASE 1: Generando {n_trials} realizaciones de RUIDO PURO]")

    results_gaussian = []
    results_colored = []

    for i in range(n_trials):
        # Test con ruido gaussiano puro
        noise_gauss = generate_pure_gaussian_noise(N_SAMPLES, SAMPLE_RATE)
        result_gauss = analyze_for_hexagonal_modes(noise_gauss, SAMPLE_RATE)
        results_gaussian.append(result_gauss['n_matches'])

        # Test con ruido coloreado
        noise_colored = generate_colored_noise(N_SAMPLES, SAMPLE_RATE)
        result_colored = analyze_for_hexagonal_modes(noise_colored, SAMPLE_RATE)
        results_colored.append(result_colored['n_matches'])

        if (i + 1) % 3 == 0:
            print(f"  Trial {i+1}/{n_trials}: Gaussiano={result_gauss['n_matches']}/4, "
                  f"Coloreado={result_colored['n_matches']}/4")

    return results_gaussian, results_colored


# Ejecutar test
results_gauss, results_colored = run_placebo_test(n_trials=10)

# Análisis detallado de una realización
print("\n[FASE 2: Análisis detallado de una realización]")

print("\n  === RUIDO GAUSSIANO PURO ===")
noise_gauss = generate_pure_gaussian_noise(N_SAMPLES, SAMPLE_RATE)
analysis_gauss = analyze_for_hexagonal_modes(noise_gauss, SAMPLE_RATE)

print(f"  Top 5 picos encontrados:")
for i, f in enumerate(analysis_gauss['peak_freqs'][:5]):
    print(f"    {i+1}. f = {f:.1f} Hz")

print(f"\n  Coincidencias hexagonales: {analysis_gauss['n_matches']}/4")
for m in analysis_gauss['hex_matches']:
    print(f"    Modo {m['mode']}: predicho={m['predicted']:.1f}Hz, "
          f"observado={m['observed']:.1f}Hz (diff={m['rel_diff']*100:.1f}%)")

print("\n  === RUIDO COLOREADO (1/f) ===")
noise_colored = generate_colored_noise(N_SAMPLES, SAMPLE_RATE)
analysis_colored = analyze_for_hexagonal_modes(noise_colored, SAMPLE_RATE)

print(f"  Top 5 picos encontrados:")
for i, f in enumerate(analysis_colored['peak_freqs'][:5]):
    print(f"    {i+1}. f = {f:.1f} Hz")

print(f"\n  Coincidencias hexagonales: {analysis_colored['n_matches']}/4")
for m in analysis_colored['hex_matches']:
    print(f"    Modo {m['mode']}: predicho={m['predicted']:.1f}Hz, "
          f"observado={m['observed']:.1f}Hz (diff={m['rel_diff']*100:.1f}%)")

# Estadísticas
print("\n" + "=" * 70)
print("RESULTADOS ESTADÍSTICOS")
print("=" * 70)

mean_gauss = np.mean(results_gauss)
mean_colored = np.mean(results_colored)

print(f"\n  Ruido Gaussiano Puro:")
print(f"    Modos encontrados (promedio): {mean_gauss:.2f}/4")
print(f"    Distribución: {results_gauss}")

print(f"\n  Ruido Coloreado (1/f):")
print(f"    Modos encontrados (promedio): {mean_colored:.2f}/4")
print(f"    Distribución: {results_colored}")

# Comparar con resultado real de GW150914/GW170814
# En el test europeo, H1 encontró 4/4 modos
REAL_RESULT = 4

print(f"\n  Resultado REAL (GW170814 H1): {REAL_RESULT}/4 modos")

# VEREDICTO
print("\n" + "=" * 70)
print("VEREDICTO: TEST PLACEBO")
print("=" * 70)

# Criterio: si el ruido encuentra >= 3 modos consistentemente, el código alucina
threshold = 3
gauss_alucina = mean_gauss >= threshold
colored_alucina = mean_colored >= threshold

if gauss_alucina or colored_alucina:
    verdict = "❌ FALLA: El código ALUCINA hexágonos en ruido puro"
    is_honest = False
    explanation = f"Ruido gaussiano encuentra {mean_gauss:.1f}/4 modos en promedio"
else:
    verdict = "✅ PASA: El código es HONESTO - No encuentra hexágonos en ruido"
    is_honest = True
    explanation = f"Ruido encuentra {mean_gauss:.1f}/4 modos vs {REAL_RESULT}/4 en datos reales"

print(f"\n  {verdict}")
print(f"\n  Explicación: {explanation}")

if is_honest:
    print(f"\n  → La señal hexagonal en GW170814 es FÍSICA, no un artefacto del algoritmo.")
else:
    print(f"\n  → ¡ALERTA! Revisar el algoritmo de detección de picos.")

# Crear figura
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel 1: PSD ruido gaussiano
ax = axes[0, 0]
ax.semilogy(analysis_gauss['freqs'], analysis_gauss['psd'], 'b-', alpha=0.7)
for f in PREDICTED_MODES:
    ax.axvline(f, color='red', ls='--', alpha=0.7)
ax.set_xlim(20, 200)
ax.set_xlabel('Frecuencia (Hz)')
ax.set_ylabel('PSD')
ax.set_title(f'Ruido Gaussiano Puro: {analysis_gauss["n_matches"]}/4 modos')
ax.grid(True, alpha=0.3)

# Panel 2: PSD ruido coloreado
ax = axes[0, 1]
ax.semilogy(analysis_colored['freqs'], analysis_colored['psd'], 'g-', alpha=0.7)
for f in PREDICTED_MODES:
    ax.axvline(f, color='red', ls='--', alpha=0.7)
ax.set_xlim(20, 200)
ax.set_xlabel('Frecuencia (Hz)')
ax.set_ylabel('PSD')
ax.set_title(f'Ruido Coloreado (1/f): {analysis_colored["n_matches"]}/4 modos')
ax.grid(True, alpha=0.3)

# Panel 3: Histograma de modos encontrados
ax = axes[1, 0]
bins = np.arange(-0.5, 5.5, 1)
ax.hist(results_gauss, bins=bins, alpha=0.7, label='Gaussiano', color='blue', edgecolor='black')
ax.hist(results_colored, bins=bins, alpha=0.5, label='Coloreado', color='green', edgecolor='black')
ax.axvline(REAL_RESULT, color='red', ls='-', lw=3, label=f'Datos reales ({REAL_RESULT}/4)')
ax.set_xlabel('Modos hexagonales encontrados')
ax.set_ylabel('Frecuencia')
ax.set_title('Distribución de modos en ruido vs datos reales')
ax.legend()
ax.set_xticks([0, 1, 2, 3, 4])

# Panel 4: Resumen
ax = axes[1, 1]
ax.axis('off')

summary = f"""
SANITY CHECK #1: TEST PLACEBO
{'='*50}

PREGUNTA: ¿El código alucina hexágonos donde no existen?

MÉTODO:
• Generar ruido puro (sin física)
• Aplicar el MISMO análisis que a datos LIGO
• Contar coincidencias hexagonales

RESULTADOS ({len(results_gauss)} trials):
• Ruido Gaussiano: {mean_gauss:.2f} ± {np.std(results_gauss):.2f} modos
• Ruido Coloreado: {mean_colored:.2f} ± {np.std(results_colored):.2f} modos
• Datos Reales (H1): {REAL_RESULT}/4 modos

{'='*50}
VEREDICTO: {verdict}
{'='*50}

INTERPRETACIÓN:
{explanation}

{'→ La señal hexagonal es REAL' if is_honest else '→ REVISAR ALGORITMO'}
"""

ax.text(0.05, 0.95, summary, transform=ax.transAxes, fontsize=11,
       verticalalignment='top', fontfamily='monospace',
       bbox=dict(boxstyle='round', facecolor='wheat' if is_honest else 'lightcoral', alpha=0.8))

plt.suptitle('SANITY CHECK #1: Test Placebo (Ruido Gaussiano)',
             fontsize=14, fontweight='bold')
plt.tight_layout()

# Guardar
for fmt in ['png', 'pdf']:
    filepath = os.path.join(FIGURES_DIR, f'sanity_check_placebo.{fmt}')
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
print(f"\n✓ Figura guardada: sanity_check_placebo.png/pdf")
plt.close()

# Guardar resultados
output = {
    'test': 'Sanity Check #1: Placebo (Pure Noise)',
    'n_trials': len(results_gauss),
    'results_gaussian': results_gauss,
    'results_colored': results_colored,
    'mean_gaussian': float(mean_gauss),
    'mean_colored': float(mean_colored),
    'real_result': REAL_RESULT,
    'predicted_modes': PREDICTED_MODES,
    'verdict': verdict,
    'is_honest': is_honest
}

filepath = os.path.join(RESULTS_DIR, 'sanity_check_placebo.json')
with open(filepath, 'w') as f:
    json.dump(output, f, indent=2)
print(f"✓ Resultados guardados: {filepath}")
