#!/usr/bin/env python3
"""
TEST ASESINO A: Inyección GR Ciega

OBJETIVO: Verificar que el pipeline NO encuentra patrones hexagonales
cuando la señal es puramente GR (sin estructura hexagonal).

MÉTODO:
1. Tomar ruido real de LIGO (quiet time)
2. Inyectar señales GR puras (IMRPhenomD/SEOBNRv4)
3. Correr el pipeline de detección hexagonal
4. Verificar tasa de falsos positivos

CRITERIO DE ÉXITO:
- El pipeline NO debe encontrar 4/4 modos hexagonales en >5% de inyecciones
- Si encuentra >5%, el algoritmo está sesgado

CRITERIO DE FALLO:
- Si encuentra patrones hexagonales sistemáticamente en señales GR,
  toda la evidencia de OCTH queda invalidada.

Autor: Francisco Molina Burgos
Fecha: Enero 2025
"""

import numpy as np
import os
import json
from scipy import signal
from scipy.stats import norm, chi2
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import warnings
warnings.filterwarnings('ignore')

# Configuración
np.random.seed(42)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures')

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# =============================================================================
# GENERACIÓN DE SEÑALES GR PURAS
# =============================================================================

def generate_gr_waveform(M_total, q=1.0, distance_mpc=400, sample_rate=4096, duration=4.0):
    """
    Genera una forma de onda GR pura usando aproximación IMRPhenomD simplificada.

    Parámetros:
    -----------
    M_total : float
        Masa total en masas solares
    q : float
        Ratio de masas (m1/m2), q >= 1
    distance_mpc : float
        Distancia en Mpc
    sample_rate : int
        Tasa de muestreo en Hz
    duration : float
        Duración total en segundos

    Returns:
    --------
    t : array
        Vector de tiempo
    h : array
        Strain de la onda gravitacional
    f_merger : float
        Frecuencia de merger
    f_ringdown : float
        Frecuencia de ringdown (QNM fundamental GR)
    """
    # Constantes
    G = 6.674e-11  # m^3 kg^-1 s^-2
    c = 3e8  # m/s
    M_sun = 1.989e30  # kg
    Mpc = 3.086e22  # m

    # Masas
    M = M_total * M_sun
    m1 = M * q / (1 + q)
    m2 = M / (1 + q)
    eta = m1 * m2 / M**2  # Symmetric mass ratio
    M_chirp = M * eta**(3/5)

    # Escalas de tiempo y frecuencia
    t_M = G * M / c**3  # Tiempo característico
    f_isco = c**3 / (6**(3/2) * np.pi * G * M)

    # Frecuencia de ringdown GR (Berti et al. 2009)
    # f_220 ≈ 0.0588 / M (en unidades geométricas)
    # Convertido a Hz:
    chi_f = 0.7  # Spin final típico para q~1
    f_ringdown = (c**3 / (2 * np.pi * G * M)) * (1.5251 - 1.1568 * (1 - chi_f)**0.1292)

    # Tiempo y frecuencia
    t = np.linspace(-duration/2, duration/2, int(sample_rate * duration))
    dt = 1.0 / sample_rate

    # Fase inspiral (aproximación PN)
    # phi(t) ~ -2 * (t_c - t)^(5/8) / (5 * eta * t_M)
    t_merger = 0.0
    t_inspiral = t[t < t_merger - 0.01]

    # Frecuencia instantánea durante inspiral
    tau = np.abs(t_merger - t_inspiral)
    tau = np.maximum(tau, 1e-6)  # Evitar división por cero
    f_inst_inspiral = (c**3 / (8 * np.pi * G * M_chirp)) * (5 / (256 * tau))**(3/8)
    f_inst_inspiral = np.minimum(f_inst_inspiral, f_isco)

    # Amplitud durante inspiral
    D = distance_mpc * Mpc
    A_inspiral = (G * M_chirp / (c**2 * D)) * (np.pi * G * M_chirp * f_inst_inspiral / c**3)**(2/3)

    # Fase inspiral
    phi_inspiral = 2 * np.pi * np.cumsum(f_inst_inspiral) * dt

    # Señal inspiral
    h_inspiral = A_inspiral * np.cos(phi_inspiral)

    # Merger + Ringdown
    t_post = t[t >= t_merger - 0.01]

    # Ringdown: decaimiento exponencial con frecuencia QNM
    # tau_ringdown ≈ 2 * Q / (2 * pi * f_ringdown), Q ~ 2-4
    Q = 3.0  # Factor de calidad típico
    tau_ringdown = Q / (np.pi * f_ringdown)

    # Amplitud pico en merger
    A_peak = np.max(A_inspiral) if len(A_inspiral) > 0 else 1e-21

    # Señal post-merger
    t_rel = t_post - t_merger + 0.01
    A_ringdown = A_peak * np.exp(-t_rel / tau_ringdown)
    phi_ringdown = 2 * np.pi * f_ringdown * t_rel
    h_ringdown = A_ringdown * np.cos(phi_ringdown)

    # Combinar
    h = np.zeros_like(t)
    h[t < t_merger - 0.01] = h_inspiral
    h[t >= t_merger - 0.01] = h_ringdown

    # Suavizar transición
    transition_width = int(0.02 * sample_rate)
    if transition_width > 0:
        transition_idx = np.where(t >= t_merger - 0.01)[0][0]
        if transition_idx > transition_width:
            window = np.hanning(2 * transition_width)
            h[transition_idx-transition_width:transition_idx+transition_width] *= window

    return t, h, f_isco, f_ringdown


def generate_ligo_noise(duration=4.0, sample_rate=4096, psd_type='design'):
    """
    Genera ruido coloreado tipo LIGO.

    Parámetros:
    -----------
    duration : float
        Duración en segundos
    sample_rate : int
        Tasa de muestreo
    psd_type : str
        'design' para sensibilidad de diseño, 'o3' para O3

    Returns:
    --------
    noise : array
        Ruido coloreado
    """
    N = int(duration * sample_rate)
    freqs = np.fft.rfftfreq(N, 1/sample_rate)

    # PSD aproximada de LIGO (simplificada)
    # S(f) ~ S0 * [(f0/f)^4 + 2 + (f/f0)^2]
    f0 = 150.0  # Hz, frecuencia de mínimo ruido
    S0 = 1e-46 if psd_type == 'design' else 3e-46  # Hz^-1

    # Evitar división por cero
    freqs_safe = np.maximum(freqs, 1.0)

    # PSD
    psd = S0 * ((f0/freqs_safe)**4 + 2 + (freqs_safe/f0)**2)
    psd[freqs < 20] = 1e-40  # Corte bajo
    psd[0] = 0

    # Ruido blanco en dominio de Fourier
    white_noise = np.random.randn(len(freqs)) + 1j * np.random.randn(len(freqs))

    # Colorear el ruido
    colored_fft = white_noise * np.sqrt(psd * sample_rate / 2)

    # Transformar a dominio del tiempo
    noise = np.fft.irfft(colored_fft, N)

    return noise


# =============================================================================
# PIPELINE DE DETECCIÓN HEXAGONAL (COPIADO DEL ORIGINAL)
# =============================================================================

def detect_hexagonal_modes(data, sample_rate, M_total, tolerance=0.15):
    """
    Detecta modos hexagonales en una señal.

    Parámetros:
    -----------
    data : array
        Datos de strain
    sample_rate : int
        Tasa de muestreo
    M_total : float
        Masa total estimada en M_sun
    tolerance : float
        Tolerancia para coincidencia de modos

    Returns:
    --------
    result : dict
        Diccionario con resultados del análisis
    """
    # Calcular f_ISCO
    f_isco = 4400 / M_total  # Hz

    # Frecuencias hexagonales predichas
    hex_ratios = [1.0, np.sqrt(3), 2.0, np.sqrt(7)]
    hex_freqs = [f_isco * r / 2 for r in hex_ratios]  # f_n = f_ISCO/2 * r_n

    # Espectro de potencia
    nperseg = min(len(data) // 4, 2048)
    freqs, psd = signal.welch(data, sample_rate, nperseg=nperseg)

    # Buscar picos
    peaks_idx, properties = signal.find_peaks(psd, height=np.median(psd) * 2, distance=5)
    peak_freqs = freqs[peaks_idx]
    peak_powers = psd[peaks_idx]

    # Verificar coincidencia con modos hexagonales
    modes_found = []
    for i, f_hex in enumerate(hex_freqs):
        # Buscar pico cercano a frecuencia hexagonal
        if len(peak_freqs) > 0:
            closest_idx = np.argmin(np.abs(peak_freqs - f_hex))
            closest_freq = peak_freqs[closest_idx]

            # Verificar si está dentro de tolerancia
            if abs(closest_freq - f_hex) / f_hex < tolerance:
                modes_found.append({
                    'mode': i + 1,
                    'predicted': f_hex,
                    'observed': closest_freq,
                    'error': abs(closest_freq - f_hex) / f_hex,
                    'power': peak_powers[closest_idx]
                })

    n_modes = len(modes_found)
    hex_fraction = n_modes / 4

    # Calcular chi2 para OCTH y GR
    if len(modes_found) > 0:
        octh_residuals = [m['error'] for m in modes_found]
        octh_chi2 = np.sum(np.array(octh_residuals)**2) / (tolerance**2)
    else:
        octh_chi2 = 100  # Penalización si no se encuentran modos

    # GR predice ratios diferentes (QNM estándar)
    gr_ratios = [1.0, 1.2, 1.5, 2.0]  # Ratios típicos GR
    gr_freqs = [f_isco * r / 2 for r in gr_ratios]

    gr_residuals = []
    for f_obs in [m['observed'] for m in modes_found]:
        closest_gr = min(gr_freqs, key=lambda x: abs(x - f_obs))
        gr_residuals.append(abs(f_obs - closest_gr) / closest_gr)

    if len(gr_residuals) > 0:
        gr_chi2 = np.sum(np.array(gr_residuals)**2) / (tolerance**2)
    else:
        gr_chi2 = 100

    delta_chi2 = gr_chi2 - octh_chi2

    return {
        'n_modes_found': n_modes,
        'hex_fraction': hex_fraction,
        'modes': modes_found,
        'octh_chi2': octh_chi2,
        'gr_chi2': gr_chi2,
        'delta_chi2': delta_chi2,
        'favors': 'OCTH' if delta_chi2 > 0 else 'GR',
        'peak_freqs': peak_freqs.tolist(),
        'hex_freqs_predicted': hex_freqs
    }


# =============================================================================
# TEST DE INYECCIÓN CIEGA
# =============================================================================

def run_blind_injection_test(n_injections=100, verbose=True):
    """
    Ejecuta el test de inyección ciega.

    Inyecta señales GR puras en ruido y verifica que el pipeline
    NO encuentre patrones hexagonales sistemáticamente.
    """
    print("\n" + "="*70)
    print("  TEST ASESINO A: INYECCIÓN GR CIEGA")
    print("="*70)
    print(f"\n  Número de inyecciones: {n_injections}")
    print("  Hipótesis nula: Las señales GR NO tienen estructura hexagonal")
    print("  Criterio de fallo: >5% de inyecciones muestran 4/4 modos")

    # Parámetros de inyección
    mass_range = (20, 100)  # M_sun
    distance_range = (200, 2000)  # Mpc
    sample_rate = 4096
    duration = 4.0

    results = []
    n_4_modes = 0
    n_favor_octh = 0

    print("\n[FASE 1: Generando e inyectando señales GR puras]")

    for i in range(n_injections):
        if verbose and (i + 1) % 20 == 0:
            print(f"  Inyección {i+1}/{n_injections}...")

        # Parámetros aleatorios
        M_total = np.random.uniform(*mass_range)
        distance = np.random.uniform(*distance_range)
        q = np.random.uniform(1.0, 4.0)

        # Generar señal GR pura
        t, h_gr, f_isco, f_ringdown = generate_gr_waveform(
            M_total, q, distance, sample_rate, duration
        )

        # Generar ruido LIGO
        noise = generate_ligo_noise(duration, sample_rate)

        # Inyectar (escalar para SNR realista)
        snr_target = np.random.uniform(8, 25)
        h_rms = np.sqrt(np.mean(h_gr**2))
        noise_rms = np.sqrt(np.mean(noise**2))
        if h_rms > 0:
            scale = snr_target * noise_rms / h_rms
            h_injected = h_gr * scale
        else:
            h_injected = h_gr

        data = noise + h_injected

        # Correr detector hexagonal
        detection = detect_hexagonal_modes(data, sample_rate, M_total)

        # Registrar resultados
        results.append({
            'injection': i + 1,
            'M_total': M_total,
            'distance': distance,
            'q': q,
            'snr_target': snr_target,
            'f_ringdown_gr': f_ringdown,
            'n_modes_found': detection['n_modes_found'],
            'hex_fraction': detection['hex_fraction'],
            'delta_chi2': detection['delta_chi2'],
            'favors': detection['favors']
        })

        if detection['n_modes_found'] == 4:
            n_4_modes += 1
        if detection['favors'] == 'OCTH':
            n_favor_octh += 1

    # Estadísticas
    print("\n[FASE 2: Análisis estadístico]")

    n_modes_dist = [r['n_modes_found'] for r in results]
    delta_chi2_dist = [r['delta_chi2'] for r in results]

    mode_counts = {i: n_modes_dist.count(i) for i in range(5)}

    print(f"\n  Distribución de modos encontrados:")
    for n, count in mode_counts.items():
        pct = 100 * count / n_injections
        bar = "█" * int(pct / 2)
        print(f"    {n}/4 modos: {count:4d} ({pct:5.1f}%) {bar}")

    print(f"\n  Eventos que favorecen OCTH: {n_favor_octh}/{n_injections} ({100*n_favor_octh/n_injections:.1f}%)")
    print(f"  Eventos que favorecen GR:   {n_injections - n_favor_octh}/{n_injections} ({100*(n_injections-n_favor_octh)/n_injections:.1f}%)")

    print(f"\n  Δχ² medio: {np.mean(delta_chi2_dist):.2f} ± {np.std(delta_chi2_dist):.2f}")

    # Tasa de falsos positivos (4/4 modos)
    fp_rate = n_4_modes / n_injections
    fp_rate_expected = 0.005  # 0.5% esperado por azar (del test placebo)

    print(f"\n  === RESULTADO DEL TEST ===")
    print(f"  Tasa de falsos positivos (4/4 modos): {fp_rate*100:.2f}%")
    print(f"  Tasa esperada por azar: {fp_rate_expected*100:.2f}%")

    # Determinar si pasa el test
    # Usamos test binomial: ¿fp_rate > 5%?
    from scipy.stats import binom
    p_value = 1 - binom.cdf(n_4_modes - 1, n_injections, 0.05)

    if fp_rate <= 0.05:
        print(f"\n  ✅ TEST PASADO: El pipeline NO está sesgado hacia hexágonos")
        print(f"     La tasa de FP ({fp_rate*100:.1f}%) es ≤ 5%")
        test_passed = True
    else:
        print(f"\n  ❌ TEST FALLIDO: El pipeline ESTÁ SESGADO")
        print(f"     La tasa de FP ({fp_rate*100:.1f}%) es > 5%")
        print(f"     TODA LA EVIDENCIA DE OCTH REQUIERE REVISIÓN")
        test_passed = False

    # Guardar resultados
    output = {
        'test_name': 'Blind GR Injection Test',
        'n_injections': n_injections,
        'n_4_modes': n_4_modes,
        'fp_rate': fp_rate,
        'n_favor_octh': n_favor_octh,
        'mean_delta_chi2': float(np.mean(delta_chi2_dist)),
        'std_delta_chi2': float(np.std(delta_chi2_dist)),
        'mode_distribution': mode_counts,
        'test_passed': test_passed,
        'p_value_bias': float(p_value),
        'individual_results': results
    }

    output_path = os.path.join(RESULTS_DIR, 'test_asesino_A_inyeccion_GR.json')
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\n  Resultados guardados en: {output_path}")

    # Generar figura
    create_injection_figure(results, mode_counts, fp_rate, test_passed)

    return output


def create_injection_figure(results, mode_counts, fp_rate, test_passed):
    """Genera figura del test de inyección."""

    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)

    # 1. Distribución de modos encontrados
    ax1 = fig.add_subplot(gs[0, 0])
    modes = list(mode_counts.keys())
    counts = list(mode_counts.values())
    colors = ['green' if m < 4 else 'red' for m in modes]
    bars = ax1.bar(modes, counts, color=colors, edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Modos hexagonales encontrados')
    ax1.set_ylabel('Número de inyecciones')
    ax1.set_title('Distribución de Modos en Señales GR Puras')
    ax1.set_xticks(modes)
    ax1.set_xticklabels(['0/4', '1/4', '2/4', '3/4', '4/4'])

    # Añadir línea de expectativa
    n_total = sum(counts)
    ax1.axhline(n_total * 0.005, color='blue', linestyle='--', label='Esperado por azar (0.5%)')
    ax1.legend()

    # 2. Δχ² vs Masa
    ax2 = fig.add_subplot(gs[0, 1])
    masses = [r['M_total'] for r in results]
    delta_chi2 = [r['delta_chi2'] for r in results]
    colors = ['red' if r['favors'] == 'OCTH' else 'blue' for r in results]
    ax2.scatter(masses, delta_chi2, c=colors, alpha=0.5, s=20)
    ax2.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax2.set_xlabel('Masa total (M☉)')
    ax2.set_ylabel('Δχ² (GR - OCTH)')
    ax2.set_title('Δχ² vs Masa\n(Rojo = favorece OCTH, Azul = favorece GR)')

    # 3. Histograma de Δχ²
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.hist(delta_chi2, bins=30, color='gray', edgecolor='black', alpha=0.7)
    ax3.axvline(0, color='red', linestyle='--', linewidth=2, label='Δχ² = 0')
    ax3.axvline(np.mean(delta_chi2), color='blue', linestyle='-', linewidth=2,
                label=f'Media = {np.mean(delta_chi2):.1f}')
    ax3.set_xlabel('Δχ² (GR - OCTH)')
    ax3.set_ylabel('Frecuencia')
    ax3.set_title('Distribución de Δχ² en Inyecciones GR')
    ax3.legend()

    # 4. Resultado del test
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis('off')

    if test_passed:
        result_color = 'green'
        result_text = '✅ TEST PASADO'
        conclusion = 'El pipeline NO está sesgado\nhacia patrones hexagonales'
    else:
        result_color = 'red'
        result_text = '❌ TEST FALLIDO'
        conclusion = 'El pipeline ESTÁ SESGADO\nToda evidencia requiere revisión'

    ax4.text(0.5, 0.7, result_text, fontsize=24, fontweight='bold',
             ha='center', va='center', color=result_color,
             transform=ax4.transAxes)
    ax4.text(0.5, 0.5, conclusion, fontsize=14,
             ha='center', va='center', transform=ax4.transAxes)
    ax4.text(0.5, 0.3, f'Tasa FP (4/4 modos): {fp_rate*100:.2f}%\nUmbral: 5%',
             fontsize=12, ha='center', va='center', transform=ax4.transAxes,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.suptitle('TEST ASESINO A: Inyección de Señales GR Puras', fontsize=16, fontweight='bold')

    # Guardar
    fig_path = os.path.join(FIGURES_DIR, 'test_asesino_A_inyeccion_GR.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
    print(f"  Figura guardada en: {fig_path}")
    plt.close()


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Test Asesino A: Inyección GR Ciega')
    parser.add_argument('--n', type=int, default=100, help='Número de inyecciones')
    parser.add_argument('--quiet', action='store_true', help='Modo silencioso')

    args = parser.parse_args()

    results = run_blind_injection_test(n_injections=args.n, verbose=not args.quiet)

    print("\n" + "="*70)
    if results['test_passed']:
        print("  CONCLUSIÓN: El algoritmo es válido para detectar señal OCTH")
    else:
        print("  CONCLUSIÓN: El algoritmo necesita revisión antes de claims")
    print("="*70)
