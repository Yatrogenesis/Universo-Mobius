#!/usr/bin/env python3
"""
TEST EUROPEO: VIRGO vs LIGO
============================

Este test es DEFINITIVO para eliminar el argumento del ruido de 60Hz.

VIRGO (Italia): Red eléctrica 50Hz
LIGO (USA): Red eléctrica 60Hz

Si los modos hexagonales aparecen en VIRGO también → NO ES RUIDO
Si solo aparecen en LIGO → Podría ser artefacto de 60Hz

Evento: GW170814 - El PRIMER evento triple detector (H1+L1+V1)
Este es el evento perfecto porque tenemos datos de los 3 detectores.

Autor: Francisco Molina Burgos & Claude
Fecha: 2026-01-09
"""

import numpy as np
import h5py
from scipy import signal
from scipy.fft import fft, fftfreq
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
import json

# Configuración
np.random.seed(42)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, 'data', 'ligo')
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures')

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# Parámetros del evento GW170814
EVENT_GPS = 1186741861.5
M1 = 30.6  # M_sun
M2 = 25.2  # M_sun
M_TOTAL = M1 + M2  # 55.8 M_sun
M_FINAL = 53.2  # M_sun (después de radiación GW)
F_ISCO = 4400 / M_TOTAL  # ~79 Hz

# Archivos
FILES = {
    'V1': 'V-V1_GWOSC_4KHZ_R1-1186741846-32.hdf5',
    'H1': 'H-H1_GWOSC_4KHZ_R1-1186741846-32.hdf5',
    'L1': 'L-L1_GWOSC_4KHZ_R1-1186741846-32.hdf5'
}

# Info de red eléctrica
POWER_LINE = {
    'V1': 50,   # Italia - Europa 50Hz
    'H1': 60,   # USA - 60Hz
    'L1': 60    # USA - 60Hz
}

def load_strain_data(detector):
    """Carga datos de strain de un detector."""
    filepath = os.path.join(DATA_DIR, FILES[detector])

    if not os.path.exists(filepath):
        print(f"  ⚠ Archivo no encontrado: {filepath}")
        return None, None, None

    with h5py.File(filepath, 'r') as f:
        # Estructura típica GWOSC
        strain = f['strain/Strain'][:]

        # Metadatos
        meta = f['meta']
        gps_start = meta['GPSstart'][()]
        duration = meta['Duration'][()]

        # Calcular sample rate desde duración y número de muestras
        sample_rate = len(strain) // duration

    dt = 1.0 / sample_rate
    time = np.arange(len(strain)) * dt

    return strain, time, sample_rate


def bandpass_filter(data, lowcut, highcut, fs, order=4):
    """Filtro pasabanda Butterworth."""
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq

    # Asegurar que los valores están en rango válido
    low = max(0.001, min(low, 0.99))
    high = max(low + 0.01, min(high, 0.999))

    b, a = signal.butter(order, [low, high], btype='band')
    return signal.filtfilt(b, a, data)


def compute_spectrogram(strain, sample_rate, nperseg=256):
    """Calcula espectrograma."""
    f, t, Sxx = signal.spectrogram(strain, fs=sample_rate, nperseg=nperseg,
                                   noverlap=nperseg//2, scaling='density')
    return f, t, Sxx


def find_peaks_in_spectrum(psd, freqs, f_min=20, f_max=500, n_peaks=10):
    """Encuentra picos en el espectro de potencia."""
    # Filtrar rango de frecuencias
    mask = (freqs >= f_min) & (freqs <= f_max)
    psd_masked = psd[mask]
    freqs_masked = freqs[mask]

    # Encontrar picos
    peaks_idx, properties = signal.find_peaks(psd_masked,
                                               height=np.median(psd_masked),
                                               distance=5)

    # Ordenar por altura
    if len(peaks_idx) > 0:
        heights = properties['peak_heights']
        sorted_idx = np.argsort(heights)[::-1][:n_peaks]
        peak_freqs = freqs_masked[peaks_idx[sorted_idx]]
        peak_heights = heights[sorted_idx]
        return peak_freqs, peak_heights

    return np.array([]), np.array([])


def check_hexagonal_modes(peak_freqs, f_isco, tolerance=0.15):
    """
    Verifica si los picos coinciden con modos hexagonales.

    Modos predichos OCTH:
    - f1 = f_ISCO / 2 (fundamental)
    - f2 = f1 * sqrt(3) ≈ f1 * 1.732
    - f3 = f1 * 2
    - f4 = f1 * sqrt(7) ≈ f1 * 2.646
    """
    f1 = f_isco / 2
    hex_ratios = [1.0, np.sqrt(3), 2.0, np.sqrt(7)]
    predicted_modes = [f1 * r for r in hex_ratios]

    matches = []
    for i, f_pred in enumerate(predicted_modes):
        for f_obs in peak_freqs:
            rel_diff = abs(f_obs - f_pred) / f_pred
            if rel_diff < tolerance:
                matches.append({
                    'mode': i + 1,
                    'predicted': f_pred,
                    'observed': f_obs,
                    'ratio': hex_ratios[i],
                    'rel_diff': rel_diff
                })
                break

    return matches, predicted_modes


def analyze_detector(detector, verbose=True):
    """Analiza los datos de un detector."""
    if verbose:
        print(f"\n{'='*60}")
        print(f"  ANALIZANDO {detector} (Red: {POWER_LINE[detector]}Hz)")
        print(f"{'='*60}")

    # Cargar datos
    strain, time, sample_rate = load_strain_data(detector)
    if strain is None:
        return None

    if verbose:
        print(f"  Datos cargados: {len(strain)} muestras @ {sample_rate}Hz")
        print(f"  Duración: {len(strain)/sample_rate:.1f}s")

    # Filtrar datos (20-500 Hz para evitar ruido de baja frecuencia)
    strain_filtered = bandpass_filter(strain, 20, 500, sample_rate)

    # Calcular PSD
    nperseg = min(4096, len(strain_filtered)//4)
    freqs, psd = signal.welch(strain_filtered, fs=sample_rate, nperseg=nperseg)

    # Encontrar picos
    peak_freqs, peak_heights = find_peaks_in_spectrum(psd, freqs, f_min=30, f_max=300)

    if verbose and len(peak_freqs) > 0:
        print(f"\n  Top 5 picos de frecuencia:")
        for i, (f, h) in enumerate(zip(peak_freqs[:5], peak_heights[:5])):
            print(f"    {i+1}. f = {f:.1f} Hz (altura: {h:.2e})")

    # Verificar modos hexagonales
    matches, predicted = check_hexagonal_modes(peak_freqs, F_ISCO)

    if verbose:
        print(f"\n  Frecuencias predichas OCTH (f_ISCO = {F_ISCO:.1f} Hz):")
        for i, (ratio, f) in enumerate(zip([1, np.sqrt(3), 2, np.sqrt(7)], predicted)):
            label = ['f₁', 'f₂=√3·f₁', 'f₃=2·f₁', 'f₄=√7·f₁'][i]
            print(f"    {label} = {f:.1f} Hz")

        print(f"\n  Coincidencias encontradas: {len(matches)}/4")
        for m in matches:
            print(f"    Modo {m['mode']}: predicho={m['predicted']:.1f}Hz, "
                  f"observado={m['observed']:.1f}Hz (diff={m['rel_diff']*100:.1f}%)")

    # Verificar ruido de línea eléctrica
    power_line_freq = POWER_LINE[detector]
    harmonics = [power_line_freq * n for n in range(1, 6)]

    power_line_contamination = []
    for harm in harmonics:
        # Buscar picos cerca de armónicos de la línea
        for f in peak_freqs:
            if abs(f - harm) < 2:  # Dentro de 2Hz
                power_line_contamination.append({
                    'harmonic': harm,
                    'peak': f
                })

    if verbose:
        print(f"\n  Contaminación de línea eléctrica ({power_line_freq}Hz):")
        if power_line_contamination:
            for c in power_line_contamination:
                print(f"    ⚠ Pico en {c['peak']:.1f}Hz cerca de armónico {c['harmonic']}Hz")
        else:
            print(f"    ✓ No se detecta contaminación significativa")

    # Calcular espectrograma para visualización
    f_spec, t_spec, Sxx = compute_spectrogram(strain_filtered, sample_rate, nperseg=256)

    return {
        'detector': detector,
        'power_line': POWER_LINE[detector],
        'sample_rate': sample_rate,
        'strain': strain_filtered,
        'time': time,
        'freqs': freqs,
        'psd': psd,
        'peak_freqs': peak_freqs.tolist() if len(peak_freqs) > 0 else [],
        'peak_heights': peak_heights.tolist() if len(peak_heights) > 0 else [],
        'hex_matches': matches,
        'predicted_modes': predicted,
        'power_line_contamination': power_line_contamination,
        'spectrogram': {'f': f_spec, 't': t_spec, 'Sxx': Sxx}
    }


def create_comparison_figure(results):
    """Crea figura comparativa de los 3 detectores."""
    print("\n[Generando figura comparativa]")

    fig = plt.figure(figsize=(16, 14))
    gs = GridSpec(4, 3, figure=fig, hspace=0.35, wspace=0.3)

    detectors = ['V1', 'H1', 'L1']
    colors = {'V1': 'green', 'H1': 'blue', 'L1': 'orange'}

    # Fila 1: Espectrogramas
    for i, det in enumerate(detectors):
        if det not in results or results[det] is None:
            continue

        ax = fig.add_subplot(gs[0, i])
        r = results[det]
        spec = r['spectrogram']

        # Limitar a frecuencias de interés
        f_mask = spec['f'] < 300

        im = ax.pcolormesh(spec['t'], spec['f'][f_mask],
                          10*np.log10(spec['Sxx'][f_mask] + 1e-50),
                          cmap='viridis', shading='auto')

        # Marcar frecuencias predichas
        for f_pred in r['predicted_modes']:
            if f_pred < 300:
                ax.axhline(f_pred, color='red', ls='--', alpha=0.7, lw=1)

        ax.set_ylabel('Frecuencia (Hz)')
        ax.set_xlabel('Tiempo (s)')
        ax.set_title(f'{det} ({POWER_LINE[det]}Hz) - Espectrograma', fontweight='bold')
        ax.set_ylim(20, 200)

    # Fila 2: PSD
    for i, det in enumerate(detectors):
        if det not in results or results[det] is None:
            continue

        ax = fig.add_subplot(gs[1, i])
        r = results[det]

        ax.semilogy(r['freqs'], r['psd'], color=colors[det], lw=0.8, alpha=0.8)

        # Marcar picos encontrados
        for f in r['peak_freqs'][:5]:
            ax.axvline(f, color='gray', ls=':', alpha=0.5)

        # Marcar modos predichos
        for j, f_pred in enumerate(r['predicted_modes']):
            color = 'red' if any(m['mode']==j+1 for m in r['hex_matches']) else 'pink'
            ax.axvline(f_pred, color=color, ls='--', lw=2, alpha=0.8)

        # Marcar armónicos de línea
        for harm in range(1, 5):
            ax.axvline(POWER_LINE[det] * harm, color='black', ls=':', alpha=0.3)

        ax.set_xlim(20, 300)
        ax.set_xlabel('Frecuencia (Hz)')
        ax.set_ylabel('PSD')
        ax.set_title(f'{det} - Densidad Espectral de Potencia')
        ax.grid(True, alpha=0.3)

    # Fila 3: Comparación directa de los 3 detectores
    ax = fig.add_subplot(gs[2, :])

    for det in detectors:
        if det not in results or results[det] is None:
            continue
        r = results[det]

        # Normalizar PSD para comparación
        psd_norm = r['psd'] / np.max(r['psd'])
        label = f"{det} ({POWER_LINE[det]}Hz)"
        ax.semilogy(r['freqs'], psd_norm, color=colors[det], lw=1.2, alpha=0.8, label=label)

    # Marcar modos hexagonales predichos
    f1 = F_ISCO / 2
    hex_modes = [f1, f1*np.sqrt(3), f1*2, f1*np.sqrt(7)]
    for j, f in enumerate(hex_modes):
        label = ['f₁', '√3·f₁', '2·f₁', '√7·f₁'][j]
        ax.axvline(f, color='red', ls='--', lw=2, alpha=0.7,
                  label=f'{label}={f:.0f}Hz' if j < 2 else '')

    ax.set_xlim(20, 250)
    ax.set_xlabel('Frecuencia (Hz)', fontsize=12)
    ax.set_ylabel('PSD normalizado', fontsize=12)
    ax.set_title('COMPARACIÓN: VIRGO (50Hz) vs LIGO (60Hz) - ¿Mismos modos?',
                fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # Fila 4: Resumen
    ax = fig.add_subplot(gs[3, :])
    ax.axis('off')

    # Contar coincidencias por detector
    summary_data = []
    for det in detectors:
        if det not in results or results[det] is None:
            continue
        r = results[det]
        n_matches = len(r['hex_matches'])
        contaminated = len(r['power_line_contamination']) > 0
        summary_data.append({
            'det': det,
            'power': POWER_LINE[det],
            'matches': n_matches,
            'contaminated': contaminated
        })

    # Verificar si VIRGO tiene los mismos modos que LIGO
    v1_modes = set(m['mode'] for m in results.get('V1', {}).get('hex_matches', []))
    h1_modes = set(m['mode'] for m in results.get('H1', {}).get('hex_matches', []))
    l1_modes = set(m['mode'] for m in results.get('L1', {}).get('hex_matches', []))

    common_modes = v1_modes & h1_modes & l1_modes if v1_modes else set()

    # Determinar veredicto
    if len(v1_modes) >= 2 and v1_modes == h1_modes:
        verdict = "✅ VERDE: VIRGO muestra los MISMOS modos que LIGO\n   → NO ES RUIDO DE 60Hz"
        verdict_color = 'green'
    elif len(v1_modes) >= 1:
        verdict = "🟡 AMARILLO: VIRGO detecta algunos modos hexagonales\n   → Evidencia parcial"
        verdict_color = 'gold'
    else:
        verdict = "⚪ BLANCO: No hay suficientes coincidencias\n   → Datos insuficientes"
        verdict_color = 'gray'

    summary_text = f"""
    {'='*70}
    TEST EUROPEO: VIRGO (Italia, 50Hz) vs LIGO (USA, 60Hz)
    {'='*70}

    EVENTO: GW170814 (Primer evento triple detector)
    MASA TOTAL: {M_TOTAL:.1f} M☉  |  f_ISCO: {F_ISCO:.1f} Hz

    RESULTADOS POR DETECTOR:
    """

    for s in summary_data:
        status = "⚠ contaminado" if s['contaminated'] else "✓ limpio"
        summary_text += f"\n    {s['det']} ({s['power']}Hz): {s['matches']}/4 modos coincidentes [{status}]"

    summary_text += f"""

    MODOS COMUNES (V1 ∩ H1 ∩ L1): {common_modes if common_modes else 'Ninguno'}

    {'='*70}
    VEREDICTO: {verdict}
    {'='*70}

    INTERPRETACIÓN:
    Si los modos aparecen en VIRGO (red 50Hz) y en LIGO (red 60Hz),
    NO pueden ser artefactos de la red eléctrica.
    → La señal hexagonal es de origen ASTROFÍSICO.
    """

    ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, fontsize=11,
           verticalalignment='top', fontfamily='monospace',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.suptitle('TEST EUROPEO: Eliminando el argumento del ruido de 60Hz',
                fontsize=16, fontweight='bold', y=0.98)

    # Guardar
    for fmt in ['png', 'pdf']:
        filepath = os.path.join(FIGURES_DIR, f'fig_european_test_virgo.{fmt}')
        plt.savefig(filepath, dpi=300, bbox_inches='tight')

    print(f"  ✓ Figura guardada: fig_european_test_virgo.png/pdf")
    plt.close()

    return {
        'v1_modes': list(v1_modes),
        'h1_modes': list(h1_modes),
        'l1_modes': list(l1_modes),
        'common_modes': list(common_modes),
        'verdict': verdict
    }


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("TEST EUROPEO: VIRGO vs LIGO")
    print("¿Los modos hexagonales son ruido de 60Hz o señal astrofísica?")
    print("=" * 70)
    print(f"\nEvento: GW170814 (Primer evento triple detector)")
    print(f"Masas: {M1} + {M2} = {M_TOTAL} M☉")
    print(f"f_ISCO predicha: {F_ISCO:.1f} Hz")

    # Analizar cada detector
    results = {}
    for detector in ['V1', 'H1', 'L1']:
        results[detector] = analyze_detector(detector, verbose=True)

    # Crear figura comparativa
    comparison = create_comparison_figure(results)

    # Guardar resultados
    output = {
        'event': 'GW170814',
        'masses': {'m1': M1, 'm2': M2, 'total': M_TOTAL, 'final': M_FINAL},
        'f_isco': F_ISCO,
        'detectors': {}
    }

    for det, r in results.items():
        if r is not None:
            output['detectors'][det] = {
                'power_line_hz': POWER_LINE[det],
                'peak_freqs': r['peak_freqs'][:10],
                'hex_matches': r['hex_matches'],
                'predicted_modes': r['predicted_modes'],
                'n_matches': len(r['hex_matches']),
                'power_line_contamination': r['power_line_contamination']
            }

    output['comparison'] = comparison

    filepath = os.path.join(RESULTS_DIR, 'test_european_virgo.json')
    with open(filepath, 'w') as f:
        json.dump(output, f, indent=2, default=float)
    print(f"\n✓ Resultados guardados: {filepath}")

    # Veredicto final
    print("\n" + "=" * 70)
    print("VEREDICTO FINAL - TEST EUROPEO")
    print("=" * 70)
    print(f"\n{comparison['verdict']}")
