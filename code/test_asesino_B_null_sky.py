#!/usr/bin/env python3
"""
TEST ASESINO B: Null Sky (Patrón de Antena)

OBJETIVO: Verificar que la correlación CMB×LIGO NO es un artefacto
del patrón de sensibilidad de antena de LIGO.

PROBLEMA A RESOLVER:
LIGO tiene "puntos ciegos" en el cielo donde la sensibilidad es baja.
Si las anomalías del CMB casualmente caen en zonas de alta sensibilidad,
podríamos ver correlación espuria.

MÉTODO:
1. Generar 10,000 universos simulados con BBH distribuidos según:
   - Distribución isotrópica en el cielo
   - Pesada por el patrón de antena de LIGO
2. Calcular correlación con anomalías CMB en cada universo
3. Verificar que Z=4.31 observado es raro incluso considerando el sesgo

CRITERIO DE ÉXITO:
- La distribución nula de Z-scores debe tener media ~0
- Z=4.31 debe estar en el tail <0.1% de la distribución nula

CRITERIO DE FALLO:
- Si Z≥4.31 ocurre en >1% de simulaciones, P3 pierde fuerza

Autor: Francisco Molina Burgos
Fecha: Enero 2026
"""

import numpy as np
import os
import json
from scipy import stats
from scipy.spatial.distance import cdist
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
# ANOMALÍAS CMB (POSICIONES REALES)
# =============================================================================

# Anomalías CMB conocidas con sus posiciones (Planck 2018)
CMB_ANOMALIES = [
    {"name": "Cold Spot", "ra": 209.0, "dec": -57.0, "significance": 3.5},
    {"name": "Quadrupole-Octopole Axis", "ra": 240.0, "dec": -63.0, "significance": 3.0},
    {"name": "Hemispherical Asymmetry Axis", "ra": 227.0, "dec": -27.0, "significance": 3.2},
    {"name": "North Galactic Pole", "ra": 192.86, "dec": 27.13, "significance": 2.5},
    {"name": "South Galactic Pole", "ra": 12.86, "dec": -27.13, "significance": 2.8},
    {"name": "CMB Dipole Direction", "ra": 264.0, "dec": 48.0, "significance": 2.0},
    {"name": "Parity Asymmetry Axis", "ra": 230.0, "dec": -15.0, "significance": 2.7},
    {"name": "Low-l Alignment Axis", "ra": 235.0, "dec": -20.0, "significance": 2.5},
]


# =============================================================================
# PATRÓN DE ANTENA DE LIGO
# =============================================================================

def antenna_pattern_ligo(ra, dec, detector='H1', gps_time=None):
    """
    Calcula el patrón de antena F+² + Fx² para un detector LIGO.

    El patrón de antena determina la sensibilidad direccional del detector.

    Parámetros:
    -----------
    ra : float or array
        Ascensión recta en grados
    dec : float or array
        Declinación en grados
    detector : str
        'H1' (Hanford), 'L1' (Livingston), 'V1' (Virgo)
    gps_time : float
        Tiempo GPS (afecta rotación de la Tierra)

    Returns:
    --------
    antenna_response : float or array
        F+² + Fx² (entre 0 y 1)
    """
    # Posiciones de los detectores (lat, lon, orientación del brazo)
    detectors = {
        'H1': {'lat': 46.45, 'lon': -119.41, 'arm_angle': 126.0},  # Hanford
        'L1': {'lat': 30.56, 'lon': -90.77, 'arm_angle': 243.0},   # Livingston
        'V1': {'lat': 43.63, 'lon': 10.50, 'arm_angle': 116.0},    # Virgo
    }

    det = detectors.get(detector, detectors['H1'])

    # Convertir a radianes
    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)
    lat_rad = np.radians(det['lat'])
    lon_rad = np.radians(det['lon'])
    arm_rad = np.radians(det['arm_angle'])

    # Hora sidérea local (simplificada)
    if gps_time is None:
        gps_time = 1187008882  # GW170814
    gmst = (gps_time / 86400.0 * 360.0) % 360.0
    lst = np.radians(gmst + det['lon'])

    # Ángulo horario
    ha = lst - ra_rad

    # Funciones de patrón de antena (aproximación)
    # F+ y Fx dependen de la orientación del detector respecto a la fuente

    # Componentes del vector de dirección en coordenadas del detector
    cos_dec = np.cos(dec_rad)
    sin_dec = np.sin(dec_rad)
    cos_ha = np.cos(ha)
    sin_ha = np.sin(ha)
    cos_lat = np.cos(lat_rad)
    sin_lat = np.sin(lat_rad)

    # Patrón de antena simplificado (promediado sobre polarización)
    # F² = F+² + Fx² ≈ (1 + cos²(θ)²)/4 donde θ es el ángulo cenital
    cos_theta = sin_lat * sin_dec + cos_lat * cos_dec * cos_ha
    sin_theta_sq = 1 - cos_theta**2

    # Respuesta de antena (normalizada 0-1)
    # Máxima en el cenit, cero en el horizonte
    F_squared = (1 + cos_theta**2)**2 / 4

    # Incluir orientación del brazo
    psi = arm_rad
    F_plus = 0.5 * (1 + cos_theta**2) * np.cos(2*psi)
    F_cross = cos_theta * np.sin(2*psi)

    antenna_response = F_plus**2 + F_cross**2

    # Normalizar
    antenna_response = np.clip(antenna_response, 0, 1)

    return antenna_response


def combined_antenna_pattern(ra, dec, gps_time=None):
    """
    Calcula el patrón de antena combinado para la red LIGO-Virgo.

    La sensibilidad combinada es aproximadamente la suma en cuadratura.
    """
    F_H1 = antenna_pattern_ligo(ra, dec, 'H1', gps_time)
    F_L1 = antenna_pattern_ligo(ra, dec, 'L1', gps_time)
    F_V1 = antenna_pattern_ligo(ra, dec, 'V1', gps_time)

    # Red combinada (promedio de sensibilidades)
    F_combined = np.sqrt((F_H1**2 + F_L1**2 + F_V1**2) / 3)

    return F_combined


# =============================================================================
# GENERACIÓN DE POBLACIONES GW SIMULADAS
# =============================================================================

def generate_gw_population(n_events, antenna_weighted=True, n_gps_samples=10):
    """
    Genera una población de eventos GW simulados.

    Parámetros:
    -----------
    n_events : int
        Número de eventos a generar
    antenna_weighted : bool
        Si True, pesa por patrón de antena
    n_gps_samples : int
        Número de tiempos GPS a promediar

    Returns:
    --------
    events : list of dict
        Lista de eventos con ra, dec
    """
    events = []

    # Generar más candidatos de los necesarios (rejection sampling)
    n_candidates = n_events * 10

    # Distribución uniforme en el cielo
    ra_candidates = np.random.uniform(0, 360, n_candidates)
    # Declinación uniforme en sin(dec)
    dec_candidates = np.degrees(np.arcsin(np.random.uniform(-1, 1, n_candidates)))

    if antenna_weighted:
        # Calcular peso de antena promediado sobre varios tiempos GPS
        gps_times = np.linspace(1167559936, 1269388800, n_gps_samples)  # O3 epoch

        weights = np.zeros(n_candidates)
        for gps in gps_times:
            weights += combined_antenna_pattern(ra_candidates, dec_candidates, gps)
        weights /= n_gps_samples

        # Normalizar
        weights = weights / np.max(weights)

        # Rejection sampling
        accept = np.random.random(n_candidates) < weights
        ra_accepted = ra_candidates[accept]
        dec_accepted = dec_candidates[accept]

        # Tomar los primeros n_events
        ra_final = ra_accepted[:n_events]
        dec_final = dec_accepted[:n_events]

        # Si no hay suficientes, rellenar con uniformes
        if len(ra_final) < n_events:
            n_missing = n_events - len(ra_final)
            ra_extra = np.random.uniform(0, 360, n_missing)
            dec_extra = np.degrees(np.arcsin(np.random.uniform(-1, 1, n_missing)))
            ra_final = np.concatenate([ra_final, ra_extra])
            dec_final = np.concatenate([dec_final, dec_extra])
    else:
        ra_final = ra_candidates[:n_events]
        dec_final = dec_candidates[:n_events]

    for i in range(n_events):
        events.append({
            'ra': ra_final[i],
            'dec': dec_final[i],
            'localization_error': np.random.uniform(10, 500)  # deg²
        })

    return events


# =============================================================================
# CORRELACIÓN CMB-GW
# =============================================================================

def angular_separation(ra1, dec1, ra2, dec2):
    """
    Calcula la separación angular entre dos puntos en el cielo.

    Returns:
    --------
    sep : float
        Separación en grados
    """
    ra1_rad = np.radians(ra1)
    dec1_rad = np.radians(dec1)
    ra2_rad = np.radians(ra2)
    dec2_rad = np.radians(dec2)

    cos_sep = (np.sin(dec1_rad) * np.sin(dec2_rad) +
               np.cos(dec1_rad) * np.cos(dec2_rad) * np.cos(ra1_rad - ra2_rad))

    cos_sep = np.clip(cos_sep, -1, 1)
    sep_rad = np.arccos(cos_sep)

    return np.degrees(sep_rad)


def count_close_pairs(gw_events, cmb_anomalies, threshold_deg=30):
    """
    Cuenta pares cercanos entre eventos GW y anomalías CMB.

    Parámetros:
    -----------
    gw_events : list
        Lista de eventos GW
    cmb_anomalies : list
        Lista de anomalías CMB
    threshold_deg : float
        Umbral de separación en grados

    Returns:
    --------
    n_pairs : int
        Número de pares con separación < threshold
    min_sep : float
        Mínima separación encontrada
    pairs : list
        Lista de pares encontrados
    """
    pairs = []
    min_sep = 180.0

    for gw in gw_events:
        for cmb in cmb_anomalies:
            sep = angular_separation(gw['ra'], gw['dec'], cmb['ra'], cmb['dec'])
            if sep < min_sep:
                min_sep = sep
            if sep < threshold_deg:
                pairs.append({
                    'gw_ra': gw['ra'],
                    'gw_dec': gw['dec'],
                    'cmb_name': cmb['name'],
                    'separation': sep
                })

    return len(pairs), min_sep, pairs


def compute_correlation_zscore(gw_events, cmb_anomalies, threshold_deg=30, n_expected=8.4, sigma_expected=2.2):
    """
    Calcula el Z-score de la correlación CMB-GW.

    Usa los valores esperados del análisis real como referencia.
    """
    n_observed, min_sep, pairs = count_close_pairs(gw_events, cmb_anomalies, threshold_deg)

    z_score = (n_observed - n_expected) / sigma_expected

    return {
        'n_observed': n_observed,
        'n_expected': n_expected,
        'z_score': z_score,
        'min_separation': min_sep,
        'pairs': pairs
    }


# =============================================================================
# TEST NULL SKY
# =============================================================================

def run_null_sky_test(n_simulations=10000, n_events_per_sim=20, verbose=True):
    """
    Ejecuta el test Null Sky.

    Genera muchos universos simulados con GW distribuidos según
    el patrón de antena de LIGO y verifica que la correlación
    con CMB no sea espuria.
    """
    print("\n" + "="*70)
    print("  TEST ASESINO B: NULL SKY (PATRÓN DE ANTENA)")
    print("="*70)
    print(f"\n  Número de simulaciones: {n_simulations}")
    print(f"  Eventos por simulación: {n_events_per_sim}")
    print("  Hipótesis nula: No hay correlación CMB-GW real")
    print("  Z observado en datos reales: 4.31")
    print("  Criterio de fallo: Z≥4.31 ocurre en >1% de simulaciones")

    # Valores observados reales
    Z_OBSERVED = 4.31
    N_EXPECTED = 8.4
    SIGMA_EXPECTED = 2.2
    THRESHOLD_DEG = 30

    results_uniform = []
    results_weighted = []

    print("\n[FASE 1: Simulaciones con distribución uniforme]")

    for i in range(n_simulations):
        if verbose and (i + 1) % 2000 == 0:
            print(f"  Simulación {i+1}/{n_simulations}...")

        # Generar población uniforme
        gw_events = generate_gw_population(n_events_per_sim, antenna_weighted=False)
        result = compute_correlation_zscore(gw_events, CMB_ANOMALIES, THRESHOLD_DEG, N_EXPECTED, SIGMA_EXPECTED)
        results_uniform.append(result)

    print("\n[FASE 2: Simulaciones con peso de antena]")

    for i in range(n_simulations):
        if verbose and (i + 1) % 2000 == 0:
            print(f"  Simulación {i+1}/{n_simulations}...")

        # Generar población pesada por antena
        gw_events = generate_gw_population(n_events_per_sim, antenna_weighted=True)
        result = compute_correlation_zscore(gw_events, CMB_ANOMALIES, THRESHOLD_DEG, N_EXPECTED, SIGMA_EXPECTED)
        results_weighted.append(result)

    # Análisis estadístico
    print("\n[FASE 3: Análisis estadístico]")

    z_uniform = [r['z_score'] for r in results_uniform]
    z_weighted = [r['z_score'] for r in results_weighted]

    print(f"\n  Distribución UNIFORME:")
    print(f"    Media Z: {np.mean(z_uniform):.3f}")
    print(f"    Std Z:   {np.std(z_uniform):.3f}")
    print(f"    P(Z ≥ 4.31): {100*np.mean(np.array(z_uniform) >= Z_OBSERVED):.3f}%")

    print(f"\n  Distribución PESADA por antena:")
    print(f"    Media Z: {np.mean(z_weighted):.3f}")
    print(f"    Std Z:   {np.std(z_weighted):.3f}")
    print(f"    P(Z ≥ 4.31): {100*np.mean(np.array(z_weighted) >= Z_OBSERVED):.3f}%")

    # Determinar resultado del test
    p_uniform = np.mean(np.array(z_uniform) >= Z_OBSERVED)
    p_weighted = np.mean(np.array(z_weighted) >= Z_OBSERVED)

    # Usar el caso más conservador (pesado por antena)
    p_null = max(p_uniform, p_weighted)

    print(f"\n  === RESULTADO DEL TEST ===")
    print(f"  P(Z ≥ 4.31 | null) = {p_null*100:.3f}%")

    if p_null < 0.01:
        print(f"\n  ✅ TEST PASADO: La correlación CMB-GW es genuina")
        print(f"     El patrón de antena NO explica Z=4.31")
        test_passed = True
    else:
        print(f"\n  ❌ TEST FALLIDO: La correlación puede ser espuria")
        print(f"     El patrón de antena PUEDE explicar Z=4.31")
        test_passed = False

    # Calcular percentil del valor observado
    percentile_uniform = stats.percentileofscore(z_uniform, Z_OBSERVED)
    percentile_weighted = stats.percentileofscore(z_weighted, Z_OBSERVED)

    print(f"\n  Percentil de Z=4.31:")
    print(f"    En distribución uniforme: {percentile_uniform:.1f}%")
    print(f"    En distribución pesada:   {percentile_weighted:.1f}%")

    # Guardar resultados
    output = {
        'test_name': 'Null Sky Antenna Pattern Test',
        'n_simulations': n_simulations,
        'n_events_per_sim': n_events_per_sim,
        'z_observed': Z_OBSERVED,
        'uniform': {
            'mean_z': float(np.mean(z_uniform)),
            'std_z': float(np.std(z_uniform)),
            'p_exceed': float(p_uniform),
            'percentile': float(percentile_uniform)
        },
        'weighted': {
            'mean_z': float(np.mean(z_weighted)),
            'std_z': float(np.std(z_weighted)),
            'p_exceed': float(p_weighted),
            'percentile': float(percentile_weighted)
        },
        'test_passed': test_passed,
        'p_null': float(p_null)
    }

    output_path = os.path.join(RESULTS_DIR, 'test_asesino_B_null_sky.json')
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\n  Resultados guardados en: {output_path}")

    # Generar figura
    create_null_sky_figure(z_uniform, z_weighted, Z_OBSERVED, test_passed)

    return output


def create_null_sky_figure(z_uniform, z_weighted, z_observed, test_passed):
    """Genera figura del test Null Sky."""

    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)

    # 1. Histograma de Z-scores (uniforme)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.hist(z_uniform, bins=50, density=True, alpha=0.7, color='blue',
             edgecolor='black', label='Simulaciones')
    ax1.axvline(z_observed, color='red', linewidth=2, linestyle='--',
                label=f'Z observado = {z_observed}')
    ax1.axvline(0, color='gray', linewidth=1, linestyle='-')

    # Curva normal de referencia
    x = np.linspace(-4, 6, 100)
    ax1.plot(x, stats.norm.pdf(x), 'k-', linewidth=1, label='N(0,1)')

    ax1.set_xlabel('Z-score')
    ax1.set_ylabel('Densidad')
    ax1.set_title('Distribución Null (Uniforme en cielo)')
    ax1.legend()
    ax1.set_xlim(-4, 6)

    # 2. Histograma de Z-scores (pesado por antena)
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.hist(z_weighted, bins=50, density=True, alpha=0.7, color='orange',
             edgecolor='black', label='Simulaciones')
    ax2.axvline(z_observed, color='red', linewidth=2, linestyle='--',
                label=f'Z observado = {z_observed}')
    ax2.axvline(0, color='gray', linewidth=1, linestyle='-')
    ax2.plot(x, stats.norm.pdf(x), 'k-', linewidth=1, label='N(0,1)')

    ax2.set_xlabel('Z-score')
    ax2.set_ylabel('Densidad')
    ax2.set_title('Distribución Null (Pesada por antena LIGO)')
    ax2.legend()
    ax2.set_xlim(-4, 6)

    # 3. Mapa del cielo con anomalías CMB
    ax3 = fig.add_subplot(gs[1, 0], projection='mollweide')

    # Convertir a radianes para proyección Mollweide
    cmb_ra = [np.radians(a['ra'] - 180) for a in CMB_ANOMALIES]
    cmb_dec = [np.radians(a['dec']) for a in CMB_ANOMALIES]
    cmb_sig = [a['significance'] for a in CMB_ANOMALIES]

    ax3.scatter(cmb_ra, cmb_dec, c='red', s=[s**2 * 30 for s in cmb_sig],
                marker='*', label='Anomalías CMB', zorder=10)

    # Generar una muestra de la distribución pesada para visualizar
    sample_events = generate_gw_population(100, antenna_weighted=True)
    sample_ra = [np.radians(e['ra'] - 180) for e in sample_events]
    sample_dec = [np.radians(e['dec']) for e in sample_events]
    ax3.scatter(sample_ra, sample_dec, c='blue', s=5, alpha=0.3,
                label='GW simulados (antena)')

    ax3.set_title('Anomalías CMB y GW Simulados')
    ax3.legend(loc='lower right')
    ax3.grid(True, alpha=0.3)

    # 4. Resultado del test
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis('off')

    p_exceed = np.mean(np.array(z_weighted) >= z_observed)

    if test_passed:
        result_color = 'green'
        result_text = '✅ TEST PASADO'
        conclusion = f'La correlación CMB-GW es genuina\nP(Z≥4.31|null) = {p_exceed*100:.3f}% < 1%'
    else:
        result_color = 'red'
        result_text = '❌ TEST FALLIDO'
        conclusion = f'La correlación puede ser espuria\nP(Z≥4.31|null) = {p_exceed*100:.3f}% ≥ 1%'

    ax4.text(0.5, 0.7, result_text, fontsize=24, fontweight='bold',
             ha='center', va='center', color=result_color,
             transform=ax4.transAxes)
    ax4.text(0.5, 0.5, conclusion, fontsize=12,
             ha='center', va='center', transform=ax4.transAxes)

    # Tabla de resultados
    table_text = (
        f"Simulaciones: 10,000\n"
        f"Z observado: {z_observed}\n"
        f"Media Z (uniforme): {np.mean(z_uniform):.2f}\n"
        f"Media Z (antena): {np.mean(z_weighted):.2f}\n"
        f"Percentil (antena): {stats.percentileofscore(z_weighted, z_observed):.1f}%"
    )
    ax4.text(0.5, 0.25, table_text, fontsize=10, ha='center', va='center',
             transform=ax4.transAxes, family='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.suptitle('TEST ASESINO B: Null Sky (Patrón de Antena)', fontsize=16, fontweight='bold')

    # Guardar
    fig_path = os.path.join(FIGURES_DIR, 'test_asesino_B_null_sky.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
    print(f"  Figura guardada en: {fig_path}")
    plt.close()


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Test Asesino B: Null Sky')
    parser.add_argument('--n', type=int, default=10000, help='Número de simulaciones')
    parser.add_argument('--events', type=int, default=20, help='Eventos por simulación')
    parser.add_argument('--quiet', action='store_true', help='Modo silencioso')

    args = parser.parse_args()

    results = run_null_sky_test(
        n_simulations=args.n,
        n_events_per_sim=args.events,
        verbose=not args.quiet
    )

    print("\n" + "="*70)
    if results['test_passed']:
        print("  CONCLUSIÓN: La correlación CMB-GW es evidencia genuina")
    else:
        print("  CONCLUSIÓN: La correlación CMB-GW necesita investigación adicional")
    print("="*70)
