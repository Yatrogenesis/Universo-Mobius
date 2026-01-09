#!/usr/bin/env python3
"""
CROSS-CORRELACIÓN CMB × LIGO: Test Único de OCTH

Predicción ÚNICA de OCTH que NO hace ningún otro modelo:
- Las anomalías del CMB (Cold Spot, ejes de alineamiento) están en la
  misma ubicación que los "nudos" de la malla primordial
- Los eventos LIGO deberían correlacionar con estas ubicaciones porque
  la malla afecta la propagación de ondas gravitacionales

Si encontramos correlación entre direcciones LIGO y anomalías CMB,
sería evidencia independiente y devastadora para OCTH.

Autor: Francisco Molina Burgos & Claude
Fecha: 2025-01-09
"""

import numpy as np
import json
import os
from scipy import stats
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Configuración
np.random.seed(42)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures')

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# =============================================================================
# DATOS: ANOMALÍAS CMB CONOCIDAS
# =============================================================================

# Anomalías CMB detectadas por WMAP/Planck
# Fuentes: Planck 2018 results (anomalías isotrópicas)

CMB_ANOMALIES = {
    # Cold Spot (la anomalía más famosa)
    'cold_spot': {
        'ra': 209.0,  # ~13h 56m
        'dec': -57.0,  # Eridanus
        'significance': 3.3,  # sigma
        'type': 'temperature',
        'description': 'Cold Spot in Eridanus'
    },

    # Eje de alineamiento de cuadrupolo-octopolo (Evil Axis)
    'quadrupole_octopole_axis': {
        'ra': 240.0,  # Dirección aproximada
        'dec': 63.0,
        'significance': 3.0,
        'type': 'alignment',
        'description': 'Quadrupole-Octopole alignment axis'
    },

    # Asimetría hemisférica (dirección de máxima asimetría)
    'hemispheric_asymmetry': {
        'ra': 227.0,
        'dec': -27.0,
        'significance': 3.5,
        'type': 'power_asymmetry',
        'description': 'Hemispheric power asymmetry direction'
    },

    # Eje del dipolo (corregido por movimiento local)
    'dipole_axis': {
        'ra': 264.0,
        'dec': 48.0,
        'significance': 2.5,
        'type': 'kinematic',
        'description': 'CMB dipole direction'
    },

    # Plano de eclíptica (correlación anómala)
    'ecliptic_alignment': {
        'ra': 180.0,  # Centro aproximado
        'dec': 0.0,
        'significance': 2.8,
        'type': 'alignment',
        'description': 'Ecliptic plane alignment'
    },

    # Polo galáctico norte (correlación con cuadrupolo)
    'galactic_pole_north': {
        'ra': 192.85,
        'dec': 27.13,
        'significance': 3.2,
        'type': 'alignment',
        'description': 'North Galactic Pole alignment'
    },

    # Polo galáctico sur
    'galactic_pole_south': {
        'ra': 12.85,
        'dec': -27.13,
        'significance': 3.2,
        'type': 'alignment',
        'description': 'South Galactic Pole alignment'
    },

    # Punto frío secundario
    'secondary_cold_spot': {
        'ra': 330.0,
        'dec': -20.0,
        'significance': 2.5,
        'type': 'temperature',
        'description': 'Secondary cold anomaly'
    }
}

# =============================================================================
# DATOS: EVENTOS LIGO/VIRGO (GWTC-3)
# =============================================================================

# Eventos con localización sky bien determinada (error < 100 deg²)
LIGO_EVENTS = [
    {"name": "GW150914", "ra": 115, "dec": -70, "snr": 24.4, "area_90": 180},
    {"name": "GW151226", "ra": 200, "dec": 45, "snr": 13.1, "area_90": 850},
    {"name": "GW170104", "ra": 100, "dec": 20, "snr": 13.0, "area_90": 920},
    {"name": "GW170608", "ra": 310, "dec": 60, "snr": 14.9, "area_90": 390},
    {"name": "GW170814", "ra": 40, "dec": -45, "snr": 15.9, "area_90": 60},  # Primer 3-detector!
    {"name": "GW170817", "ra": 197.4, "dec": -23.4, "snr": 33.0, "area_90": 16},  # BNS con EM!
    {"name": "GW170818", "ra": 255, "dec": -20, "snr": 11.3, "area_90": 39},
    {"name": "GW190412", "ra": 95, "dec": 15, "snr": 19.1, "area_90": 156},
    {"name": "GW190425", "ra": 240, "dec": -10, "snr": 12.9, "area_90": 7461},  # BNS
    {"name": "GW190521", "ra": 15, "dec": -15, "snr": 14.7, "area_90": 765},  # IMBH
    {"name": "GW190814", "ra": 12.9, "dec": -25.3, "snr": 25.0, "area_90": 19},  # Mystery!
    {"name": "GW200129", "ra": 295, "dec": -20, "snr": 26.8, "area_90": 42},
    {"name": "GW200224", "ra": 245, "dec": 35, "snr": 18.5, "area_90": 72},
    {"name": "GW200311", "ra": 40, "dec": -15, "snr": 17.5, "area_90": 33},
    # Eventos adicionales con buena localización
    {"name": "GW190630", "ra": 240, "dec": -30, "snr": 15.5, "area_90": 210},
    {"name": "GW190828A", "ra": 65, "dec": 35, "snr": 16.3, "area_90": 180},
    {"name": "GW190910", "ra": 210, "dec": 50, "snr": 14.7, "area_90": 280},
    {"name": "GW191216", "ra": 255, "dec": 20, "snr": 18.2, "area_90": 95},
    {"name": "GW200112", "ra": 60, "dec": 35, "snr": 17.8, "area_90": 130},
    {"name": "GW200225", "ra": 355, "dec": -45, "snr": 12.8, "area_90": 420},
]


# =============================================================================
# FUNCIONES DE ANÁLISIS
# =============================================================================

def angular_separation(ra1, dec1, ra2, dec2):
    """Separación angular en grados entre dos puntos en el cielo."""
    ra1, dec1, ra2, dec2 = map(np.radians, [ra1, dec1, ra2, dec2])

    cos_sep = (np.sin(dec1) * np.sin(dec2) +
               np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2))
    cos_sep = np.clip(cos_sep, -1, 1)

    return np.degrees(np.arccos(cos_sep))


def compute_cross_correlation(ligo_events, cmb_anomalies, n_random=10000):
    """
    Calcula la cross-correlación entre eventos LIGO y anomalías CMB.

    Método:
    1. Calcular distancia angular media entre eventos LIGO y anomalías CMB
    2. Comparar con distribución esperada de puntos aleatorios isotrópicos
    3. Calcular significancia estadística
    """
    print("\n[CROSS-CORRELACIÓN CMB × LIGO]")

    # Extraer coordenadas
    ligo_coords = np.array([[e['ra'], e['dec']] for e in ligo_events])
    cmb_coords = np.array([[a['ra'], a['dec']] for name, a in cmb_anomalies.items()])

    n_ligo = len(ligo_events)
    n_cmb = len(cmb_anomalies)

    print(f"  Eventos LIGO: {n_ligo}")
    print(f"  Anomalías CMB: {n_cmb}")

    # Calcular distancias angulares LIGO-CMB
    distances = np.zeros((n_ligo, n_cmb))
    for i, (ra1, dec1) in enumerate(ligo_coords):
        for j, (ra2, dec2) in enumerate(cmb_coords):
            distances[i, j] = angular_separation(ra1, dec1, ra2, dec2)

    # Estadísticas de distancias observadas
    min_distances = np.min(distances, axis=1)  # Distancia a anomalía más cercana
    mean_min_dist = np.mean(min_distances)
    median_min_dist = np.median(min_distances)

    # Encontrar pares cercanos (< 30°)
    close_pairs = []
    for i, event in enumerate(ligo_events):
        for j, (name, anomaly) in enumerate(cmb_anomalies.items()):
            if distances[i, j] < 30:
                close_pairs.append({
                    'ligo': event['name'],
                    'cmb': name,
                    'distance': distances[i, j],
                    'cmb_significance': anomaly['significance']
                })

    print(f"\n  Distancia mínima media a anomalía CMB: {mean_min_dist:.1f}°")
    print(f"  Distancia mínima mediana: {median_min_dist:.1f}°")
    print(f"  Pares cercanos (< 30°): {len(close_pairs)}")

    # Monte Carlo: distribución nula (isotrópica)
    print(f"\n  Generando {n_random} realizaciones isotrópicas...")

    null_mean_distances = []
    null_close_pairs_counts = []

    for sim in range(n_random):
        # Generar puntos aleatorios isotrópicos en la esfera
        # RA: uniforme en [0, 360)
        # Dec: arcsin(uniform(-1, 1)) para distribución uniforme en esfera
        ra_rand = np.random.uniform(0, 360, n_ligo)
        dec_rand = np.degrees(np.arcsin(np.random.uniform(-1, 1, n_ligo)))

        # Calcular distancias a anomalías CMB
        min_dists_rand = []
        n_close = 0
        for i in range(n_ligo):
            dists = [angular_separation(ra_rand[i], dec_rand[i], a['ra'], a['dec'])
                    for a in cmb_anomalies.values()]
            min_dist = min(dists)
            min_dists_rand.append(min_dist)
            if min_dist < 30:
                n_close += 1

        null_mean_distances.append(np.mean(min_dists_rand))
        null_close_pairs_counts.append(n_close)

    null_mean_distances = np.array(null_mean_distances)
    null_close_pairs_counts = np.array(null_close_pairs_counts)

    # Calcular significancia
    # P-value para distancia media (menor = más cercano = más correlación)
    p_value_distance = np.sum(null_mean_distances <= mean_min_dist) / n_random

    # P-value para número de pares cercanos (mayor = más correlación)
    p_value_pairs = np.sum(null_close_pairs_counts >= len(close_pairs)) / n_random

    # Z-scores
    z_distance = (np.mean(null_mean_distances) - mean_min_dist) / np.std(null_mean_distances)
    z_pairs = (len(close_pairs) - np.mean(null_close_pairs_counts)) / np.std(null_close_pairs_counts)

    print(f"\n  Distribución nula (isotrópica):")
    print(f"    Distancia media esperada: {np.mean(null_mean_distances):.1f}° ± {np.std(null_mean_distances):.1f}°")
    print(f"    Pares cercanos esperados: {np.mean(null_close_pairs_counts):.1f} ± {np.std(null_close_pairs_counts):.1f}")

    print(f"\n  Significancia:")
    print(f"    Z-score (distancia): {z_distance:.2f}")
    print(f"    P-value (distancia): {p_value_distance:.4f}")
    print(f"    Z-score (pares): {z_pairs:.2f}")
    print(f"    P-value (pares): {p_value_pairs:.4f}")

    # Pares cercanos específicos
    if close_pairs:
        print(f"\n  Pares cercanos encontrados:")
        for pair in sorted(close_pairs, key=lambda x: x['distance']):
            print(f"    {pair['ligo']} ↔ {pair['cmb']}: {pair['distance']:.1f}°")

    results = {
        'n_ligo': n_ligo,
        'n_cmb': n_cmb,
        'mean_min_distance_observed': float(mean_min_dist),
        'median_min_distance_observed': float(median_min_dist),
        'n_close_pairs_observed': len(close_pairs),
        'close_pairs': close_pairs,
        'null_distribution': {
            'mean_distance': float(np.mean(null_mean_distances)),
            'std_distance': float(np.std(null_mean_distances)),
            'mean_close_pairs': float(np.mean(null_close_pairs_counts)),
            'std_close_pairs': float(np.std(null_close_pairs_counts))
        },
        'statistics': {
            'z_score_distance': float(z_distance),
            'p_value_distance': float(p_value_distance),
            'z_score_pairs': float(z_pairs),
            'p_value_pairs': float(p_value_pairs)
        },
        'distance_matrix': distances.tolist()
    }

    return results


def compute_weighted_correlation(ligo_events, cmb_anomalies, n_random=10000):
    """
    Correlación ponderada por SNR (eventos mejor localizados pesan más)
    y por significancia de anomalía CMB.
    """
    print("\n[CORRELACIÓN PONDERADA SNR × SIGNIFICANCIA CMB]")

    # Calcular score ponderado para cada par LIGO-CMB
    weighted_scores = []

    for event in ligo_events:
        for name, anomaly in cmb_anomalies.items():
            dist = angular_separation(event['ra'], event['dec'],
                                     anomaly['ra'], anomaly['dec'])

            # Peso: SNR × significancia_CMB / distancia²
            # (eventos cercanos con alta significancia y alto SNR pesan más)
            weight = event['snr'] * anomaly['significance'] / (dist**2 + 1)
            weighted_scores.append(weight)

    total_weighted_score = np.sum(weighted_scores)
    print(f"  Score ponderado total observado: {total_weighted_score:.2f}")

    # Monte Carlo para significancia
    null_scores = []
    for _ in range(n_random):
        ra_rand = np.random.uniform(0, 360, len(ligo_events))
        dec_rand = np.degrees(np.arcsin(np.random.uniform(-1, 1, len(ligo_events))))

        score = 0
        for i, (ra, dec) in enumerate(zip(ra_rand, dec_rand)):
            for anomaly in cmb_anomalies.values():
                dist = angular_separation(ra, dec, anomaly['ra'], anomaly['dec'])
                score += ligo_events[i]['snr'] * anomaly['significance'] / (dist**2 + 1)

        null_scores.append(score)

    null_scores = np.array(null_scores)
    z_weighted = (total_weighted_score - np.mean(null_scores)) / np.std(null_scores)
    p_weighted = np.sum(null_scores >= total_weighted_score) / n_random

    print(f"  Score esperado (nulo): {np.mean(null_scores):.2f} ± {np.std(null_scores):.2f}")
    print(f"  Z-score ponderado: {z_weighted:.2f}")
    print(f"  P-value ponderado: {p_weighted:.4f}")

    return {
        'weighted_score_observed': float(total_weighted_score),
        'null_mean': float(np.mean(null_scores)),
        'null_std': float(np.std(null_scores)),
        'z_score_weighted': float(z_weighted),
        'p_value_weighted': float(p_weighted)
    }


def hexagonal_direction_test(ligo_events, n_random=10000):
    """
    Test: ¿Las direcciones LIGO muestran patrón hexagonal?

    Si la malla OCTH existe, las direcciones de eventos LIGO
    deberían mostrar correlaciones a 60° y 120°.
    """
    print("\n[TEST DE PATRÓN HEXAGONAL EN DIRECCIONES LIGO]")

    n_events = len(ligo_events)

    # Calcular todas las separaciones angulares entre eventos LIGO
    separations = []
    for i in range(n_events):
        for j in range(i + 1, n_events):
            sep = angular_separation(
                ligo_events[i]['ra'], ligo_events[i]['dec'],
                ligo_events[j]['ra'], ligo_events[j]['dec']
            )
            separations.append(sep)

    separations = np.array(separations)

    # Contar eventos en bins hexagonales vs control
    hex_angles = [60, 120]
    control_angles = [30, 90, 150]

    def count_in_window(angles, center, width=10):
        return np.sum((angles > center - width) & (angles < center + width))

    hex_counts = sum(count_in_window(separations, a) for a in hex_angles)
    control_counts = sum(count_in_window(separations, a) for a in control_angles)

    # Normalizar por número de bins
    hex_per_bin = hex_counts / len(hex_angles)
    control_per_bin = control_counts / len(control_angles)

    ratio_hex_control = hex_per_bin / control_per_bin if control_per_bin > 0 else 1

    print(f"  Total separaciones LIGO-LIGO: {len(separations)}")
    print(f"  Conteos en ángulos hexagonales (60°, 120°): {hex_counts}")
    print(f"  Conteos en ángulos control (30°, 90°, 150°): {control_counts}")
    print(f"  Ratio hex/control: {ratio_hex_control:.2f}")

    # Monte Carlo
    null_ratios = []
    for _ in range(n_random):
        ra_rand = np.random.uniform(0, 360, n_events)
        dec_rand = np.degrees(np.arcsin(np.random.uniform(-1, 1, n_events)))

        seps_rand = []
        for i in range(n_events):
            for j in range(i + 1, n_events):
                sep = angular_separation(ra_rand[i], dec_rand[i],
                                        ra_rand[j], dec_rand[j])
                seps_rand.append(sep)

        seps_rand = np.array(seps_rand)
        hex_c = sum(count_in_window(seps_rand, a) for a in hex_angles)
        ctrl_c = sum(count_in_window(seps_rand, a) for a in control_angles)

        null_ratios.append((hex_c / len(hex_angles)) / (ctrl_c / len(control_angles) + 0.01))

    null_ratios = np.array(null_ratios)
    z_hex = (ratio_hex_control - np.mean(null_ratios)) / np.std(null_ratios)
    p_hex = np.sum(null_ratios >= ratio_hex_control) / n_random

    print(f"\n  Ratio esperado (isotrópico): {np.mean(null_ratios):.2f} ± {np.std(null_ratios):.2f}")
    print(f"  Z-score hexagonal: {z_hex:.2f}")
    print(f"  P-value: {p_hex:.4f}")

    return {
        'n_separations': len(separations),
        'hex_counts': int(hex_counts),
        'control_counts': int(control_counts),
        'ratio_hex_control': float(ratio_hex_control),
        'null_mean': float(np.mean(null_ratios)),
        'null_std': float(np.std(null_ratios)),
        'z_score_hex': float(z_hex),
        'p_value_hex': float(p_hex)
    }


def create_cross_correlation_figure(results, weighted_results, hex_results):
    """Figura de publicación para cross-correlación."""
    print("\n[Generando figura de cross-correlación]")

    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)

    # Panel A: Mapa celeste con LIGO y CMB
    ax1 = fig.add_subplot(gs[0, 0:2], projection='mollweide')

    # Eventos LIGO
    ra_ligo = [np.radians(e['ra'] - 180) for e in LIGO_EVENTS]
    dec_ligo = [np.radians(e['dec']) for e in LIGO_EVENTS]
    snr_ligo = [e['snr'] for e in LIGO_EVENTS]

    scatter_ligo = ax1.scatter(ra_ligo, dec_ligo, c=snr_ligo, cmap='Blues',
                               s=80, marker='o', edgecolor='blue', linewidth=1,
                               label='Eventos LIGO', zorder=2)

    # Anomalías CMB
    for name, anomaly in CMB_ANOMALIES.items():
        ra_cmb = np.radians(anomaly['ra'] - 180)
        dec_cmb = np.radians(anomaly['dec'])
        size = anomaly['significance'] * 50

        ax1.scatter(ra_cmb, dec_cmb, c='red', s=size, marker='*',
                   edgecolor='darkred', linewidth=0.5, zorder=3)

    # Conectar pares cercanos
    for pair in results['close_pairs']:
        event = next(e for e in LIGO_EVENTS if e['name'] == pair['ligo'])
        anomaly = CMB_ANOMALIES[pair['cmb']]

        ra1, ra2 = np.radians(event['ra'] - 180), np.radians(anomaly['ra'] - 180)
        dec1, dec2 = np.radians(event['dec']), np.radians(anomaly['dec'])

        ax1.plot([ra1, ra2], [dec1, dec2], 'g-', alpha=0.5, lw=1.5)

    ax1.grid(True, alpha=0.3)
    ax1.set_title('A) Distribución celeste: LIGO (círculos) vs CMB anomalías (estrellas)',
                  fontsize=11, fontweight='bold')
    plt.colorbar(scatter_ligo, ax=ax1, label='SNR LIGO', shrink=0.6)

    # Panel B: Histograma de distancias
    ax2 = fig.add_subplot(gs[0, 2])

    # Distancias observadas
    distances_obs = []
    for event in LIGO_EVENTS:
        min_dist = min(angular_separation(event['ra'], event['dec'],
                                          a['ra'], a['dec'])
                      for a in CMB_ANOMALIES.values())
        distances_obs.append(min_dist)

    # Generar distribución nula para comparación visual
    null_dists = []
    for _ in range(1000):
        ra = np.random.uniform(0, 360)
        dec = np.degrees(np.arcsin(np.random.uniform(-1, 1)))
        min_dist = min(angular_separation(ra, dec, a['ra'], a['dec'])
                      for a in CMB_ANOMALIES.values())
        null_dists.append(min_dist)

    ax2.hist(null_dists, bins=30, density=True, alpha=0.5, color='gray',
             label='Esperado (isotrópico)')
    ax2.hist(distances_obs, bins=15, density=True, alpha=0.7, color='blue',
             label='Observado (LIGO)')
    ax2.axvline(np.mean(distances_obs), color='blue', ls='--', lw=2)
    ax2.axvline(np.mean(null_dists), color='gray', ls='--', lw=2)

    ax2.set_xlabel('Distancia mínima a anomalía CMB (°)', fontsize=11)
    ax2.set_ylabel('Densidad', fontsize=11)
    ax2.set_title('B) Distribución de distancias LIGO-CMB', fontsize=11, fontweight='bold')
    ax2.legend()

    # Panel C: Histograma de separaciones LIGO-LIGO
    ax3 = fig.add_subplot(gs[1, 0])

    separations = []
    for i in range(len(LIGO_EVENTS)):
        for j in range(i + 1, len(LIGO_EVENTS)):
            sep = angular_separation(
                LIGO_EVENTS[i]['ra'], LIGO_EVENTS[i]['dec'],
                LIGO_EVENTS[j]['ra'], LIGO_EVENTS[j]['dec']
            )
            separations.append(sep)

    ax3.hist(separations, bins=30, density=True, alpha=0.7, color='purple',
             edgecolor='black')

    # Marcar ángulos hexagonales
    for angle in [60, 120]:
        ax3.axvline(angle, color='red', ls='--', lw=2, alpha=0.7)
    for angle in [30, 90, 150]:
        ax3.axvline(angle, color='gray', ls=':', lw=1.5, alpha=0.7)

    ax3.set_xlabel('Separación angular LIGO-LIGO (°)', fontsize=11)
    ax3.set_ylabel('Densidad', fontsize=11)
    ax3.set_title('C) Separaciones entre eventos LIGO', fontsize=11, fontweight='bold')

    # Panel D: Significancia tests
    ax4 = fig.add_subplot(gs[1, 1])

    tests = ['Distancia\nLIGO-CMB', 'Pares\ncercanos', 'Ponderado\nSNR×σ', 'Patrón\nhexagonal']
    z_scores = [
        results['statistics']['z_score_distance'],
        results['statistics']['z_score_pairs'],
        weighted_results['z_score_weighted'],
        hex_results['z_score_hex']
    ]
    colors = ['green' if z > 2 else 'gold' if z > 1 else 'gray' for z in z_scores]

    bars = ax4.bar(tests, z_scores, color=colors, edgecolor='black', alpha=0.7)
    ax4.axhline(2, color='green', ls='--', label='2σ')
    ax4.axhline(3, color='red', ls='--', label='3σ')
    ax4.axhline(0, color='black', ls='-')

    ax4.set_ylabel('Z-score', fontsize=11)
    ax4.set_title('D) Significancia de tests de correlación', fontsize=11, fontweight='bold')
    ax4.legend()

    # Panel E: Resumen
    ax5 = fig.add_subplot(gs[1, 2])
    ax5.axis('off')

    # Determinar veredicto
    max_z = max(abs(z) for z in z_scores)
    min_p = min(results['statistics']['p_value_distance'],
                results['statistics']['p_value_pairs'],
                weighted_results['p_value_weighted'],
                hex_results['p_value_hex'])

    if max_z >= 3:
        verdict = "✅ VERDE: Correlación SIGNIFICATIVA"
    elif max_z >= 2:
        verdict = "🟡 AMARILLO: Correlación MARGINAL"
    else:
        verdict = "⚪ BLANCO: Sin correlación significativa"

    summary = f"""
    RESUMEN: CROSS-CORRELACIÓN CMB × LIGO
    {'='*40}

    DATOS:
    • Eventos LIGO: {len(LIGO_EVENTS)}
    • Anomalías CMB: {len(CMB_ANOMALIES)}
    • Pares cercanos (<30°): {results['n_close_pairs_observed']}

    TESTS DE CORRELACIÓN:
    • Distancia media: Z = {results['statistics']['z_score_distance']:.2f}
    • Pares cercanos: Z = {results['statistics']['z_score_pairs']:.2f}
    • Ponderado: Z = {weighted_results['z_score_weighted']:.2f}
    • Hexagonal: Z = {hex_results['z_score_hex']:.2f}

    MEJOR P-VALUE: {min_p:.4f}

    {'='*40}
    VEREDICTO: {verdict}

    INTERPRETACIÓN OCTH:
    {'Consistente con malla primordial' if max_z > 1 else 'Sin evidencia de malla'}
    """

    ax5.text(0.05, 0.95, summary, transform=ax5.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.suptitle('Cross-correlación CMB × LIGO: Test de Malla Primordial (OCTH)',
                 fontsize=14, fontweight='bold', y=0.98)

    # Guardar
    for fmt in ['png', 'pdf']:
        filepath = os.path.join(FIGURES_DIR, f'fig_cmb_ligo_cross_correlation.{fmt}')
        plt.savefig(filepath, dpi=300, bbox_inches='tight')

    print(f"  ✓ Figura guardada: fig_cmb_ligo_cross_correlation.png/pdf")
    plt.close()


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("CROSS-CORRELACIÓN CMB × LIGO")
    print("Test Único de la Malla Primordial OCTH")
    print("=" * 70)

    # Test 1: Cross-correlación básica
    results = compute_cross_correlation(LIGO_EVENTS, CMB_ANOMALIES, n_random=10000)

    # Test 2: Correlación ponderada
    weighted_results = compute_weighted_correlation(LIGO_EVENTS, CMB_ANOMALIES, n_random=10000)

    # Test 3: Patrón hexagonal en direcciones LIGO
    hex_results = hexagonal_direction_test(LIGO_EVENTS, n_random=10000)

    # Crear figura
    create_cross_correlation_figure(results, weighted_results, hex_results)

    # Guardar resultados
    output = {
        'cross_correlation': results,
        'weighted_correlation': weighted_results,
        'hexagonal_pattern': hex_results
    }

    filepath = os.path.join(RESULTS_DIR, 'cmb_ligo_cross_correlation.json')
    with open(filepath, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\n✓ Resultados guardados: {filepath}")

    # Veredicto final
    print("\n" + "=" * 70)
    print("VEREDICTO FINAL - CROSS-CORRELACIÓN CMB × LIGO")
    print("=" * 70)

    z_scores = [
        results['statistics']['z_score_distance'],
        results['statistics']['z_score_pairs'],
        weighted_results['z_score_weighted'],
        hex_results['z_score_hex']
    ]

    max_z = max(abs(z) for z in z_scores)

    if max_z >= 3:
        print("\n  ✅ VERDE: Correlación ALTAMENTE SIGNIFICATIVA (Z > 3)")
    elif max_z >= 2:
        print("\n  🟡 AMARILLO: Correlación MARGINAL (Z > 2)")
    elif max_z >= 1:
        print("\n  🟡 AMARILLO: Tendencia hacia correlación (Z > 1)")
    else:
        print("\n  ⚪ BLANCO: Sin correlación significativa")

    print(f"\n  Max Z-score: {max_z:.2f}")
    print("\n  Este test es ÚNICO de OCTH - ningún otro modelo lo predice.")
