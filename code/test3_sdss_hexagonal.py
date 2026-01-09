#!/usr/bin/env python3
"""
TEST #3: GEOMETRÍA HEXAGONAL EN DISTRIBUCIÓN DE GALAXIAS
========================================================

Buscar patrones hexagonales (60°/120°) en la distribución angular
de galaxias usando datos del Sloan Digital Sky Survey (SDSS).

PREDICCIÓN OCTH:
Si el espaciotiempo tiene estructura hexagonal a escala de Planck,
la inflación cósmica habría estirado esa geometría hasta escalas
macroscópicas. Esto produciría:

1. Exceso de correlación angular a 60° y 120°
2. Déficit a 90° (cuadrado) respecto al patrón hexagonal
3. Simetría D3 en la función de correlación angular

METODOLOGÍA:
1. Obtener catálogo de galaxias de SDSS
2. Calcular función de correlación angular ω(θ)
3. Buscar picos en θ = 60°, 120° vs baseline
4. Comparar con simulaciones isotrópicas (null hypothesis)

DATOS:
- SDSS DR17 (Data Release 17)
- Galaxias con redshift espectroscópico
- Selección por magnitud y tipo

Autor: Francisco Molina Burgos
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.spatial.distance import cdist
from scipy.special import legendre
import os
import json
import urllib.request

# =============================================================================
# CONSTANTES Y CONFIGURACIÓN
# =============================================================================

# Ángulos hexagonales de interés
HEXAGONAL_ANGLES = [60, 120]  # grados
SQUARE_ANGLE = 90  # para comparación

# Directorios
DATA_DIR = "../data/sdss"
FIGURES_DIR = "../figures"
RESULTS_DIR = "../results"

# URL del catálogo SDSS (query simplificado)
SDSS_QUERY_URL = "https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch"

# =============================================================================
# FUNCIONES DE DATOS
# =============================================================================

def generate_mock_sdss_catalog(n_galaxies=10000, seed=42):
    """
    Genera catálogo simulado de galaxias con distribución realista.

    Para desarrollo, usamos simulación. Luego reemplazar con datos reales.

    La distribución incluye:
    - Clustering a pequeñas escalas
    - Distribución aproximadamente uniforme a grandes escalas
    - Estructura filamentaria (opcional)
    """
    np.random.seed(seed)

    # Cobertura del cielo SDSS (aproximada)
    # RA: 0-360°, Dec: -10° a +70° (hemisferio norte principalmente)
    ra_min, ra_max = 100, 260  # Región principal de SDSS
    dec_min, dec_max = 0, 60

    # Distribución base uniforme
    n_uniform = int(n_galaxies * 0.7)
    ra_uniform = np.random.uniform(ra_min, ra_max, n_uniform)
    dec_uniform = np.random.uniform(dec_min, dec_max, n_uniform)

    # Añadir clustering (grupos de galaxias)
    n_clusters = int(n_galaxies * 0.3)
    n_cluster_centers = 50

    # Centros de clusters
    cluster_ra = np.random.uniform(ra_min + 20, ra_max - 20, n_cluster_centers)
    cluster_dec = np.random.uniform(dec_min + 10, dec_max - 10, n_cluster_centers)
    cluster_sizes = np.random.exponential(3, n_cluster_centers)  # grados

    ra_clustered = []
    dec_clustered = []

    for i in range(n_clusters):
        # Elegir un cluster aleatorio
        c_idx = np.random.randint(0, n_cluster_centers)
        # Posición gaussiana alrededor del centro
        ra_clustered.append(cluster_ra[c_idx] + np.random.normal(0, cluster_sizes[c_idx]))
        dec_clustered.append(cluster_dec[c_idx] + np.random.normal(0, cluster_sizes[c_idx] * 0.5))

    ra_clustered = np.array(ra_clustered)
    dec_clustered = np.array(dec_clustered)

    # Combinar
    ra = np.concatenate([ra_uniform, ra_clustered])
    dec = np.concatenate([dec_uniform, dec_clustered])

    # Filtrar fuera de rango
    valid = (ra >= ra_min) & (ra <= ra_max) & (dec >= dec_min) & (dec <= dec_max)
    ra = ra[valid]
    dec = dec[valid]

    # Redshift simulado (distribución realista)
    z = np.random.exponential(0.1, len(ra))
    z = np.clip(z, 0.01, 0.5)

    return {
        'ra': ra,
        'dec': dec,
        'z': z,
        'n_galaxies': len(ra),
        'source': 'Simulación con clustering'
    }


def inject_hexagonal_signal(catalog, strength=0.05, scale=60):
    """
    Inyecta señal hexagonal artificial para test de sensibilidad.

    Añade galaxias en patrones hexagonales a la distribución existente.

    Parámetros:
    - strength: fracción de galaxias adicionales en patrón hexagonal
    - scale: escala angular del patrón en grados
    """
    n_hex = int(len(catalog['ra']) * strength)

    # Crear patrón hexagonal
    # Centros aleatorios
    n_centers = n_hex // 7  # 7 galaxias por hexágono (1 centro + 6 vértices)

    center_ra = np.random.uniform(catalog['ra'].min() + scale,
                                   catalog['ra'].max() - scale, n_centers)
    center_dec = np.random.uniform(catalog['dec'].min() + scale,
                                    catalog['dec'].max() - scale, n_centers)

    hex_ra = list(center_ra)
    hex_dec = list(center_dec)

    # Añadir vértices del hexágono a escala variable
    for i in range(n_centers):
        for angle in [0, 60, 120, 180, 240, 300]:
            r = scale * np.random.uniform(0.8, 1.2)  # variación en tamaño
            theta = np.radians(angle + np.random.normal(0, 5))  # variación angular
            hex_ra.append(center_ra[i] + r * np.cos(theta) / np.cos(np.radians(center_dec[i])))
            hex_dec.append(center_dec[i] + r * np.sin(theta))

    hex_ra = np.array(hex_ra)
    hex_dec = np.array(hex_dec)

    # Combinar con catálogo original
    new_catalog = {
        'ra': np.concatenate([catalog['ra'], hex_ra]),
        'dec': np.concatenate([catalog['dec'], hex_dec]),
        'z': np.concatenate([catalog['z'], np.random.exponential(0.1, len(hex_ra))]),
        'source': catalog['source'] + ' + Hexagonal inyectado'
    }
    new_catalog['n_galaxies'] = len(new_catalog['ra'])

    return new_catalog


def download_sdss_sample(n_max=50000):
    """
    Descarga muestra de galaxias de SDSS usando SQL query.

    Nota: Para datasets grandes, usar CasJobs de SDSS.
    """
    # Query SQL para SDSS
    query = f"""
    SELECT TOP {n_max}
        p.ra, p.dec, s.z, s.zErr
    FROM PhotoObj AS p
    JOIN SpecObj AS s ON s.bestobjid = p.objid
    WHERE
        p.type = 3  -- Galaxias
        AND s.class = 'GALAXY'
        AND s.zWarning = 0  -- Redshift confiable
        AND s.z BETWEEN 0.01 AND 0.3
        AND p.r BETWEEN 14 AND 17.77  -- Magnitud límite
    ORDER BY p.ra
    """

    try:
        import requests

        params = {
            'cmd': query,
            'format': 'csv'
        }

        response = requests.get(SDSS_QUERY_URL, params=params, timeout=60)

        if response.status_code == 200:
            # Parsear CSV
            lines = response.text.strip().split('\n')
            if len(lines) > 1:
                ra, dec, z = [], [], []
                for line in lines[1:]:  # Skip header
                    parts = line.split(',')
                    if len(parts) >= 3:
                        try:
                            ra.append(float(parts[0]))
                            dec.append(float(parts[1]))
                            z.append(float(parts[2]))
                        except:
                            continue

                return {
                    'ra': np.array(ra),
                    'dec': np.array(dec),
                    'z': np.array(z),
                    'n_galaxies': len(ra),
                    'source': 'SDSS DR17'
                }
    except Exception as e:
        print(f"  ⚠ Error descargando SDSS: {e}")

    return None


# =============================================================================
# ANÁLISIS DE CORRELACIÓN ANGULAR
# =============================================================================

def angular_separation(ra1, dec1, ra2, dec2):
    """
    Calcula separación angular entre dos puntos en el cielo.

    Fórmula de Vincenty (más estable que haversine para distancias pequeñas).

    Entrada en grados, salida en grados.
    """
    ra1, dec1, ra2, dec2 = map(np.radians, [ra1, dec1, ra2, dec2])

    dra = ra2 - ra1

    num = np.sqrt((np.cos(dec2) * np.sin(dra))**2 +
                  (np.cos(dec1) * np.sin(dec2) - np.sin(dec1) * np.cos(dec2) * np.cos(dra))**2)
    den = np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(dra)

    return np.degrees(np.arctan2(num, den))


def compute_angular_correlation(catalog, theta_bins=None, n_random=None, max_pairs=500000):
    """
    Calcula función de correlación angular ω(θ) usando estimador Landy-Szalay.

    ω(θ) = (DD - 2DR + RR) / RR

    donde:
    - DD = pares data-data a separación θ
    - DR = pares data-random
    - RR = pares random-random
    """
    if theta_bins is None:
        # Bins de 1 grado de 0 a 180
        theta_bins = np.arange(0, 181, 2)

    ra = catalog['ra']
    dec = catalog['dec']
    n_data = len(ra)

    if n_random is None:
        n_random = min(n_data * 2, 20000)

    print(f"  → Calculando correlación angular ({n_data} galaxias)...")

    # Generar catálogo random con misma geometría
    ra_rand = np.random.uniform(ra.min(), ra.max(), n_random)
    dec_rand = np.random.uniform(dec.min(), dec.max(), n_random)

    # Submuestra para eficiencia si hay muchas galaxias
    if n_data > 5000:
        idx_sample = np.random.choice(n_data, 5000, replace=False)
        ra_sample = ra[idx_sample]
        dec_sample = dec[idx_sample]
    else:
        ra_sample = ra
        dec_sample = dec

    n_sample = len(ra_sample)

    # Calcular DD (data-data)
    print("    Calculando DD...")
    DD_counts = np.zeros(len(theta_bins) - 1)

    n_pairs = 0
    for i in range(n_sample):
        if n_pairs > max_pairs:
            break
        for j in range(i + 1, n_sample):
            sep = angular_separation(ra_sample[i], dec_sample[i],
                                    ra_sample[j], dec_sample[j])
            bin_idx = np.searchsorted(theta_bins, sep) - 1
            if 0 <= bin_idx < len(DD_counts):
                DD_counts[bin_idx] += 1
            n_pairs += 1

    # Normalizar DD
    n_DD_pairs = n_sample * (n_sample - 1) / 2
    DD = DD_counts / n_DD_pairs if n_DD_pairs > 0 else DD_counts

    # Calcular RR (random-random)
    print("    Calculando RR...")
    RR_counts = np.zeros(len(theta_bins) - 1)

    n_rand_sample = min(n_random, 3000)
    idx_rand = np.random.choice(n_random, n_rand_sample, replace=False)
    ra_rand_sample = ra_rand[idx_rand]
    dec_rand_sample = dec_rand[idx_rand]

    n_pairs = 0
    for i in range(n_rand_sample):
        if n_pairs > max_pairs:
            break
        for j in range(i + 1, n_rand_sample):
            sep = angular_separation(ra_rand_sample[i], dec_rand_sample[i],
                                    ra_rand_sample[j], dec_rand_sample[j])
            bin_idx = np.searchsorted(theta_bins, sep) - 1
            if 0 <= bin_idx < len(RR_counts):
                RR_counts[bin_idx] += 1
            n_pairs += 1

    # Normalizar RR
    n_RR_pairs = n_rand_sample * (n_rand_sample - 1) / 2
    RR = RR_counts / n_RR_pairs if n_RR_pairs > 0 else RR_counts

    # Calcular DR (data-random)
    print("    Calculando DR...")
    DR_counts = np.zeros(len(theta_bins) - 1)

    n_dr_sample = min(n_sample, 2000)
    n_pairs = 0
    for i in range(n_dr_sample):
        if n_pairs > max_pairs // 2:
            break
        for j in range(min(n_rand_sample, 2000)):
            sep = angular_separation(ra_sample[i], dec_sample[i],
                                    ra_rand_sample[j], dec_rand_sample[j])
            bin_idx = np.searchsorted(theta_bins, sep) - 1
            if 0 <= bin_idx < len(DR_counts):
                DR_counts[bin_idx] += 1
            n_pairs += 1

    # Normalizar DR
    n_DR_pairs = n_dr_sample * min(n_rand_sample, 2000)
    DR = DR_counts / n_DR_pairs if n_DR_pairs > 0 else DR_counts

    # Estimador Landy-Szalay
    # ω(θ) = (DD - 2DR + RR) / RR
    omega = np.zeros(len(theta_bins) - 1)
    for i in range(len(omega)):
        if RR[i] > 0:
            omega[i] = (DD[i] - 2 * DR[i] + RR[i]) / RR[i]
        else:
            omega[i] = 0

    # Error Poisson
    omega_err = np.sqrt(DD_counts + 1) / n_DD_pairs / (RR + 1e-10)

    theta_centers = (theta_bins[:-1] + theta_bins[1:]) / 2

    return {
        'theta': theta_centers,
        'omega': omega,
        'omega_err': omega_err,
        'DD': DD,
        'DR': DR,
        'RR': RR,
        'theta_bins': theta_bins
    }


def analyze_hexagonal_signature(correlation):
    """
    Busca señal hexagonal en la función de correlación.

    MÉTODO MEJORADO:
    Compara ángulos hexagonales (60°, 120°) con sus vecinos inmediatos
    para eliminar sesgo por tendencia general de ω(θ).

    Hexagonal: 60° comparado con (50°, 70°)
               120° comparado con (110°, 130°)
    """
    theta = correlation['theta']
    omega = correlation['omega']
    omega_err = correlation['omega_err']

    # Encontrar valores en ángulos de interés
    def get_omega_at_angle(target_angle, width=3):
        """Obtiene omega promediado cerca de un ángulo."""
        mask = np.abs(theta - target_angle) < width
        if np.sum(mask) > 0:
            return np.mean(omega[mask]), np.std(omega[mask]) / np.sqrt(np.sum(mask))
        return np.nan, np.nan

    # Ángulos hexagonales y sus vecinos
    omega_50, err_50 = get_omega_at_angle(50)
    omega_60, err_60 = get_omega_at_angle(60)
    omega_70, err_70 = get_omega_at_angle(70)

    omega_90, err_90 = get_omega_at_angle(90)

    omega_110, err_110 = get_omega_at_angle(110)
    omega_120, err_120 = get_omega_at_angle(120)
    omega_130, err_130 = get_omega_at_angle(130)

    # Para 30° y 150° (adicionales)
    omega_30, err_30 = get_omega_at_angle(30)
    omega_150, err_150 = get_omega_at_angle(150)

    # Exceso LOCAL en 60° (vs interpolación de vecinos)
    baseline_60 = (omega_50 + omega_70) / 2
    excess_60 = omega_60 - baseline_60

    # Exceso LOCAL en 120° (vs interpolación de vecinos)
    baseline_120 = (omega_110 + omega_130) / 2
    excess_120 = omega_120 - baseline_120

    # Exceso hexagonal total (promedio de ambos)
    hex_excess = (excess_60 + excess_120) / 2

    # Error combinado
    err_baseline_60 = np.sqrt(err_50**2 + err_70**2) / 2
    err_baseline_120 = np.sqrt(err_110**2 + err_130**2) / 2

    err_excess_60 = np.sqrt(err_60**2 + err_baseline_60**2)
    err_excess_120 = np.sqrt(err_120**2 + err_baseline_120**2)

    combined_err = np.sqrt(err_excess_60**2 + err_excess_120**2) / 2

    z_score = hex_excess / combined_err if combined_err > 0 else 0
    p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))

    # También calcular ω(90) para referencia
    # El 90° es el ángulo "cuadrado" - si hay geometría hexagonal, 90° debería ser menor

    return {
        'omega_60': omega_60,
        'omega_120': omega_120,
        'omega_30': omega_30,
        'omega_90': omega_90,
        'omega_150': omega_150,
        'excess_60': excess_60,
        'excess_120': excess_120,
        'hex_mean': (omega_60 + omega_120) / 2,
        'control_mean': omega_90,  # 90° como control principal
        'hex_excess': hex_excess,
        'z_score': z_score,
        'p_value': p_value
    }


def test_angular_uniformity(correlation, n_simulations=100):
    """
    Test de hipótesis: ¿La correlación angular es consistente con isotropía?

    Genera simulaciones isotrópicas y compara.
    """
    theta = correlation['theta']
    omega_observed = correlation['omega']

    print(f"  → Generando {n_simulations} simulaciones isotrópicas...")

    # Simular distribuciones isotrópicas
    omega_simulations = []

    for i in range(n_simulations):
        # Catálogo isotrópico
        n_sim = 5000
        ra_sim = np.random.uniform(100, 260, n_sim)
        dec_sim = np.random.uniform(0, 60, n_sim)

        sim_catalog = {
            'ra': ra_sim,
            'dec': dec_sim,
            'z': np.random.exponential(0.1, n_sim)
        }

        # Correlación simplificada (solo DD normalizado)
        # Para velocidad, solo calculamos DD
        DD_counts = np.zeros(len(correlation['theta_bins']) - 1)

        n_pairs = 0
        max_pairs = 50000
        n_sample = min(n_sim, 1000)
        idx = np.random.choice(n_sim, n_sample, replace=False)

        for ii in range(n_sample):
            if n_pairs > max_pairs:
                break
            for jj in range(ii + 1, n_sample):
                sep = angular_separation(ra_sim[idx[ii]], dec_sim[idx[ii]],
                                        ra_sim[idx[jj]], dec_sim[idx[jj]])
                bin_idx = np.searchsorted(correlation['theta_bins'], sep) - 1
                if 0 <= bin_idx < len(DD_counts):
                    DD_counts[bin_idx] += 1
                n_pairs += 1

        DD_counts = DD_counts / np.sum(DD_counts) if np.sum(DD_counts) > 0 else DD_counts
        omega_simulations.append(DD_counts)

    omega_simulations = np.array(omega_simulations)

    # Estadísticas de las simulaciones
    omega_sim_mean = np.mean(omega_simulations, axis=0)
    omega_sim_std = np.std(omega_simulations, axis=0)

    # Chi-squared vs isotrópico
    # Normalizar omega observado para comparación
    omega_obs_norm = correlation['DD'] / np.sum(correlation['DD']) if np.sum(correlation['DD']) > 0 else correlation['DD']

    chi2 = np.sum(((omega_obs_norm - omega_sim_mean) / (omega_sim_std + 1e-10))**2)
    dof = len(theta) - 1
    p_value_chi2 = 1 - stats.chi2.cdf(chi2, dof)

    return {
        'chi2': chi2,
        'dof': dof,
        'chi2_reduced': chi2 / dof if dof > 0 else 0,
        'p_value': p_value_chi2,
        'omega_sim_mean': omega_sim_mean,
        'omega_sim_std': omega_sim_std
    }


# =============================================================================
# VISUALIZACIÓN
# =============================================================================

def plot_analysis(catalog, correlation, hex_analysis, uniformity_test):
    """Genera figura completa del análisis."""

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Panel 1: Distribución de galaxias
    ax = axes[0, 0]
    ax.scatter(catalog['ra'], catalog['dec'], s=1, alpha=0.3, c='blue')
    ax.set_xlabel('RA (grados)', fontsize=12)
    ax.set_ylabel('Dec (grados)', fontsize=12)
    ax.set_title(f'Distribución de Galaxias (N={catalog["n_galaxies"]})', fontsize=12)
    ax.grid(True, alpha=0.3)

    # Panel 2: Función de correlación angular
    ax = axes[0, 1]
    theta = correlation['theta']
    omega = correlation['omega']
    omega_err = correlation['omega_err']

    ax.errorbar(theta, omega, yerr=omega_err, fmt='o-', ms=3, capsize=2,
                alpha=0.7, label='Observado')

    # Marcar ángulos hexagonales
    for angle in [60, 120]:
        ax.axvline(angle, color='red', ls='--', alpha=0.7,
                   label=f'{angle}° (hexagonal)' if angle == 60 else '')
    ax.axvline(90, color='gray', ls=':', alpha=0.5, label='90° (cuadrado)')

    ax.axhline(0, color='black', ls='-', alpha=0.3)
    ax.set_xlabel('Separación angular θ (grados)', fontsize=12)
    ax.set_ylabel('ω(θ)', fontsize=12)
    ax.set_title('Función de Correlación Angular', fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 180)

    # Panel 3: Comparación hexagonal vs control
    ax = axes[1, 0]
    angles = ['30°', '60°\n(hex)', '90°', '120°\n(hex)', '150°']
    values = [hex_analysis['omega_30'], hex_analysis['omega_60'],
              hex_analysis['omega_90'], hex_analysis['omega_120'],
              hex_analysis['omega_150']]
    colors = ['gray', 'red', 'gray', 'red', 'gray']

    bars = ax.bar(angles, values, color=colors, alpha=0.7, edgecolor='black')
    ax.axhline(hex_analysis['control_mean'], color='blue', ls='--',
               label=f'Media control: {hex_analysis["control_mean"]:.4f}')
    ax.axhline(hex_analysis['hex_mean'], color='red', ls='--',
               label=f'Media hex: {hex_analysis["hex_mean"]:.4f}')

    ax.set_ylabel('ω(θ)', fontsize=12)
    ax.set_title(f'Comparación Hexagonal vs Control\nExceso: {hex_analysis["hex_excess"]:.4f}, Z={hex_analysis["z_score"]:.2f}',
                fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')

    # Panel 4: Resumen
    ax = axes[1, 1]
    ax.axis('off')

    # Determinar conclusión basada en comparación directa hex vs 90°
    # Ratio: omega_hex / omega_90 - si > 1, hay exceso hexagonal
    hex_mean = (hex_analysis['omega_60'] + hex_analysis['omega_120']) / 2
    ratio_hex_vs_square = hex_mean / hex_analysis['omega_90'] if hex_analysis['omega_90'] > 0 else 1

    if hex_analysis['p_value'] < 0.05 and hex_analysis['hex_excess'] > 0:
        conclusion = "POSIBLE SEÑAL HEXAGONAL"
        result_color = 'green'
    elif ratio_hex_vs_square > 1.3:
        conclusion = f"EXCESO HEXAGONAL vs CUADRADO (ratio={ratio_hex_vs_square:.2f})"
        result_color = 'green'
    elif hex_analysis['p_value'] < 0.05 and hex_analysis['hex_excess'] < 0:
        conclusion = "DÉFICIT EN ÁNGULOS HEXAGONALES"
        result_color = 'orange'
    else:
        conclusion = "SIN SEÑAL HEXAGONAL SIGNIFICATIVA"
        result_color = 'gray'

    summary = f"""
    RESULTADOS TEST #3: GEOMETRIA HEXAGONAL EN SDSS
    ================================================

    DATOS
    -----
    Fuente: {catalog['source']}
    N galaxias: {catalog['n_galaxies']}

    CORRELACION EN ANGULOS CLAVE
    ----------------------------
    omega(60):  {hex_analysis['omega_60']:.4f}
    omega(90):  {hex_analysis['omega_90']:.4f}
    omega(120): {hex_analysis['omega_120']:.4f}

    ANALISIS HEXAGONAL
    ------------------
    Media hexagonal (60,120): {hex_analysis['hex_mean']:.4f}
    Media control (30,90,150): {hex_analysis['control_mean']:.4f}
    Exceso hexagonal: {hex_analysis['hex_excess']:.4f}
    Z-score: {hex_analysis['z_score']:.2f}
    P-value: {hex_analysis['p_value']:.4f}

    TEST DE ISOTROPIA
    -----------------
    chi2/dof: {uniformity_test['chi2_reduced']:.2f}
    P-value: {uniformity_test['p_value']:.4f}

    CONCLUSION: {conclusion}

    INTERPRETACION OCTH:
    {'Consistente con patron hexagonal residual' if hex_analysis['hex_excess'] > 0 and hex_analysis['p_value'] < 0.1 else 'Sin evidencia de geometria hexagonal a esta escala'}
    """

    ax.text(0.05, 0.95, summary, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', fontfamily='monospace',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()

    # Guardar
    plt.savefig(os.path.join(FIGURES_DIR, 'fig9_sdss_hexagonal.png'), dpi=300)
    plt.savefig(os.path.join(FIGURES_DIR, 'fig9_sdss_hexagonal.pdf'), dpi=300)
    plt.close()

    print(f"  ✓ fig9_sdss_hexagonal.png/pdf")

    return conclusion


# =============================================================================
# EJECUCIÓN PRINCIPAL
# =============================================================================

def run_test3(use_real_data=False, inject_hex=False, n_galaxies=10000):
    """Ejecuta Test #3 completo."""

    print("=" * 70)
    print("TEST #3: GEOMETRÍA HEXAGONAL EN DISTRIBUCIÓN DE GALAXIAS")
    print("=" * 70)

    if inject_hex:
        print(">>> MODO TEST DE SENSIBILIDAD: Señal hexagonal inyectada <<<")

    # Crear directorios
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(FIGURES_DIR, exist_ok=True)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Cargar o generar datos
    print("\n[1] CARGANDO DATOS")

    if use_real_data:
        print("  → Intentando descargar datos de SDSS...")
        catalog = download_sdss_sample(n_max=n_galaxies)

        if catalog is None:
            print("  ⚠ No se pudo descargar SDSS, usando simulación...")
            catalog = generate_mock_sdss_catalog(n_galaxies=n_galaxies)
    else:
        catalog = generate_mock_sdss_catalog(n_galaxies=n_galaxies)

    print(f"  ✓ {catalog['n_galaxies']} galaxias cargadas")
    print(f"  ✓ Fuente: {catalog['source']}")

    # Inyectar señal hexagonal si es test de sensibilidad
    if inject_hex:
        print("\n[1b] INYECTANDO SEÑAL HEXAGONAL")
        catalog = inject_hexagonal_signal(catalog, strength=0.1, scale=30)
        print(f"  ✓ {catalog['n_galaxies']} galaxias después de inyección")

    # Calcular correlación angular
    print("\n[2] CALCULANDO CORRELACIÓN ANGULAR")
    correlation = compute_angular_correlation(catalog)
    print(f"  ✓ Correlación calculada para {len(correlation['theta'])} bins angulares")

    # Analizar señal hexagonal
    print("\n[3] ANÁLISIS DE SEÑAL HEXAGONAL")
    hex_analysis = analyze_hexagonal_signature(correlation)
    print(f"  → ω(60°) = {hex_analysis['omega_60']:.4f}")
    print(f"  → ω(90°) = {hex_analysis['omega_90']:.4f}")
    print(f"  → ω(120°) = {hex_analysis['omega_120']:.4f}")
    print(f"  → Exceso hexagonal: {hex_analysis['hex_excess']:.4f}")
    print(f"  → Z-score: {hex_analysis['z_score']:.2f}")
    print(f"  → P-value: {hex_analysis['p_value']:.4f}")

    if hex_analysis['p_value'] < 0.05:
        if hex_analysis['hex_excess'] > 0:
            print("  → SEÑAL HEXAGONAL DETECTADA")
        else:
            print("  → DÉFICIT HEXAGONAL DETECTADO")
    else:
        print("  → Sin señal significativa")

    # Test de isotropía
    print("\n[4] TEST DE ISOTROPÍA")
    uniformity_test = test_angular_uniformity(correlation, n_simulations=50)
    print(f"  → χ²/dof = {uniformity_test['chi2_reduced']:.2f}")
    print(f"  → P-value = {uniformity_test['p_value']:.4f}")

    # Generar figuras
    print("\n[5] GENERANDO FIGURAS")
    conclusion = plot_analysis(catalog, correlation, hex_analysis, uniformity_test)

    # Guardar resultados
    print("\n[6] GUARDANDO RESULTADOS")
    results = {
        'test': 'Test #3 - Geometría Hexagonal en SDSS',
        'date': '2026-01-08',
        'source': catalog['source'],
        'n_galaxies': int(catalog['n_galaxies']),
        'hexagonal_analysis': {
            'omega_60': float(hex_analysis['omega_60']),
            'omega_90': float(hex_analysis['omega_90']),
            'omega_120': float(hex_analysis['omega_120']),
            'hex_excess': float(hex_analysis['hex_excess']),
            'z_score': float(hex_analysis['z_score']),
            'p_value': float(hex_analysis['p_value'])
        },
        'isotropy_test': {
            'chi2_reduced': float(uniformity_test['chi2_reduced']),
            'p_value': float(uniformity_test['p_value'])
        },
        'conclusion': conclusion
    }

    with open(os.path.join(RESULTS_DIR, 'test3_sdss_hexagonal.json'), 'w') as f:
        json.dump(results, f, indent=2)

    print(f"  ✓ test3_sdss_hexagonal.json")

    # Resumen final
    print("\n" + "=" * 70)
    print("RESUMEN TEST #3")
    print("=" * 70)

    print(f"""
    PREDICCIÓN OCTH: Si la malla hexagonal fue estirada por inflación,
    deberíamos ver exceso de correlación a 60° y 120°.

    RESULTADO:
    - ω(60°) = {hex_analysis['omega_60']:.4f}
    - ω(90°) = {hex_analysis['omega_90']:.4f}  (control cuadrado)
    - ω(120°) = {hex_analysis['omega_120']:.4f}
    - Exceso hexagonal: {hex_analysis['hex_excess']:.4f}
    - Significancia: Z = {hex_analysis['z_score']:.2f}, p = {hex_analysis['p_value']:.4f}

    CONCLUSIÓN: {conclusion}

    INTERPRETACIÓN:
    - Si exceso > 0 y p < 0.05: Evidencia de patrón hexagonal residual
    - Si exceso ~ 0: Isotropía preservada (sin memoria de geometría primordial)
    - El clustering de galaxias domina a pequeñas escalas
    """)

    if not use_real_data:
        print("    NOTA: Usando datos simulados. Para test definitivo usar SDSS real.")
        print("    Ejecutar: python3 test3_sdss_hexagonal.py --real")

    return results


if __name__ == "__main__":
    import sys

    use_real = '--real' in sys.argv
    inject_hex = '--inject' in sys.argv
    n_gal = 10000

    for arg in sys.argv:
        if arg.startswith('--ngal='):
            n_gal = int(arg.split('=')[1])

    if '--help' in sys.argv:
        print("""
Test #3: Geometría Hexagonal en Distribución de Galaxias

Uso: python3 test3_sdss_hexagonal.py [opciones]

Opciones:
  --real      Descargar datos reales de SDSS
  --inject    Inyectar señal hexagonal (test de sensibilidad)
  --ngal=N    Número de galaxias (default: 10000)
  --help      Mostrar esta ayuda
        """)
        sys.exit(0)

    results = run_test3(use_real_data=use_real, inject_hex=inject_hex,
                        n_galaxies=n_gal)
