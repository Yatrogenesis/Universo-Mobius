"""
JWST Little Red Dots - OCTH Antipodal Correlation Analysis
===========================================================

Fecha: 2026-01-12
Autor: Francisco Molina-Burgos
Afiliacion: Avermex Research Division, Merida, Yucatan, Mexico

Este script analiza el catalogo ALT DR1 de JWST buscando correlaciones
antipodrales predichas por la teoria OCTH (estructura de Mobius).

Hipotesis OCTH:
- Si los LRDs son defectos topologicos primordiales en la banda de Mobius,
  deberian mostrar anti-correlacion en direcciones antipodrales (180 grados)
- Similar a la señal detectada en CMB (Z = -6.08 sigma)
"""

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from pathlib import Path
import json
from datetime import datetime
from scipy import stats

# Configuracion
DATA_DIR = Path(r"H:\Claude dev\Universo-Mobius\data\jwst_lrd")
RESULTS_DIR = Path(r"H:\Claude dev\Universo-Mobius\results\jwst_lrd")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def load_alt_catalog():
    """Cargar catalogo ALT DR1 de JWST"""
    fits_path = DATA_DIR / "ALT_DR1_public.fits"
    with fits.open(fits_path) as hdul:
        data = hdul[1].data
        catalog = {
            'id': data['id'],
            'ra': data['ra'],
            'dec': data['dec'],
            'z': data['z_ALT'],
            'n_lines': data['n_lines_detected']
        }
    return catalog


def compute_angular_separations(ra, dec):
    """Calcular separaciones angulares entre todos los pares de fuentes"""
    coords = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    n = len(coords)

    # Calcular separaciones (solo triangulo superior)
    separations = []
    for i in range(n):
        for j in range(i+1, n):
            sep = coords[i].separation(coords[j]).degree
            separations.append(sep)

    return np.array(separations)


def compute_antipodal_correlation(ra, dec, n_bins=36, n_bootstrap=1000):
    """
    Calcular correlacion antipodral similar al analisis CMB

    En OCTH, esperamos anti-correlacion cerca de 180 grados
    (fuentes en lados opuestos del universo Mobius)
    """
    coords = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    n = len(coords)

    # Crear histograma de separaciones angulares
    bins = np.linspace(0, 180, n_bins + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    separations = compute_angular_separations(ra, dec)
    hist, _ = np.histogram(separations, bins=bins)

    # Normalizar por area esferica (correccion geometrica)
    # El numero esperado de pares en un bin depende del area del casquete esferico
    bin_areas = np.abs(np.cos(np.radians(bins[:-1])) - np.cos(np.radians(bins[1:])))
    expected = hist.sum() * bin_areas / bin_areas.sum()

    # Calcular exceso/deficit relativo al esperado
    with np.errstate(divide='ignore', invalid='ignore'):
        excess = (hist - expected) / np.sqrt(expected)
        excess = np.nan_to_num(excess, nan=0.0, posinf=0.0, neginf=0.0)

    # Bootstrap para estimar errores
    bootstrap_excess = np.zeros((n_bootstrap, n_bins))
    for b in range(n_bootstrap):
        idx = np.random.choice(n, size=n, replace=True)
        boot_sep = compute_angular_separations(ra[idx], dec[idx])
        boot_hist, _ = np.histogram(boot_sep, bins=bins)
        boot_expected = boot_hist.sum() * bin_areas / bin_areas.sum()
        with np.errstate(divide='ignore', invalid='ignore'):
            boot_excess = (boot_hist - boot_expected) / np.sqrt(boot_expected)
            boot_excess = np.nan_to_num(boot_excess, nan=0.0, posinf=0.0, neginf=0.0)
        bootstrap_excess[b] = boot_excess

    excess_std = np.std(bootstrap_excess, axis=0)

    # Significancia en el bin antipodral (170-180 grados)
    antipodal_idx = np.where(bin_centers >= 170)[0]
    if len(antipodal_idx) > 0:
        antipodal_signal = np.mean(excess[antipodal_idx])
        antipodal_std = np.mean(excess_std[antipodal_idx])
        z_score_antipodal = antipodal_signal / antipodal_std if antipodal_std > 0 else 0
    else:
        antipodal_signal = 0
        z_score_antipodal = 0

    return {
        'bin_centers': bin_centers.tolist(),
        'excess': excess.tolist(),
        'excess_std': excess_std.tolist(),
        'antipodal_signal': float(antipodal_signal),
        'z_score_antipodal': float(z_score_antipodal),
        'n_sources': n,
        'n_pairs': len(separations)
    }


def analyze_redshift_distribution(z):
    """Analizar distribucion de redshift"""
    valid_z = z[z > 0]
    return {
        'n_valid': len(valid_z),
        'z_min': float(np.min(valid_z)),
        'z_max': float(np.max(valid_z)),
        'z_mean': float(np.mean(valid_z)),
        'z_median': float(np.median(valid_z)),
        'z_std': float(np.std(valid_z))
    }


def analyze_hexagonal_pattern(ra, dec, n_bins=12):
    """
    Buscar patron hexagonal en distribucion angular

    OCTH predice exceso de pares separados por multiplos de 60 grados
    debido a la estructura hexagonal del reticulado de Planck
    """
    separations = compute_angular_separations(ra, dec)

    # Histograma modulo 60 grados
    sep_mod60 = separations % 60
    bins = np.linspace(0, 60, n_bins + 1)
    hist, _ = np.histogram(sep_mod60, bins=bins)

    # Test chi-cuadrado contra distribucion uniforme
    expected = np.mean(hist)
    chi2, p_value = stats.chisquare(hist)

    # Buscar picos cerca de 0 y 60 (vertices del hexagono)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    peak_idx = np.argmax(hist)
    peak_angle = bin_centers[peak_idx]

    return {
        'bin_centers': bin_centers.tolist(),
        'counts': hist.tolist(),
        'chi2': float(chi2),
        'p_value': float(p_value),
        'peak_angle': float(peak_angle),
        'is_hexagonal': p_value < 0.05 and (peak_angle < 10 or peak_angle > 50)
    }


def main():
    print("=" * 70)
    print("JWST Little Red Dots - OCTH Analysis")
    print("=" * 70)
    print(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    # Cargar datos
    print("Cargando catalogo ALT DR1...")
    catalog = load_alt_catalog()
    ra = np.array(catalog['ra'])
    dec = np.array(catalog['dec'])
    z = np.array(catalog['z'])

    print(f"  Fuentes totales: {len(ra)}")
    print(f"  RA range: [{ra.min():.2f}, {ra.max():.2f}] deg")
    print(f"  Dec range: [{dec.min():.2f}, {dec.max():.2f}] deg")
    print()

    # Filtrar por redshift alto (z > 4) para enfocarnos en LRDs
    high_z_mask = z > 4
    ra_hz = ra[high_z_mask]
    dec_hz = dec[high_z_mask]
    z_hz = z[high_z_mask]
    print(f"Fuentes con z > 4 (candidatos LRD): {len(ra_hz)}")
    print()

    # Analisis de redshift
    print("Analizando distribucion de redshift...")
    z_stats = analyze_redshift_distribution(z)
    print(f"  z range: [{z_stats['z_min']:.2f}, {z_stats['z_max']:.2f}]")
    print(f"  z mean: {z_stats['z_mean']:.2f}")
    print(f"  z median: {z_stats['z_median']:.2f}")
    print()

    # Correlacion antipodral (todas las fuentes)
    print("Calculando correlacion antipodral (todas las fuentes)...")
    print("  (Esto puede tomar unos minutos con bootstrap)")
    antipodal_all = compute_antipodal_correlation(ra, dec, n_bootstrap=100)
    print(f"  Señal antipodal: {antipodal_all['antipodal_signal']:.3f}")
    print(f"  Z-score: {antipodal_all['z_score_antipodal']:.2f} sigma")
    print()

    # Correlacion antipodral (solo alto z)
    if len(ra_hz) > 10:
        print("Calculando correlacion antipodral (z > 4)...")
        antipodal_hz = compute_antipodal_correlation(ra_hz, dec_hz, n_bootstrap=100)
        print(f"  Señal antipodal: {antipodal_hz['antipodal_signal']:.3f}")
        print(f"  Z-score: {antipodal_hz['z_score_antipodal']:.2f} sigma")
        print()
    else:
        antipodal_hz = None
        print("Insuficientes fuentes con z > 4 para analisis")
        print()

    # Patron hexagonal
    print("Buscando patron hexagonal...")
    hexagonal = analyze_hexagonal_pattern(ra, dec)
    print(f"  Chi2: {hexagonal['chi2']:.2f}")
    print(f"  p-value: {hexagonal['p_value']:.4f}")
    print(f"  Pico en: {hexagonal['peak_angle']:.1f} grados")
    print(f"  Patron hexagonal detectado: {hexagonal['is_hexagonal']}")
    print()

    # Guardar resultados
    results = {
        'analysis_date': datetime.now().isoformat(),
        'catalog': 'ALT_DR1_public.fits',
        'n_sources_total': len(ra),
        'n_sources_high_z': len(ra_hz),
        'redshift_stats': z_stats,
        'antipodal_all': antipodal_all,
        'antipodal_high_z': antipodal_hz,
        'hexagonal_pattern': hexagonal,
        'octh_interpretation': {
            'hypothesis': 'LRDs as primordial topological defects in Mobius structure',
            'prediction_antipodal': 'Anti-correlation at 180 degrees (like CMB)',
            'prediction_hexagonal': 'Excess at 0/60/120/180 degree separations',
            'cmb_reference_z_score': -6.08
        }
    }

    results_path = RESULTS_DIR / "JWST_LRD_OCTH_analysis.json"
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Resultados guardados en: {results_path}")

    # Resumen OCTH
    print()
    print("=" * 70)
    print("INTERPRETACION OCTH")
    print("=" * 70)

    z_antipodal = antipodal_all['z_score_antipodal']
    if z_antipodal < -2:
        print(f"ANTI-CORRELACION DETECTADA: Z = {z_antipodal:.2f} sigma")
        print("Consistente con prediccion OCTH de estructura Mobius")
    elif z_antipodal > 2:
        print(f"CORRELACION POSITIVA: Z = {z_antipodal:.2f} sigma")
        print("Inconsistente con OCTH - requiere investigacion adicional")
    else:
        print(f"Sin señal significativa: Z = {z_antipodal:.2f} sigma")
        print("Datos insuficientes o cobertura del cielo limitada")

    if hexagonal['is_hexagonal']:
        print()
        print("PATRON HEXAGONAL DETECTADO")
        print("Consistente con reticulo hexagonal de Planck (OCTH)")

    print()
    print("=" * 70)

    return results


if __name__ == "__main__":
    results = main()
