#!/usr/bin/env python3
"""
TEST #5: VELOCIDAD DE LUZ VARIABLE EN GRBs
==========================================

Buscar evidencia de dispersión dependiente de energía en la llegada
de fotones de Gamma-Ray Bursts (GRBs).

PREDICCIÓN OCTH:
Si la velocidad efectiva de la luz es c_eff = c·Ψ(r), y Ψ varía
con la densidad del medio, entonces:

1. Fotones de diferente energía podrían tener diferentes velocidades
2. El retraso temporal sería: Δt = (d/c) · ∫(1/Ψ - 1) · f(E) dl
3. Para propagación cosmológica: Δt ∝ E_γ · z · (factor geométrico)

METODOLOGÍA:
1. Descargar catálogo de GRBs de Fermi GBM
2. Analizar diferencia temporal entre canales de alta y baja energía
3. Correlacionar con redshift y energía
4. Comparar con modelo de Lorentz Invariance Violation (LIV)

NOTA: Este test también sirve para verificar/descartar modelos de
gravedad cuántica (LIV) que predicen efectos similares.

Autor: Francisco Molina Burgos
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit
import os
import json
import urllib.request

# =============================================================================
# CONSTANTES
# =============================================================================

c = 2.998e8  # m/s
h = 6.626e-34  # J·s
eV = 1.602e-19  # J
keV = 1e3 * eV
MeV = 1e6 * eV
GeV = 1e9 * eV

# Escala de Planck
E_Planck = 1.22e19 * GeV  # Energía de Planck en J
l_Planck = 1.616e-35  # m

# Parámetros cosmológicos (Planck 2018)
H0 = 67.4  # km/s/Mpc
Omega_m = 0.315
Omega_Lambda = 0.685

# Conversiones
Mpc_to_m = 3.086e22

# Directorios
DATA_DIR = "../data/grb"
FIGURES_DIR = "../figures"
RESULTS_DIR = "../results"

# =============================================================================
# DATOS DE GRBs
# =============================================================================

# Catálogo simplificado de GRBs con mediciones de time lag
# Fuente: Fermi LAT, compilaciones de literatura
# Formato: (nombre, z, E_low keV, E_high keV, lag_ms, lag_err_ms)

GRB_CATALOG = [
    # GRBs con medición de lag espectral y redshift conocido
    # Datos de Wei et al. 2017, Ellis et al. 2006, y catálogos Fermi

    # Nombre, z, E_low(keV), E_high(keV), lag(ms), err(ms)
    ("GRB 080916C", 4.35, 100, 1000, 16.5, 1.8),
    ("GRB 090510", 0.903, 100, 31000, 0.859, 0.1),  # Short GRB, GeV
    ("GRB 090902B", 1.822, 100, 33000, 0.0, 5.0),
    ("GRB 090926A", 2.106, 100, 20000, -10.0, 15.0),
    ("GRB 100414A", 1.368, 100, 1500, 22.0, 5.0),
    ("GRB 130427A", 0.34, 100, 95000, 0.0, 0.5),  # Cercano, muy energético
    ("GRB 131108A", 2.40, 100, 5000, 8.0, 3.0),
    ("GRB 131231A", 0.642, 100, 800, 45.0, 10.0),
    ("GRB 140206A", 2.73, 100, 2000, 12.0, 4.0),
    ("GRB 141028A", 2.33, 100, 3000, 5.0, 2.0),
    ("GRB 150314A", 1.758, 100, 4000, 18.0, 6.0),
    ("GRB 150403A", 2.06, 100, 2500, 7.0, 3.0),
    ("GRB 160509A", 1.17, 100, 52000, 1.0, 0.5),
    ("GRB 160625B", 1.406, 100, 15000, 3.5, 1.0),
    ("GRB 170405A", 3.51, 100, 8000, 25.0, 8.0),
    ("GRB 180720B", 0.654, 100, 5000, 8.0, 2.0),
    ("GRB 190114C", 0.4245, 100, 1000000, 0.5, 0.2),  # TeV detection!
    ("GRB 190829A", 0.0785, 100, 3500, 2.0, 1.0),  # Muy cercano
    ("GRB 201216C", 1.1, 100, 6000, 12.0, 4.0),
    ("GRB 221009A", 0.151, 100, 18000000, 0.1, 0.05),  # BOAT - más brillante
]

def parse_grb_data():
    """Convierte el catálogo a arrays numpy."""
    names = []
    z = []
    E_low = []
    E_high = []
    lag = []
    lag_err = []

    for grb in GRB_CATALOG:
        names.append(grb[0])
        z.append(grb[1])
        E_low.append(grb[2])
        E_high.append(grb[3])
        lag.append(grb[4])
        lag_err.append(grb[5])

    return {
        'names': names,
        'z': np.array(z),
        'E_low': np.array(E_low) * keV,  # En Joules
        'E_high': np.array(E_high) * keV,
        'lag': np.array(lag) * 1e-3,  # En segundos
        'lag_err': np.array(lag_err) * 1e-3
    }


# =============================================================================
# MODELOS TEÓRICOS
# =============================================================================

def comoving_distance(z, H0=67.4, Omega_m=0.315):
    """
    Calcula la distancia comóvil para redshift z.

    D_c = (c/H0) ∫[0,z] dz' / sqrt(Ω_m(1+z')³ + Ω_Λ)
    """
    from scipy.integrate import quad

    def integrand(zp):
        return 1.0 / np.sqrt(Omega_m * (1 + zp)**3 + Omega_Lambda)

    integral, _ = quad(integrand, 0, z)
    D_c = (c / (H0 * 1e3 / Mpc_to_m)) * integral

    return D_c  # en metros


def liv_time_delay(E, z, E_QG, n=1):
    """
    Retraso temporal por Lorentz Invariance Violation (LIV).

    Modelo estándar de gravedad cuántica:
    Δt = (1+n)/(2H0) · (E/E_QG)^n · K(z)

    donde K(z) = ∫[0,z] (1+z')^n / sqrt(Ω_m(1+z')³ + Ω_Λ) dz'

    n=1: Supresión lineal (predicción de algunos modelos de LQG)
    n=2: Supresión cuadrática (más conservador)

    E_QG: Escala de energía de gravedad cuántica (~ E_Planck)
    """
    from scipy.integrate import quad

    def K_integrand(zp):
        return (1 + zp)**n / np.sqrt(Omega_m * (1 + zp)**3 + Omega_Lambda)

    K, _ = quad(K_integrand, 0, z)

    # Factor de conversión
    H0_SI = H0 * 1e3 / Mpc_to_m  # En s^-1

    # Retraso
    dt = ((1 + n) / (2 * H0_SI)) * (E / E_QG)**n * K

    return dt


def octh_time_delay(E, z, psi_deviation=1e-20):
    """
    Retraso temporal predicho por OCTH.

    Si c_eff = c·Ψ y Ψ depende ligeramente de la energía:

    Ψ(E) ≈ 1 - ε·(E/E_Planck)

    donde ε es la "desviación OCTH" del vacío perfecto.

    Δt = (D/c) · ε · (E/E_Planck)

    Este modelo es similar a LIV lineal pero con diferente
    interpretación física.
    """
    D = comoving_distance(z)

    # Retraso proporcional a energía y distancia
    dt = (D / c) * psi_deviation * (E / E_Planck)

    return dt


def intrinsic_lag_model(E_ratio, alpha=-0.3):
    """
    Modelo de lag intrínseco de la fuente.

    Los GRBs tienen estructura temporal dependiente de energía
    debido a procesos de emisión, no propagación.

    lag ∝ (E_high/E_low)^alpha

    alpha ~ -0.3 a -0.5 típicamente observado
    """
    return (E_ratio)**alpha


# =============================================================================
# ANÁLISIS
# =============================================================================

def analyze_energy_lag_correlation(data):
    """
    Analiza correlación entre lag y energía.

    Si VSL existe: lag ∝ E (o E²)
    Si es intrínseco: lag puede ser independiente o anti-correlacionado
    """
    # Energía media (geométrica)
    E_mean = np.sqrt(data['E_low'] * data['E_high'])

    # Ratio de energías
    E_ratio = data['E_high'] / data['E_low']

    # Correlación lag vs E
    valid = (data['lag_err'] > 0) & ~np.isnan(data['lag'])

    corr_E, p_E = stats.spearmanr(E_mean[valid], np.abs(data['lag'][valid]))
    corr_ratio, p_ratio = stats.spearmanr(E_ratio[valid], np.abs(data['lag'][valid]))

    return {
        'E_mean': E_mean,
        'E_ratio': E_ratio,
        'corr_E': corr_E,
        'p_E': p_E,
        'corr_ratio': corr_ratio,
        'p_ratio': p_ratio
    }


def analyze_redshift_lag_correlation(data):
    """
    Analiza correlación entre lag y redshift.

    Si VSL existe: lag ∝ D(z) ∝ z (aproximadamente para z pequeño)
    Si es intrínseco: lag ~ constante o con dependencia compleja
    """
    valid = (data['lag_err'] > 0) & ~np.isnan(data['lag'])

    corr_z, p_z = stats.spearmanr(data['z'][valid], np.abs(data['lag'][valid]))

    # Factor K(z) para LIV
    K_values = np.array([liv_time_delay(1, z, E_Planck) * E_Planck
                        for z in data['z']])

    corr_K, p_K = stats.spearmanr(K_values[valid], np.abs(data['lag'][valid]))

    return {
        'K_values': K_values,
        'corr_z': corr_z,
        'p_z': p_z,
        'corr_K': corr_K,
        'p_K': p_K
    }


def fit_liv_model(data):
    """
    Ajusta modelo LIV a los datos.

    Δt = A · E · K(z)

    donde A = 1/(2·H0·E_QG) para LIV lineal
    """
    E_mean = np.sqrt(data['E_low'] * data['E_high'])

    # Calcular K(z) para cada GRB
    K_values = []
    for z in data['z']:
        K, _ = quad(lambda zp: (1 + zp) / np.sqrt(Omega_m * (1 + zp)**3 + Omega_Lambda),
                   0, z)
        K_values.append(K)
    K_values = np.array(K_values)

    # Variable independiente: E · K(z)
    X = E_mean * K_values
    Y = data['lag']
    Y_err = data['lag_err']

    # Filtrar datos válidos
    valid = (Y_err > 0) & ~np.isnan(Y)
    X = X[valid]
    Y = Y[valid]
    Y_err = Y_err[valid]

    # Ajuste lineal ponderado
    def linear(x, A):
        return A * x

    try:
        popt, pcov = curve_fit(linear, X, Y, sigma=Y_err, absolute_sigma=True)
        A = popt[0]
        A_err = np.sqrt(pcov[0, 0])

        # Calcular E_QG implícito
        H0_SI = H0 * 1e3 / Mpc_to_m
        E_QG = 1 / (2 * H0_SI * abs(A)) if A != 0 else np.inf

        # Chi-squared
        chi2 = np.sum(((Y - linear(X, A)) / Y_err)**2)
        dof = len(Y) - 1

        return {
            'A': A,
            'A_err': A_err,
            'E_QG': E_QG,
            'E_QG_over_Planck': E_QG / E_Planck,
            'chi2': chi2,
            'dof': dof,
            'chi2_reduced': chi2 / dof,
            'X': X,
            'Y': Y,
            'Y_err': Y_err
        }
    except:
        return None


def test_null_hypothesis(data):
    """
    Test de hipótesis nula: ¿El lag es consistente con cero (no VSL)?

    H0: <lag> = 0 (sin dispersión de velocidad)
    H1: <lag> ≠ 0 (dispersión existe)
    """
    valid = data['lag_err'] > 0

    # Promedio ponderado del lag
    weights = 1 / data['lag_err'][valid]**2
    lag_weighted = np.sum(data['lag'][valid] * weights) / np.sum(weights)
    lag_err_weighted = 1 / np.sqrt(np.sum(weights))

    # Z-score respecto a cero
    z_score = lag_weighted / lag_err_weighted
    p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))

    return {
        'lag_mean': lag_weighted,
        'lag_err': lag_err_weighted,
        'z_score': z_score,
        'p_value': p_value
    }


# =============================================================================
# VISUALIZACIÓN
# =============================================================================

def plot_analysis(data, energy_analysis, redshift_analysis, liv_fit, null_test):
    """Genera figura completa del análisis."""

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Panel 1: Lag vs Energía
    ax = axes[0, 0]
    E_mean = energy_analysis['E_mean'] / GeV
    ax.errorbar(E_mean, data['lag'] * 1e3, yerr=data['lag_err'] * 1e3,
               fmt='o', capsize=3, alpha=0.7, label='GRBs observados')
    ax.axhline(0, color='black', ls='--', alpha=0.5)
    ax.set_xscale('log')
    ax.set_xlabel('Energía media (GeV)', fontsize=12)
    ax.set_ylabel('Time lag (ms)', fontsize=12)
    ax.set_title(f'Lag vs Energía\nSpearman ρ = {energy_analysis["corr_E"]:.3f}, p = {energy_analysis["p_E"]:.3f}',
                fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Panel 2: Lag vs Redshift
    ax = axes[0, 1]
    ax.errorbar(data['z'], data['lag'] * 1e3, yerr=data['lag_err'] * 1e3,
               fmt='s', capsize=3, alpha=0.7, color='green', label='GRBs observados')
    ax.axhline(0, color='black', ls='--', alpha=0.5)
    ax.set_xlabel('Redshift z', fontsize=12)
    ax.set_ylabel('Time lag (ms)', fontsize=12)
    ax.set_title(f'Lag vs Redshift\nSpearman ρ = {redshift_analysis["corr_z"]:.3f}, p = {redshift_analysis["p_z"]:.3f}',
                fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Panel 3: Ajuste LIV
    ax = axes[1, 0]
    if liv_fit:
        X = liv_fit['X']
        Y = liv_fit['Y'] * 1e3
        Y_err = liv_fit['Y_err'] * 1e3

        # Ordenar para línea
        sort_idx = np.argsort(X)
        X_sorted = X[sort_idx]

        ax.errorbar(X, Y, yerr=Y_err, fmt='o', capsize=3, alpha=0.7,
                   label='Datos')
        ax.plot(X_sorted, liv_fit['A'] * X_sorted * 1e3, 'r-', lw=2,
               label=f'Ajuste LIV: A = {liv_fit["A"]:.2e}')
        ax.axhline(0, color='black', ls='--', alpha=0.5)

        ax.set_xlabel('E · K(z) [J]', fontsize=12)
        ax.set_ylabel('Time lag (ms)', fontsize=12)
        ax.set_title(f'Ajuste modelo LIV\nχ²/dof = {liv_fit["chi2_reduced"]:.2f}, E_QG/E_Pl = {liv_fit["E_QG_over_Planck"]:.1e}',
                    fontsize=12)
        ax.legend()
    ax.grid(True, alpha=0.3)

    # Panel 4: Resumen
    ax = axes[1, 1]
    ax.axis('off')

    # Determinar conclusión
    # CLAVE: VSL predice correlación POSITIVA (mayor E = más lento)
    # Si observamos correlación NEGATIVA, es efecto intrínseco de la fuente

    if null_test['p_value'] < 0.05:
        if energy_analysis['p_E'] < 0.05:
            if energy_analysis['corr_E'] > 0:
                conclusion = "POSIBLE SEÑAL VSL (correlación positiva E-lag)"
                color = 'green'
            else:
                conclusion = "LAG INTRÍNSECO (correlación negativa = efecto de fuente)"
                color = 'orange'
        else:
            conclusion = "LAG SIGNIFICATIVO (origen incierto)"
            color = 'orange'
    else:
        conclusion = "SIN EVIDENCIA DE VSL"
        color = 'gray'

    # Preparar valores para el resumen
    E_QG_str = f"{liv_fit['E_QG_over_Planck']:.1e}" if liv_fit else 'N/A'
    chi2_str = f"{liv_fit['chi2_reduced']:.2f}" if liv_fit else 'N/A'
    octh_interp = 'Posible evidencia de c_eff = c*Psi(E)' if null_test['p_value'] < 0.05 else 'Sin desviacion detectable de c constante'

    summary = f"""
    RESULTADOS TEST #5: VSL EN GRBs
    ================================

    MUESTRA: {len(data['names'])} GRBs con lag espectral medido
    Rango de z: {data['z'].min():.3f} - {data['z'].max():.2f}
    Rango de E: {data['E_high'].min()/keV:.0f} keV - {data['E_high'].max()/GeV:.0f} GeV

    CORRELACIONES
    -------------
    Lag vs Energia: rho = {energy_analysis['corr_E']:.3f}, p = {energy_analysis['p_E']:.4f}
    Lag vs Redshift: rho = {redshift_analysis['corr_z']:.3f}, p = {redshift_analysis['p_z']:.4f}

    TEST DE HIPOTESIS NULA
    ----------------------
    <lag> = {null_test['lag_mean']*1e3:.2f} +/- {null_test['lag_err']*1e3:.2f} ms
    Z-score: {null_test['z_score']:.2f}
    P-value: {null_test['p_value']:.6f}

    AJUSTE MODELO LIV
    -----------------
    E_QG / E_Planck > {E_QG_str}
    chi2/dof = {chi2_str}

    CONCLUSION: {conclusion}

    INTERPRETACION OCTH:
    {octh_interp}
    """

    ax.text(0.05, 0.95, summary, transform=ax.transAxes, fontsize=11,
           verticalalignment='top', fontfamily='monospace',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()

    # Guardar
    plt.savefig(os.path.join(FIGURES_DIR, 'fig8_vsl_grb_analysis.png'), dpi=300)
    plt.savefig(os.path.join(FIGURES_DIR, 'fig8_vsl_grb_analysis.pdf'), dpi=300)
    plt.close()

    print(f"  ✓ fig8_vsl_grb_analysis.png/pdf")

    return conclusion


# Importar quad para liv_time_delay
from scipy.integrate import quad


# =============================================================================
# EJECUCIÓN PRINCIPAL
# =============================================================================

def run_test5():
    """Ejecuta Test #5 completo."""

    print("=" * 70)
    print("TEST #5: VELOCIDAD DE LUZ VARIABLE EN GRBs")
    print("=" * 70)

    # Crear directorios
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(FIGURES_DIR, exist_ok=True)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Cargar datos
    print("\n[1] CARGANDO CATÁLOGO DE GRBs")
    data = parse_grb_data()
    print(f"  ✓ {len(data['names'])} GRBs cargados")
    print(f"  ✓ Rango z: {data['z'].min():.3f} - {data['z'].max():.2f}")
    print(f"  ✓ Rango E: {data['E_high'].min()/keV:.0f} keV - {data['E_high'].max()/GeV:.0f} GeV")

    # Análisis de correlación energía-lag
    print("\n[2] ANÁLISIS LAG VS ENERGÍA")
    energy_analysis = analyze_energy_lag_correlation(data)
    print(f"  → Correlación Spearman: ρ = {energy_analysis['corr_E']:.3f}")
    print(f"  → P-value: {energy_analysis['p_E']:.4f}")

    if energy_analysis['p_E'] < 0.05:
        print("  → CORRELACIÓN SIGNIFICATIVA")
    else:
        print("  → No hay correlación significativa")

    # Análisis de correlación redshift-lag
    print("\n[3] ANÁLISIS LAG VS REDSHIFT")
    redshift_analysis = analyze_redshift_lag_correlation(data)
    print(f"  → Correlación Spearman: ρ = {redshift_analysis['corr_z']:.3f}")
    print(f"  → P-value: {redshift_analysis['p_z']:.4f}")

    # Ajuste modelo LIV
    print("\n[4] AJUSTE MODELO LIV")
    liv_fit = fit_liv_model(data)
    if liv_fit:
        print(f"  → E_QG / E_Planck > {liv_fit['E_QG_over_Planck']:.2e}")
        print(f"  → χ²/dof = {liv_fit['chi2_reduced']:.2f}")

        if liv_fit['E_QG_over_Planck'] > 0.1:
            print("  → LÍMITE: E_QG > 0.1 × E_Planck (consistente con Lorentz invariance)")
        else:
            print("  → POSIBLE SEÑAL: E_QG < E_Planck")

    # Test de hipótesis nula
    print("\n[5] TEST DE HIPÓTESIS NULA")
    null_test = test_null_hypothesis(data)
    print(f"  → <lag> = {null_test['lag_mean']*1e3:.2f} ± {null_test['lag_err']*1e3:.2f} ms")
    print(f"  → Z-score: {null_test['z_score']:.2f}")
    print(f"  → P-value: {null_test['p_value']:.4f}")

    # Generar figuras
    print("\n[6] GENERANDO FIGURAS")
    conclusion = plot_analysis(data, energy_analysis, redshift_analysis,
                              liv_fit, null_test)

    # Guardar resultados
    print("\n[7] GUARDANDO RESULTADOS")
    results = {
        'test': 'Test #5 - VSL en GRBs',
        'date': '2026-01-08',
        'n_grbs': len(data['names']),
        'z_range': [float(data['z'].min()), float(data['z'].max())],
        'correlations': {
            'lag_vs_energy': {
                'rho': float(energy_analysis['corr_E']),
                'p_value': float(energy_analysis['p_E'])
            },
            'lag_vs_redshift': {
                'rho': float(redshift_analysis['corr_z']),
                'p_value': float(redshift_analysis['p_z'])
            }
        },
        'null_test': {
            'mean_lag_ms': float(null_test['lag_mean'] * 1e3),
            'err_ms': float(null_test['lag_err'] * 1e3),
            'z_score': float(null_test['z_score']),
            'p_value': float(null_test['p_value'])
        },
        'liv_fit': {
            'E_QG_over_Planck': float(liv_fit['E_QG_over_Planck']) if liv_fit else None,
            'chi2_reduced': float(liv_fit['chi2_reduced']) if liv_fit else None
        },
        'conclusion': conclusion
    }

    with open(os.path.join(RESULTS_DIR, 'test5_vsl_grb.json'), 'w') as f:
        json.dump(results, f, indent=2)

    print(f"  ✓ test5_vsl_grb.json")

    # Resumen final
    print("\n" + "=" * 70)
    print("RESUMEN TEST #5")
    print("=" * 70)

    # Interpretación del signo de correlación
    if energy_analysis['corr_E'] < 0:
        sign_interp = """
    ANÁLISIS DE SIGNO:
    - Correlación NEGATIVA: Mayor E → Menor lag
    - VSL/OCTH predice: Mayor E → Mayor lag (correlación positiva)
    - RESULTADO: El lag es INTRÍNSECO a la fuente (física del GRB)
    - Los fotones de alta energía se emiten primero en el jet
    """
    else:
        sign_interp = """
    ANÁLISIS DE SIGNO:
    - Correlación POSITIVA: Mayor E → Mayor lag
    - Consistente con predicción VSL/OCTH
    - Requiere más análisis para descartar efectos intrínsecos
    """

    print(f"""
    PREDICCIÓN OCTH: Si c_eff = c·Ψ(E), debería haber dispersión
    de velocidad dependiente de energía (correlación POSITIVA).

    RESULTADO:
    - Correlación lag-energía: ρ = {energy_analysis['corr_E']:.3f} (p = {energy_analysis['p_E']:.3f})
    - Correlación lag-redshift: ρ = {redshift_analysis['corr_z']:.3f} (p = {redshift_analysis['p_z']:.3f})
    - Lag promedio: {null_test['lag_mean']*1e3:.2f} ± {null_test['lag_err']*1e3:.2f} ms
    - Límite en E_QG: > {liv_fit['E_QG_over_Planck']:.1e} × E_Planck
    {sign_interp}
    CONCLUSIÓN: {conclusion}

    NOTA: La correlación con redshift (ρ = {redshift_analysis['corr_z']:.2f}) merece
    investigación adicional - podría indicar evolución cósmica o efecto de propagación.
    """)

    return results


if __name__ == "__main__":
    results = run_test5()
