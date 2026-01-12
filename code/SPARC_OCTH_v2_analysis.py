"""
SPARC Rotation Curves - OCTH v2 Analysis
=========================================

Fecha: 2026-01-12
Autor: Francisco Molina-Burgos
Afiliacion: Avermex Research Division, Merida, Yucatan, Mexico

Version 2: Modelo OCTH mejorado con permeabilidad fisicamente motivada.

El modelo v1 usaba transicion exponencial que no captura la fisica.
El modelo v2 usa una permeabilidad que naturalmente produce curvas planas.

Clave: Para v^2 = GM/r * (1/Psi) = constante a gran r
       Necesitamos Psi(r) ~ r/r_0 a grandes distancias

Nuevo modelo de permeabilidad:
Psi(r) = Psi_core + (r/r_0) * f(r/r_0)

donde f es una funcion de interpolacion que asegura Psi > 0 siempre.
"""

import numpy as np
from scipy.optimize import minimize
from pathlib import Path
import json
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Configuracion
DATA_DIR = Path(r"H:\Claude dev\Universo-Mobius\data\sparc\rotation_curves")
RESULTS_DIR = Path(r"H:\Claude dev\Universo-Mobius\results\sparc_octh_v2")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Constantes
G = 4.302e-6  # kpc * (km/s)^2 / M_sun


def load_rotation_curve(filepath):
    """Cargar curva de rotacion SPARC"""
    data = {'r': [], 'v_obs': [], 'v_err': [], 'v_gas': [], 'v_disk': [], 'v_bul': []}
    distance = None

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('# Distance'):
                distance = float(line.split('=')[1].replace('Mpc', '').strip())
            elif line.startswith('#') or not line:
                continue
            else:
                parts = line.split()
                if len(parts) >= 6:
                    data['r'].append(float(parts[0]))
                    data['v_obs'].append(float(parts[1]))
                    data['v_err'].append(float(parts[2]))
                    data['v_gas'].append(float(parts[3]))
                    data['v_disk'].append(float(parts[4]))
                    data['v_bul'].append(float(parts[5]))

    for key in data:
        data[key] = np.array(data[key])

    data['distance'] = distance
    data['v_bar'] = np.sqrt(data['v_gas']**2 + data['v_disk']**2 + data['v_bul']**2)

    return data


def v_octh_v2(r, v_bar, r_0, psi_core, alpha):
    """
    OCTH v2: Modelo de permeabilidad mejorado

    Psi(r) = psi_core + (r/r_0)^alpha

    Parametros:
    - r_0: escala de transicion (kpc)
    - psi_core: permeabilidad central (0 < psi_core < 1)
    - alpha: exponente (0 < alpha < 2, tipicamente ~1)

    Para alpha = 1: v -> constante a gran r (curvas planas)
    Para alpha < 1: v sube a gran r
    Para alpha > 1: v baja a gran r
    """
    x = r / r_0
    psi = psi_core + np.power(x, alpha)
    psi = np.maximum(psi, 0.01)  # Evitar division por cero
    v_squared = v_bar**2 / psi
    return np.sqrt(np.maximum(v_squared, 0))


def v_octh_v2b(r, v_bar, r_0, v_flat):
    """
    OCTH v2b: Modelo simplificado con v_flat explicito

    Definimos Psi tal que v -> v_flat a grandes distancias:

    v^2 = v_bar^2 / Psi = v_flat^2 a gran r
    => Psi = v_bar^2 / v_flat^2 a gran r

    Modelo de interpolacion:
    Psi(r) = Psi_newton * (1 - f) + Psi_flat * f

    donde f = tanh(r/r_0) es funcion de transicion
    y Psi_newton = 1 (comportamiento Newtoniano cerca)
    """
    f = np.tanh(r / r_0)

    # Psi que da comportamiento Newtoniano
    psi_newton = 1.0

    # Psi que da curva plana con v_flat
    # Cuidado: v_bar puede ser cero, evitar division
    v_bar_safe = np.maximum(v_bar, 1.0)
    psi_flat = (v_bar_safe / v_flat)**2
    psi_flat = np.maximum(psi_flat, 0.01)

    # Interpolacion
    psi = psi_newton * (1 - f) + psi_flat * f
    psi = np.maximum(psi, 0.01)

    v_squared = v_bar**2 / psi
    return np.sqrt(np.maximum(v_squared, 0))


def v_octh_v2c(r, v_bar, r_0, beta):
    """
    OCTH v2c: Modelo tipo MOND pero con base metrica

    En MOND: a = sqrt(a_N * a_0) para a_N << a_0
    En OCTH: v^2 = v_bar^2 / Psi donde Psi depende de aceleracion local

    Psi = 1 / (1 + (r_0/r)^beta)

    A pequeño r: Psi -> 0, v >> v_bar
    A grande r: Psi -> 1, v -> v_bar

    Pero esto es al reves... necesitamos:
    A pequeño r: Newtoniano
    A grande r: Curva plana (v > v_bar_Keplerian)

    Mejor:
    Psi = 1 / (1 + (r/r_0)^beta)
    """
    x = r / r_0
    psi = 1.0 / (1.0 + np.power(x, beta))
    psi = np.maximum(psi, 0.01)
    v_squared = v_bar**2 / psi
    return np.sqrt(np.maximum(v_squared, 0))


def v_nfw(r, v_bar, r_s, rho_0):
    """NFW dark matter halo model"""
    x = r / r_s
    m_nfw = 4 * np.pi * rho_0 * r_s**3 * (np.log(1 + x) - x / (1 + x))
    v_dm_sq = G * m_nfw / r
    v_dm_sq = np.maximum(v_dm_sq, 0)
    v_total_sq = v_bar**2 + v_dm_sq
    return np.sqrt(np.maximum(v_total_sq, 0))


def v_mond(r, v_bar, a_0=1.2e-10):
    """MOND model"""
    a_bar = v_bar**2 / r
    a_bar = np.maximum(a_bar, 1e-15)
    x = a_bar / a_0
    nu = 1 / (1 - np.exp(-np.sqrt(x)))
    nu = np.where(np.isfinite(nu), nu, 1.0)
    v_mond_sq = v_bar**2 * nu
    return np.sqrt(np.maximum(v_mond_sq, 0))


def chi_squared(v_model, v_obs, v_err, n_params):
    """Chi-cuadrado reducido"""
    chi2 = np.sum(((v_obs - v_model) / v_err)**2)
    dof = len(v_obs) - n_params
    if dof <= 0:
        dof = 1
    return chi2 / dof


def fit_octh_v2c(data):
    """Ajustar modelo OCTH v2c"""
    r = data['r']
    v_obs = data['v_obs']
    v_err = data['v_err']
    v_bar = data['v_bar']

    def residuals(params):
        r_0, beta = params
        if r_0 <= 0 or beta <= 0 or beta > 3:
            return 1e10
        v_model = v_octh_v2c(r, v_bar, r_0, beta)
        return np.sum(((v_obs - v_model) / v_err)**2)

    # Probar varios puntos iniciales
    best_result = None
    best_chi2 = 1e10

    for r0_init in [1.0, 5.0, 10.0, 20.0]:
        for beta_init in [0.5, 1.0, 1.5, 2.0]:
            try:
                result = minimize(residuals, x0=[r0_init, beta_init],
                                bounds=[(0.1, 100), (0.1, 3.0)],
                                method='L-BFGS-B')
                if result.fun < best_chi2:
                    best_chi2 = result.fun
                    best_result = result
            except:
                continue

    if best_result is None:
        return {'chi2_red': 1e10, 'success': False}

    r_0, beta = best_result.x
    v_model = v_octh_v2c(r, v_bar, r_0, beta)
    chi2_red = chi_squared(v_model, v_obs, v_err, 2)

    return {
        'r_0': float(r_0),
        'beta': float(beta),
        'chi2_red': float(chi2_red),
        'v_model': v_model.tolist(),
        'success': best_result.success
    }


def fit_octh_v2b(data):
    """Ajustar modelo OCTH v2b (con v_flat explicito)"""
    r = data['r']
    v_obs = data['v_obs']
    v_err = data['v_err']
    v_bar = data['v_bar']

    # Estimar v_flat de datos
    v_flat_est = np.mean(v_obs[-3:]) if len(v_obs) >= 3 else v_obs[-1]

    def residuals(params):
        r_0, v_flat = params
        if r_0 <= 0 or v_flat <= 0:
            return 1e10
        v_model = v_octh_v2b(r, v_bar, r_0, v_flat)
        return np.sum(((v_obs - v_model) / v_err)**2)

    result = minimize(residuals, x0=[5.0, v_flat_est],
                     bounds=[(0.1, 100), (10, 500)],
                     method='L-BFGS-B')

    r_0, v_flat = result.x
    v_model = v_octh_v2b(r, v_bar, r_0, v_flat)
    chi2_red = chi_squared(v_model, v_obs, v_err, 2)

    return {
        'r_0': float(r_0),
        'v_flat': float(v_flat),
        'chi2_red': float(chi2_red),
        'v_model': v_model.tolist(),
        'success': result.success
    }


def fit_nfw(data):
    """Ajustar modelo NFW"""
    r = data['r']
    v_obs = data['v_obs']
    v_err = data['v_err']
    v_bar = data['v_bar']

    def residuals(params):
        r_s, log_rho_0 = params
        if r_s <= 0:
            return 1e10
        rho_0 = 10**log_rho_0
        v_model = v_nfw(r, v_bar, r_s, rho_0)
        return np.sum(((v_obs - v_model) / v_err)**2)

    result = minimize(residuals, x0=[10.0, 7.0],
                     bounds=[(0.1, 100), (4, 10)],
                     method='L-BFGS-B')

    r_s, log_rho_0 = result.x
    rho_0 = 10**log_rho_0
    v_model = v_nfw(r, v_bar, r_s, rho_0)
    chi2_red = chi_squared(v_model, v_obs, v_err, 2)

    return {
        'r_s': float(r_s),
        'rho_0': float(rho_0),
        'chi2_red': float(chi2_red),
        'v_model': v_model.tolist(),
        'success': result.success
    }


def fit_mond(data):
    """MOND (sin parametros libres)"""
    r = data['r']
    v_obs = data['v_obs']
    v_err = data['v_err']
    v_bar = data['v_bar']

    v_model = v_mond(r, v_bar)
    chi2_red = chi_squared(v_model, v_obs, v_err, 0)

    return {
        'a_0': 1.2e-10,
        'chi2_red': float(chi2_red),
        'v_model': v_model.tolist()
    }


def analyze_galaxy(filepath):
    """Analizar una galaxia"""
    name = filepath.stem.replace('_rotmod', '')
    print(f"  {name}...", end=' ')

    try:
        data = load_rotation_curve(filepath)

        if len(data['r']) < 5:
            print("SKIP")
            return None

        # Ajustar modelos
        octh_v2c = fit_octh_v2c(data)
        octh_v2b = fit_octh_v2b(data)
        nfw_fit = fit_nfw(data)
        mond_fit = fit_mond(data)

        # Mejor OCTH
        if octh_v2c['chi2_red'] < octh_v2b['chi2_red']:
            octh_best = octh_v2c
            octh_best['model'] = 'v2c'
        else:
            octh_best = octh_v2b
            octh_best['model'] = 'v2b'

        # Determinar mejor modelo global
        chi2_values = {
            'OCTH': octh_best['chi2_red'],
            'NFW': nfw_fit['chi2_red'],
            'MOND': mond_fit['chi2_red']
        }
        best_model = min(chi2_values, key=chi2_values.get)

        result = {
            'name': name,
            'distance_mpc': data['distance'],
            'n_points': len(data['r']),
            'r_max_kpc': float(data['r'].max()),
            'v_max_obs': float(data['v_obs'].max()),
            'octh_v2c': octh_v2c,
            'octh_v2b': octh_v2b,
            'octh_best': octh_best,
            'nfw': nfw_fit,
            'mond': mond_fit,
            'best_model': best_model,
            'octh_vs_nfw': octh_best['chi2_red'] / nfw_fit['chi2_red'] if nfw_fit['chi2_red'] > 0 else 999
        }

        status = "OCTH!" if best_model == 'OCTH' else best_model
        print(f"chi2: OCTH={octh_best['chi2_red']:.1f}, NFW={nfw_fit['chi2_red']:.1f} -> {status}")

        return result

    except Exception as e:
        print(f"ERROR: {e}")
        return None


def main():
    print("=" * 70)
    print("SPARC Rotation Curves - OCTH v2 Analysis")
    print("Improved permeability models")
    print("=" * 70)
    print(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    files = list(DATA_DIR.glob("*_rotmod.dat"))
    print(f"Galaxias: {len(files)}")
    print()

    results = []
    octh_wins = 0
    nfw_wins = 0
    mond_wins = 0

    for filepath in files:
        result = analyze_galaxy(filepath)
        if result:
            results.append(result)
            if result['best_model'] == 'OCTH':
                octh_wins += 1
            elif result['best_model'] == 'NFW':
                nfw_wins += 1
            else:
                mond_wins += 1

    # Estadisticas
    print()
    print("=" * 70)
    print("RESULTADOS OCTH v2")
    print("=" * 70)
    print(f"Galaxias analizadas: {len(results)}")
    print()
    print(f"Mejor modelo:")
    print(f"  OCTH: {octh_wins} ({100*octh_wins/len(results):.1f}%)")
    print(f"  NFW:  {nfw_wins} ({100*nfw_wins/len(results):.1f}%)")
    print(f"  MOND: {mond_wins} ({100*mond_wins/len(results):.1f}%)")
    print()

    # Chi2 promedio
    chi2_octh = np.mean([r['octh_best']['chi2_red'] for r in results])
    chi2_nfw = np.mean([r['nfw']['chi2_red'] for r in results])
    chi2_mond = np.mean([r['mond']['chi2_red'] for r in results])

    print(f"Chi2 reducido promedio:")
    print(f"  OCTH: {chi2_octh:.2f}")
    print(f"  NFW:  {chi2_nfw:.2f}")
    print(f"  MOND: {chi2_mond:.2f}")
    print(f"  Ratio OCTH/NFW: {chi2_octh/chi2_nfw:.2f}")
    print()

    # Galaxias donde OCTH gana claramente
    octh_clear = [r for r in results if r['octh_vs_nfw'] < 0.9]
    print(f"Galaxias donde OCTH < 90% chi2 de NFW: {len(octh_clear)}")
    for r in sorted(octh_clear, key=lambda x: x['octh_vs_nfw'])[:15]:
        print(f"  {r['name']}: ratio={r['octh_vs_nfw']:.2f}, model={r['octh_best'].get('model', 'v2c')}")

    # Guardar
    output = {
        'analysis_date': datetime.now().isoformat(),
        'n_galaxies': len(results),
        'summary': {
            'octh_wins': octh_wins,
            'nfw_wins': nfw_wins,
            'mond_wins': mond_wins,
            'chi2_octh_mean': float(chi2_octh),
            'chi2_nfw_mean': float(chi2_nfw),
            'chi2_mond_mean': float(chi2_mond),
            'ratio_octh_nfw': float(chi2_octh/chi2_nfw)
        },
        'galaxies': results
    }

    output_path = RESULTS_DIR / "SPARC_OCTH_v2_analysis.json"
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nGuardado: {output_path}")

    return output


if __name__ == "__main__":
    results = main()
