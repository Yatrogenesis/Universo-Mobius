"""
SPARC Rotation Curves - OCTH vs Dark Matter Analysis
=====================================================

Fecha: 2026-01-12
Autor: Francisco Molina-Burgos
Afiliacion: Avermex Research Division, Merida, Yucatan, Mexico

Este script ajusta curvas de rotacion de galaxias SPARC usando:
1. Modelo OCTH (permeabilidad Ψ)
2. Modelo NFW (halo de materia oscura)
3. Modelo MOND (gravedad modificada)

Objetivo: Demostrar que OCTH explica curvas planas sin materia oscura.
"""

import numpy as np
from scipy.optimize import curve_fit, minimize
from pathlib import Path
import json
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Configuracion
DATA_DIR = Path(r"H:\Claude dev\Universo-Mobius\data\sparc\rotation_curves")
RESULTS_DIR = Path(r"H:\Claude dev\Universo-Mobius\results\sparc_octh")
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


def v_octh(r, v_bar, r_psi, psi_0):
    """
    Velocidad OCTH: v^2 = v_bar^2 / Ψ(r)

    Modelo de permeabilidad:
    Ψ(r) = psi_0 + (1 - psi_0) * (1 - exp(-r/r_psi))

    - psi_0: permeabilidad en el centro (< 1)
    - r_psi: escala de transicion

    Cuando r >> r_psi: Ψ → 1 (campo plano)
    Cuando r → 0: Ψ → psi_0 (baja permeabilidad)
    """
    psi = psi_0 + (1 - psi_0) * (1 - np.exp(-r / r_psi))
    # Evitar division por cero
    psi = np.maximum(psi, 0.01)
    v_squared = v_bar**2 / psi
    return np.sqrt(np.maximum(v_squared, 0))


def v_nfw(r, v_bar, r_s, rho_0):
    """
    Velocidad con halo NFW (materia oscura)

    M_dm(r) = 4π ρ_0 r_s^3 [ln(1 + r/r_s) - r/(r+r_s)]
    v_dm^2 = G M_dm / r
    v_total^2 = v_bar^2 + v_dm^2
    """
    x = r / r_s
    # Masa encerrada NFW
    m_nfw = 4 * np.pi * rho_0 * r_s**3 * (np.log(1 + x) - x / (1 + x))
    v_dm_sq = G * m_nfw / r
    v_dm_sq = np.maximum(v_dm_sq, 0)
    v_total_sq = v_bar**2 + v_dm_sq
    return np.sqrt(np.maximum(v_total_sq, 0))


def v_mond(r, v_bar, a_0=1.2e-10):
    """
    Velocidad MOND (Modified Newtonian Dynamics)

    En regimen deep-MOND: a = sqrt(a_N * a_0)
    v^4 = G M a_0 = v_bar^4 * (a_0 / a_N)

    Para simplificar, usamos interpolacion estandar.
    """
    # Aceleracion Newtoniana de barionica
    a_bar = v_bar**2 / r
    # Evitar division por cero
    a_bar = np.maximum(a_bar, 1e-15)

    # Funcion de interpolacion simple (standard)
    x = a_bar / a_0
    # nu(x) = 1 para x >> 1, nu(x) = 1/sqrt(x) para x << 1
    nu = 1 / (1 - np.exp(-np.sqrt(x)))
    nu = np.where(np.isfinite(nu), nu, 1.0)

    v_mond_sq = v_bar**2 * nu
    return np.sqrt(np.maximum(v_mond_sq, 0))


def chi_squared(v_model, v_obs, v_err):
    """Calcular chi-cuadrado reducido"""
    chi2 = np.sum(((v_obs - v_model) / v_err)**2)
    return chi2 / (len(v_obs) - 2)  # -2 por parametros


def fit_octh(data):
    """Ajustar modelo OCTH a curva de rotacion"""
    r = data['r']
    v_obs = data['v_obs']
    v_err = data['v_err']
    v_bar = data['v_bar']

    def residuals(params):
        r_psi, psi_0 = params
        if r_psi <= 0 or psi_0 <= 0 or psi_0 >= 1:
            return 1e10
        v_model = v_octh(r, v_bar, r_psi, psi_0)
        return np.sum(((v_obs - v_model) / v_err)**2)

    # Optimizacion
    result = minimize(residuals, x0=[5.0, 0.5],
                      bounds=[(0.1, 100), (0.01, 0.99)],
                      method='L-BFGS-B')

    r_psi, psi_0 = result.x
    v_model = v_octh(r, v_bar, r_psi, psi_0)
    chi2_red = chi_squared(v_model, v_obs, v_err)

    return {
        'r_psi': float(r_psi),
        'psi_0': float(psi_0),
        'chi2_red': float(chi2_red),
        'v_model': v_model.tolist(),
        'success': result.success
    }


def fit_nfw(data):
    """Ajustar modelo NFW a curva de rotacion"""
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

    # Optimizacion
    result = minimize(residuals, x0=[10.0, 7.0],
                      bounds=[(0.1, 100), (4, 10)],
                      method='L-BFGS-B')

    r_s, log_rho_0 = result.x
    rho_0 = 10**log_rho_0
    v_model = v_nfw(r, v_bar, r_s, rho_0)
    chi2_red = chi_squared(v_model, v_obs, v_err)

    return {
        'r_s': float(r_s),
        'rho_0': float(rho_0),
        'chi2_red': float(chi2_red),
        'v_model': v_model.tolist(),
        'success': result.success
    }


def fit_mond(data):
    """Calcular prediccion MOND (sin parametros libres adicionales)"""
    r = data['r']
    v_obs = data['v_obs']
    v_err = data['v_err']
    v_bar = data['v_bar']

    v_model = v_mond(r, v_bar)
    chi2_red = chi_squared(v_model, v_obs, v_err)

    return {
        'a_0': 1.2e-10,
        'chi2_red': float(chi2_red),
        'v_model': v_model.tolist()
    }


def analyze_galaxy(filepath):
    """Analizar una galaxia completa"""
    name = filepath.stem.replace('_rotmod', '')
    print(f"  Analizando {name}...", end=' ')

    try:
        data = load_rotation_curve(filepath)

        if len(data['r']) < 5:
            print("SKIP (pocos datos)")
            return None

        # Ajustar modelos
        octh_fit = fit_octh(data)
        nfw_fit = fit_nfw(data)
        mond_fit = fit_mond(data)

        # Determinar mejor modelo
        chi2_values = {
            'OCTH': octh_fit['chi2_red'],
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
            'octh': octh_fit,
            'nfw': nfw_fit,
            'mond': mond_fit,
            'best_model': best_model,
            'octh_vs_nfw': octh_fit['chi2_red'] / nfw_fit['chi2_red'] if nfw_fit['chi2_red'] > 0 else 0
        }

        status = "OCTH!" if best_model == 'OCTH' else best_model
        print(f"chi2: OCTH={octh_fit['chi2_red']:.2f}, NFW={nfw_fit['chi2_red']:.2f}, MOND={mond_fit['chi2_red']:.2f} -> {status}")

        return result

    except Exception as e:
        print(f"ERROR: {e}")
        return None


def main():
    print("=" * 70)
    print("SPARC Rotation Curves - OCTH vs Dark Matter Analysis")
    print("=" * 70)
    print(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    # Listar archivos
    files = list(DATA_DIR.glob("*_rotmod.dat"))
    print(f"Galaxias encontradas: {len(files)}")
    print()

    # Analizar todas las galaxias
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
    print("RESULTADOS")
    print("=" * 70)
    print(f"Galaxias analizadas: {len(results)}")
    print()
    print(f"Mejor modelo por galaxia:")
    print(f"  OCTH: {octh_wins} ({100*octh_wins/len(results):.1f}%)")
    print(f"  NFW:  {nfw_wins} ({100*nfw_wins/len(results):.1f}%)")
    print(f"  MOND: {mond_wins} ({100*mond_wins/len(results):.1f}%)")
    print()

    # Chi2 promedio
    chi2_octh = np.mean([r['octh']['chi2_red'] for r in results])
    chi2_nfw = np.mean([r['nfw']['chi2_red'] for r in results])
    chi2_mond = np.mean([r['mond']['chi2_red'] for r in results])

    print(f"Chi2 reducido promedio:")
    print(f"  OCTH: {chi2_octh:.2f}")
    print(f"  NFW:  {chi2_nfw:.2f}")
    print(f"  MOND: {chi2_mond:.2f}")
    print()

    # Parametros OCTH tipicos
    psi_0_values = [r['octh']['psi_0'] for r in results if r['octh']['success']]
    r_psi_values = [r['octh']['r_psi'] for r in results if r['octh']['success']]

    print(f"Parametros OCTH tipicos:")
    print(f"  psi_0 (permeabilidad central): {np.mean(psi_0_values):.3f} +/- {np.std(psi_0_values):.3f}")
    print(f"  r_psi (escala kpc): {np.mean(r_psi_values):.2f} +/- {np.std(r_psi_values):.2f}")
    print()

    # Galaxias donde OCTH es claramente mejor
    octh_clear_wins = [r for r in results if r['octh_vs_nfw'] < 0.8]
    print(f"Galaxias donde OCTH es >20% mejor que NFW: {len(octh_clear_wins)}")
    for r in sorted(octh_clear_wins, key=lambda x: x['octh_vs_nfw'])[:10]:
        print(f"  {r['name']}: OCTH/NFW = {r['octh_vs_nfw']:.2f}")

    # Guardar resultados
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
            'psi_0_mean': float(np.mean(psi_0_values)),
            'psi_0_std': float(np.std(psi_0_values)),
            'r_psi_mean': float(np.mean(r_psi_values)),
            'r_psi_std': float(np.std(r_psi_values))
        },
        'galaxies': results
    }

    output_path = RESULTS_DIR / "SPARC_OCTH_analysis.json"
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nResultados guardados en: {output_path}")

    # Interpretacion
    print()
    print("=" * 70)
    print("INTERPRETACION OCTH")
    print("=" * 70)
    if chi2_octh <= chi2_nfw:
        print("OCTH ajusta las curvas de rotacion TAN BIEN O MEJOR que materia oscura NFW")
        print("SIN necesidad de halo de materia oscura.")
        print()
        print(f"Permeabilidad central tipica: Ψ_0 = {np.mean(psi_0_values):.2f}")
        print("Esto significa que cerca del centro galactico, el tiempo 'fluye'")
        print(f"a ~{np.mean(psi_0_values)*100:.0f}% de la tasa en el vacio.")
        print()
        print("La escala de transicion tipica es ~{:.0f} kpc, comparable".format(np.mean(r_psi_values)))
        print("al tamano tipico de las galaxias.")
    else:
        print("NFW ajusta mejor en promedio, pero OCTH es competitivo")
        print(f"con ratio chi2 = {chi2_octh/chi2_nfw:.2f}")

    print()
    print("=" * 70)

    return output


if __name__ == "__main__":
    results = main()
