"""
SPARC - OCTH con CERO parametros libres
=======================================

Usando la ecuacion derivada:
   Psi(r) = sqrt(a_bar(r) / a_0)

donde a_0 = 1.2e-10 m/s^2 (constante universal)

Este modelo tiene CERO parametros libres por galaxia.
Solo depende de la fisica local (v_bar, r).

Autor: Francisco Molina-Burgos
Fecha: 2026-01-12
"""

import numpy as np
from pathlib import Path
import json
from datetime import datetime

DATA_DIR = Path(r"H:\Claude dev\Universo-Mobius\data\sparc\rotation_curves")
RESULTS_DIR = Path(r"H:\Claude dev\Universo-Mobius\results\sparc_octh_zero_params")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Constante fundamental
# a_0 en (km/s)^2 / kpc
# a_0 = 1.2e-10 m/s^2
# 1 (km/s)^2/kpc = 3.24e-14 m/s^2
# a_0 = 1.2e-10 / 3.24e-14 = 3704 (km/s)^2/kpc
A_0 = 3704.0  # (km/s)^2 / kpc


def load_rotation_curve(filepath):
    data = {'r': [], 'v_obs': [], 'v_err': [], 'v_bar': []}
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
                    r = float(parts[0])
                    v_obs = float(parts[1])
                    v_err = float(parts[2])
                    v_gas = float(parts[3])
                    v_disk = float(parts[4])
                    v_bul = float(parts[5])

                    v_bar = np.sqrt(v_gas**2 + v_disk**2 + v_bul**2)

                    data['r'].append(r)
                    data['v_obs'].append(v_obs)
                    data['v_err'].append(v_err)
                    data['v_bar'].append(v_bar)

    for key in data:
        data[key] = np.array(data[key])

    data['distance'] = distance
    return data


def v_octh_zero_params(r, v_bar):
    """
    OCTH con cero parametros libres.

    Psi(r) = sqrt(a_bar(r) / a_0)
    v^2 = v_bar^2 / Psi = v_bar^2 / sqrt(a_bar/a_0)
        = v_bar^2 * sqrt(a_0/a_bar)
        = v_bar^2 * sqrt(a_0 * r / v_bar^2)
        = v_bar * sqrt(a_0 * r)
    v = sqrt(v_bar * sqrt(a_0 * r))
    v = (v_bar^2 * a_0 * r)^(1/4)
    """
    # Evitar r = 0
    r_safe = np.maximum(r, 0.01)
    v_bar_safe = np.maximum(v_bar, 0.1)

    # Aceleracion barionica
    a_bar = v_bar_safe**2 / r_safe

    # Permeabilidad (limitar para evitar problemas numericos)
    psi = np.sqrt(a_bar / A_0)
    psi = np.clip(psi, 0.01, 10.0)

    # Velocidad OCTH
    v_octh = v_bar_safe / np.sqrt(psi)

    return v_octh, psi


def v_octh_interpolated(r, v_bar):
    """
    OCTH con funcion de interpolacion tipo MOND.

    En lugar de Psi = sqrt(a_bar/a_0) puro, usamos:
    Psi = sqrt(a_bar/a_0) * mu(a_bar/a_0)

    donde mu es funcion de interpolacion:
    mu(x) = x / (1 + x) para transicion suave

    Esto da:
    - a_bar >> a_0: Psi -> sqrt(a_bar/a_0) * 1 = sqrt(a_bar/a_0) >> 1 -> Newtoniano
    - a_bar << a_0: Psi -> sqrt(a_bar/a_0) * (a_bar/a_0) = (a_bar/a_0)^1.5 << 1
    """
    r_safe = np.maximum(r, 0.01)
    v_bar_safe = np.maximum(v_bar, 0.1)

    a_bar = v_bar_safe**2 / r_safe
    x = a_bar / A_0

    # Funcion de interpolacion simple
    # mu(x) -> 1 para x >> 1, mu(x) -> x para x << 1
    mu = x / (1 + x)

    # Permeabilidad con interpolacion
    # En regimen Newtoniano (x>>1): psi -> sqrt(x) * 1 = sqrt(a_bar/a_0)
    # En regimen MOND (x<<1): psi -> sqrt(x) * x = x^1.5

    # Para que v^2 = v_bar^2/psi de resultados MOND:
    # Necesitamos v^4 = v_bar^4 * a_0/a_bar en regimen MOND
    # Esto requiere psi = sqrt(a_bar/a_0) en deep MOND

    # Interpolacion suave
    psi_newton = 1.0  # Newtoniano puro
    psi_mond = np.sqrt(x)  # Deep MOND

    # Funcion de transicion
    nu = 1.0 / np.sqrt(1 + 1/x)  # nu -> 1 para x>>1, nu -> sqrt(x) para x<<1

    psi = psi_newton * (1 - nu) + psi_mond * nu
    psi = np.clip(psi, 0.01, 10.0)

    v_octh = v_bar_safe / np.sqrt(psi)

    return v_octh, psi


def v_mond_standard(r, v_bar):
    """MOND estandar para comparacion"""
    r_safe = np.maximum(r, 0.01)
    v_bar_safe = np.maximum(v_bar, 0.1)

    a_bar = v_bar_safe**2 / r_safe
    x = a_bar / A_0

    # Funcion de interpolacion estandar
    nu = 1.0 / (1.0 - np.exp(-np.sqrt(x)))
    nu = np.where(np.isfinite(nu), nu, 1.0)

    v_mond = v_bar_safe * np.sqrt(nu)

    return v_mond


def chi_squared(v_model, v_obs, v_err, n_params=0):
    """Chi^2 reducido"""
    chi2 = np.sum(((v_obs - v_model) / v_err)**2)
    dof = len(v_obs) - n_params
    if dof <= 0:
        dof = 1
    return chi2 / dof


def analyze_galaxy(filepath):
    name = filepath.stem.replace('_rotmod', '')
    print(f"  {name}...", end=' ')

    try:
        data = load_rotation_curve(filepath)
        if len(data['r']) < 5:
            print("SKIP")
            return None

        r = data['r']
        v_obs = data['v_obs']
        v_err = data['v_err']
        v_bar = data['v_bar']

        # OCTH cero parametros
        v_octh, psi_octh = v_octh_zero_params(r, v_bar)
        chi2_octh = chi_squared(v_octh, v_obs, v_err, n_params=0)

        # OCTH interpolado
        v_octh_int, psi_int = v_octh_interpolated(r, v_bar)
        chi2_octh_int = chi_squared(v_octh_int, v_obs, v_err, n_params=0)

        # MOND para comparacion
        v_mond = v_mond_standard(r, v_bar)
        chi2_mond = chi_squared(v_mond, v_obs, v_err, n_params=0)

        # Mejor modelo
        models = {
            'OCTH_pure': chi2_octh,
            'OCTH_interp': chi2_octh_int,
            'MOND': chi2_mond
        }
        best = min(models, key=models.get)

        result = {
            'name': name,
            'n_points': len(r),
            'r_max': float(r.max()),
            'v_max': float(v_obs.max()),
            'chi2_octh_pure': float(chi2_octh),
            'chi2_octh_interp': float(chi2_octh_int),
            'chi2_mond': float(chi2_mond),
            'psi_mean': float(np.mean(psi_octh)),
            'psi_outer': float(np.mean(psi_octh[-3:])) if len(psi_octh) >= 3 else float(psi_octh[-1]),
            'best_model': best
        }

        print(f"chi2: OCTH={chi2_octh:.1f}, MOND={chi2_mond:.1f} -> {best}")

        return result

    except Exception as e:
        print(f"ERROR: {e}")
        return None


def main():
    print("=" * 70)
    print("SPARC - OCTH con CERO PARAMETROS LIBRES")
    print("=" * 70)
    print(f"Fecha: {datetime.now()}")
    print()
    print("Ecuacion: Psi(r) = sqrt(a_bar(r) / a_0)")
    print(f"a_0 = 1.2e-10 m/s^2 = {A_0} (km/s)^2/kpc")
    print("CERO parametros libres por galaxia")
    print()

    files = list(DATA_DIR.glob("*_rotmod.dat"))
    print(f"Galaxias: {len(files)}")
    print()

    results = []
    octh_pure_wins = 0
    octh_int_wins = 0
    mond_wins = 0

    for f in files:
        r = analyze_galaxy(f)
        if r:
            results.append(r)
            if r['best_model'] == 'OCTH_pure':
                octh_pure_wins += 1
            elif r['best_model'] == 'OCTH_interp':
                octh_int_wins += 1
            else:
                mond_wins += 1

    # Estadisticas
    print()
    print("=" * 70)
    print("RESULTADOS - CERO PARAMETROS")
    print("=" * 70)
    print(f"Analizadas: {len(results)}")
    print(f"OCTH puro gana:    {octh_pure_wins} ({100*octh_pure_wins/len(results):.1f}%)")
    print(f"OCTH interp gana:  {octh_int_wins} ({100*octh_int_wins/len(results):.1f}%)")
    print(f"MOND gana:         {mond_wins} ({100*mond_wins/len(results):.1f}%)")
    print()

    chi2_octh = np.mean([r['chi2_octh_pure'] for r in results])
    chi2_int = np.mean([r['chi2_octh_interp'] for r in results])
    chi2_mond = np.mean([r['chi2_mond'] for r in results])

    print(f"Chi2 promedio (CERO params):")
    print(f"  OCTH puro:   {chi2_octh:.2f}")
    print(f"  OCTH interp: {chi2_int:.2f}")
    print(f"  MOND:        {chi2_mond:.2f}")

    # Mediana (menos sensible a outliers)
    chi2_octh_med = np.median([r['chi2_octh_pure'] for r in results])
    chi2_int_med = np.median([r['chi2_octh_interp'] for r in results])
    chi2_mond_med = np.median([r['chi2_mond'] for r in results])

    print()
    print(f"Chi2 mediana:")
    print(f"  OCTH puro:   {chi2_octh_med:.2f}")
    print(f"  OCTH interp: {chi2_int_med:.2f}")
    print(f"  MOND:        {chi2_mond_med:.2f}")

    # Guardar
    output = {
        'date': datetime.now().isoformat(),
        'n_galaxies': len(results),
        'n_params': 0,
        'a_0': A_0,
        'summary': {
            'octh_pure_wins': octh_pure_wins,
            'octh_int_wins': octh_int_wins,
            'mond_wins': mond_wins,
            'chi2_octh_mean': float(chi2_octh),
            'chi2_octh_interp_mean': float(chi2_int),
            'chi2_mond_mean': float(chi2_mond),
            'chi2_octh_median': float(chi2_octh_med),
            'chi2_mond_median': float(chi2_mond_med)
        },
        'galaxies': results
    }

    out_path = RESULTS_DIR / "SPARC_OCTH_zero_params.json"
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\nGuardado: {out_path}")

    return output


if __name__ == "__main__":
    main()
