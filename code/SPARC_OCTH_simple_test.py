"""
SPARC - Test Simple: Ψ constante por galaxia
============================================

Hipotesis: En OCTH, la permeabilidad NO varia con el radio
dentro de una galaxia, sino que es determinada por la topologia
global del universo.

Si Ψ = constante < 1, entonces:
v² = v_bar² / Ψ = v_bar² × (1/Ψ)

Esto equivale a multiplicar la masa efectiva por (1/Ψ).

Comparacion:
- NFW: M_total = M_bar + M_dark (halo)
- OCTH simple: M_eff = M_bar / Ψ = M_bar × (1/Ψ)
"""

import numpy as np
from scipy.optimize import minimize_scalar
from pathlib import Path
import json
from datetime import datetime

DATA_DIR = Path(r"H:\Claude dev\Universo-Mobius\data\sparc\rotation_curves")
RESULTS_DIR = Path(r"H:\Claude dev\Universo-Mobius\results\sparc_octh_simple")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def load_rotation_curve(filepath):
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


def fit_psi_constant(data):
    """
    Ajustar Ψ constante: v² = v_bar² / Ψ

    Solo UN parametro libre: Ψ
    """
    v_obs = data['v_obs']
    v_err = data['v_err']
    v_bar = data['v_bar']

    # Evitar divisiones por cero
    v_bar_safe = np.maximum(v_bar, 1.0)

    def chi2(psi):
        if psi <= 0.01 or psi > 2.0:
            return 1e10
        v_model = v_bar_safe / np.sqrt(psi)
        return np.sum(((v_obs - v_model) / v_err)**2)

    # Buscar mejor Ψ
    result = minimize_scalar(chi2, bounds=(0.01, 2.0), method='bounded')
    psi_best = result.x

    v_model = v_bar_safe / np.sqrt(psi_best)
    chi2_red = np.sum(((v_obs - v_model) / v_err)**2) / (len(v_obs) - 1)

    # Factor de amplificacion de masa
    mass_factor = 1.0 / psi_best

    return {
        'psi': float(psi_best),
        'mass_factor': float(mass_factor),
        'chi2_red': float(chi2_red),
        'v_model': v_model.tolist(),
        'n_params': 1
    }


def fit_psi_with_transition(data):
    """
    Modelo con transicion: Ψ decrece con r

    Ψ(r) = Ψ_inf + (Ψ_0 - Ψ_inf) × exp(-r/r_t)

    - Ψ_0: permeabilidad en el centro (cerca de 1, Newtoniano)
    - Ψ_inf: permeabilidad a grandes distancias (< 1, curvas planas)
    - r_t: radio de transicion
    """
    r = data['r']
    v_obs = data['v_obs']
    v_err = data['v_err']
    v_bar = data['v_bar']

    v_bar_safe = np.maximum(v_bar, 1.0)

    def chi2(params):
        psi_0, psi_inf, r_t = params
        if psi_0 <= 0 or psi_inf <= 0 or psi_0 > 2 or psi_inf > psi_0 or r_t <= 0:
            return 1e10

        psi = psi_inf + (psi_0 - psi_inf) * np.exp(-r / r_t)
        psi = np.maximum(psi, 0.01)
        v_model = v_bar_safe / np.sqrt(psi)
        return np.sum(((v_obs - v_model) / v_err)**2)

    # Optimizar
    from scipy.optimize import minimize

    best_result = None
    best_chi2 = 1e10

    for psi0 in [0.8, 1.0, 1.2]:
        for psi_inf in [0.1, 0.2, 0.3, 0.5]:
            for rt in [2.0, 5.0, 10.0]:
                if psi_inf >= psi0:
                    continue
                try:
                    result = minimize(chi2, x0=[psi0, psi_inf, rt],
                                    bounds=[(0.5, 1.5), (0.01, 0.99), (0.1, 50)],
                                    method='L-BFGS-B')
                    if result.fun < best_chi2:
                        best_chi2 = result.fun
                        best_result = result
                except:
                    continue

    if best_result is None:
        return {'chi2_red': 1e10, 'success': False, 'n_params': 3}

    psi_0, psi_inf, r_t = best_result.x
    psi = psi_inf + (psi_0 - psi_inf) * np.exp(-r / r_t)
    v_model = v_bar_safe / np.sqrt(np.maximum(psi, 0.01))
    chi2_red = best_chi2 / (len(v_obs) - 3)

    return {
        'psi_0': float(psi_0),
        'psi_inf': float(psi_inf),
        'r_t': float(r_t),
        'chi2_red': float(chi2_red),
        'v_model': v_model.tolist(),
        'success': best_result.success,
        'n_params': 3
    }


def fit_nfw(data):
    """NFW para comparacion"""
    from scipy.optimize import minimize

    G = 4.302e-6
    r = data['r']
    v_obs = data['v_obs']
    v_err = data['v_err']
    v_bar = data['v_bar']

    def residuals(params):
        r_s, log_rho = params
        if r_s <= 0:
            return 1e10
        rho_0 = 10**log_rho
        x = r / r_s
        m_nfw = 4 * np.pi * rho_0 * r_s**3 * (np.log(1 + x) - x / (1 + x))
        v_dm_sq = G * m_nfw / r
        v_dm_sq = np.maximum(v_dm_sq, 0)
        v_model = np.sqrt(v_bar**2 + v_dm_sq)
        return np.sum(((v_obs - v_model) / v_err)**2)

    result = minimize(residuals, x0=[10.0, 7.0],
                     bounds=[(0.1, 100), (4, 10)],
                     method='L-BFGS-B')

    r_s, log_rho = result.x
    rho_0 = 10**log_rho
    x = r / r_s
    m_nfw = 4 * np.pi * rho_0 * r_s**3 * (np.log(1 + x) - x / (1 + x))
    v_dm_sq = G * m_nfw / r
    v_model = np.sqrt(v_bar**2 + np.maximum(v_dm_sq, 0))
    chi2_red = result.fun / (len(v_obs) - 2)

    return {
        'r_s': float(r_s),
        'rho_0': float(rho_0),
        'chi2_red': float(chi2_red),
        'v_model': v_model.tolist(),
        'n_params': 2
    }


def analyze_galaxy(filepath):
    name = filepath.stem.replace('_rotmod', '')
    print(f"  {name}...", end=' ')

    try:
        data = load_rotation_curve(filepath)
        if len(data['r']) < 5:
            print("SKIP")
            return None

        # Ajustar modelos
        octh_const = fit_psi_constant(data)
        octh_trans = fit_psi_with_transition(data)
        nfw = fit_nfw(data)

        # Mejor OCTH
        if octh_trans.get('success', False) and octh_trans['chi2_red'] < octh_const['chi2_red']:
            octh_best = octh_trans
            octh_best['type'] = 'transition'
        else:
            octh_best = octh_const
            octh_best['type'] = 'constant'

        # Comparar
        chi2_octh = octh_best['chi2_red']
        chi2_nfw = nfw['chi2_red']

        ratio = chi2_octh / chi2_nfw if chi2_nfw > 0 else 999
        winner = "OCTH" if ratio < 1.0 else "NFW"

        result = {
            'name': name,
            'n_points': len(data['r']),
            'r_max': float(data['r'].max()),
            'v_max': float(data['v_obs'].max()),
            'octh_constant': octh_const,
            'octh_transition': octh_trans,
            'octh_best': octh_best,
            'nfw': nfw,
            'ratio_octh_nfw': float(ratio),
            'winner': winner
        }

        psi_str = f"Psi={octh_const['psi']:.3f}" if octh_best['type'] == 'constant' else f"Psi_inf={octh_trans.get('psi_inf', 0):.3f}"
        print(f"chi2: OCTH={chi2_octh:.1f} NFW={chi2_nfw:.1f} -> {winner} ({psi_str})")

        return result

    except Exception as e:
        print(f"ERROR: {e}")
        return None


def main():
    print("=" * 70)
    print("SPARC - OCTH Modelo Simple (Psi constante o con transicion)")
    print("=" * 70)
    print(f"Fecha: {datetime.now()}")
    print()
    print("Modelo: v^2 = v_bar^2 / Psi")
    print("  - Psi < 1: masa efectiva mayor que barionica")
    print("  - Psi = 1: dinamica Newtoniana pura")
    print()

    files = list(DATA_DIR.glob("*_rotmod.dat"))
    print(f"Galaxias: {len(files)}")
    print()

    results = []
    octh_wins = 0
    nfw_wins = 0

    for f in files:
        r = analyze_galaxy(f)
        if r:
            results.append(r)
            if r['winner'] == 'OCTH':
                octh_wins += 1
            else:
                nfw_wins += 1

    # Estadisticas
    print()
    print("=" * 70)
    print("RESULTADOS")
    print("=" * 70)
    print(f"Analizadas: {len(results)}")
    print(f"OCTH gana: {octh_wins} ({100*octh_wins/len(results):.1f}%)")
    print(f"NFW gana:  {nfw_wins} ({100*nfw_wins/len(results):.1f}%)")
    print()

    # Chi2 promedio
    chi2_octh = np.mean([r['octh_best']['chi2_red'] for r in results])
    chi2_nfw = np.mean([r['nfw']['chi2_red'] for r in results])
    print(f"Chi2 promedio:")
    print(f"  OCTH: {chi2_octh:.2f}")
    print(f"  NFW:  {chi2_nfw:.2f}")
    print(f"  Ratio: {chi2_octh/chi2_nfw:.2f}")
    print()

    # Distribucion de Psi
    psi_values = [r['octh_constant']['psi'] for r in results]
    print(f"Distribucion de Psi (modelo constante):")
    print(f"  Media: {np.mean(psi_values):.3f}")
    print(f"  Std:   {np.std(psi_values):.3f}")
    print(f"  Min:   {np.min(psi_values):.3f}")
    print(f"  Max:   {np.max(psi_values):.3f}")
    print()

    # Donde OCTH es claramente mejor
    octh_clear = [r for r in results if r['ratio_octh_nfw'] < 0.8]
    print(f"Galaxias donde OCTH es >20% mejor que NFW: {len(octh_clear)}")
    for r in sorted(octh_clear, key=lambda x: x['ratio_octh_nfw'])[:10]:
        print(f"  {r['name']}: ratio={r['ratio_octh_nfw']:.2f}, Psi={r['octh_constant']['psi']:.3f}")

    # Guardar
    output = {
        'date': datetime.now().isoformat(),
        'n_galaxies': len(results),
        'summary': {
            'octh_wins': octh_wins,
            'nfw_wins': nfw_wins,
            'chi2_octh_mean': float(chi2_octh),
            'chi2_nfw_mean': float(chi2_nfw),
            'psi_mean': float(np.mean(psi_values)),
            'psi_std': float(np.std(psi_values))
        },
        'galaxies': results
    }

    out_path = RESULTS_DIR / "SPARC_OCTH_simple.json"
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nGuardado: {out_path}")

    return output


if __name__ == "__main__":
    main()
