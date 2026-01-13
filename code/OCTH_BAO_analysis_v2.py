#!/usr/bin/env python3
"""
OCTH BAO ANALYSIS v2.0 - Analisis corregido de Oscilaciones Acusticas Barionicas
=================================================================================

Este script realiza el analisis BAO corregido para OCTH.

CORRECCIONES vs v1.0:
1. Comparacion via theta_s (observable directo) en lugar de r_s
2. c_s NO se modifica en OCTH (solo H(z) cambia)
3. Uso correcto de distancias comoving

Author: F. Molina-Burgos
Date: 13 January 2026
Version: 2.0
"""

import numpy as np
import json
import os
from datetime import datetime

# Importar modulo core
from OCTH_cosmology_core import (
    PLANCK2018, OCTH_DEFAULT,
    CosmologyParams, OCTHParams,
    H_LCDM, H_OCTH, Psi,
    r_s, r_s_LCDM, r_s_OCTH,
    D_V, D_M, D_H, D_C,
    DV_rd, DM_rd, DH_rd,
    theta_s,
    z_drag_EH98, z_star_EH98,
    compute_all_LCDM, compute_all_OCTH
)

# =============================================================================
# DATOS OBSERVACIONALES BAO
# =============================================================================

# Mediciones BAO isotropicas (D_V/r_d)
BAO_ISOTROPIC = [
    {'name': '6dFGS', 'z': 0.106, 'DV_rd': 2.976, 'error': 0.133, 'ref': 'Beutler+2011'},
    {'name': 'SDSS MGS', 'z': 0.15, 'DV_rd': 4.466, 'error': 0.168, 'ref': 'Ross+2015'},
    {'name': 'BOSS DR12 lowz', 'z': 0.32, 'DV_rd': 8.467, 'error': 0.167, 'ref': 'Alam+2017'},
    {'name': 'BOSS DR12 CMASS', 'z': 0.57, 'DV_rd': 13.773, 'error': 0.134, 'ref': 'Alam+2017'},
]

# Mediciones BAO anisotropicas (D_M/r_d y D_H/r_d)
BAO_ANISOTROPIC = [
    {'name': 'BOSS DR12 z=0.38', 'z': 0.38, 'DM_rd': 10.27, 'DM_err': 0.15,
     'DH_rd': 25.00, 'DH_err': 0.76, 'ref': 'Alam+2017'},
    {'name': 'BOSS DR12 z=0.51', 'z': 0.51, 'DM_rd': 13.38, 'DM_err': 0.18,
     'DH_rd': 22.33, 'DH_err': 0.58, 'ref': 'Alam+2017'},
    {'name': 'BOSS DR12 z=0.61', 'z': 0.61, 'DM_rd': 15.45, 'DM_err': 0.20,
     'DH_rd': 20.75, 'DH_err': 0.60, 'ref': 'Alam+2017'},
    {'name': 'eBOSS DR16 LRG', 'z': 0.70, 'DM_rd': 17.65, 'DM_err': 0.30,
     'DH_rd': 19.78, 'DH_err': 0.46, 'ref': 'eBOSS+2020'},
    {'name': 'eBOSS DR16 QSO', 'z': 1.48, 'DM_rd': 30.21, 'DM_err': 0.79,
     'DH_rd': 13.23, 'DH_err': 0.47, 'ref': 'eBOSS+2020'},
    {'name': 'eBOSS DR16 Lya', 'z': 2.33, 'DM_rd': 37.41, 'DM_err': 1.86,
     'DH_rd': 8.93, 'DH_err': 0.28, 'ref': 'eBOSS+2020'},
]

# DESI Year 1 (2024)
DESI_Y1 = [
    {'name': 'DESI BGS', 'z': 0.295, 'DV_rd': 7.93, 'error': 0.15, 'ref': 'DESI+2024'},
    {'name': 'DESI LRG1', 'z': 0.510, 'DV_rd': 13.62, 'error': 0.25, 'ref': 'DESI+2024'},
    {'name': 'DESI LRG2', 'z': 0.706, 'DV_rd': 16.85, 'error': 0.32, 'ref': 'DESI+2024'},
    {'name': 'DESI LRG3+ELG1', 'z': 0.930, 'DV_rd': 21.71, 'error': 0.28, 'ref': 'DESI+2024'},
    {'name': 'DESI ELG2', 'z': 1.317, 'DV_rd': 27.79, 'error': 0.69, 'ref': 'DESI+2024'},
    {'name': 'DESI QSO', 'z': 1.491, 'DV_rd': 26.07, 'error': 0.67, 'ref': 'DESI+2024'},
]

# CMB (Planck 2018)
CMB_DATA = {
    'theta_s': 0.0104110,  # rad - observable directo
    'theta_s_error': 0.0000031,  # ~0.03%
    'r_d_derived': 147.09,  # Mpc - derivado asumiendo LCDM
    'r_d_error': 0.26,
}


# =============================================================================
# FUNCIONES DE ANALISIS
# =============================================================================

def analyze_bao_isotropic(bao_data: list,
                          H_func,
                          r_d: float,
                          model_name: str) -> dict:
    """
    Analiza datos BAO isotropicos.

    Args:
        bao_data: Lista de mediciones BAO
        H_func: Funcion H(z)
        r_d: Sound horizon en Mpc
        model_name: Nombre del modelo

    Returns:
        Diccionario con resultados
    """
    results = []
    chi2 = 0.0

    for bao in bao_data:
        z = bao['z']
        obs = bao['DV_rd']
        err = bao['error']

        # Prediccion del modelo
        pred = DV_rd(z, H_func, r_d)

        # Tension
        tension = (pred - obs) / err
        chi2 += tension**2

        results.append({
            'name': bao['name'],
            'z': z,
            'observed': obs,
            'predicted': pred,
            'error': err,
            'tension_sigma': tension,
            'ref': bao.get('ref', '')
        })

    return {
        'model': model_name,
        'r_d_used': r_d,
        'measurements': results,
        'chi2': chi2,
        'n_points': len(bao_data),
        'chi2_reduced': chi2 / len(bao_data) if bao_data else 0,
    }


def analyze_bao_anisotropic(bao_data: list,
                            H_func,
                            r_d: float,
                            model_name: str) -> dict:
    """
    Analiza datos BAO anisotropicos (D_M/r_d y D_H/r_d).

    Args:
        bao_data: Lista de mediciones
        H_func: Funcion H(z)
        r_d: Sound horizon
        model_name: Nombre del modelo

    Returns:
        Diccionario con resultados
    """
    results = []
    chi2 = 0.0

    for bao in bao_data:
        z = bao['z']

        # Predicciones
        DM_pred = DM_rd(z, H_func, r_d)
        DH_pred = DH_rd(z, H_func, r_d)

        # Tensiones
        tension_DM = (DM_pred - bao['DM_rd']) / bao['DM_err']
        tension_DH = (DH_pred - bao['DH_rd']) / bao['DH_err']

        chi2 += tension_DM**2 + tension_DH**2

        results.append({
            'name': bao['name'],
            'z': z,
            'DM_rd_obs': bao['DM_rd'],
            'DM_rd_pred': DM_pred,
            'DM_err': bao['DM_err'],
            'DM_tension': tension_DM,
            'DH_rd_obs': bao['DH_rd'],
            'DH_rd_pred': DH_pred,
            'DH_err': bao['DH_err'],
            'DH_tension': tension_DH,
        })

    n_points = len(bao_data) * 2  # 2 observables por punto

    return {
        'model': model_name,
        'r_d_used': r_d,
        'measurements': results,
        'chi2': chi2,
        'n_points': n_points,
        'chi2_reduced': chi2 / n_points if n_points else 0,
    }


def analyze_cmb_theta(H_func,
                      r_d: float,
                      model_name: str,
                      cosmo: CosmologyParams = PLANCK2018) -> dict:
    """
    Analiza el observable CMB theta_s.

    Este es el TEST CORRECTO para comparar con CMB, no r_s.

    Args:
        H_func: Funcion H(z)
        r_d: Sound horizon
        model_name: Nombre del modelo
        cosmo: Parametros cosmologicos

    Returns:
        Diccionario con resultados
    """
    # Calcular theta_s del modelo
    z_star = z_star_EH98(cosmo)
    D_C_star = D_C(z_star, H_func)
    theta_pred = r_d / D_C_star

    # Comparar con Planck
    theta_obs = CMB_DATA['theta_s']
    theta_err = CMB_DATA['theta_s_error']

    tension = (theta_pred - theta_obs) / theta_err

    return {
        'model': model_name,
        'z_star': z_star,
        'r_d': r_d,
        'D_C_star': D_C_star,
        'theta_s_predicted': theta_pred,
        'theta_s_observed': theta_obs,
        'theta_s_error': theta_err,
        'tension_sigma': tension,
        'theta_s_predicted_arcmin': np.degrees(theta_pred) * 60,
        'theta_s_observed_arcmin': np.degrees(theta_obs) * 60,
    }


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def run_full_analysis():
    """Ejecuta el analisis BAO completo para LCDM y OCTH."""

    print("=" * 80)
    print("OCTH BAO ANALYSIS v2.0")
    print("=" * 80)
    print(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M')}")

    # Calcular cantidades base
    print("\n" + "-" * 80)
    print("1. CANTIDADES BASE")
    print("-" * 80)

    cosmo = PLANCK2018
    octh = OCTH_DEFAULT

    z_d = z_drag_EH98(cosmo)
    z_s = z_star_EH98(cosmo)

    # LCDM
    H_lcdm = lambda z: H_LCDM(z, cosmo)
    r_d_lcdm = r_s_LCDM(cosmo)

    # OCTH
    H_octh = lambda z: H_OCTH(z, cosmo, octh)
    r_d_octh = r_s_OCTH(cosmo, octh)

    print(f"\n{'Cantidad':<25} {'LCDM':>15} {'OCTH':>15} {'Diff':>10}")
    print("-" * 70)
    print(f"{'z_drag':<25} {z_d:>15.2f} {z_d:>15.2f} {'-':>10}")
    print(f"{'z_star':<25} {z_s:>15.2f} {z_s:>15.2f} {'-':>10}")
    print(f"{'r_d (Mpc)':<25} {r_d_lcdm:>15.2f} {r_d_octh:>15.2f} {100*(r_d_octh/r_d_lcdm-1):>+9.1f}%")
    print(f"{'Psi(z_drag)':<25} {1.0:>15.4f} {Psi(z_d, octh):>15.4f} {'-':>10}")

    # ==========================================================================
    # CMB ANALYSIS (theta_s)
    # ==========================================================================
    print("\n" + "-" * 80)
    print("2. ANALISIS CMB (theta_s) - COMPARACION CORRECTA")
    print("-" * 80)

    cmb_lcdm = analyze_cmb_theta(H_lcdm, r_d_lcdm, 'LCDM', cosmo)
    cmb_octh = analyze_cmb_theta(H_octh, r_d_octh, 'OCTH', cosmo)

    print(f"\nObservable CMB: theta_s (angulo subtendido por horizonte acustico)")
    print(f"Planck observado: {cmb_lcdm['theta_s_observed']:.6f} rad = {cmb_lcdm['theta_s_observed_arcmin']:.2f} arcmin")
    print(f"\n{'Modelo':<10} {'theta_s (rad)':>15} {'arcmin':>10} {'Tension':>12}")
    print("-" * 50)
    print(f"{'LCDM':<10} {cmb_lcdm['theta_s_predicted']:>15.6f} {cmb_lcdm['theta_s_predicted_arcmin']:>10.2f} {cmb_lcdm['tension_sigma']:>+11.1f} sigma")
    print(f"{'OCTH':<10} {cmb_octh['theta_s_predicted']:>15.6f} {cmb_octh['theta_s_predicted_arcmin']:>10.2f} {cmb_octh['tension_sigma']:>+11.1f} sigma")

    print(f"\nNOTA: Las formulas analiticas de E&H98 no son exactas.")
    print(f"      Para precision necesaria usar CLASS/CAMB.")
    print(f"      La comparacion relativa LCDM vs OCTH es valida.")

    # ==========================================================================
    # BAO ISOTROPICO (BOSS/eBOSS tradicional)
    # ==========================================================================
    print("\n" + "-" * 80)
    print("3. BAO ISOTROPICO (D_V/r_d)")
    print("-" * 80)

    iso_lcdm = analyze_bao_isotropic(BAO_ISOTROPIC, H_lcdm, r_d_lcdm, 'LCDM')
    iso_octh = analyze_bao_isotropic(BAO_ISOTROPIC, H_octh, r_d_octh, 'OCTH')

    print(f"\n{'Survey':<20} {'z':>6} {'Obs':>8} {'LCDM':>8} {'OCTH':>8} {'T_L':>8} {'T_O':>8}")
    print("-" * 75)
    for i, bao in enumerate(BAO_ISOTROPIC):
        l = iso_lcdm['measurements'][i]
        o = iso_octh['measurements'][i]
        print(f"{bao['name']:<20} {bao['z']:>6.3f} {bao['DV_rd']:>8.2f} {l['predicted']:>8.2f} {o['predicted']:>8.2f} {l['tension_sigma']:>+7.1f}s {o['tension_sigma']:>+7.1f}s")

    print(f"\nchi^2/dof:  LCDM = {iso_lcdm['chi2']:.1f}/{iso_lcdm['n_points']} = {iso_lcdm['chi2_reduced']:.2f}")
    print(f"            OCTH = {iso_octh['chi2']:.1f}/{iso_octh['n_points']} = {iso_octh['chi2_reduced']:.2f}")

    # ==========================================================================
    # DESI Y1 (2024)
    # ==========================================================================
    print("\n" + "-" * 80)
    print("4. DESI Year 1 (2024)")
    print("-" * 80)

    desi_lcdm = analyze_bao_isotropic(DESI_Y1, H_lcdm, r_d_lcdm, 'LCDM')
    desi_octh = analyze_bao_isotropic(DESI_Y1, H_octh, r_d_octh, 'OCTH')

    print(f"\n{'Survey':<20} {'z':>6} {'Obs':>8} {'LCDM':>8} {'OCTH':>8} {'T_L':>8} {'T_O':>8}")
    print("-" * 75)
    for i, bao in enumerate(DESI_Y1):
        l = desi_lcdm['measurements'][i]
        o = desi_octh['measurements'][i]
        print(f"{bao['name']:<20} {bao['z']:>6.3f} {bao['DV_rd']:>8.2f} {l['predicted']:>8.2f} {o['predicted']:>8.2f} {l['tension_sigma']:>+7.1f}s {o['tension_sigma']:>+7.1f}s")

    print(f"\nchi^2/dof:  LCDM = {desi_lcdm['chi2']:.1f}/{desi_lcdm['n_points']} = {desi_lcdm['chi2_reduced']:.2f}")
    print(f"            OCTH = {desi_octh['chi2']:.1f}/{desi_octh['n_points']} = {desi_octh['chi2_reduced']:.2f}")

    # ==========================================================================
    # BAO ANISOTROPICO
    # ==========================================================================
    print("\n" + "-" * 80)
    print("5. BAO ANISOTROPICO (D_M/r_d y D_H/r_d)")
    print("-" * 80)

    aniso_lcdm = analyze_bao_anisotropic(BAO_ANISOTROPIC, H_lcdm, r_d_lcdm, 'LCDM')
    aniso_octh = analyze_bao_anisotropic(BAO_ANISOTROPIC, H_octh, r_d_octh, 'OCTH')

    print(f"\n{'Survey':<22} {'z':>5} {'DM obs':>8} {'DM pred':>8} {'T_DM':>7} {'DH obs':>8} {'DH pred':>8} {'T_DH':>7}")
    print("-" * 90)
    for i, result in enumerate(aniso_octh['measurements']):
        l = aniso_lcdm['measurements'][i]
        o = result
        print(f"{o['name']:<22} {o['z']:>5.2f} {o['DM_rd_obs']:>8.2f} {o['DM_rd_pred']:>8.2f} {o['DM_tension']:>+6.1f}s {o['DH_rd_obs']:>8.2f} {o['DH_rd_pred']:>8.2f} {o['DH_tension']:>+6.1f}s")

    print(f"\nchi^2/dof:  LCDM = {aniso_lcdm['chi2']:.1f}/{aniso_lcdm['n_points']} = {aniso_lcdm['chi2_reduced']:.2f}")
    print(f"            OCTH = {aniso_octh['chi2']:.1f}/{aniso_octh['n_points']} = {aniso_octh['chi2_reduced']:.2f}")

    # ==========================================================================
    # RESUMEN
    # ==========================================================================
    print("\n" + "=" * 80)
    print("RESUMEN DEL ANALISIS BAO")
    print("=" * 80)

    # Total chi2
    total_lcdm = iso_lcdm['chi2'] + desi_lcdm['chi2'] + aniso_lcdm['chi2']
    total_octh = iso_octh['chi2'] + desi_octh['chi2'] + aniso_octh['chi2']
    total_dof = iso_lcdm['n_points'] + desi_lcdm['n_points'] + aniso_lcdm['n_points']

    print(f"\nSound horizon at drag epoch:")
    print(f"  LCDM: r_d = {r_d_lcdm:.2f} Mpc")
    print(f"  OCTH: r_d = {r_d_octh:.2f} Mpc ({100*(r_d_octh/r_d_lcdm-1):+.1f}%)")

    print(f"\nChi-squared total (BAO isotropico + DESI + anisotropico):")
    print(f"  LCDM: chi^2 = {total_lcdm:.1f} / {total_dof} dof = {total_lcdm/total_dof:.2f}")
    print(f"  OCTH: chi^2 = {total_octh:.1f} / {total_dof} dof = {total_octh/total_dof:.2f}")

    print(f"\nInterpretacion:")
    if total_octh < total_lcdm:
        print(f"  OCTH tiene MEJOR chi^2 que LCDM (Delta chi^2 = {total_octh - total_lcdm:.1f})")
    else:
        print(f"  LCDM tiene mejor chi^2 (Delta chi^2 = {total_octh - total_lcdm:.1f})")

    print(f"""
NOTA IMPORTANTE:
----------------
Este analisis usa formulas analiticas simplificadas (Eisenstein & Hu 1998).
Para resultados de precision se requiere:
1. Implementar OCTH en CLASS o CAMB
2. Hacer MCMC completo de parametros cosmologicos
3. Incluir covarianzas entre mediciones BAO

La comparacion RELATIVA LCDM vs OCTH es indicativa, pero las tensiones
absolutas con datos pueden tener errores sistematicos de ~5%.
""")

    # ==========================================================================
    # GUARDAR RESULTADOS
    # ==========================================================================
    output = {
        'date': datetime.now().isoformat(),
        'version': '2.0',
        'parameters': {
            'cosmology': {
                'H_0': cosmo.H_0,
                'Omega_m': cosmo.Omega_m,
                'Omega_b': cosmo.Omega_b,
                'Omega_Lambda': cosmo.Omega_Lambda,
            },
            'octh': {
                'epsilon': octh.epsilon,
                'sigma': octh.sigma,
                'z_topo': octh.z_topo,
            }
        },
        'base_quantities': {
            'z_drag': z_d,
            'z_star': z_s,
            'r_d_lcdm': r_d_lcdm,
            'r_d_octh': r_d_octh,
            'r_d_ratio': r_d_octh / r_d_lcdm,
        },
        'cmb_analysis': {
            'lcdm': cmb_lcdm,
            'octh': cmb_octh,
        },
        'bao_isotropic': {
            'lcdm': iso_lcdm,
            'octh': iso_octh,
        },
        'desi_y1': {
            'lcdm': desi_lcdm,
            'octh': desi_octh,
        },
        'bao_anisotropic': {
            'lcdm': aniso_lcdm,
            'octh': aniso_octh,
        },
        'total_chi2': {
            'lcdm': total_lcdm,
            'octh': total_octh,
            'dof': total_dof,
            'delta_chi2': total_octh - total_lcdm,
        }
    }

    os.makedirs('../results', exist_ok=True)
    output_file = '../results/OCTH_BAO_analysis_v2.json'
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2, default=float)

    print(f"\nResultados guardados en: {output_file}")
    print("\n" + "=" * 80)

    return output


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    results = run_full_analysis()
