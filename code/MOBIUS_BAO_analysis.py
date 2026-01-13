#!/usr/bin/env python3
"""
MOBIUS GRAVITY: Baryon Acoustic Oscillations (BAO) Analysis

This script computes the BAO scale (sound horizon at drag epoch) in MOBIUS
and compares with observations from SDSS, BOSS, and DESI.

The key modification: during the drag epoch, Psi affects the baryon-photon
coupling and thus the sound horizon.

Author: F. Molina-Burgos
Date: 12 January 2026
"""

import numpy as np
from scipy.integrate import quad, odeint
from scipy.interpolate import interp1d
import json
import os

print("=" * 70)
print("MOBIUS GRAVITY: BAO ANALYSIS")
print("=" * 70)

# =============================================================================
# COSMOLOGICAL PARAMETERS (Planck 2018)
# =============================================================================

c = 2.998e5           # km/s (for cosmological calculations)
c_mps = 2.998e8       # m/s (for MOBIUS scales)
H_0 = 67.4            # km/s/Mpc
h = H_0 / 100

Omega_b = 0.0493      # Baryon density
Omega_c = 0.264       # CDM density (in standard model)
Omega_m = 0.315       # Total matter
Omega_r = 9.0e-5      # Radiation
Omega_Lambda = 0.685  # Dark energy

# CMB temperature
T_CMB = 2.7255        # K

# MOBIUS parameters
a_0 = 1.2e-10         # m/s^2
# H_0 in SI: 67.4 km/s/Mpc = 67.4 / (3.086e19 km/Mpc) s^-1 = 2.18e-18 s^-1
H0_SI = H_0 / 3.086e19  # s^-1
a_H = c_mps * H0_SI     # m/s^2 = 6.55e-10
epsilon = a_0 / a_H     # = 0.183
sigma = 0.6
z_topo = 1089         # Peak of topological effect

print(f"\nCosmological parameters:")
print(f"  H_0 = {H_0} km/s/Mpc")
print(f"  Omega_b = {Omega_b}")
print(f"  Omega_m = {Omega_m}")
print(f"  Omega_Lambda = {Omega_Lambda}")
print(f"\nMOBIUS parameters:")
print(f"  epsilon = {epsilon:.4f}")
print(f"  sigma = {sigma}")

# =============================================================================
# PSI FUNCTION
# =============================================================================

def Psi_MOBIUS(z):
    """
    MOBIUS Psi function - topological modification to gravity.
    Psi < 1 near recombination enhances gravitational effects.
    """
    ln_a = -np.log(1 + z)
    ln_a_topo = -np.log(1 + z_topo)

    # Gaussian profile centered at recombination
    f_topo = np.exp(-0.5 * ((ln_a - ln_a_topo) / sigma)**2)

    return 1.0 - epsilon * f_topo

# =============================================================================
# HUBBLE PARAMETER
# =============================================================================

def E_squared(z, use_mobius=True):
    """
    E(z)^2 = H(z)^2 / H_0^2
    In MOBIUS, matter density is effectively enhanced by 1/Psi^2
    """
    a = 1.0 / (1 + z)

    if use_mobius:
        Psi = Psi_MOBIUS(z)
        # Effective matter density enhanced
        Omega_m_eff = Omega_m / Psi**2
    else:
        Omega_m_eff = Omega_m

    return Omega_r * (1+z)**4 + Omega_m_eff * (1+z)**3 + Omega_Lambda

def H(z, use_mobius=True):
    """Hubble parameter in km/s/Mpc"""
    return H_0 * np.sqrt(E_squared(z, use_mobius))

# =============================================================================
# SOUND SPEED AND DRAG EPOCH
# =============================================================================

def R_baryon(z):
    """
    Baryon-to-photon momentum density ratio
    R = 3 * rho_b / (4 * rho_gamma)
    """
    # rho_b ~ Omega_b * (1+z)^3
    # rho_gamma ~ Omega_gamma * (1+z)^4
    # R ~ (3/4) * (Omega_b / Omega_gamma) * (1+z)^-1

    Omega_gamma = 2.47e-5 / h**2  # Photon density parameter
    return 0.75 * (Omega_b / Omega_gamma) / (1 + z)

def c_s(z, use_mobius=True):
    """
    Sound speed in baryon-photon fluid
    c_s = c / sqrt(3 * (1 + R))

    In MOBIUS, the effective baryon loading is modified.
    """
    R = R_baryon(z)

    if use_mobius:
        Psi = Psi_MOBIUS(z)
        # Deeper potential wells -> more baryon loading
        # Effective R is enhanced by 1/Psi^2
        R_eff = R / Psi**2
    else:
        R_eff = R

    return c / np.sqrt(3.0 * (1.0 + R_eff))

def z_drag_fitting(Omega_b_h2, Omega_m_h2):
    """
    Fitting formula for drag epoch (Eisenstein & Hu 1998)
    """
    b1 = 0.313 * Omega_m_h2**(-0.419) * (1 + 0.607 * Omega_m_h2**0.674)
    b2 = 0.238 * Omega_m_h2**0.223
    z_d = 1291 * Omega_m_h2**0.251 / (1 + 0.659 * Omega_m_h2**0.828) * (1 + b1 * Omega_b_h2**b2)
    return z_d

# =============================================================================
# SOUND HORIZON CALCULATION
# =============================================================================

def sound_horizon(z_final, use_mobius=True, n_points=1000):
    """
    Calculate sound horizon r_s at redshift z_final

    r_s = integral_z^inf c_s(z') / H(z') dz'
    """
    z_max = 1e6  # Start from very early times

    def integrand(z):
        return c_s(z, use_mobius) / H(z, use_mobius)

    # Numerical integration
    result, error = quad(integrand, z_final, z_max, limit=200)

    return result  # in Mpc

print("\n" + "=" * 70)
print("SOUND HORIZON CALCULATIONS")
print("=" * 70)

# Drag epoch
Omega_b_h2 = Omega_b * h**2
Omega_m_h2 = Omega_m * h**2
z_drag = z_drag_fitting(Omega_b_h2, Omega_m_h2)

print(f"\nDrag epoch:")
print(f"  z_drag = {z_drag:.1f}")

# Sound horizon at drag epoch
r_s_standard = sound_horizon(z_drag, use_mobius=False)
r_s_mobius = sound_horizon(z_drag, use_mobius=True)

print(f"\nSound horizon at drag epoch:")
print(f"  r_s (standard) = {r_s_standard:.2f} Mpc")
print(f"  r_s (MOBIUS)   = {r_s_mobius:.2f} Mpc")
print(f"  Ratio = {r_s_mobius/r_s_standard:.4f}")

# =============================================================================
# COMPARISON WITH OBSERVATIONS
# =============================================================================

print("\n" + "=" * 70)
print("COMPARISON WITH BAO OBSERVATIONS")
print("=" * 70)

# Observed BAO measurements (r_d values)
observations = {
    "Planck 2018 (CMB)": {"r_d": 147.09, "error": 0.26, "z_eff": 1089},
    "SDSS DR12 (2017)": {"r_d": 147.78, "error": 0.74, "z_eff": 0.51},
    "BOSS DR14 (2018)": {"r_d": 147.33, "error": 0.49, "z_eff": 0.38},
    "eBOSS DR16 (2020)": {"r_d": 147.21, "error": 0.48, "z_eff": 0.70},
    "DESI 2024": {"r_d": 142.1, "error": 2.0, "z_eff": 0.51},  # Preliminary
}

print(f"\n{'Survey':<25} {'r_d obs':<12} {'r_d MOBIUS':<12} {'Tension':>10}")
print("-" * 60)

results = {}
for name, obs in observations.items():
    r_d_obs = obs["r_d"]
    r_d_err = obs["error"]

    # For CMB, use z_drag; for low-z, r_d is measured via D_V/r_d
    if obs["z_eff"] > 100:
        r_d_pred = r_s_mobius
    else:
        # At low-z, MOBIUS predicts slightly different expansion
        # But r_d is anchored at drag epoch
        r_d_pred = r_s_mobius

    tension = (r_d_pred - r_d_obs) / r_d_err

    results[name] = {
        "r_d_observed": r_d_obs,
        "r_d_predicted": r_d_pred,
        "tension_sigma": tension
    }

    print(f"{name:<25} {r_d_obs:.2f}+-{r_d_err:.2f}  {r_d_pred:.2f} Mpc   {tension:+.2f} sigma")

# =============================================================================
# D_V / r_d MEASUREMENTS
# =============================================================================

print("\n" + "=" * 70)
print("D_V/r_d MEASUREMENTS")
print("=" * 70)

def D_V(z, use_mobius=True):
    """
    Volume-averaged distance
    D_V(z) = [c*z/H(z) * D_A(z)^2]^(1/3)
    """
    # Comoving distance
    def integrand(z_prime):
        return c / H(z_prime, use_mobius)

    D_C, _ = quad(integrand, 0, z, limit=100)
    D_A = D_C / (1 + z)  # Angular diameter distance

    # Volume distance
    return (c * z / H(z, use_mobius) * D_A**2)**(1.0/3.0)

# BAO measurements as D_V/r_d
bao_measurements = [
    {"z": 0.15, "DV_rd": 4.47, "error": 0.17, "survey": "6dFGS"},
    {"z": 0.32, "DV_rd": 8.47, "error": 0.17, "survey": "BOSS lowz"},
    {"z": 0.51, "DV_rd": 13.38, "error": 0.18, "survey": "BOSS CMASS"},
    {"z": 0.70, "DV_rd": 17.86, "error": 0.33, "survey": "eBOSS LRG"},
    {"z": 1.48, "DV_rd": 30.21, "error": 0.79, "survey": "eBOSS QSO"},
    {"z": 2.33, "DV_rd": 37.6, "error": 1.9, "survey": "eBOSS Lya"},
]

print(f"\n{'Survey':<15} {'z':>6} {'(D_V/r_d)_obs':>14} {'(D_V/r_d)_MOB':>14} {'Tension':>10}")
print("-" * 65)

chi2_total = 0
dof = 0

for m in bao_measurements:
    z = m["z"]
    DV_rd_obs = m["DV_rd"]
    err = m["error"]

    # MOBIUS prediction
    D_V_mobius = D_V(z, use_mobius=True)
    DV_rd_mobius = D_V_mobius / r_s_mobius

    tension = (DV_rd_mobius - DV_rd_obs) / err
    chi2_total += tension**2
    dof += 1

    results[m["survey"]] = {
        "z": z,
        "DV_rd_observed": DV_rd_obs,
        "DV_rd_predicted": DV_rd_mobius,
        "tension_sigma": tension
    }

    print(f"{m['survey']:<15} {z:>6.2f} {DV_rd_obs:>10.2f}+-{err:.2f} {DV_rd_mobius:>14.2f} {tension:>+10.2f} sigma")

print(f"\nTotal chi^2 = {chi2_total:.2f} for {dof} data points")
print(f"Reduced chi^2 = {chi2_total/dof:.2f}")

# =============================================================================
# D_A AND D_H MEASUREMENTS (ANISOTROPIC BAO)
# =============================================================================

print("\n" + "=" * 70)
print("ANISOTROPIC BAO: D_A/r_d AND D_H/r_d")
print("=" * 70)

def D_A_comoving(z, use_mobius=True):
    """Comoving angular diameter distance"""
    def integrand(z_prime):
        return c / H(z_prime, use_mobius)
    D_C, _ = quad(integrand, 0, z, limit=100)
    return D_C

def D_H(z, use_mobius=True):
    """Hubble distance D_H = c/H(z)"""
    return c / H(z, use_mobius)

# BOSS/eBOSS anisotropic measurements
aniso_bao = [
    {"z": 0.38, "DA_rd": 10.27, "DA_err": 0.15, "DH_rd": 25.00, "DH_err": 0.76, "survey": "BOSS lowz"},
    {"z": 0.51, "DA_rd": 13.38, "DA_err": 0.17, "DH_rd": 22.33, "DH_err": 0.58, "survey": "BOSS CMASS"},
    {"z": 0.70, "DA_rd": 17.65, "DA_err": 0.30, "DH_rd": 19.78, "DH_err": 0.46, "survey": "eBOSS LRG"},
]

print(f"\n{'Survey':<12} {'z':>5} {'DA/rd obs':>10} {'DA/rd MOB':>10} {'DH/rd obs':>10} {'DH/rd MOB':>10}")
print("-" * 70)

chi2_aniso = 0
for m in aniso_bao:
    z = m["z"]

    # MOBIUS predictions
    DA_mobius = D_A_comoving(z, use_mobius=True) / r_s_mobius
    DH_mobius = D_H(z, use_mobius=True) / r_s_mobius

    tension_DA = (DA_mobius - m["DA_rd"]) / m["DA_err"]
    tension_DH = (DH_mobius - m["DH_rd"]) / m["DH_err"]

    chi2_aniso += tension_DA**2 + tension_DH**2

    print(f"{m['survey']:<12} {z:>5.2f} {m['DA_rd']:>10.2f} {DA_mobius:>10.2f} {m['DH_rd']:>10.2f} {DH_mobius:>10.2f}")

print(f"\nAnisotropic chi^2 = {chi2_aniso:.2f}")

# =============================================================================
# DESI 2024 SPECIFIC ANALYSIS
# =============================================================================

print("\n" + "=" * 70)
print("DESI 2024 PRELIMINARY RESULTS")
print("=" * 70)

# DESI found tension with Planck - let's see if MOBIUS helps
print("""
DESI 2024 Key Finding:
- D_V(z=0.51)/r_d = 13.62 +/- 0.25
- This implies r_d ~ 137-142 Mpc (if using Planck D_V)
- Tension with Planck r_d = 147.09 Mpc

MOBIUS Analysis:
""")

DESI_DV_rd = 13.62
DESI_err = 0.25

# Standard cosmology prediction
D_V_standard = D_V(0.51, use_mobius=False)
D_V_mobius_051 = D_V(0.51, use_mobius=True)

print(f"At z = 0.51:")
print(f"  D_V (standard) = {D_V_standard:.2f} Mpc")
print(f"  D_V (MOBIUS)   = {D_V_mobius_051:.2f} Mpc")
print(f"\n  r_s (standard) = {r_s_standard:.2f} Mpc")
print(f"  r_s (MOBIUS)   = {r_s_mobius:.2f} Mpc")

DV_rd_standard = D_V_standard / r_s_standard
DV_rd_mobius = D_V_mobius_051 / r_s_mobius

print(f"\n  D_V/r_d (standard) = {DV_rd_standard:.2f}")
print(f"  D_V/r_d (MOBIUS)   = {DV_rd_mobius:.2f}")
print(f"  D_V/r_d (DESI)     = {DESI_DV_rd:.2f} +/- {DESI_err}")

tension_standard = (DV_rd_standard - DESI_DV_rd) / DESI_err
tension_mobius = (DV_rd_mobius - DESI_DV_rd) / DESI_err

print(f"\n  Tension with DESI:")
print(f"    Standard: {tension_standard:+.2f} sigma")
print(f"    MOBIUS:   {tension_mobius:+.2f} sigma")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("BAO ANALYSIS SUMMARY")
print("=" * 70)

summary = {
    "r_s_standard_Mpc": r_s_standard,
    "r_s_mobius_Mpc": r_s_mobius,
    "z_drag": z_drag,
    "chi2_isotropic": chi2_total,
    "chi2_anisotropic": chi2_aniso,
    "dof": dof,
    "tension_with_DESI_sigma": tension_mobius
}

print(f"""
RESULTS:
--------
Sound horizon at drag epoch:
  Standard: r_s = {r_s_standard:.2f} Mpc
  MOBIUS:   r_s = {r_s_mobius:.2f} Mpc
  Ratio:    {r_s_mobius/r_s_standard:.4f}

Chi-squared (isotropic BAO):
  chi^2 = {chi2_total:.2f} / {dof} dof
  Reduced chi^2 = {chi2_total/dof:.2f}

DESI Tension:
  Standard: {tension_standard:+.2f} sigma
  MOBIUS:   {tension_mobius:+.2f} sigma

INTERPRETATION:
---------------
MOBIUS predicts a SMALLER sound horizon due to enhanced baryon loading
from deeper potential wells (1/Psi^2 effect).

This goes in the SAME DIRECTION as the DESI discrepancy with Planck!

The effect is modest ({100*(1-r_s_mobius/r_s_standard):.1f}% reduction) but
could help explain part of the DESI tension.
""")

# =============================================================================
# SAVE RESULTS
# =============================================================================

os.makedirs('../results', exist_ok=True)
output_data = {
    "summary": summary,
    "detailed_results": results,
    "parameters": {
        "epsilon": epsilon,
        "sigma": sigma,
        "z_topo": z_topo,
        "H_0": H_0,
        "Omega_m": Omega_m,
        "Omega_b": Omega_b
    }
}

with open('../results/MOBIUS_BAO_analysis.json', 'w') as f:
    json.dump(output_data, f, indent=2)

print(f"Results saved to: ../results/MOBIUS_BAO_analysis.json")

print("\n" + "=" * 70)
print("BAO ANALYSIS COMPLETE")
print("=" * 70)
