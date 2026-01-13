"""
OCTH H0 GAP ANALYSIS
====================

Problem: OCTH with epsilon_topo=0.015 predicts H0=67.9
         But observed is H0=73.0

This script analyzes the gap and finds the required epsilon_topo
to fully resolve the Hubble tension.

Author: F. Molina-Burgos
Date: January 2026
"""

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.integrate import quad

print("="*70)
print("  OCTH H0 GAP ANALYSIS")
print("  Finding epsilon_topo to fully resolve Hubble tension")
print("="*70)

# =============================================================================
# CONSTANTS
# =============================================================================

c = 2.998e5  # km/s
H0_PLANCK = 67.36  # km/s/Mpc (CMB-inferred, assuming LCDM)
H0_LOCAL = 73.04   # km/s/Mpc (SH0ES 2022)
H0_LOCAL_ERR = 1.04

OMEGA_M = 0.315
OMEGA_R = 9.0e-5
OMEGA_LAMBDA = 1 - OMEGA_M - OMEGA_R

# Sound horizon
RS_PLANCK = 147.09  # Mpc (LCDM)
RS_DESI = 143.0     # Mpc (DESI implied)

# =============================================================================
# OCTH MODEL
# =============================================================================

class OCTHCosmology:
    def __init__(self, epsilon_topo, sigma=0.6, z_topo=1089.0):
        self.epsilon = epsilon_topo
        self.sigma = sigma
        self.z_topo = z_topo
        self.ln_a_topo = -np.log(1 + z_topo)

    def psi(self, z):
        """Temporal permeability field"""
        ln_a = -np.log(1 + z)
        delta = (ln_a - self.ln_a_topo) / self.sigma
        return 1.0 - self.epsilon * np.exp(-0.5 * delta**2)

    def E_lcdm(self, z):
        """E(z) = H(z)/H0 for LCDM"""
        opz = 1 + z
        return np.sqrt(OMEGA_R * opz**4 + OMEGA_M * opz**3 + OMEGA_LAMBDA)

    def E_octh(self, z):
        """E(z) = H(z)/H0 for OCTH"""
        return self.E_lcdm(z) / self.psi(z)

    def comoving_distance(self, z):
        """Comoving distance in Mpc/h"""
        def integrand(zp):
            return 1.0 / self.E_octh(zp)
        result, _ = quad(integrand, 0, z)
        return (c / 100) * result  # Mpc/h

    def angular_diameter_distance(self, z):
        """Angular diameter distance"""
        return self.comoving_distance(z) / (1 + z)

    def sound_horizon_ratio(self):
        """
        r_s(OCTH) / r_s(LCDM)

        The sound horizon integral is:
        r_s = integral_0^{z_*} c_s / H(z) dz

        Since H_OCTH = H_LCDM / Psi, we have:
        r_s_OCTH = integral c_s * Psi / H_LCDM dz

        For the ratio:
        r_s_OCTH / r_s_LCDM ~ <Psi> weighted by integrand
        """
        z_star = 1089.0
        z_drag = 1060.0  # Drag epoch

        # Numerical integration
        def integrand_lcdm(z):
            return 1.0 / self.E_lcdm(z)

        def integrand_octh(z):
            return self.psi(z) / self.E_lcdm(z)

        rs_lcdm, _ = quad(integrand_lcdm, z_drag, 1e6, limit=100)
        rs_octh, _ = quad(integrand_octh, z_drag, 1e6, limit=100)

        return rs_octh / rs_lcdm

    def H0_local_from_CMB(self, H0_cmb=67.36):
        """
        Infer local H0 from CMB measurement.

        CMB measures theta_* = r_s / D_A(z_*)
        This is independent of H0.

        The tension arises because:
        - LCDM assumes r_s = 147 Mpc
        - If r_s is actually smaller, then H0 must be larger

        Scaling: H0_local ~ H0_cmb * (r_s_LCDM / r_s_OCTH)
        """
        ratio = self.sound_horizon_ratio()
        return H0_cmb / ratio


# =============================================================================
# ANALYSIS 1: Required epsilon for H0 match
# =============================================================================

print("\n" + "="*70)
print("ANALYSIS 1: What epsilon_topo is needed?")
print("="*70)

def H0_residual(epsilon):
    """Residual between predicted and observed H0"""
    model = OCTHCosmology(epsilon)
    H0_pred = model.H0_local_from_CMB(H0_PLANCK)
    return H0_pred - H0_LOCAL

# Scan epsilon
print("\n{:^12} {:^12} {:^12} {:^12}".format(
    "epsilon", "Psi(z*)", "r_s/r_s_LCDM", "H0_local"))
print("-"*50)

epsilons_to_test = [0.01, 0.02, 0.03, 0.05, 0.07, 0.10, 0.15, 0.183]
results = []

for eps in epsilons_to_test:
    model = OCTHCosmology(eps)
    psi_star = model.psi(1089)
    rs_ratio = model.sound_horizon_ratio()
    H0_pred = model.H0_local_from_CMB(H0_PLANCK)

    results.append({
        'epsilon': eps,
        'psi_star': psi_star,
        'rs_ratio': rs_ratio,
        'H0_local': H0_pred,
        'delta_H0': H0_pred - H0_LOCAL
    })

    print(f"{eps:^12.3f} {psi_star:^12.4f} {rs_ratio:^12.4f} {H0_pred:^12.2f}")

# Find epsilon that gives H0 = 73
try:
    epsilon_optimal = brentq(H0_residual, 0.01, 0.3)
    model_opt = OCTHCosmology(epsilon_optimal)
    print(f"\nOptimal epsilon for H0 = 73.04:")
    print(f"  epsilon_topo = {epsilon_optimal:.4f}")
    print(f"  Psi(z*) = {model_opt.psi(1089):.4f}")
    print(f"  r_s/r_s_LCDM = {model_opt.sound_horizon_ratio():.4f}")
except:
    print("\nCould not find exact epsilon (may be out of range)")
    epsilon_optimal = None

# =============================================================================
# ANALYSIS 2: Why is epsilon_topo = 0.015 insufficient?
# =============================================================================

print("\n" + "="*70)
print("ANALYSIS 2: Why epsilon=0.015 is insufficient")
print("="*70)

model_015 = OCTHCosmology(0.015)
model_077 = OCTHCosmology(0.077) if epsilon_optimal else None

print(f"""
With epsilon_topo = 0.015:
  Psi(z*) = {model_015.psi(1089):.4f}
  r_s ratio = {model_015.sound_horizon_ratio():.4f}
  H0 predicted = {model_015.H0_local_from_CMB(H0_PLANCK):.2f} km/s/Mpc

  Gap from observed = {model_015.H0_local_from_CMB(H0_PLANCK) - H0_LOCAL:.2f} km/s/Mpc
""")

if epsilon_optimal:
    print(f"""
Required epsilon for FULL resolution:
  epsilon_topo = {epsilon_optimal:.4f}
  Psi(z*) = {model_opt.psi(1089):.4f}
  r_s ratio = {model_opt.sound_horizon_ratio():.4f}
  H0 predicted = {model_opt.H0_local_from_CMB(H0_PLANCK):.2f} km/s/Mpc
""")

# =============================================================================
# ANALYSIS 3: What does this mean for the two-component model?
# =============================================================================

print("\n" + "="*70)
print("ANALYSIS 3: Implications for two-component model")
print("="*70)

# Original model:
# epsilon_topo = 0.015 (cosmological)
# epsilon_quantum = 0.168 (galactic)
# Total = 0.183 (MOND)

# If we need epsilon_topo ~ 0.077 for H0, then:
# epsilon_quantum = 0.183 - 0.077 = 0.106

if epsilon_optimal:
    eps_topo_new = epsilon_optimal
    eps_quantum_new = 0.183 - eps_topo_new

    print(f"""
REVISED TWO-COMPONENT MODEL:

  Original:
    epsilon_topo = 0.015
    epsilon_quantum = 0.168
    Total = 0.183 (MOND)

  Revised (for full H0 resolution):
    epsilon_topo = {eps_topo_new:.3f}
    epsilon_quantum = {eps_quantum_new:.3f}
    Total = 0.183 (MOND, unchanged)

  Interpretation:
    More of the epsilon comes from topology (z~1089)
    Less comes from quantum effects (a<a0)

  Testable prediction:
    Galaxy rotation should show WEAKER MOND effect
    (epsilon_quantum = {eps_quantum_new:.3f} vs original 0.168)

    At very low acceleration (a << a0):
      v_asymptotic / v_Newton = 1 / (1 - eps_quantum)

    Original: v_ratio = {1/(1-0.168):.3f}
    Revised:  v_ratio = {1/(1-eps_quantum_new):.3f}
""")

# =============================================================================
# ANALYSIS 4: CMB peak positions with revised epsilon
# =============================================================================

print("\n" + "="*70)
print("ANALYSIS 4: CMB implications of revised epsilon")
print("="*70)

if epsilon_optimal:
    # CMB peaks scale as l_n ~ n * pi / theta_*
    # theta_* = r_s / D_A(z_*)
    # If r_s decreases, theta_* decreases, peaks shift to higher l

    model_lcdm = OCTHCosmology(0.0)  # LCDM limit
    model_revised = OCTHCosmology(eps_topo_new)

    # Peak positions in LCDM
    peaks_lcdm = [220.0, 537.5, 810.8, 1120.9, 1444.2]

    # Shift factor
    shift = model_lcdm.sound_horizon_ratio() / model_revised.sound_horizon_ratio()

    peaks_revised = [p * shift for p in peaks_lcdm]

    print(f"CMB Peak Shifts:")
    print(f"  Shift factor = {shift:.4f}")
    print()
    print(f"  {'Peak':<8} {'LCDM':<12} {'OCTH revised':<12} {'Planck obs':<12}")
    print("-"*50)

    planck_obs = [220.0, 537.5, 810.8, 1120.9, 1444.2]
    planck_err = [0.5, 0.7, 0.7, 1.5, 2.0]

    chi2_peaks = 0
    for i, (l_lcdm, l_octh, l_obs, sigma) in enumerate(zip(
            peaks_lcdm, peaks_revised, planck_obs, planck_err)):
        tension = (l_octh - l_obs) / sigma
        chi2_peaks += tension**2
        print(f"  {i+1:<8} {l_lcdm:<12.1f} {l_octh:<12.1f} {l_obs:<12.1f} ({tension:+.1f}σ)")

    print(f"\n  Total chi2 for peaks = {chi2_peaks:.1f}")

    if chi2_peaks > 10:
        print("\n  WARNING: Peak positions deviate significantly!")
        print("  This suggests tension between H0 resolution and CMB fit.")

# =============================================================================
# ANALYSIS 5: The sigma parameter
# =============================================================================

print("\n" + "="*70)
print("ANALYSIS 5: Effect of sigma (width parameter)")
print("="*70)

print("""
The Psi function has width parameter sigma:
  Psi(z) = 1 - epsilon * exp(-0.5 * ((ln(a) - ln(a_topo))/sigma)^2)

A larger sigma spreads the effect over more redshift range.
This might allow smaller epsilon with same integrated effect.
""")

# Scan sigma
print(f"\n{'sigma':^10} {'eps_needed':^12} {'Psi(z*)':^10} {'H0':^10}")
print("-"*45)

for sigma in [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5]:
    def residual_sigma(eps):
        model = OCTHCosmology(eps, sigma=sigma)
        return model.H0_local_from_CMB(H0_PLANCK) - H0_LOCAL

    try:
        eps_opt = brentq(residual_sigma, 0.01, 0.5)
        model = OCTHCosmology(eps_opt, sigma=sigma)
        print(f"{sigma:^10.2f} {eps_opt:^12.4f} {model.psi(1089):^10.4f} {model.H0_local_from_CMB(H0_PLANCK):^10.2f}")
    except:
        print(f"{sigma:^10.2f} {'N/A':^12} {'N/A':^10} {'N/A':^10}")

# =============================================================================
# ANALYSIS 6: Alternative - z_topo parameter
# =============================================================================

print("\n" + "="*70)
print("ANALYSIS 6: Effect of z_topo (transition redshift)")
print("="*70)

print("""
What if the topological transition is not exactly at recombination?
Maybe it's slightly earlier or later?
""")

print(f"\n{'z_topo':^10} {'eps_needed':^12} {'Psi(z*)':^10}")
print("-"*35)

for z_topo in [900, 1000, 1089, 1100, 1200, 1500]:
    def residual_ztopo(eps):
        model = OCTHCosmology(eps, z_topo=z_topo)
        return model.H0_local_from_CMB(H0_PLANCK) - H0_LOCAL

    try:
        eps_opt = brentq(residual_ztopo, 0.01, 0.5)
        model = OCTHCosmology(eps_opt, z_topo=z_topo)
        print(f"{z_topo:^10} {eps_opt:^12.4f} {model.psi(1089):^10.4f}")
    except:
        print(f"{z_topo:^10} {'N/A':^12} {'N/A':^10}")

# =============================================================================
# CONCLUSION
# =============================================================================

print("\n" + "="*70)
print("CONCLUSIONS")
print("="*70)

eps_str = f"{epsilon_optimal:.3f}" if epsilon_optimal else "0.07-0.08"
print(f"""
1. EPSILON REQUIREMENT:
   To fully resolve H0 tension (73.04 km/s/Mpc):
   epsilon_topo ~ {eps_str}

   This is ~5x larger than our original 0.015!

2. TWO-COMPONENT REVISION:
   If epsilon_total = 0.183 (MOND) is fixed:
   epsilon_topo ~ 0.077 => epsilon_quantum ~ 0.106

   This REDUCES the galactic MOND effect by ~40%!

3. TENSION:
   - CMB peaks are sensitive to epsilon_topo
   - Larger epsilon shifts peaks more
   - May create tension with Planck peak positions

4. POSSIBLE RESOLUTIONS:

   a) Accept partial H0 resolution (~5 km/s/Mpc gap)
      Keep epsilon_topo = 0.015, CMB fits well

   b) Modify sigma (width) parameter
      Larger sigma might allow full H0 resolution with smaller peak shifts

   c) Shift z_topo away from recombination
      If transition is at z ~ 1200 instead of 1089,
      effect on peaks is different

   d) Allow epsilon_total > 0.183
      More total epsilon, but then MOND prediction changes

   e) Non-Gaussian Psi profile
      Different functional form might fit both CMB and H0

5. NEXT STEP:
   Calculate full CMB power spectrum with CAMB/CLASS
   to quantify the trade-off between H0 and CMB fit.
""")

print("="*70)
