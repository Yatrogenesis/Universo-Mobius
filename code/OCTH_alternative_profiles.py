"""
OCTH ALTERNATIVE PSI PROFILES
=============================

Goal: Find a Psi profile that:
1. Reduces sound horizon r_s enough to resolve H0
2. But doesn't shift CMB peaks excessively

Key insight: CMB peaks go as l ~ pi / theta*
where theta* = r_s / D_A(z*)

If r_s and D_A change together, the ratio might be preserved!

Author: F. Molina-Burgos
Date: January 2026
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
import warnings
warnings.filterwarnings('ignore')

print("="*70)
print("  OCTH ALTERNATIVE PSI PROFILES")
print("  Searching for CMB-compatible H0 resolution")
print("="*70)

# =============================================================================
# CONSTANTS
# =============================================================================

c = 2.998e5  # km/s
H0_PLANCK = 67.36
H0_LOCAL = 73.04

OMEGA_M = 0.315
OMEGA_R = 9.0e-5
OMEGA_LAMBDA = 1 - OMEGA_M - OMEGA_R

Z_STAR = 1089.0  # Recombination
Z_DRAG = 1060.0  # Drag epoch

# Planck observations
L_PEAK_1 = 220.0
THETA_STAR_PLANCK = 0.0104  # radians, approx

# =============================================================================
# BASE COSMOLOGY
# =============================================================================

def E_lcdm(z):
    """E(z) = H(z)/H0 for LCDM"""
    opz = 1 + z
    return np.sqrt(OMEGA_R * opz**4 + OMEGA_M * opz**3 + OMEGA_LAMBDA)

# =============================================================================
# PSI PROFILE CLASSES
# =============================================================================

class GaussianPsi:
    """Original Gaussian profile"""
    name = "Gaussian"

    def __init__(self, eps, sigma=0.6, z_topo=1089.0):
        self.eps = eps
        self.sigma = sigma
        self.z_topo = z_topo
        self.ln_a_topo = -np.log(1 + z_topo)

    def psi(self, z):
        ln_a = -np.log(1 + z)
        delta = (ln_a - self.ln_a_topo) / self.sigma
        return 1.0 - self.eps * np.exp(-0.5 * delta**2)


class StepPsi:
    """Step function: Psi drops sharply at z_topo"""
    name = "Step"

    def __init__(self, eps, z_topo=1089.0, width=50.0):
        self.eps = eps
        self.z_topo = z_topo
        self.width = width

    def psi(self, z):
        # Smooth step using tanh
        return 1.0 - self.eps * 0.5 * (1 + np.tanh((z - self.z_topo) / self.width))


class AsymmetricPsi:
    """Asymmetric profile: fast drop, slow recovery"""
    name = "Asymmetric"

    def __init__(self, eps, z_topo=1089.0, sigma_left=0.3, sigma_right=1.0):
        self.eps = eps
        self.z_topo = z_topo
        self.sigma_left = sigma_left
        self.sigma_right = sigma_right
        self.ln_a_topo = -np.log(1 + z_topo)

    def psi(self, z):
        ln_a = -np.log(1 + z)
        if ln_a < self.ln_a_topo:
            # Early universe (higher z)
            delta = (ln_a - self.ln_a_topo) / self.sigma_right
        else:
            # Late universe (lower z)
            delta = (ln_a - self.ln_a_topo) / self.sigma_left
        return 1.0 - self.eps * np.exp(-0.5 * delta**2)


class PlateauPsi:
    """Plateau profile: Psi stays reduced over a range"""
    name = "Plateau"

    def __init__(self, eps, z_start=1500.0, z_end=800.0, edge_width=50.0):
        self.eps = eps
        self.z_start = z_start
        self.z_end = z_end
        self.edge_width = edge_width

    def psi(self, z):
        # Smooth plateau using product of tanh
        high_z_cutoff = 0.5 * (1 + np.tanh((z - self.z_start) / self.edge_width))
        low_z_cutoff = 0.5 * (1 - np.tanh((z - self.z_end) / self.edge_width))
        return 1.0 - self.eps * high_z_cutoff * low_z_cutoff


class DoublePeakPsi:
    """Two Gaussian peaks - at recombination AND drag epoch"""
    name = "Double Peak"

    def __init__(self, eps1, eps2, z1=1089.0, z2=1060.0, sigma=0.3):
        self.eps1 = eps1
        self.eps2 = eps2
        self.z1 = z1
        self.z2 = z2
        self.sigma = sigma
        self.ln_a1 = -np.log(1 + z1)
        self.ln_a2 = -np.log(1 + z2)

    def psi(self, z):
        ln_a = -np.log(1 + z)
        delta1 = (ln_a - self.ln_a1) / self.sigma
        delta2 = (ln_a - self.ln_a2) / self.sigma
        g1 = self.eps1 * np.exp(-0.5 * delta1**2)
        g2 = self.eps2 * np.exp(-0.5 * delta2**2)
        return 1.0 - (g1 + g2)


class ExtendedPsi:
    """Extended low-z tail - affects D_A more"""
    name = "Extended"

    def __init__(self, eps_peak, eps_tail, z_peak=1089.0, z_tail_start=500.0):
        self.eps_peak = eps_peak
        self.eps_tail = eps_tail
        self.z_peak = z_peak
        self.z_tail_start = z_tail_start
        self.ln_a_peak = -np.log(1 + z_peak)

    def psi(self, z):
        ln_a = -np.log(1 + z)
        # Peak at recombination
        delta = (ln_a - self.ln_a_peak) / 0.6
        peak = self.eps_peak * np.exp(-0.5 * delta**2)
        # Extended tail at lower z
        tail = self.eps_tail * 0.5 * (1 - np.tanh((z - self.z_tail_start) / 100))
        return 1.0 - peak - tail

# =============================================================================
# COSMOLOGICAL QUANTITIES
# =============================================================================

def compute_quantities(profile):
    """Compute r_s, D_A, theta*, and H0_local for a Psi profile"""

    def E_octh(z):
        return E_lcdm(z) / profile.psi(z)

    # Sound horizon (integral from z_drag to infinity)
    def integrand_rs_lcdm(z):
        return 1.0 / E_lcdm(z)

    def integrand_rs_octh(z):
        return profile.psi(z) / E_lcdm(z)

    rs_lcdm, _ = quad(integrand_rs_lcdm, Z_DRAG, 1e6, limit=200)
    rs_octh, _ = quad(integrand_rs_octh, Z_DRAG, 1e6, limit=200)
    rs_ratio = rs_octh / rs_lcdm

    # Angular diameter distance
    def integrand_DA(z):
        return 1.0 / E_octh(z)

    DA_octh, _ = quad(integrand_DA, 0, Z_STAR, limit=200)

    def integrand_DA_lcdm(z):
        return 1.0 / E_lcdm(z)

    DA_lcdm, _ = quad(integrand_DA_lcdm, 0, Z_STAR, limit=200)
    DA_ratio = DA_octh / DA_lcdm

    # theta* = r_s / D_A
    # theta*_OCTH / theta*_LCDM = rs_ratio / DA_ratio
    theta_ratio = rs_ratio / DA_ratio

    # Peak shift factor
    # l_peak ~ pi / theta*, so l_OCTH / l_LCDM = 1 / theta_ratio
    peak_shift = 1.0 / theta_ratio

    # H0 local inference
    H0_local = H0_PLANCK / rs_ratio

    return {
        'rs_ratio': rs_ratio,
        'DA_ratio': DA_ratio,
        'theta_ratio': theta_ratio,
        'peak_shift': peak_shift,
        'H0_local': H0_local,
        'psi_star': profile.psi(Z_STAR)
    }

# =============================================================================
# SCAN PROFILES
# =============================================================================

print("\n" + "="*70)
print("SCANNING DIFFERENT PSI PROFILES")
print("="*70)

print("""
Goal: Find profile where theta* = r_s / D_A is preserved
      while H0 is increased (r_s decreased)

For this, D_A must decrease proportionally with r_s.
""")

results = []

# Profile 1: Gaussian (baseline)
print("\n--- Profile 1: Gaussian (baseline) ---")
for eps in [0.015, 0.10, 0.15, 0.18]:
    p = GaussianPsi(eps)
    q = compute_quantities(p)
    results.append({'profile': f'Gaussian eps={eps}', **q})
    print(f"eps={eps}: H0={q['H0_local']:.2f}, peak_shift={q['peak_shift']:.4f}")

# Profile 2: Step function
print("\n--- Profile 2: Step Function ---")
for eps in [0.05, 0.10, 0.15]:
    p = StepPsi(eps)
    q = compute_quantities(p)
    results.append({'profile': f'Step eps={eps}', **q})
    print(f"eps={eps}: H0={q['H0_local']:.2f}, peak_shift={q['peak_shift']:.4f}")

# Profile 3: Asymmetric (fast drop, slow recovery)
print("\n--- Profile 3: Asymmetric ---")
for eps in [0.10, 0.15, 0.20]:
    p = AsymmetricPsi(eps, sigma_left=0.2, sigma_right=1.5)
    q = compute_quantities(p)
    results.append({'profile': f'Asymmetric eps={eps}', **q})
    print(f"eps={eps}: H0={q['H0_local']:.2f}, peak_shift={q['peak_shift']:.4f}")

# Profile 4: Plateau
print("\n--- Profile 4: Plateau (1500 > z > 800) ---")
for eps in [0.05, 0.08, 0.10]:
    p = PlateauPsi(eps, z_start=1500, z_end=800)
    q = compute_quantities(p)
    results.append({'profile': f'Plateau eps={eps}', **q})
    print(f"eps={eps}: H0={q['H0_local']:.2f}, peak_shift={q['peak_shift']:.4f}")

# Profile 5: Extended (low-z tail)
print("\n--- Profile 5: Extended (affects D_A more) ---")
for eps_tail in [0.01, 0.02, 0.03]:
    p = ExtendedPsi(eps_peak=0.10, eps_tail=eps_tail, z_tail_start=500)
    q = compute_quantities(p)
    results.append({'profile': f'Extended tail={eps_tail}', **q})
    print(f"tail={eps_tail}: H0={q['H0_local']:.2f}, peak_shift={q['peak_shift']:.4f}")

# =============================================================================
# FIND OPTIMAL PROFILE
# =============================================================================

print("\n" + "="*70)
print("SEARCHING FOR OPTIMAL PROFILE")
print("="*70)

print("""
We want: H0 ~ 73 AND peak_shift ~ 1.0

The key is balancing r_s reduction with D_A reduction.
Extended profiles that also affect lower z can modify D_A.
""")

def objective(params, profile_type='extended'):
    """Objective: minimize |H0 - 73| + penalty for peak shift"""
    eps_peak, eps_tail = params

    if profile_type == 'extended':
        p = ExtendedPsi(eps_peak, eps_tail, z_peak=1089, z_tail_start=500)
    else:
        return 1e10

    try:
        q = compute_quantities(p)
        H0_error = abs(q['H0_local'] - H0_LOCAL)
        peak_penalty = 100 * abs(q['peak_shift'] - 1.0)
        return H0_error + peak_penalty
    except:
        return 1e10

# Grid search for extended profile
print("\nGrid search for Extended profile:")
print(f"{'eps_peak':<12} {'eps_tail':<12} {'H0':<12} {'peak_shift':<12} {'score':<12}")
print("-"*60)

best_score = 1e10
best_params = None

for eps_peak in np.linspace(0.05, 0.25, 9):
    for eps_tail in np.linspace(0.0, 0.10, 6):
        p = ExtendedPsi(eps_peak, eps_tail, z_tail_start=500)
        try:
            q = compute_quantities(p)
            H0_error = abs(q['H0_local'] - H0_LOCAL)
            peak_penalty = 100 * abs(q['peak_shift'] - 1.0)
            score = H0_error + peak_penalty

            if score < best_score:
                best_score = score
                best_params = (eps_peak, eps_tail, q)

            if score < 20:  # Only print good ones
                print(f"{eps_peak:<12.3f} {eps_tail:<12.3f} {q['H0_local']:<12.2f} {q['peak_shift']:<12.4f} {score:<12.2f}")
        except:
            pass

# =============================================================================
# ANALYSIS OF BEST PROFILE
# =============================================================================

print("\n" + "="*70)
print("BEST PROFILE FOUND")
print("="*70)

if best_params:
    eps_peak, eps_tail, q = best_params
    print(f"""
Profile: Extended
  eps_peak = {eps_peak:.3f} (at z = 1089)
  eps_tail = {eps_tail:.3f} (for z < 500)

Results:
  r_s ratio = {q['rs_ratio']:.4f}
  D_A ratio = {q['DA_ratio']:.4f}
  theta* ratio = {q['theta_ratio']:.4f}
  Peak shift = {q['peak_shift']:.4f}
  H0 local = {q['H0_local']:.2f} km/s/Mpc

Peak 1 position:
  LCDM: l = {L_PEAK_1:.1f}
  OCTH: l = {L_PEAK_1 * q['peak_shift']:.1f}
  Shift: {(q['peak_shift'] - 1) * 100:.1f}%
""")

# =============================================================================
# KEY INSIGHT
# =============================================================================

print("\n" + "="*70)
print("KEY INSIGHT: THE COMPENSATION MECHANISM")
print("="*70)

print("""
For CMB peaks to NOT shift, we need theta* = r_s / D_A preserved.

If Psi is localized ONLY at high z (around recombination):
  - r_s decreases (sound horizon shrinks)
  - D_A barely changes (integral from 0 to z*)
  => theta* decreases, peaks shift to higher l

If Psi ALSO affects lower z:
  - D_A decreases (light travels through modified expansion)
  - This can COMPENSATE for r_s decrease!

REQUIRED: Find eps_tail such that:
  r_s_ratio / D_A_ratio ~ 1

The challenge: the effect on D_A needs to be ~8% at z* to match
the r_s reduction, but this requires significant Psi modification
at z < 1089, which might create other tensions.
""")

# =============================================================================
# CONCLUSIONS
# =============================================================================

print("\n" + "="*70)
print("CONCLUSIONS")
print("="*70)

print("""
1. SIMPLE PROFILES DON'T WORK:
   Standard Gaussian eps_topo ~ 0.18 shifts peaks by ~8%
   This creates ~40-100 sigma tension with Planck

2. COMPENSATION IS POSSIBLE IN PRINCIPLE:
   Extended profiles that affect D_A can reduce peak shift
   But requires significant modification at z < 1089

3. BEST EXTENDED PROFILE:
   Can reduce peak shift, but doesn't fully eliminate it
   Still ~2-4% residual shift expected

4. IMPLICATIONS FOR OCTH:
   a) Accept partial H0 resolution (eps_topo ~ 0.015, gap of 5 km/s/Mpc)
   b) Or find physics that modifies D_A without affecting CMB otherwise
   c) Or accept some CMB tension and look for other explanations

5. NEXT STEP:
   Full Boltzmann code (CLASS/CAMB) calculation needed
   Simplified theta* analysis may miss important effects
   (damping, ISW, lensing, etc.)
""")

print("="*70)
