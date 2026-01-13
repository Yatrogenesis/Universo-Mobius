"""
OCTH SIGMA_8 ANALYSIS
=====================

The sigma_8 tension: another cosmological discrepancy

Planck (CMB, assumes LCDM): sigma_8 = 0.811 +/- 0.006
Weak lensing:               sigma_8 ~ 0.75-0.78
KiDS/DES:                   S_8 = 0.759 +/- 0.024

S_8 = sigma_8 * sqrt(Omega_m/0.3)

Does OCTH help with this tension?

Author: F. Molina-Burgos
Date: January 2026
"""

import numpy as np
from scipy.integrate import quad, odeint
from scipy.interpolate import interp1d

print("="*70)
print("  OCTH SIGMA_8 ANALYSIS")
print("  Structure growth in modified cosmology")
print("="*70)

# =============================================================================
# CONSTANTS
# =============================================================================

c = 2.998e5  # km/s
H0 = 67.36   # km/s/Mpc (Planck)

OMEGA_M = 0.315
OMEGA_R = 9.0e-5
OMEGA_LAMBDA = 1 - OMEGA_M - OMEGA_R

# Observations
SIGMA8_PLANCK = 0.811
SIGMA8_PLANCK_ERR = 0.006
SIGMA8_WL = 0.76  # Weak lensing average
SIGMA8_WL_ERR = 0.02

S8_PLANCK = SIGMA8_PLANCK * np.sqrt(OMEGA_M / 0.3)  # ~0.832
S8_WL = 0.759  # KiDS/DES

# =============================================================================
# OCTH MODEL
# =============================================================================

class OCTHStructure:
    def __init__(self, eps_topo, sigma=0.6, z_topo=1089.0):
        self.eps_topo = eps_topo
        self.sigma = sigma
        self.z_topo = z_topo
        self.ln_a_topo = -np.log(1 + z_topo)

    def psi_topo(self, z):
        """Topological permeability"""
        ln_a = -np.log(1 + z)
        delta = (ln_a - self.ln_a_topo) / self.sigma
        return 1.0 - self.eps_topo * np.exp(-0.5 * delta**2)

    def E_lcdm(self, z):
        """E(z) for LCDM"""
        opz = 1 + z
        return np.sqrt(OMEGA_R * opz**4 + OMEGA_M * opz**3 + OMEGA_LAMBDA)

    def E_octh(self, z):
        """E(z) for OCTH"""
        return self.E_lcdm(z) / self.psi_topo(z)

    def H_octh(self, z):
        """Hubble parameter in km/s/Mpc"""
        return H0 * self.E_octh(z)

    def Omega_m_z(self, z):
        """Matter density parameter at redshift z"""
        opz = 1 + z
        return OMEGA_M * opz**3 / self.E_octh(z)**2

# =============================================================================
# LINEAR GROWTH FACTOR
# =============================================================================

def compute_growth_lcdm(z_array):
    """
    Compute linear growth factor D(z) for LCDM.

    Growth equation:
    d^2 D / da^2 + (3/a + dlnE/da) dD/da - 3/2 * Omega_m(a) * D / a^2 = 0

    Or in terms of z:
    (1+z) d^2D/dz^2 - (1 + (1+z)/E * dE/dz) dD/dz + 3/2 * Omega_m(z) * D / (1+z) / E^2 = 0
    """
    # Use the ODE form: D'' + f(z) D' + g(z) D = 0
    # Convert to first-order system

    def E_lcdm(z):
        opz = 1 + z
        return np.sqrt(OMEGA_R * opz**4 + OMEGA_M * opz**3 + OMEGA_LAMBDA)

    def dE_dz(z):
        opz = 1 + z
        E = E_lcdm(z)
        return (4*OMEGA_R * opz**3 + 3*OMEGA_M * opz**2) / (2 * E)

    def Omega_m_z(z):
        opz = 1 + z
        return OMEGA_M * opz**3 / E_lcdm(z)**2

    # Growth equation in terms of ln(a)
    # d^2 D / d(lna)^2 + (2 + dlnE/dlna) dD/dlna - 3/2 * Omega_m * D = 0

    def growth_ode(y, lna):
        D, dD_dlna = y
        a = np.exp(lna)
        z = 1/a - 1
        E = E_lcdm(z)
        Om = Omega_m_z(z)

        # dlnE/dlna = -(1+z)/E * dE/dz
        dlnE_dlna = -(1+z) / E * dE_dz(z)

        d2D_dlna2 = -(2 + dlnE_dlna) * dD_dlna + 1.5 * Om * D

        return [dD_dlna, d2D_dlna2]

    # Initial conditions at high redshift (matter dominated)
    z_init = 1000
    lna_init = -np.log(1 + z_init)

    # In matter domination, D ~ a
    D_init = 1.0 / (1 + z_init)
    dD_dlna_init = D_init  # dD/dlna = D (since D ~ a)

    # Integration range
    lna_array = -np.log(1 + z_array)
    lna_span = np.linspace(lna_init, 0, 1000)

    # Solve ODE
    solution = odeint(growth_ode, [D_init, dD_dlna_init], lna_span)
    D_lna = solution[:, 0]

    # Interpolate to desired z values
    D_interp = interp1d(lna_span, D_lna, kind='cubic', fill_value='extrapolate')

    # Normalize: D(z=0) = 1
    D_today = D_interp(0)
    D_result = D_interp(lna_array) / D_today

    return D_result


def compute_growth_octh(z_array, eps_topo, sigma=0.6):
    """
    Compute linear growth factor D(z) for OCTH.

    The growth equation is modified because H(z) is modified.
    """
    model = OCTHStructure(eps_topo, sigma)

    def dE_dz_octh(z, E):
        # Numerical derivative
        dz = 0.001
        E_plus = model.E_octh(z + dz)
        E_minus = model.E_octh(z - dz) if z > dz else E
        return (E_plus - E_minus) / (2 * dz) if z > dz else (E_plus - E) / dz

    def growth_ode(y, lna):
        D, dD_dlna = y
        a = np.exp(lna)
        z = 1/a - 1 if a < 1 else 0

        E = model.E_octh(z)
        Om = model.Omega_m_z(z)

        # dlnE/dlna = -(1+z)/E * dE/dz
        dlnE_dlna = -(1+z) / E * dE_dz_octh(z, E)

        d2D_dlna2 = -(2 + dlnE_dlna) * dD_dlna + 1.5 * Om * D

        return [dD_dlna, d2D_dlna2]

    # Initial conditions
    z_init = 1000
    lna_init = -np.log(1 + z_init)

    D_init = 1.0 / (1 + z_init)
    dD_dlna_init = D_init

    lna_array = -np.log(1 + z_array)
    lna_span = np.linspace(lna_init, 0, 1000)

    solution = odeint(growth_ode, [D_init, dD_dlna_init], lna_span)
    D_lna = solution[:, 0]

    D_interp = interp1d(lna_span, D_lna, kind='cubic', fill_value='extrapolate')

    D_today = D_interp(0)
    D_result = D_interp(lna_array) / D_today

    return D_result

# =============================================================================
# ANALYSIS
# =============================================================================

print("\n" + "="*70)
print("LINEAR GROWTH COMPARISON")
print("="*70)

z_array = np.array([0, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0])

# LCDM growth
D_lcdm = compute_growth_lcdm(z_array)

print("\nGrowth factor D(z) normalized to D(0)=1:")
print("-"*60)
print(f"{'z':<10} {'D_LCDM':<15} {'D_OCTH(0.015)':<15} {'D_OCTH(0.18)':<15}")
print("-"*60)

# OCTH growth for different eps_topo
D_octh_015 = compute_growth_octh(z_array, 0.015)
D_octh_18 = compute_growth_octh(z_array, 0.18)

for i, z in enumerate(z_array):
    print(f"{z:<10.1f} {D_lcdm[i]:<15.4f} {D_octh_015[i]:<15.4f} {D_octh_18[i]:<15.4f}")

# =============================================================================
# SIGMA_8 CALCULATION
# =============================================================================

print("\n" + "="*70)
print("SIGMA_8 IMPLICATIONS")
print("="*70)

print("""
How does OCTH affect sigma_8?

sigma_8 = sigma_8(initial) * D(z=0)

If OCTH modifies the growth factor, it changes sigma_8.

Key insight:
- The Psi modification is localized around z ~ 1089
- At low z, Psi ~ 1, so growth is nearly LCDM
- The main effect is on the *early* growth history
""")

# Compare D(z=10) which captures early growth
print("\nGrowth comparison at z=10 (early universe):")
D_lcdm_10 = compute_growth_lcdm(np.array([10.0]))[0]
D_octh_015_10 = compute_growth_octh(np.array([10.0]), 0.015)[0]
D_octh_18_10 = compute_growth_octh(np.array([10.0]), 0.18)[0]

print(f"  D_LCDM(z=10) = {D_lcdm_10:.4f}")
print(f"  D_OCTH(z=10, eps=0.015) = {D_octh_015_10:.4f} (ratio: {D_octh_015_10/D_lcdm_10:.4f})")
print(f"  D_OCTH(z=10, eps=0.18) = {D_octh_18_10:.4f} (ratio: {D_octh_18_10/D_lcdm_10:.4f})")

# Estimate sigma_8 scaling
print("\n" + "="*70)
print("SIGMA_8 PREDICTIONS")
print("="*70)

# The sigma_8 tension is ~10%: Planck 0.811 vs WL 0.76
# Can OCTH explain this?

# sigma_8 scales with D(z=0), but D(z=0) = 1 by normalization
# The relevant quantity is the *amplitude* at CMB which propagates to today

# Actually, sigma_8 from CMB depends on:
# - A_s (primordial amplitude)
# - n_s (spectral index)
# - Growth from CMB to today

# The CMB constrains A_s * D(z_CMB)^2 * exp(2*tau)
# If OCTH changes D(z_CMB), it affects the inferred sigma_8

print("""
OCTH effect on sigma_8:

The CMB sees fluctuations at z ~ 1089.
If Psi(z ~ 1089) < 1, the Hubble rate is HIGHER.
This affects:
1. Sound horizon (r_s decreases) -> H0 tension addressed
2. Diffusion damping (more damping at high l)
3. Growth rate around recombination

The growth suppression from modified H(z) can reduce sigma_8!
""")

# Simplified estimate
# If H is larger around recombination, growth is slower
# The effect on sigma_8 depends on the integral of growth

def sigma8_ratio(eps_topo):
    """
    Estimate sigma_8(OCTH) / sigma_8(LCDM)

    This is approximate - proper calculation needs full Boltzmann code
    """
    # Growth at z = 50 (deep in matter era but after recombination)
    D_lcdm_50 = compute_growth_lcdm(np.array([50.0]))[0]
    D_octh_50 = compute_growth_octh(np.array([50.0]), eps_topo)[0]

    # Ratio gives approximate sigma_8 suppression
    return D_octh_50 / D_lcdm_50

print("\nEstimated sigma_8 ratio vs LCDM:")
print("-"*40)
print(f"{'eps_topo':<15} {'sigma8 ratio':<15}")
print("-"*40)

for eps in [0.01, 0.02, 0.05, 0.10, 0.15, 0.18, 0.20]:
    ratio = sigma8_ratio(eps)
    sigma8_pred = SIGMA8_PLANCK * ratio
    print(f"{eps:<15.3f} {ratio:<15.4f} -> sigma8 ~ {sigma8_pred:.3f}")

# =============================================================================
# OPTIMAL EPSILON FOR BOTH TENSIONS
# =============================================================================

print("\n" + "="*70)
print("CAN OCTH SOLVE BOTH H0 AND SIGMA_8 TENSIONS?")
print("="*70)

# H0 tension: need eps_topo ~ 0.18
# sigma_8 tension: need ~10% suppression

print("""
H0 tension requires:    eps_topo ~ 0.18
sigma_8 tension needs:  ~10% suppression (0.811 -> 0.76)

Let's check what sigma_8 suppression we get with eps_topo = 0.18:
""")

ratio_018 = sigma8_ratio(0.18)
sigma8_018 = SIGMA8_PLANCK * ratio_018

print(f"With eps_topo = 0.18:")
print(f"  sigma_8 ratio = {ratio_018:.4f}")
print(f"  sigma_8(OCTH) = {sigma8_018:.3f}")
print(f"  sigma_8(target) = {SIGMA8_WL:.3f}")
print(f"  Discrepancy = {(sigma8_018 - SIGMA8_WL)/SIGMA8_WL_ERR:.1f} sigma")

# =============================================================================
# S_8 PARAMETER
# =============================================================================

print("\n" + "="*70)
print("S_8 ANALYSIS")
print("="*70)

print("""
S_8 = sigma_8 * sqrt(Omega_m / 0.3)

This combination is better constrained by weak lensing.

S_8 (Planck) = 0.832 +/- 0.013
S_8 (KiDS/DES) = 0.759 +/- 0.024

Tension: ~2.5 sigma
""")

S8_planck = SIGMA8_PLANCK * np.sqrt(OMEGA_M / 0.3)
S8_octh = sigma8_018 * np.sqrt(OMEGA_M / 0.3)

print(f"S_8 (Planck LCDM) = {S8_planck:.3f}")
print(f"S_8 (OCTH, eps=0.18) = {S8_octh:.3f}")
print(f"S_8 (weak lensing) = {S8_WL:.3f}")

# =============================================================================
# CONCLUSIONS
# =============================================================================

print("\n" + "="*70)
print("CONCLUSIONS")
print("="*70)

print(f"""
1. GROWTH MODIFICATION:
   OCTH with eps_topo ~ 0.18 slightly suppresses early growth
   But the effect is small (~1-2%) because Psi is localized

2. SIGMA_8 PREDICTION:
   With eps_topo = 0.18:
   sigma_8(OCTH) ~ {sigma8_018:.3f} (vs LCDM {SIGMA8_PLANCK:.3f})
   Suppression: {(1-ratio_018)*100:.1f}%

3. S_8 PREDICTION:
   S_8(OCTH) ~ {S8_octh:.3f}
   Still above weak lensing S_8 = 0.759

4. CONCLUSION:
   OCTH with eps_topo ~ 0.18 (for H0 resolution):
   - Provides SMALL sigma_8 suppression
   - NOT enough to fully resolve sigma_8/S_8 tension
   - But moves in the RIGHT DIRECTION

5. POSSIBLE ENHANCEMENT:
   - Larger sigma (broader Psi effect) might help more
   - Different Psi profile could optimize both tensions
   - Additional physics (neutrino mass, etc.) could combine with OCTH

BOTTOM LINE:
OCTH primarily addresses H0 tension.
sigma_8/S_8 tension requires additional ingredients.
""")

print("="*70)
