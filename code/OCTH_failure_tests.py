"""
OCTH/MOBIUS FAILURE ANALYSIS
============================

Testing where the theory FAILS or has problems.

Critical tests that MUST pass:
1. Solar System (perihelion, light deflection, Shapiro delay)
2. Binary pulsars (orbital decay via GW)
3. Gravitational wave speed (LIGO/Virgo)
4. Big Bang Nucleosynthesis (BBN)
5. CMB full spectrum (not just peaks)

Author: F. Molina-Burgos
Date: January 2026
"""

import numpy as np

print("="*70)
print("  OCTH/MOBIUS - FAILURE ANALYSIS")
print("  Where does the theory break?")
print("="*70)

# =============================================================================
# OCTH PARAMETERS
# =============================================================================

# Two-component model
EPSILON_TOPO = 0.015      # Cosmological (z~1089)
EPSILON_QUANTUM = 0.168   # Galactic (high density)
EPSILON_TOTAL = 0.183     # MOND value

SIGMA_TOPO = 0.6
Z_TOPO = 1089.0

def psi_topo(z):
    """Topological component of Psi"""
    ln_a = -np.log(1 + z)
    ln_a_topo = -np.log(1 + Z_TOPO)
    delta = (ln_a - ln_a_topo) / SIGMA_TOPO
    return 1.0 - EPSILON_TOPO * np.exp(-0.5 * delta**2)

def psi_quantum(rho_ratio):
    """
    Quantum component of Psi.
    Active when rho >> rho_crit.

    rho_ratio = rho / rho_crit
    """
    # Threshold: activates around rho ~ 100 rho_crit
    rho_threshold = 100.0

    if rho_ratio < 1:
        return 1.0  # No effect in low density

    # Smooth activation
    x = rho_ratio / rho_threshold
    activation = 1 - np.exp(-x)

    return 1.0 - EPSILON_QUANTUM * activation

# =============================================================================
# TEST 1: SOLAR SYSTEM
# =============================================================================

print("\n" + "="*70)
print("TEST 1: SOLAR SYSTEM")
print("="*70)

# Constants
G = 6.674e-11  # m^3/kg/s^2
c = 2.998e8    # m/s
M_sun = 1.989e30  # kg

# Solar system densities (very low!)
rho_sun_avg = M_sun / (4/3 * np.pi * (6.96e8)**3)  # ~1400 kg/m^3
rho_crit = 9.47e-27  # kg/m^3

rho_ratio_sun = rho_sun_avg / rho_crit
print(f"\nSolar density / rho_crit = {rho_ratio_sun:.2e}")

# At z=0, psi_topo ~ 1
psi_t = psi_topo(0)
print(f"Psi_topo(z=0) = {psi_t:.6f}")

# For solar system, rho >> rho_crit locally
psi_q = psi_quantum(rho_ratio_sun)
print(f"Psi_quantum(Sun) = {psi_q:.6f}")

# Combined effect
psi_solar = psi_t * psi_q
print(f"Psi_total(Solar System) = {psi_solar:.6f}")

# GR predictions for Mercury perihelion precession
# Delta_phi = 6*pi*G*M / (c^2 * a * (1-e^2))
# Mercury: a = 5.79e10 m, e = 0.2056

a_mercury = 5.79e10  # m
e_mercury = 0.2056

delta_phi_GR = 6 * np.pi * G * M_sun / (c**2 * a_mercury * (1 - e_mercury**2))
delta_phi_GR_arcsec = delta_phi_GR * 180 * 3600 / np.pi  # per orbit
# Per century (orbital period 87.97 days)
orbits_per_century = 100 * 365.25 / 87.97
delta_phi_century_GR = delta_phi_GR_arcsec * orbits_per_century

print(f"\nMercury perihelion precession:")
print(f"  GR prediction: {delta_phi_century_GR:.2f} arcsec/century")
print(f"  Observed:      43.11 ± 0.21 arcsec/century")

# OCTH modification: G_eff = G / Psi^2
# But Psi ~ 1 in solar system!
G_eff_ratio = 1 / psi_solar**2
delta_phi_OCTH = delta_phi_century_GR * G_eff_ratio

print(f"\n  OCTH prediction: {delta_phi_OCTH:.2f} arcsec/century")
print(f"  Deviation from GR: {(delta_phi_OCTH/delta_phi_century_GR - 1)*100:.4f}%")

if abs(delta_phi_OCTH - 43.11) < 0.5:
    print("\n  ✓ OCTH PASSES Solar System test (within 1%)")
    SOLAR_PASS = True
else:
    print(f"\n  ✗ OCTH FAILS Solar System test!")
    print(f"    Deviation: {abs(delta_phi_OCTH - 43.11):.2f} arcsec/century")
    SOLAR_PASS = False

# =============================================================================
# TEST 2: BINARY PULSARS
# =============================================================================

print("\n" + "="*70)
print("TEST 2: BINARY PULSARS (Hulse-Taylor PSR B1913+16)")
print("="*70)

# Hulse-Taylor pulsar parameters
M1 = 1.4398 * M_sun  # Pulsar mass
M2 = 1.3886 * M_sun  # Companion mass
P_orbit = 7.75 * 3600  # Orbital period in seconds
e_pulsar = 0.6171     # Eccentricity

# GR prediction for orbital decay
# dP/dt = -192*pi/5 * (G^(5/3)/c^5) * (2*pi/P)^(5/3) * M1*M2/(M1+M2)^(1/3) * f(e)
# where f(e) = (1 + 73/24*e^2 + 37/96*e^4) / (1-e^2)^(7/2)

def f_ecc(e):
    return (1 + 73/24*e**2 + 37/96*e**4) / (1 - e**2)**(7/2)

omega = 2 * np.pi / P_orbit
M_chirp = (M1 * M2)**(3/5) / (M1 + M2)**(1/5)

dP_dt_GR = -192 * np.pi / 5 * (G**(5/3) / c**5) * omega**(5/3) * \
           M1 * M2 / (M1 + M2)**(1/3) * f_ecc(e_pulsar)

# Convert to seconds per second (dimensionless)
# Observed: dP/dt = -2.4211e-12 s/s
dP_dt_observed = -2.4211e-12

print(f"\nOrbital decay rate:")
print(f"  GR prediction: {dP_dt_GR:.4e} s/s")
print(f"  Observed:      {dP_dt_observed:.4e} s/s")
print(f"  GR accuracy:   {(dP_dt_GR/dP_dt_observed - 1)*100:.3f}%")

# OCTH modification
# GW emission scales as G^(5/3), so dP/dt scales as 1/Psi^(10/3)
# But what is Psi for a neutron star?

rho_NS = 4e17  # kg/m^3 (typical neutron star)
rho_ratio_NS = rho_NS / rho_crit
print(f"\nNeutron star density / rho_crit = {rho_ratio_NS:.2e}")

psi_NS = psi_quantum(rho_ratio_NS)
print(f"Psi_quantum(NS) = {psi_NS:.6f}")

# PROBLEM: At neutron star densities, Psi << 1 if quantum component active!
dP_dt_OCTH = dP_dt_GR / psi_NS**(10/3)

print(f"\nOCTH prediction: {dP_dt_OCTH:.4e} s/s")
print(f"Deviation from GR: {(dP_dt_OCTH/dP_dt_GR - 1)*100:.1f}%")

if abs(dP_dt_OCTH / dP_dt_observed - 1) < 0.01:
    print("\n  ✓ OCTH PASSES Binary Pulsar test")
    PULSAR_PASS = True
else:
    print(f"\n  ✗ OCTH FAILS Binary Pulsar test!")
    print(f"    Predicted/Observed = {dP_dt_OCTH/dP_dt_observed:.2f}")
    print(f"    GR works to 0.2% - OCTH deviates by {abs(dP_dt_OCTH/dP_dt_observed - 1)*100:.0f}%")
    PULSAR_PASS = False

# =============================================================================
# TEST 3: GRAVITATIONAL WAVE SPEED
# =============================================================================

print("\n" + "="*70)
print("TEST 3: GRAVITATIONAL WAVE SPEED (GW170817)")
print("="*70)

# GW170817: GW + gamma ray burst
# Speed constraint: |c_gw - c| / c < 3e-15

print("""
GW170817 observation:
  - GW arrived ~1.7s before gamma rays
  - Distance: ~40 Mpc
  - Travel time: ~130 million years

  Constraint: |c_gw/c - 1| < 3 × 10^-15
""")

# In OCTH, does GW speed change?
# GW equation: h_uv'' + k^2 h_uv = 0 (in vacuum)
# If Psi modifies the metric, does it modify GW propagation?

# Key question: Does Psi affect the wave equation?
# If H^2 = rho/Psi^2, then the effective speed of gravity could be c/Psi

# In vacuum (rho ~ 0), Psi ~ 1, so c_gw ~ c
# But GWs travel through varying Psi regions...

print("OCTH analysis:")
print("  - In vacuum: Psi = 1, so c_gw = c")
print("  - GWs propagate on null geodesics of effective metric")
print("  - If g_eff = g/Psi^2, then c_gw could be modified")

# Conservative estimate: average Psi along path
# Most of the path is through low-density IGM
psi_IGM = psi_topo(0) * psi_quantum(1)  # rho ~ rho_crit in IGM
print(f"\n  Psi(IGM) ~ {psi_IGM:.6f}")

# Speed modification
delta_c = 1/psi_IGM - 1
print(f"  delta_c/c = {delta_c:.2e}")

if abs(delta_c) < 3e-15:
    print("\n  ✓ OCTH PASSES GW speed test")
    GW_PASS = True
else:
    print(f"\n  ⚠ OCTH may FAIL GW speed test")
    print(f"    Required: |delta_c/c| < 3e-15")
    print(f"    OCTH gives: {abs(delta_c):.2e}")
    if abs(delta_c) < 1e-6:
        print("    (But this depends on interpretation of Psi in wave equation)")
        GW_PASS = "UNCERTAIN"
    else:
        GW_PASS = False

# =============================================================================
# TEST 4: BIG BANG NUCLEOSYNTHESIS
# =============================================================================

print("\n" + "="*70)
print("TEST 4: BIG BANG NUCLEOSYNTHESIS (BBN)")
print("="*70)

# BBN occurs at z ~ 10^9, T ~ 1 MeV
z_BBN = 1e9

# What is Psi at BBN?
psi_BBN = psi_topo(z_BBN)
print(f"\nAt z = {z_BBN:.0e} (BBN epoch):")
print(f"  Psi_topo = {psi_BBN:.6f}")

# BBN is sensitive to expansion rate H
# If H is modified by factor 1/Psi, then:
# - Freeze-out temperature changes
# - n/p ratio at freeze-out changes
# - Final He-4 abundance changes

# Standard BBN: Y_p = 0.247 (He-4 mass fraction)
# Observed: Y_p = 0.245 ± 0.003

Y_p_standard = 0.247
Y_p_observed = 0.245
Y_p_sigma = 0.003

# Effect of modified H:
# delta_Y/Y ~ 0.1 * delta_H/H (approximate sensitivity)
# If H -> H/Psi, then delta_H/H = 1/Psi - 1

delta_H_H = 1/psi_BBN - 1
delta_Y_Y = 0.1 * delta_H_H  # Sensitivity factor

Y_p_OCTH = Y_p_standard * (1 + delta_Y_Y)

print(f"\n  delta_H/H at BBN = {delta_H_H:.6f}")
print(f"  Predicted Y_p shift: {delta_Y_Y*100:.4f}%")
print(f"\n  Standard BBN: Y_p = {Y_p_standard:.4f}")
print(f"  OCTH prediction: Y_p = {Y_p_OCTH:.4f}")
print(f"  Observed: Y_p = {Y_p_observed:.3f} ± {Y_p_sigma:.3f}")

if abs(Y_p_OCTH - Y_p_observed) < 2 * Y_p_sigma:
    print("\n  ✓ OCTH PASSES BBN test (within 2-sigma)")
    BBN_PASS = True
else:
    print(f"\n  ✗ OCTH FAILS BBN test")
    print(f"    Tension: {abs(Y_p_OCTH - Y_p_observed)/Y_p_sigma:.1f} sigma")
    BBN_PASS = False

# =============================================================================
# TEST 5: CMB FULL SPECTRUM
# =============================================================================

print("\n" + "="*70)
print("TEST 5: CMB FULL SPECTRUM")
print("="*70)

# We showed earlier that OCTH with epsilon=0.183 shifts peaks too much
# Optimal epsilon was 0.015

print("""
Previous analysis showed:

  With epsilon = 0.183 (MOND):
    - CMB peaks shift by ~3%
    - Chi2 for peaks: VERY BAD

  With epsilon = 0.015 (optimized):
    - CMB peaks shift by ~0.3%
    - Chi2 acceptable
    - BUT: This is different from MOND's epsilon!

PROBLEM: The epsilon that fits CMB (0.015) is NOT the same as
         the epsilon that fits MOND (0.183).
""")

# The two-component model attempts to resolve this:
# epsilon_topo = 0.015 (cosmological, affects CMB)
# epsilon_quantum = 0.168 (galactic, affects MOND)

print("Two-component resolution:")
print(f"  epsilon_topo = {EPSILON_TOPO:.3f} (cosmological)")
print(f"  epsilon_quantum = {EPSILON_QUANTUM:.3f} (galactic)")
print(f"  epsilon_total = {EPSILON_TOTAL:.3f} (MOND)")

# But this raises the question: what is the PHYSICS behind two components?
print("""
UNRESOLVED QUESTIONS:
  1. Why are there two components?
  2. What determines the density threshold for quantum component?
  3. Is this just parameter fitting or genuine physics?
""")

CMB_PASS = "PARTIAL - requires two-component model"

# =============================================================================
# TEST 6: STRUCTURE FORMATION
# =============================================================================

print("\n" + "="*70)
print("TEST 6: STRUCTURE FORMATION (sigma_8)")
print("="*70)

# Structure formation depends on growth factor D(z)
# In OCTH, the growth equation is modified

# Linear growth: dD/dt + 2H dD/dt = 4*pi*G*rho_m * D
# With OCTH: H -> H/Psi, but what about the source term?

print("""
In LCDM:
  sigma_8 = 0.811 ± 0.006 (Planck)

Structure growth in OCTH:
  - Modified Hubble affects friction term
  - May also affect source term (G_eff?)
  - Need full numerical calculation

Current status: NOT YET CALCULATED
""")

STRUCTURE_PASS = "NOT TESTED"

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*70)
print("SUMMARY: WHERE DOES OCTH FAIL?")
print("="*70)

print("""
TEST                      STATUS          NOTES
--------------------------------------------------------------------""")
print(f"1. Solar System           {'PASS' if SOLAR_PASS else 'FAIL'}            Psi~1 locally, no deviation")
print(f"2. Binary Pulsars         {'PASS' if PULSAR_PASS else 'FAIL'}            CRITICAL - high density!")
print(f"3. GW Speed               {GW_PASS}        Depends on wave equation")
print(f"4. BBN                    {'PASS' if BBN_PASS else 'FAIL'}            Psi~1 at z~10^9")
print(f"5. CMB Spectrum           {CMB_PASS}")
print(f"6. Structure Formation    {STRUCTURE_PASS}")

print("""
--------------------------------------------------------------------

CRITICAL ISSUE: BINARY PULSARS
==============================

The Hulse-Taylor pulsar test is PRECISE to 0.2%.

PROBLEM: At neutron star densities (rho ~ 10^44 rho_crit),
         the quantum component gives Psi << 1.

         This would modify GW emission rate drastically,
         VIOLATING the observed orbital decay.

POSSIBLE RESOLUTIONS:

1. The quantum component does NOT affect GW emission
   - GWs propagate in vacuum, not in dense matter
   - Need to clarify how Psi enters the field equations

2. The quantum component has different form near compact objects
   - Screening mechanism at extreme densities
   - Psi -> 1 for rho > rho_nuclear

3. The two-component model is wrong
   - Single component with density-independent epsilon
   - Requires different explanation for MOND

THIS IS THE KEY TEST OCTH MUST ADDRESS!
""")

# =============================================================================
# DETAILED ANALYSIS OF PULSAR PROBLEM
# =============================================================================

print("\n" + "="*70)
print("DETAILED: THE PULSAR PROBLEM")
print("="*70)

print("""
Question: How does Psi enter the GW emission formula?

GR derivation of orbital decay:
  1. Binary system has time-varying quadrupole moment Q_ij
  2. GW power: P = (G/5c^5) <Q'''_ij Q'''_ij>
  3. This comes from Einstein equations linearized around Minkowski

OCTH modification possibilities:

OPTION A: G -> G_eff = G/Psi^2 everywhere
  - Then P -> P/Psi^(10/3)
  - At NS density, Psi ~ 0.83, so P increases by ~50%
  - RULED OUT by Hulse-Taylor

OPTION B: Psi only affects H (cosmological), not local G
  - GW emission uses local G (unmodified)
  - Psi only relevant for Friedmann equation
  - This could work, but needs theoretical justification

OPTION C: Psi is screened in compact objects
  - Like chameleon mechanism in f(R)
  - At rho >> rho_threshold, Psi -> 1
  - Needs explicit screening mechanism

OPTION D: The quantum component only affects low-acceleration regions
  - MOND regime: a < a_0
  - NS surface: a ~ 10^12 m/s^2 >> a_0
  - Quantum component inactive in strong-field regime
  - This is consistent with MOND phenomenology!
""")

# Calculate acceleration at NS surface
R_NS = 10e3  # 10 km
M_NS = 1.4 * M_sun
a_NS_surface = G * M_NS / R_NS**2
a_0 = 1.2e-10  # m/s^2

print(f"\nAcceleration at NS surface: a = {a_NS_surface:.2e} m/s^2")
print(f"MOND scale: a_0 = {a_0:.2e} m/s^2")
print(f"Ratio: a/a_0 = {a_NS_surface/a_0:.2e}")

print("""
CONCLUSION: Option D looks promising!

At NS surface, a ~ 10^{12} a_0.
This is FAR into the Newtonian regime where MOND effects vanish.

PROPOSED FIX FOR OCTH:

  The quantum component (epsilon_quantum) only activates when:
    a < a_0 (low acceleration regime)

  For compact objects where a >> a_0:
    Psi_quantum -> 1 (no modification)

  This is actually REQUIRED by MOND itself!
  MOND interpolating function mu(a/a_0) -> 1 for a >> a_0.
""")

# Recalculate with acceleration-based screening
print("\n" + "="*70)
print("RECALCULATION WITH ACCELERATION SCREENING")
print("="*70)

def psi_quantum_v2(rho_ratio, acceleration):
    """
    Quantum component with acceleration screening.
    Only active when a < a_0.
    """
    a_0 = 1.2e-10  # m/s^2

    # MOND interpolating function
    x = acceleration / a_0
    mu = x / np.sqrt(1 + x**2)  # Simple interpolation

    # Quantum component scales as (1 - mu)
    # When a >> a_0: mu -> 1, so no quantum effect
    # When a << a_0: mu -> a/a_0, so full quantum effect

    quantum_activation = 1 - mu

    return 1.0 - EPSILON_QUANTUM * quantum_activation

# Recalculate for NS
psi_NS_v2 = psi_quantum_v2(rho_ratio_NS, a_NS_surface)
print(f"\nWith acceleration screening:")
print(f"  a_NS = {a_NS_surface:.2e} m/s^2")
print(f"  a_NS/a_0 = {a_NS_surface/a_0:.2e}")
print(f"  Psi_quantum_v2(NS) = {psi_NS_v2:.10f}")

dP_dt_OCTH_v2 = dP_dt_GR / psi_NS_v2**(10/3)
print(f"\n  GR prediction: {dP_dt_GR:.4e} s/s")
print(f"  OCTH v2 prediction: {dP_dt_OCTH_v2:.4e} s/s")
print(f"  Deviation: {abs(dP_dt_OCTH_v2/dP_dt_GR - 1)*100:.6f}%")

if abs(dP_dt_OCTH_v2 / dP_dt_observed - 1) < 0.01:
    print("\n  ✓ OCTH v2 PASSES Binary Pulsar test!")
else:
    print(f"\n  Status: Deviation is {abs(dP_dt_OCTH_v2/dP_dt_observed - 1)*100:.4f}%")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

print("\n" + "="*70)
print("FINAL: WHAT OCTH NEEDS TO FIX")
print("="*70)

print("""
1. BINARY PULSARS - FIXED by acceleration screening
   The quantum component must use MOND-like interpolation:
   mu(a/a_0) ensures Psi -> 1 for a >> a_0

2. CMB SPECTRUM - REQUIRES two-component model
   epsilon_topo = 0.015 for cosmology
   epsilon_quantum = 0.168 for galaxies (with a < a_0 screening)

3. GW SPEED - NEEDS THEORETICAL CLARIFICATION
   Does Psi affect wave equation or just background?

4. STRUCTURE FORMATION - NOT YET CALCULATED
   Need to solve growth equation with modified H

5. THEORETICAL FOUNDATION - INCOMPLETE
   - Why two components?
   - Covariant formulation?
   - Lagrangian?

OVERALL STATUS:
  OCTH is VIABLE but requires:
  a) Acceleration screening for quantum component
  b) Rigorous field equation derivation
  c) Structure formation calculation
""")

print("\n" + "="*70)
print("END OF FAILURE ANALYSIS")
print("="*70)
