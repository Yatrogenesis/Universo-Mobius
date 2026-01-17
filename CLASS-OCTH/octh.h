/**
 * @file octh.h
 *
 * OCTH (Observational Cosmological Topological Hypothesis) modification
 * Based on PIT (Primacia Informacional Topologica) framework
 * Part of Universo-Mobius project
 *
 * Author: Francisco Molina-Burgos
 * Avermex Research Division
 * January 2026
 *
 * PIT Framework:
 *   - Topology is fundamental, metric is emergent
 *   - Topological transitions PRECEDE metric changes
 *   - a0 = c * H0 / (2*pi) derived from topology
 *   - epsilon_geometric = 1/(2*pi) ~ 0.159
 *
 * The OCTH hypothesis modifies Hubble expansion:
 *   H_OCTH(z) = H_LCDM(z) / Psi(z)
 *
 * PIT-derived profile:
 *   Psi(z) = sqrt(1 - epsilon(z))
 *   epsilon = epsilon_topo for z < z_rec (post-recombination)
 *   epsilon = transition for z ~ z_rec
 *   epsilon = 0 for z > z_rec (pre-recombination, preserves BBN)
 */

#ifndef __OCTH__
#define __OCTH__

#include <math.h>

/* Mathematical constants */
#define OCTH_PI 3.14159265358979323846
#define OCTH_TWO_PI (2.0 * OCTH_PI)

/* PIT-derived fundamental parameter */
#define OCTH_EPSILON_GEOMETRIC (1.0 / OCTH_TWO_PI)  /* ~ 0.159, derived from topology */

/* OCTH default parameters */
#define OCTH_DEFAULT_AMPLITUDE   0.0023
#define OCTH_DEFAULT_Z_TOPO      0.35
#define OCTH_DEFAULT_SIGMA       0.12

/* PIT profile parameters */
#define OCTH_PIT_EPSILON_TOPO    0.14     /* Adjusted for H0 = 72.7 km/s/Mpc */
#define OCTH_PIT_Z_REC           1089.0   /* Recombination redshift */
#define OCTH_PIT_SIGMA_REC       50.0     /* Width of transition at recombination */

/* Extended OCTH profile parameters (legacy) */
#define OCTH_EXTENDED_EPS_PEAK      0.175
#define OCTH_EXTENDED_EPS_TAIL      0.080
#define OCTH_EXTENDED_Z_TAIL_START  500.0
#define OCTH_EXTENDED_Z_TAIL_WIDTH  200.0

/**
 * OCTH profile types
 */
typedef enum {
  octh_off = 0,       /* Standard LCDM (no OCTH modification) */
  octh_gaussian,      /* Simple Gaussian profile (legacy) */
  octh_extended,      /* Extended profile with CMB-era tail (legacy) */
  octh_pit            /* PIT-derived profile (recommended) */
} octh_profile_type;

/**
 * OCTH parameters structure
 */
struct octh_parameters {
  int is_enabled;            /* 0 = LCDM, 1 = OCTH active */
  octh_profile_type profile; /* Which profile to use */

  /* Legacy parameters */
  double amplitude;          /* A parameter for gaussian */
  double z_topo;             /* Topological redshift for gaussian */
  double sigma;              /* Width of Gaussian */

  /* Extended profile (legacy) */
  double eps_peak;           /* Peak amplitude for extended */
  double eps_tail;           /* Tail amplitude for CMB era */
  double z_tail_start;       /* Where tail starts */
  double z_tail_width;       /* Width of tail transition */

  /* PIT-derived profile parameters */
  double epsilon_topo;       /* Topological epsilon (can adjust from 1/2pi) */
  double z_rec;              /* Recombination redshift */
  double sigma_rec;          /* Width of recombination transition */
};

/**
 * Initialize OCTH parameters with defaults
 */
static inline void octh_init_defaults(struct octh_parameters * pocth) {
  pocth->is_enabled = 0;  /* LCDM by default */
  pocth->profile = octh_off;

  /* Legacy parameters */
  pocth->amplitude = OCTH_DEFAULT_AMPLITUDE;
  pocth->z_topo = OCTH_DEFAULT_Z_TOPO;
  pocth->sigma = OCTH_DEFAULT_SIGMA;

  /* Extended profile */
  pocth->eps_peak = OCTH_EXTENDED_EPS_PEAK;
  pocth->eps_tail = OCTH_EXTENDED_EPS_TAIL;
  pocth->z_tail_start = OCTH_EXTENDED_Z_TAIL_START;
  pocth->z_tail_width = OCTH_EXTENDED_Z_TAIL_WIDTH;

  /* PIT-derived defaults */
  pocth->epsilon_topo = OCTH_PIT_EPSILON_TOPO;
  pocth->z_rec = OCTH_PIT_Z_REC;
  pocth->sigma_rec = OCTH_PIT_SIGMA_REC;
}

/**
 * Compute Psi(z) - the OCTH metric modification function
 *
 * For LCDM: Psi(z) = 1 (no modification)
 * For OCTH: Psi(z) < 1, causing H_OCTH > H_LCDM
 *
 * PIT Profile Physics:
 * -------------------
 * The topological transition occurs at z_rec (recombination).
 * BEFORE recombination (z > z_rec): Universe is in "pre-transition" state, Psi = 1
 * AFTER recombination (z < z_rec): Topological effect propagates, Psi < 1
 *
 * This is consistent with PIT axiom: topology changes BEFORE metric changes.
 * The recombination marks when matter topology "freezes out".
 *
 * @param pocth  Pointer to OCTH parameters
 * @param z      Redshift
 * @return       Psi(z) value (always between 0 and 1)
 */
static inline double octh_psi(struct octh_parameters * pocth, double z) {

  double psi;
  double epsilon;

  /* If OCTH is disabled, return 1 (standard LCDM) */
  if (!pocth->is_enabled || pocth->profile == octh_off) {
    return 1.0;
  }

  if (pocth->profile == octh_gaussian) {
    /* Simple Gaussian profile (legacy) */
    double dz = z - pocth->z_topo;
    double gaussian = exp(-dz * dz / (2.0 * pocth->sigma * pocth->sigma));
    psi = sqrt(1.0 - pocth->amplitude * gaussian);
  }
  else if (pocth->profile == octh_extended) {
    /* Extended profile with CMB-era tail (legacy) */
    double dz = z - pocth->z_topo;

    /* Low-z Gaussian peak */
    double peak = pocth->eps_peak * exp(-dz * dz / (2.0 * pocth->sigma * pocth->sigma));

    /* High-z tail that smoothly activates around z_tail_start */
    double tail_transition = 0.5 * (1.0 + tanh((z - pocth->z_tail_start) / pocth->z_tail_width));
    double tail = pocth->eps_tail * tail_transition;

    /* Combine: use maximum of peak and tail */
    epsilon = (peak > tail) ? peak : tail;

    /* Ensure stability */
    if (epsilon >= 1.0) epsilon = 0.999;
    if (epsilon < 0.0) epsilon = 0.0;

    psi = sqrt(1.0 - epsilon);
  }
  else if (pocth->profile == octh_pit) {
    /**
     * PIT-derived profile
     *
     * Physical interpretation:
     * - z > z_rec: Before topological transition, epsilon = 0
     * - z ~ z_rec: Transition occurs (smooth via tanh)
     * - z < z_rec: After transition, epsilon = epsilon_topo (constant)
     *
     * The effect at z=0 is a CONSEQUENCE of the transition at z_rec,
     * propagated via the modified expansion history.
     */

    /* Transition function: 0 for z >> z_rec, 1 for z << z_rec */
    /* Note: tanh((z_rec - z)/sigma) goes from -1 to +1 as z goes from high to low */
    double transition = 0.5 * (1.0 + tanh((pocth->z_rec - z) / pocth->sigma_rec));

    /* epsilon(z) = epsilon_topo * transition */
    /* This gives epsilon = 0 at z >> z_rec, and epsilon = epsilon_topo at z << z_rec */
    epsilon = pocth->epsilon_topo * transition;

    /* Ensure stability */
    if (epsilon >= 1.0) epsilon = 0.999;
    if (epsilon < 0.0) epsilon = 0.0;

    psi = sqrt(1.0 - epsilon);
  }
  else {
    psi = 1.0;
  }

  return psi;
}

/**
 * Compute dPsi/dz - derivative of Psi with respect to z
 *
 * Needed for correct H' calculation
 *
 * @param pocth  Pointer to OCTH parameters
 * @param z      Redshift
 * @return       dPsi/dz
 */
static inline double octh_dpsi_dz(struct octh_parameters * pocth, double z) {

  if (!pocth->is_enabled || pocth->profile == octh_off) {
    return 0.0;
  }

  double psi = octh_psi(pocth, z);
  double depsilon_dz;

  if (pocth->profile == octh_gaussian) {
    double dz = z - pocth->z_topo;
    double gaussian = exp(-dz * dz / (2.0 * pocth->sigma * pocth->sigma));
    /* d(epsilon)/dz = -A * gaussian * (z - z_topo) / sigma^2 */
    depsilon_dz = -pocth->amplitude * gaussian * dz / (pocth->sigma * pocth->sigma);
  }
  else if (pocth->profile == octh_pit) {
    /* epsilon = epsilon_topo * 0.5 * (1 + tanh((z_rec - z)/sigma)) */
    /* d(epsilon)/dz = epsilon_topo * 0.5 * (-1/sigma) * sech^2((z_rec - z)/sigma) */
    double arg = (pocth->z_rec - z) / pocth->sigma_rec;
    double sech_arg = 1.0 / cosh(arg);
    depsilon_dz = -pocth->epsilon_topo * 0.5 * sech_arg * sech_arg / pocth->sigma_rec;
  }
  else {
    depsilon_dz = 0.0;
  }

  /* dPsi/dz = d(sqrt(1-epsilon))/dz = -depsilon_dz / (2 * sqrt(1-epsilon)) */
  /* = -depsilon_dz / (2 * psi) */
  if (psi > 0.001) {
    return -depsilon_dz / (2.0 * psi);
  }
  return 0.0;
}

/**
 * Compute H_OCTH from H_LCDM
 *
 * H_OCTH = H_LCDM / Psi(z)
 *
 * Since Psi < 1, this gives H_OCTH > H_LCDM
 */
static inline double octh_modify_H(struct octh_parameters * pocth, double H_lcdm, double z) {
  double psi = octh_psi(pocth, z);
  return H_lcdm / psi;
}

/**
 * Compute H'_OCTH from H_LCDM and H'_LCDM
 *
 * H' = dH/d(conformal time) = dH/dz * dz/dtau
 *
 * With OCTH: H_octh = H_lcdm / psi
 * dH_octh/dz = (1/psi) * dH_lcdm/dz - (H_lcdm/psi^2) * dpsi/dz
 *
 * And dz/dtau = -a * H * (1+z) = -(1+z)^2 * H
 *
 * So H' = dH/dz * [-(1+z)^2 * H]
 */
static inline double octh_modify_H_prime(struct octh_parameters * pocth,
                                          double H_lcdm, double H_prime_lcdm,
                                          double H_octh, double z) {
  if (!pocth->is_enabled || pocth->profile == octh_off) {
    return H_prime_lcdm;
  }

  double psi = octh_psi(pocth, z);
  double dpsi_dz = octh_dpsi_dz(pocth, z);

  /* H'_lcdm comes from CLASS, it's d(H_lcdm)/d(conf_time) */
  /* We need to convert to dH/dz then apply OCTH modification */

  /* For now, approximate: H'_octh / H_octh ~ H'_lcdm / H_lcdm */
  /* This is valid when dpsi/dz is small */
  double H_prime_octh = H_prime_lcdm / psi;

  /* Correction term from psi variation */
  /* dH_octh/dtau = d(H_lcdm/psi)/dtau = H'_lcdm/psi - H_lcdm * dpsi/dtau / psi^2 */
  /* dpsi/dtau = dpsi/dz * dz/dtau = dpsi/dz * [-(1+z)^2 * H] */
  double a = 1.0 / (1.0 + z);
  double dz_dtau = -(1.0 + z) * (1.0 + z) * H_octh;
  double dpsi_dtau = dpsi_dz * dz_dtau;

  H_prime_octh = H_prime_lcdm / psi - H_lcdm * dpsi_dtau / (psi * psi);

  return H_prime_octh;
}

/**
 * Calculate theoretical epsilon from PIT
 *
 * epsilon_geometric = 1/(2*pi) ~ 0.159
 *
 * This is the fundamental value derived from the topological
 * structure of the cosmological horizon as the H1 cycle.
 */
static inline double octh_pit_epsilon_geometric(void) {
  return OCTH_EPSILON_GEOMETRIC;
}

/**
 * Calculate H0 boost factor from epsilon_topo
 *
 * H0_OCTH = H0_LCDM / Psi(z=0)
 *
 * At z=0, if z << z_rec, then Psi(0) ~ sqrt(1 - epsilon_topo)
 * So H0_OCTH/H0_LCDM = 1/sqrt(1 - epsilon_topo)
 *
 * For epsilon_topo = 0.14: factor = 1.079, so H0 = 72.7 km/s/Mpc
 * For epsilon_topo = 0.159 (geometric): factor = 1.091, so H0 = 73.5 km/s/Mpc
 */
static inline double octh_h0_boost_factor(double epsilon_topo) {
  if (epsilon_topo >= 1.0) return 100.0;  /* Unphysical */
  if (epsilon_topo <= 0.0) return 1.0;    /* No effect */
  return 1.0 / sqrt(1.0 - epsilon_topo);
}

#endif /* __OCTH__ */
