/**
 * @file octh.h
 *
 * OCTH (Observational Cosmological Topological Hypothesis) modification
 * Part of Universo-Mobius project
 *
 * Author: Francisco Molina-Burgos
 * Avermex Research Division
 * January 2026
 *
 * Extended profile matches CLASS-RS implementation exactly
 */

#ifndef __OCTH__
#define __OCTH__

#include <math.h>

#define OCTH_PI 3.14159265358979323846
#define OCTH_TWO_PI (2.0 * OCTH_PI)

/* Extended profile parameters (CLASS-RS compatible) */
#define OCTH_EXT_EPS_PEAK      0.175
#define OCTH_EXT_EPS_TAIL      0.080
#define OCTH_EXT_Z_TOPO        1089.0
#define OCTH_EXT_Z_TAIL_START  500.0
#define OCTH_EXT_TAIL_WIDTH    100.0
#define OCTH_EXT_SIGMA         0.6

/* Legacy defaults */
#define OCTH_DEFAULT_AMPLITUDE 0.0023
#define OCTH_DEFAULT_Z_TOPO    0.35
#define OCTH_DEFAULT_SIGMA     0.12

#define OCTH_PIT_EPSILON_TOPO  0.14
#define OCTH_PIT_Z_REC         1089.0
#define OCTH_PIT_SIGMA_REC     50.0

typedef enum {
  octh_off = 0,
  octh_gaussian,
  octh_extended,
  octh_pit
} octh_profile_type;

struct octh_parameters {
  int is_enabled;
  octh_profile_type profile;

  /* Legacy Gaussian */
  double amplitude;
  double z_topo;
  double sigma;

  /* Extended (CLASS-RS) */
  double eps_peak;
  double eps_tail;
  double z_tail_start;
  double tail_width;

  /* PIT */
  double epsilon_topo;
  double z_rec;
  double sigma_rec;
};

static inline void octh_init_defaults(struct octh_parameters * pocth) {
  pocth->is_enabled = 0;
  pocth->profile = octh_off;

  pocth->amplitude = OCTH_DEFAULT_AMPLITUDE;
  pocth->z_topo = OCTH_EXT_Z_TOPO;
  pocth->sigma = OCTH_EXT_SIGMA;

  pocth->eps_peak = OCTH_EXT_EPS_PEAK;
  pocth->eps_tail = OCTH_EXT_EPS_TAIL;
  pocth->z_tail_start = OCTH_EXT_Z_TAIL_START;
  pocth->tail_width = OCTH_EXT_TAIL_WIDTH;

  pocth->epsilon_topo = OCTH_PIT_EPSILON_TOPO;
  pocth->z_rec = OCTH_PIT_Z_REC;
  pocth->sigma_rec = OCTH_PIT_SIGMA_REC;
}

/**
 * Compute Psi(z)
 *
 * Extended profile (CLASS-RS compatible):
 *   Psi(z) = 1 - peak - tail
 *   peak = eps_peak * exp(-0.5 * ((ln(a) - ln(a_topo))/sigma)^2)
 *   tail = eps_tail * 0.5 * (1 - tanh((z - z_tail)/width))
 *
 * At z=0: Psi ~ 0.92 -> H0 = 73.3 km/s/Mpc
 */
static inline double octh_psi(struct octh_parameters * pocth, double z) {

  if (!pocth->is_enabled || pocth->profile == octh_off) {
    return 1.0;
  }

  double psi = 1.0;

  if (pocth->profile == octh_gaussian) {
    /* Legacy Gaussian */
    double dz = z - pocth->z_topo;
    double gaussian = exp(-dz * dz / (2.0 * pocth->sigma * pocth->sigma));
    psi = sqrt(1.0 - pocth->amplitude * gaussian);
  }
  else if (pocth->profile == octh_extended) {
    /* CLASS-RS compatible Extended profile */
    double a = 1.0 / (1.0 + z);
    double ln_a = log(a);
    double a_topo = 1.0 / (1.0 + pocth->z_topo);
    double ln_a_topo = log(a_topo);

    /* Peak in ln(a) space */
    double delta = (ln_a - ln_a_topo) / pocth->sigma;
    double peak = pocth->eps_peak * exp(-0.5 * delta * delta);

    /* Tail active for z < z_tail_start */
    double tail_arg = (z - pocth->z_tail_start) / pocth->tail_width;
    double tail = pocth->eps_tail * 0.5 * (1.0 - tanh(tail_arg));

    /* Additive combination */
    psi = 1.0 - peak - tail;

    if (psi < 0.01) psi = 0.01;
    if (psi > 1.0) psi = 1.0;
  }
  else if (pocth->profile == octh_pit) {
    /* PIT profile */
    double transition = 0.5 * (1.0 + tanh((pocth->z_rec - z) / pocth->sigma_rec));
    double epsilon = pocth->epsilon_topo * transition;
    if (epsilon >= 0.999) epsilon = 0.999;
    if (epsilon < 0.0) epsilon = 0.0;
    psi = sqrt(1.0 - epsilon);
  }

  return psi;
}

static inline double octh_dpsi_dz(struct octh_parameters * pocth, double z) {
  if (!pocth->is_enabled || pocth->profile == octh_off) return 0.0;
  double dz = 0.01;
  return (octh_psi(pocth, z + dz) - octh_psi(pocth, z - dz)) / (2.0 * dz);
}

static inline double octh_modify_H(struct octh_parameters * pocth, double H_lcdm, double z) {
  return H_lcdm / octh_psi(pocth, z);
}

static inline double octh_modify_H_prime(struct octh_parameters * pocth,
                                          double H_lcdm, double H_prime_lcdm,
                                          double H_octh, double z) {
  if (!pocth->is_enabled || pocth->profile == octh_off) return H_prime_lcdm;
  double psi = octh_psi(pocth, z);
  double dpsi_dz = octh_dpsi_dz(pocth, z);
  double dz_dtau = -(1.0 + z) * (1.0 + z) * H_octh;
  double dpsi_dtau = dpsi_dz * dz_dtau;
  return H_prime_lcdm / psi - H_lcdm * dpsi_dtau / (psi * psi);
}

#endif /* __OCTH__ */
