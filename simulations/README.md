# MÖBIUS Cosmological Simulations

This directory contains documentation and configuration files for reproducing
the MÖBIUS cosmological simulations.

## Overview

The simulations validate the MÖBIUS cosmological modification by:

1. Running the CLASS Boltzmann solver with OCTH modifications
2. Generating initial conditions for N-body simulations
3. Running MP-GADGET with modified expansion history
4. Performing MCMC parameter estimation with Planck data

## Directory Structure

```
simulations/
├── class/
│   ├── octh.h              # OCTH header file
│   ├── background_mod.c    # Modified background.c section
│   ├── input_mod.c         # Modified input.c section
│   └── mobius_test.ini     # Configuration file
├── gadget/
│   ├── mobius_gadget.param # Simulation parameters
│   └── cosmology_mod.c     # Hubble modification
├── cobaya/
│   └── planck_full.yaml    # MCMC configuration
└── results/
    ├── pk_z0.txt           # Power spectrum at z=0
    └── mcmc_summary.txt    # MCMC results summary
```

## Installation

### Requirements

- Ubuntu 20.04+ (or WSL2 on Windows)
- GCC 9+
- MPI (OpenMPI or MPICH)
- GSL (GNU Scientific Library)
- FFTW3
- Python 3.12+
- pip packages: cobaya, classy, numpy, scipy

### Step 1: Install CLASS with OCTH

```bash
# Clone CLASS v3.2.3
git clone -b v3.2.3 https://github.com/lesgourg/class_public.git class
cd class

# Copy OCTH modifications
cp ../simulations/class/octh.h include/
# Apply modifications to background.c and input.c (see code_mods/)

# Compile
make -j4
```

### Step 2: Install MP-GADGET

```bash
# Clone MP-GADGET
git clone https://github.com/MP-Gadget/MP-Gadget.git gadget
cd gadget

# Build
make -j4

# The MÖBIUS modification is in the expansion rate calculation
# See gadget/cosmology_mod.c for the Ψ(z) implementation
```

### Step 3: Install Cobaya and classy

```bash
python -m venv venv
source venv/bin/activate
pip install cobaya classy numpy scipy

# Install Planck likelihoods
cobaya-install planck_2018_lowl.TT
cobaya-install planck_2018_highl_plik.TTTEEE_lite
```

## Running Simulations

### CLASS Power Spectrum

```bash
cd class
./class simulations/class/mobius_test.ini
```

Output: `output/mobius_test_pk.dat`

### N-body Simulation

```bash
cd gadget

# Generate initial conditions
mpirun -np 2 ./genic/MP-GenIC genic_params.param

# Run simulation
mpirun -np 2 ./gadget/MP-Gadget mobius_gadget.param
```

Output: `output/powerspectrum-*.txt`, `output/PART_*/`, `output/PIG_*/`

### MCMC Parameter Estimation

```bash
source venv/bin/activate
cobaya-run simulations/cobaya/planck_full.yaml
```

Output: `chains/mobius_planck_full.*`

## Results Summary

### N-body Simulation (January 17, 2026)

- Box: 50 Mpc/h
- Particles: 64³ = 262,144
- Redshift range: z=99 to z=0
- Power spectra: 157 files
- Snapshots: 41 files
- FOF catalogs: 69 files
- Total output: 113 MB

### MCMC Results (Planck 2018)

- Samples: 500 accepted
- Acceptance rate: 50.6%
- H₀ = 69.36 ± 0.09 km/s/Mpc (ΛCDM baseline)
- With MÖBIUS: H₀ → 73.1 km/s/Mpc

## Code Modifications

### CLASS: octh.h

The temporal permeability function:

```c
static inline double octh_psi(struct octh_parameters * pocth, double z) {
  if (!pocth->is_enabled) return 1.0;

  double tail_arg = (z - pocth->z_tail_start) / pocth->tail_width;
  double tail = pocth->eps_tail * 0.5 * (1.0 - tanh(tail_arg));
  double psi = 1.0 - tail;

  return psi;
}
```

### CLASS: background.c

The Hubble parameter modification:

```c
double octh_psi_value = 1.0;
if (pba->octh.is_enabled) {
  double z = 1.0/a - 1.0;
  octh_psi_value = octh_psi(&(pba->octh), z);
}
pvecback[pba->index_bg_H] = sqrt(rho_tot - pba->K/a/a) / octh_psi_value;
```

## Reproducibility

All results can be reproduced using the configuration files in this directory.
The random seed for initial conditions is fixed for reproducibility.

## Author

Francisco Molina-Burgos
Avermex Research Division
Mérida, Yucatán, México

January 2026
