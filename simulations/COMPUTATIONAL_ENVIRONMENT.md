# Computational Environment

This document describes the hardware and software environment where the
MÖBIUS cosmological simulations were performed.

## Hardware Specifications

| Component | Specification |
|-----------|---------------|
| **CPU** | 12th Gen Intel Core i7-12650H |
| **Cores** | 10 cores (6P + 4E), 16 threads |
| **Base Clock** | 2.3 GHz (boost to 4.7 GHz) |
| **RAM** | 8 GB DDR5 |
| **Storage** | 1 TB NVMe SSD |
| **GPU** | Intel Iris Xe (integrated) |

## Software Environment

| Software | Version |
|----------|---------|
| **Host OS** | Windows 11 |
| **WSL** | WSL2 (Windows Subsystem for Linux) |
| **Linux Distribution** | Ubuntu 24.04.3 LTS |
| **Kernel** | Linux (WSL2 kernel) |
| **Python** | 3.12.3 |
| **GCC** | 13.2.0 |
| **MPI** | OpenMPI 4.1.6 |

## Cosmological Codes

| Code | Version | Purpose |
|------|---------|---------|
| **CLASS** | 3.2.3 + OCTH | Boltzmann solver |
| **MP-GADGET** | 5.0.1 | N-body simulation |
| **Cobaya** | 3.6 | MCMC sampler |
| **classy** | 3.3.4.0 | CLASS Python wrapper |
| **CAMB** | 1.5.8 | Alternative Boltzmann |

## Libraries

| Library | Version | Purpose |
|---------|---------|---------|
| **GSL** | 2.7.1 | Scientific computing |
| **FFTW3** | 3.3.10 | Fast Fourier transforms |
| **HDF5** | 1.14.3 | Data storage |
| **NumPy** | 1.26.4 | Numerical arrays |
| **SciPy** | 1.12.0 | Scientific computing |

## Simulation Performance

### N-body (GADGET-MÖBIUS)

- Configuration: 50 Mpc/h box, 64³ particles
- MPI processes: 2
- OpenMP threads: 2 per process
- Wall time: ~5 minutes (z=99 to z=0)
- Output: 113 MB

### MCMC (Cobaya)

- Configuration: 500 samples, Planck likelihoods
- CPU time: ~180 minutes
- Acceptance rate: 50.6%

## Location

The simulations were performed at:

**Avermex Research Division**
Mérida, Yucatán, México

## Date

January 17, 2026

## Notes

1. All simulations were run on a consumer-grade laptop, demonstrating
   that the MÖBIUS validation can be reproduced without specialized
   HPC resources.

2. The N-body simulation used a modest resolution (64³) for proof of
   concept. Production runs would require higher resolution and larger
   boxes.

3. WSL2 provides near-native Linux performance, making it suitable for
   scientific computing on Windows systems.

## Reproducibility

To reproduce these results:

1. Install WSL2 with Ubuntu 24.04
2. Install required packages: `apt install build-essential gsl-bin libgsl-dev libfftw3-dev openmpi-bin libhdf5-dev`
3. Create Python virtual environment with required packages
4. Follow instructions in `simulations/README.md`

## Author

Francisco Molina-Burgos
fmolina@avermex.com
