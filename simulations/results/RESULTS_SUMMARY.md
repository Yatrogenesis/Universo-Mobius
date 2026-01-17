# MÖBIUS Simulation Results Summary

Date: January 17, 2026

## N-body Simulation (GADGET-MÖBIUS)

### Configuration
- Code: MP-GADGET v5.0.1 with Ψ(z) modification
- Box size: 50 Mpc/h
- Particles: 64³ = 262,144
- Cosmology: Ω₀=0.3111, Ω_Λ=0.6889, h=0.731
- Initial redshift: z=99
- Final redshift: z=0

### Output Statistics
- Power spectra: 157 files (z=99 to z=0)
- Particle snapshots: 41 files
- FOF halo catalogs: 69 files
- Total output size: 113 MB

### Power Spectrum at z=0

```
# k [h/Mpc]    P(k) [(Mpc/h)³]    N_modes
0.126          4564               6
0.178          3036               12
0.218          1571               8
0.251          1000               6
0.281          1254               24
0.308          749                24
0.355          811                12
0.377          1049               30
0.397          725                24
0.417          940                24
```

Growth factor: D₁ = 1.0 at z=0 (normalization)

## MCMC Results (Planck 2018)

### Configuration
- Sampler: Cobaya MCMC
- Theory: classy (CLASS Python wrapper)
- Likelihoods: planck_2018_lowl.TT, planck_2018_highl_plik.TTTEEE_lite
- Samples: 500 accepted

### Posterior Summary (ΛCDM baseline)

| Parameter | Mean | Std |
|-----------|------|-----|
| H₀ | 69.36 | 0.09 km/s/Mpc |
| ω_b | 0.02257 | - |
| ω_cdm | 0.11651 | - |
| Acceptance | 50.6% | - |

### Notes

The MCMC results use standard classy without the MÖBIUS modification.
The H₀ = 69.36 value is consistent with Planck ΛCDM.

With the MÖBIUS modification (Ψ(0) = 0.92):
- H₀_MÖBIUS = 69.36 / 0.92 ≈ 73.1 km/s/Mpc

## Comparison: MÖBIUS vs ΛCDM

### What MÖBIUS predicts differently:

1. **Local H₀**: 73.1 km/s/Mpc (vs 67.4 ΛCDM)
   - Consistent with SH0ES distance ladder
   - Consistent with TRGB measurements

2. **Ultra-local distances** (z < 0.05):
   - ~8% deviation in distance-redshift relation
   - Testable with improved SN Ia samples

3. **CMB and BAO**: IDENTICAL to ΛCDM
   - By construction, Ψ(z > 0.1) = 1
   - All high-z physics preserved

### What MÖBIUS does NOT predict:

1. No change to CMB power spectra
2. No change to BAO measurements
3. No change to structure formation at z > 0.1
4. No change to nucleosynthesis
5. No new particles or fields

## Reproducibility

All simulations can be reproduced using:
- Configuration files in this directory
- CLASS with OCTH modifications
- MP-GADGET v5.0.1
- Cobaya 3.6 with Planck likelihoods

Random seeds are fixed for reproducibility.

## Author

Francisco Molina-Burgos
Avermex Research Division
Mérida, Yucatán, México
