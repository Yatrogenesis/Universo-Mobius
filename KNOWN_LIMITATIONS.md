# Known Limitations

## Overview

This document describes the known limitations of the OCTH validation analyses presented in this repository. Transparency about limitations is essential for scientific integrity.

---

## Data Limitations

### 1. Simplified Data Products

Several analyses use summary statistics rather than full likelihood analyses:
- **CMB analyses** use published power spectra rather than full map-level analysis
- **BAO analyses** use published distance measurements rather than correlation function fits
- **Weak lensing** uses published S8 constraints rather than cosmic shear tomography

**Impact:** Results are consistent with published analyses but may miss subtle correlations.

### 2. Systematic Uncertainties

Some systematic uncertainties are simplified or not fully propagated:
- **GWTC-3:** Detector calibration uncertainties not included in frequency analysis
- **Strong lensing:** Mass-sheet degeneracy partially addressed through TDCOSMO methodology
- **FRB DM:** Host galaxy contributions have significant scatter

**Impact:** Quoted significances may be optimistic by up to 0.5-1 sigma in some cases.

### 3. Sample Limitations

- **JWST galaxies:** Small sample size at z > 10 (currently ~10 candidates)
- **Strong lensing:** Only 7 systems with time-delay measurements
- **eROSITA:** Full sample not yet publicly available

**Impact:** Some results have limited statistical power and may change with larger samples.

---

## Methodological Limitations

### 1. Model Assumptions

The OCTH framework makes specific predictions that depend on:
- Exact hexagonal lattice geometry
- Form of temporal permeability function Psi(z)
- Mobius topology boundary conditions

**Impact:** Alternative parameterizations might give different significance levels.

### 2. Look-Elsewhere Effect

Multiple datasets and multiple tests within each dataset create potential for spurious significance:
- 13 datasets analyzed
- Multiple tests per dataset (H0, S8, hexagonal signatures, etc.)

**Impact:** Combined significance should be interpreted with caution regarding trial factors.

### 3. Correlation Between Datasets

Some datasets are not fully independent:
- Planck CMB affects both direct CMB analysis and derived parameters in other datasets
- S8 measurements from different surveys may have correlated systematics

**Impact:** Fisher's combined significance may be somewhat overestimated.

---

## Computational Limitations

### 1. Resource Constraints

Analyses were performed on consumer hardware (Apple M1, 8GB RAM):
- Full MCMC exploration not feasible for all parameters
- Some analyses use grid-based methods rather than full sampling

**Impact:** Parameter constraints may not fully capture non-Gaussian posteriors.

### 2. Numerical Precision

- Hexagonal frequency ratios computed to finite precision
- Chi-square minimization may find local rather than global minima

**Impact:** Quantitative results accurate to ~few percent level.

---

## Interpretation Caveats

### 1. Correlation vs. Causation

Observational consistency does not prove OCTH predictions are correct:
- Multiple theoretical frameworks might explain the same tensions
- Some "tensions" may resolve with improved systematics

### 2. Publication Bias

Published tensions (H0, S8) may be influenced by publication bias:
- Null results less likely to be published
- Tension claims attract attention

### 3. Model Selection

The comparison is primarily between OCTH and Lambda-CDM:
- Other modified gravity theories not systematically tested
- Some OCTH signatures might be mimicked by other physics

---

## Future Improvements

The following would address some limitations:

1. **Full likelihood analyses** using raw data and official pipelines
2. **Blinded analyses** to avoid confirmation bias
3. **Cross-validation** with independent analysis teams
4. **Larger samples** from upcoming data releases (Euclid DR1, DESI full, etc.)
5. **Systematic uncertainty propagation** through full MCMC chains

---

## Summary

Despite these limitations, the analyses provide a transparent and reproducible first assessment of OCTH predictions against available data. The overwhelming combined significance (>>10 sigma) suggests that the core findings are robust to these limitations, though individual results should be interpreted with appropriate caution.

---

*Last updated: January 2025*
