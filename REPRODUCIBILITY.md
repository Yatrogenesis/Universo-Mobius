# Reproducibility Guide

## Overview

All analyses in this repository are fully reproducible using publicly available data and open-source Python libraries.

---

## Requirements

### Python Environment

```bash
python >= 3.8
```

### Dependencies

```bash
pip install numpy scipy matplotlib h5py requests astropy
```

### Optional (for extended analyses)

```bash
pip install healpy camb emcee corner
```

---

## Data Sources and DOIs

All data used in this analysis are publicly available:

| Dataset | Source | Reference |
|---------|--------|-----------|
| GWTC-3 | https://gwosc.org/ | LIGO/Virgo/KAGRA Collaboration 2023 |
| NANOGrav 15yr | https://zenodo.org/record/7967584 | DOI: 10.5281/zenodo.7967584 |
| ACT DR6 | https://act.princeton.edu/ | Qu et al. 2024 |
| SPT-3G | https://pole.uchicago.edu/ | Balkenhol et al. 2023 |
| DESI DR1 | https://data.desi.lbl.gov/ | DESI Collaboration 2024 |
| Euclid ERO | https://euclid.esac.esa.int/ | Euclid Collaboration 2024 |
| JWST COSMOS-Web | https://archive.stsci.edu/ | Casey et al. 2023 |
| Pantheon+ | https://pantheonplussh0es.github.io/ | Scolnic et al. 2022 |
| Planck 2018 | https://pla.esac.esa.int/ | Planck Collaboration 2020 |
| SDSS/BOSS | https://www.sdss.org/dr18/ | BOSS Collaboration |
| DES Y3 | https://des.ncsa.illinois.edu/ | DES Collaboration 2022 |
| eROSITA | https://erosita.mpe.mpg.de/ | eROSITA Collaboration 2024 |
| CHIME/FRB | https://www.chime-frb.ca/ | CHIME/FRB Collaboration |
| TDCOSMO | https://shsuyu.github.io/H0LiCOW/ | TDCOSMO Collaboration |

---

## Running the Analyses

### Individual Datasets

Each analysis script is self-contained:

```bash
# Example: GWTC-3 analysis
python code/GWTC3_full_raw_analysis.py

# Example: NANOGrav analysis
python code/NANOGrav_analysis.py
```

### Full Validation Suite

```bash
# Run all analyses
for script in code/*.py; do
    python "$script"
done
```

---

## Output Structure

```
results/
    full_raw_analysis/     # GWTC-3 results
    nanograv/              # NANOGrav results
    cmb_act_spt/           # ACT/SPT CMB results
    desi_dr1/              # DESI results
    euclid/                # Euclid results
    jwst/                  # JWST results
    pantheon/              # Pantheon+ results
    planck_cmb/            # Planck results
    sdss_boss/             # SDSS/BOSS results
    des_y3/                # DES results
    erosita/               # eROSITA results
    chime_frb/             # CHIME results
    strong_lensing/        # Strong lensing results

figures/
    full_raw/              # GWTC-3 figures
    nanograv/              # NANOGrav figures
    ...                    # etc.
```

---

## Statistical Methods

### Hexagonal Fit (GWTC-3)

The hexagonal frequency structure is tested using chi-square analysis comparing:
- OCTH hexagonal ratios: 1 : sqrt(3) : 2 : sqrt(7)
- Standard GR predictions

### Fisher's Method (Combined Significance)

Independent p-values are combined using Fisher's method:

```
chi2_combined = -2 * sum(log(p_i))
```

This follows a chi-square distribution with 2k degrees of freedom.

### Tension Calculations

Parameter tensions are computed as:

```
tension = |x1 - x2| / sqrt(sigma1^2 + sigma2^2)
```

---

## Verification

All results can be independently verified by:

1. Downloading raw data from the sources listed above
2. Running the analysis scripts
3. Comparing output JSON files with reported values

---

## Citation

If using this analysis methodology, please cite:

```
Molina-Burgos, F. (2026). OCTH Validation Against Astronomical Datasets.
GitHub: https://github.com/Yatrogenesis/Universo-Mobius
```

---

## Contact

Repository: https://github.com/Yatrogenesis/Universo-Mobius
