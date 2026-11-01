# Correction notice (2026-11-01)

**Author:** Francisco Molina-Burgos, Avermex Research Division

This notice supersedes the statistical claims of release v1.0.0 (DOI 10.5281/zenodo.18408178). An independent re-analysis with correct null distributions, together with a blind test against data published after the January 2026 predictions, shows that none of the significance claims hold. They should not be cited as evidence for OCTH or for a Möbius topology of the Universe. The original material remains in the repository as a historical record.

## 1. Claims of v1.0.0 and their status

| Claim (v1.0.0) | Status | Reason |
|---|---|---|
| CMB antipodal anti-correlation -0.048, Z = -6.08 (Planck SMICA) | Not significant | The value -0.0480227 reproduces exactly. The significance used Z = corr x sqrt(N_pix), which treats correlated pixels as independent, and a pixel-shuffling null. Against 2000 isotropic LCDM skies (Planck 2018 best-fit spectrum, same beam, pixel window and resolution): Z = -1.68, p = 0.041 (one-sided). With the galactic-cut scan replayed inside every simulation (look-elsewhere): p = 0.10. WMAP ILC: Z = -2.53, p = 0.0025 |
| The antipodal statistic is a Möbius signature | Withdrawn | On the full sky, the antipodal correlation equals sum (-1)^l (2l+1) C_l / sum (2l+1) C_l. That is the known low-l parity asymmetry (Planck 2018 VII). Verified: -0.052 from the pseudo-C_l against -0.048 direct (Planck); -0.093 against -0.092 (WMAP). It is not specific to any topology |
| GWTC-3 raw data: 80/81 events favour hexagonal modes, Z = 6.69 | Artefact of construction | The analysed band is [0.3 f1, 4 f1]. The GR template (32000/M x {1, 1.47, 1.87, 2.24}) lies entirely outside it, so the GR chi^2 has a data-independent lower bound of 268.35; the minimum observed over 80 events is 268.49. The hexagonal template falls inside the band and wins on pure noise as well. Recurring peaks are instrumental: H1 21.9 Hz, L1 60 Hz mains, V1 100/150 Hz, and the O1 calibration lines at H1 36.7 Hz and L1 34.7 Hz (GWOSC O1 instrumental lines; arXiv:1801.07204). The README states 33 events; the JSON contains 80 |
| CMB x LIGO cross-correlation, Z = 4.31 | Invalid method | Hand-entered, rounded positions treated as points (90 % areas up to 7461 deg²); isotropic null despite the anisotropic network sensitivity; antipodal anomaly axes counted as separate points; four statistics tried and the largest reported |
| "Galactic poles 60°/90°, Z = 5.0" | Unverified | The code that produces it was not found in the repository |
| Psi formulation reproduces Schwarzschild | Inconsistent | Psi = sqrt(1 - rho/rho_Planck) is identified with sqrt(1 - r_s/r). At the Earth's surface the two differ by 84 orders of magnitude (1e-93 against 1.4e-9). With flat g_ij the metric gives gamma_PPN = 0 (half the GR light deflection), excluded by Cassini (|gamma - 1| < 2.3e-5). `OCTH_solar_system.py` does not compute gamma |
| Hubble tension explained by a local void (results/OCTH_hubble_tension_v2.json, "100 % agreement") | Withdrawn | The code first derives Psi = 1/sqrt(1 + beta delta) (line 88), which gives Psi > 1 in a void and therefore a LOWER local H0. It then switches to Psi = (1 + delta)^alpha (line 178), which reverses the sign, and solves alpha from the very ratio it claims to explain (line 215). The "100 % agreement" is true by construction. Data test (correction/h0): with Pantheon+SH0ES distances, H0 inside z < 0.07 is 72.50 +- 0.30 and outside z > 0.1 is 72.93 +- 0.22. There is no step, while the void mechanism needs a drop of 5.64 km/s/Mpc (diagonal errors only) |

## 2. Blind test of the O4 predictions (gist of 2026-01-09)

The predictions were published on 2026-01-09 (gist 0ddbb1a15441ca2e6c998ae7b1690163, single revision). GWTC-5.0 (O4b, 161 candidates) appeared on GWOSC on 2026-05-22 and on Zenodo on 2026-05-26, after the predictions. The test was preregistered before the data were opened (correction/ligo_o4b/PREREGISTRO_LIGO_O4B.md, in Spanish).

| Prediction | O4b result | Verdict |
|---|---|---|
| At least 60 % of well-localised events (90 % area < 500 deg²) within 30° of four CMB anomaly positions | 8/40 = 20 %; isotropic expectation 25 %; Z = -0.72 | Fails. Meets the gist's own falsification criterion (Z < 2) |
| Ringdown frequencies f1 x {1, sqrt3, 2, sqrt7}, f1 = 2200/M | LVK ringdown fits (GWTC-5.0 TGR release, 5 events, 55 fits): f_220 between 110 and 210 Hz; 0/55 medians inside any OCTH band, whose highest mode lies at 60-96 Hz | Fails |
| Signal at the hexagonal frequencies in > 60 % of events | 24/35 events show excess power there, but broadband 30-300 Hz power shows excess in 33/35. The statistic detects the presence of a compact-binary signal, not hexagonal structure | Not discriminating |
| The ratios sqrt3 : 2 : sqrt7 are special (exploratory specificity test) | Against 500 random ratio combs in the same on-source spectrum: 0/35 events with p < 0.05; median p = 0.47 | No evidence |

## 3. What remains open

- Non-trivial, including non-orientable, spatial topology is a legitimate open question. The relevant framework is that of the COMPACT collaboration (Copi et al., arXiv:2510.05030; Samandar et al., arXiv:2407.09400). Matched-circle searches already exclude a nearest clone closer than about 0.985 of the diameter of the last-scattering surface, for both orientations (Vaudrevange et al. 2012, arXiv:1206.2939). The statistics of this repository do not test that question.
- The CMB parity asymmetry is real but known (about 2-3 sigma). Topology is one possible origin among others.

## 4. Reproduction

- `correction/cmb/`: the oracle for the antipodal statistic and results with 2000 LCDM simulations.
- `correction/ligo_o4b/`: preregistration, scripts and results of the O4b test.
- `correction/h0/`: the Hubble-bubble test with Pantheon+SH0ES.
- `correction/AUDIT_2026-10-31_es.md`: the full audit (in Spanish).
