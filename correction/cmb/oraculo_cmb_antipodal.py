#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Francisco Molina-Burgos, Avermex Research Division
"""Oracle for Universo-Mobius Test #1 (antipodal CMB correlation).

Reproduces the repository statistic exactly (Pearson corr of T(n) vs T(-n) on the
ud_graded NSIDE=64 map) and replaces its null (corr*sqrt(npix), pixel shuffling) with
the correct one: isotropic Gaussian LCDM skies drawn from the Planck 2018 best-fit
TT spectrum, processed identically. Also replays the galactic-cut scan with the same
"pick the most negative Z" rule inside every simulation (look-elsewhere).

Usage: oraculo_cmb_antipodal.py <map.fits> <theory_cl.txt> <nsims> <out.json> [unit_scale]
"""
import json, sys
import healpy as hp
import numpy as np

MAP, CLF, NSIM, OUT = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
SCALE = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
NSIDE, NSIDE_SIM, BEAM_ARCMIN = 64, 256, 5.0
CUTS = list(range(5, 45, 5))  # same scan as galactic_poles_test


def antipodal_index(nside):
    th, ph = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    return hp.ang2pix(nside, np.pi - th, (ph + np.pi) % (2 * np.pi))


ANTI = antipodal_index(NSIDE)
_, B_LAT = hp.pix2ang(NSIDE, np.arange(hp.nside2npix(NSIDE)), lonlat=True)


def stats(m):
    full = np.corrcoef(m, m[ANTI])[0, 1]
    scan = []
    for c in CUTS:
        ok = (np.abs(B_LAT) >= c) & (np.abs(B_LAT[ANTI]) >= c)
        if ok.sum() > 5000:
            r = np.corrcoef(m[ok], m[ANTI][ok])[0, 1]
            scan.append((c, r, r * np.sqrt(ok.sum())))
    best = min(scan, key=lambda x: x[2])
    return full, best


def parity_form(m):
    """corr predicted from the pseudo-Cl alone: sum (-1)^l (2l+1) C_l / sum (2l+1) C_l."""
    cl = hp.anafast(m - m.mean(), lmax=3 * NSIDE - 1)
    l = np.arange(len(cl))
    w = (2 * l + 1) * cl
    w[0] = 0.0
    return float(np.sum(w * (-1.0) ** l) / np.sum(w))


m = hp.read_map(MAP, field=0) * SCALE
if hp.get_nside(m) != NSIDE:
    m = hp.ud_grade(m, NSIDE)
obs_full, obs_best = stats(m)
obs_parity = parity_form(m)

d = np.loadtxt(CLF)
ell, dtt = d[:, 0].astype(int), d[:, 1]
lmax = 3 * NSIDE_SIM - 1
cl = np.zeros(lmax + 1)
sel = (ell >= 2) & (ell <= lmax)
cl[ell[sel]] = dtt[sel] * 2 * np.pi / (ell[sel] * (ell[sel] + 1))  # uK^2
rng = np.random.default_rng(20261031)
sim_full, sim_bestz = [], []
for i in range(NSIM):
    np.random.seed(int(rng.integers(2**31)))
    s = hp.synfast(cl, NSIDE_SIM, lmax=lmax, fwhm=np.radians(BEAM_ARCMIN / 60), pixwin=True)
    s = hp.ud_grade(s, NSIDE)
    f, b = stats(s)
    sim_full.append(f)
    sim_bestz.append(b[2])
sim_full, sim_bestz = np.array(sim_full), np.array(sim_bestz)

res = dict(
    map=MAP, nside=NSIDE, nsims=NSIM, seed=20261031,
    observed=dict(corr_full=obs_full, corr_from_pseudoCl_parity=obs_parity,
                  repo_style_Z_full=obs_full * np.sqrt(len(m)),
                  best_cut=dict(cut=obs_best[0], corr=obs_best[1], repo_style_Z=obs_best[2])),
    null_LCDM=dict(corr_full_mean=float(sim_full.mean()), corr_full_std=float(sim_full.std()),
                   Z_true_full=float((obs_full - sim_full.mean()) / sim_full.std()),
                   p_one_sided_full=float((np.sum(sim_full <= obs_full) + 1) / (NSIM + 1)),
                   repo_style_bestZ_mean=float(sim_bestz.mean()),
                   repo_style_bestZ_frac_le_observed=float((np.sum(sim_bestz <= obs_best[2]) + 1) / (NSIM + 1))),
)
json.dump(res, open(OUT, "w"), indent=1)
print(json.dumps(res, indent=1))
