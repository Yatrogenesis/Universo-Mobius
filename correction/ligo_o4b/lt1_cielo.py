#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Francisco Molina-Burgos, Avermex Research Division
"""LT1 of PREREGISTRO_LIGO_O4B.md: OCTH gist sky-clustering prediction against GWTC-5.0 (O4b).

Usage: lt1_cielo.py <skymap_dir> <PESummaryTable.hdf5> <out.json>
"""
import glob, json, os, re, sys
import healpy as hp
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u
from scipy.stats import binomtest

SKY, SUMMARY, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
CAP = np.radians(30.0)
gal = SkyCoord(l=220 * u.deg, b=-20 * u.deg, frame='galactic').icrs
TARGETS = {  # gist section 2.4, equatorial (ICRS) degrees
    'SGP_region': (192.75, -27.0), 'Cold_Spot': (48.75, -19.0),
    'Quadrupole_axis': (165.0, 25.0), 'Hemispherical_axis': (float(gal.ra.deg), float(gal.dec.deg))}
TV = np.array([hp.ang2vec(ra, dec, lonlat=True) for ra, dec in TARGETS.values()])

# Isotropic fraction of the sky inside the union of the four caps (NSIDE 256).
NS = 256
pv = np.array(hp.pix2vec(NS, np.arange(hp.nside2npix(NS)))).T
f0 = float(np.mean((pv @ TV.T).max(axis=1) >= np.cos(CAP)))


def moc_pixels(path):
    t = Table.read(path)
    uniq = np.asarray(t['UNIQ'], dtype=np.int64)
    order = (np.log2(uniq).astype(int) // 2) - 1
    ipix = uniq - 4 * (4 ** order.astype(np.int64))
    dens = np.asarray(t['PROBDENSITY'], dtype=float)
    vec = np.zeros((len(uniq), 3))
    area = np.zeros(len(uniq))
    for o in np.unique(order):
        s = order == o
        ns = 2 ** int(o)
        vec[s] = np.array(hp.pix2vec(ns, ipix[s], nest=True)).T
        area[s] = hp.nside2pixarea(ns)
    return vec, dens, area


events = {}
for f in sorted(glob.glob(os.path.join(SKY, '*_Skymap_PEDataRelease.fits.gz'))):
    ev = re.search(r'-(GW\d{6}_\d{6})-', os.path.basename(f)).group(1)
    events.setdefault(ev, f)  # alphabetical first map per event (preregistered rule)

rows = []
for ev, f in sorted(events.items()):
    vec, dens, area = moc_pixels(f)
    prob = dens * area
    order = np.argsort(dens)[::-1]
    cum = np.cumsum(prob[order])
    a90 = float(np.degrees(np.degrees(area[order][: np.searchsorted(cum, 0.9 * cum[-1]) + 1].sum())))
    vmax = vec[np.argmax(dens)]
    cosd = vec @ TV.T
    in_caps = cosd.max(axis=1) >= np.cos(CAP)
    rows.append(dict(event=ev, map=os.path.basename(f), area90_deg2=a90,
                     maxpix_min_sep_deg=float(np.degrees(np.arccos(np.clip((TV @ vmax).max(), -1, 1)))),
                     prob_in_caps=float(prob[in_caps].sum() / prob.sum())))

well = [r for r in rows if r['area90_deg2'] < 500]
k1 = sum(r['maxpix_min_sep_deg'] < 30 for r in well)
k2 = sum(r['prob_in_caps'] > 0.5 for r in well)
n = len(well)


def summary(k):
    frac = k / n if n else float('nan')
    z = (k - n * f0) / np.sqrt(n * f0 * (1 - f0)) if n else float('nan')
    return dict(k=k, n=n, fraction=frac, Z_vs_isotropic=float(z),
                p_binomial_greater=float(binomtest(k, n, f0, alternative='greater').pvalue) if n else None,
                gist_prediction_met_ge_0p60=bool(frac >= 0.60), gist_falsified_Z_lt_2=bool(z < 2))


res = dict(targets_icrs_deg=TARGETS, cap_deg=30, f0_isotropic=f0, n_events_with_map=len(rows),
           n_well_localized_lt_500deg2=n, primary_maxpix=summary(k1), secondary_prob_gt_half=summary(k2),
           events=rows)
json.dump(res, open(OUT, 'w'), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != 'events'}, indent=1))
