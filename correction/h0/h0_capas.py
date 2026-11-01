#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Francisco Molina-Burgos, Avermex Research Division
"""Hubble-bubble test of the OCTH local-void H0 explanation with Pantheon+SH0ES distances.

If H_obs = H_true / Psi_local inside a void of radius ~300 Mpc (z ~ 0.07), H0 inferred from
SH0ES-calibrated distances must fall from ~73 to ~67.4 beyond the void edge. Flat LCDM shape with
Omega_m = 0.334 (Pantheon+ cosmology fit); diagonal errors only (correlations ignored).
Usage: h0_capas.py <Pantheon+SH0ES.dat> <out.json>
"""
import json, sys
import numpy as np
from scipy.integrate import quad

d = np.genfromtxt(sys.argv[1], names=True, dtype=None, encoding=None)
C = 299792.458
OM = 0.334
sel = (d['IS_CALIBRATOR'] == 0) & (d['zHD'] > 0.01)
z, mu, err = d['zHD'][sel], d['MU_SH0ES'][sel], d['MU_SH0ES_ERR_DIAG'][sel]
zh = d['zHEL'][sel]
dc = np.array([quad(lambda x: 1 / np.sqrt(OM * (1 + x) ** 3 + 1 - OM), 0, zi)[0] for zi in z])
dl_unit = (1 + zh) * dc * C  # luminosity distance times H0, in km/s
dl = 10 ** ((mu - 25) / 5)  # Mpc
h0 = dl_unit / dl
sig = h0 * np.log(10) / 5 * err
edges = [0.01, 0.023, 0.04, 0.07, 0.1, 0.15, 0.25, 0.4, 0.7, 1.0, 2.3]
rows = []
for a, b in zip(edges[:-1], edges[1:]):
    m = (z >= a) & (z < b)
    if m.sum() < 5:
        continue
    w = 1 / sig[m] ** 2
    mean = float(np.sum(w * h0[m]) / np.sum(w))
    rows.append(dict(z=[a, b], n=int(m.sum()), H0=mean, err=float(1 / np.sqrt(np.sum(w)))))
inner = [r for r in rows if r['z'][1] <= 0.07]
outer = [r for r in rows if r['z'][0] >= 0.1]
wa = lambda rs: (sum(r['H0'] / r['err'] ** 2 for r in rs) / sum(1 / r['err'] ** 2 for r in rs),
                 1 / np.sqrt(sum(1 / r['err'] ** 2 for r in rs)))
hi, ei = wa(inner)
ho, eo = wa(outer)
res = dict(bins=rows, H0_inside_z_lt_0p07=[hi, ei], H0_outside_z_gt_0p1=[ho, eo],
           step=hi - ho, step_sigma=(hi - ho) / np.hypot(ei, eo),
           ocht_void_required_step=73.04 - 67.4, note='diagonal errors only; SN covariance ignored')
json.dump(res, open(sys.argv[2], 'w'), indent=1)
for r in rows:
    print(f"z {r['z'][0]:.3f}-{r['z'][1]:.3f} n={r['n']:4d} H0 = {r['H0']:.2f} +- {r['err']:.2f}")
print(f"inside z<0.07: {hi:.2f} +- {ei:.2f}; outside z>0.1: {ho:.2f} +- {eo:.2f}; step {hi-ho:.2f} ({(hi-ho)/np.hypot(ei,eo):.1f} sigma); OCTH needs 5.64")
