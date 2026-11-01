#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Francisco Molina-Burgos, Avermex Research Division
"""LT3 of PREREGISTRO_LIGO_O4B.md: LVK-measured ringdown frequencies (GWTC-5.0 TGR release, RD)
against the OCTH set f1 x {1, sqrt3, 2, sqrt7}, f1 = 2200/M (+-5 %).

Every start time in the release is reported (0-20 M after peak); no selection by outcome.
Usage: lt3_ringdown.py <ringdown_dir> <PESummaryTable.hdf5> <out.json>
"""
import glob, json, os, re, sys
import h5py
import numpy as np

RD, SUMMARY, OUT = sys.argv[1:4]
RATIOS = np.array([1.0, np.sqrt(3), 2.0, np.sqrt(7)])

t = h5py.File(SUMMARY, 'r')['summary_info'][:]
pe = {}
for r in sorted(t, key=lambda r: r['result_label']):
    pe.setdefault(r['superevent_id'].decode(), r)


def inside(f, M):
    return bool(np.any(np.abs(f / (RATIOS * 2200.0 / M) - 1) <= 0.05))


rows = []
for sdir in sorted(glob.glob(os.path.join(RD, 'S2*'))):
    sid = os.path.basename(sdir)
    r = pe.get(sid)
    if r is None:
        rows.append(dict(superevent=sid, status='not_in_PE_summary'))
        continue
    M = float(r['total_mass_source_median'])
    Mdet = M * (1 + float(r['redshift_median']))
    for run in sorted(glob.glob(os.path.join(sdir, 'production', '*'))):
        m = re.fullmatch(r'([0-9+]+)_start0M_stop20M_step2M_T[0-9.]+s_ds4', os.path.basename(run))
        if m is None:
            continue  # injections (_INJ), parametrised-deviation runs and noise-only runs are not measurements
        modes = m.group(1).split('+')
        engines = sorted(glob.glob(os.path.join(run, 'engine', '*', 'result.nc')))
        for k, eng in enumerate(engines):
            with h5py.File(eng, 'r') as f:
                fs = f['posterior/f'][:].reshape(-1, len(modes))
            for j, mode in enumerate(modes):
                q = np.percentile(fs[:, j], [5, 50, 95])
                rows.append(dict(superevent=sid, gw_name=r['gw_name'].decode(), run=os.path.basename(run),
                                 start_index=k, n_starts=len(engines), mode=mode,
                                 f_hz_5_50_95=q.tolist(), M_src=M, M_det=Mdet,
                                 ocht_src_hz=(RATIOS * 2200.0 / M).tolist(),
                                 ocht_det_hz=(RATIOS * 2200.0 / Mdet).tolist(),
                                 median_in_ocht_band_src=inside(q[1], M),
                                 median_in_ocht_band_det=inside(q[1], Mdet),
                                 frac_posterior_in_ocht_band_src=float(np.mean([inside(x, M) for x in fs[::10, j]]))))

f220 = [x for x in rows if x.get('mode') == '220' and x['run'].startswith('220_')]
res = dict(n_rows=len(rows),
           summary_220_only=dict(n_fits=len(f220),
                                 n_median_in_ocht_src=sum(x['median_in_ocht_band_src'] for x in f220),
                                 n_median_in_ocht_det=sum(x['median_in_ocht_band_det'] for x in f220)),
           rows=rows)
json.dump(res, open(OUT, 'w'), indent=1)
print(json.dumps(res['summary_220_only'], indent=1))
for x in f220:
    if x['start_index'] in (0, x['n_starts'] // 2, x['n_starts'] - 1):
        print(x['gw_name'], 'start', x['start_index'], 'f220 =', [round(v) for v in x['f_hz_5_50_95']],
              'OCTH src', [round(v) for v in x['ocht_src_hz']])
