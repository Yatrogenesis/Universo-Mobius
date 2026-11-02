#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Francisco Molina-Burgos, Avermex Research Division
"""LT4 of PREREGISTRO_LIGO_O4B.md: post-merger-only hexagonal detector with off-source and random-comb nulls.

For each selected event and detector: on-source window [tc-0.1, tc+0.4] s, PSD by Welch on
+-64 s excluding +-4 s, statistic S = mean over modes of the mean whitened power in f_n (1 +- 0.05),
null from 200 off-source windows of the same length in the same 4096-s file(s), > 8 s from tc and
inside CBC_CAT1 data. Positive control: same machinery in the GR f_220 band (1 +- 0.15).
Usage: lt2_potencia.py <PESummaryTable.hdf5> <events_GWTC-5.0.json> <strain_dir> <out_dir>
"""
import json, os, sys, urllib.request
import h5py
import numpy as np
from scipy.signal import butter, sosfiltfilt, welch
from scipy.signal.windows import tukey
from scipy.stats import combine_pvalues

SUMMARY, EVENTS, STRAIN, OUT = sys.argv[1:5]
os.makedirs(STRAIN, exist_ok=True)
os.makedirs(OUT, exist_ok=True)
RATIOS = np.array([1.0, np.sqrt(3), 2.0, np.sqrt(7)])
F220 = 1.5251 - 1.1568 * (1 - 0.69) ** 0.1292  # Berti, Cardoso & Will 2006 fit, spin fixed 0.69
TSUN_HZ = 1 / (2 * np.pi * 4.925490947e-6)  # c^3 / (2 pi G Msun) in Hz
NOFF, SEED = 200, 20261101
W0, W1 = -0.1, 0.4


def load_events():
    t = h5py.File(SUMMARY, 'r')['summary_info'][:]
    ev = {}
    for r in sorted(t, key=lambda r: r['result_label']):
        ev.setdefault(r['gw_name'].decode(), r)  # alphabetical-first PE result (preregistered)
    gps = {e['name']: e['gps'] for e in json.load(open(EVENTS))}
    sel = []
    for name, r in sorted(ev.items()):
        if r['mass_2_source_median'] > 3 and r['network_matched_filter_snr_median'] >= 12 and name in gps:
            sel.append(dict(event=name, gps=gps[name], snr=float(r['network_matched_filter_snr_median']),
                            M_src=float(r['total_mass_source_median']), z=float(r['redshift_median']),
                            Mf_src=float(r['final_mass_source_median'])))
    return sel


def fetch(det, start):
    d = start // 1048576 * 1048576
    fn = f"{det[0]}-{det}_GWOSC_O4b_4KHZ_R1-{start}-4096.hdf5"
    path = os.path.join(STRAIN, fn)
    if not os.path.exists(path):
        url = f"https://gwosc.org/archive/data/O4b_4KHZ_R1/{d}/{fn}"
        try:
            urllib.request.urlretrieve(url, path + '.part')
        except Exception:
            return None
        os.replace(path + '.part', path)
    return path


def load_strain(det, gps):
    """Strain and 1-Hz CBC_CAT1 validity covering gps +- 2100 s (one or two 4096-s files)."""
    s0 = int(gps) // 4096 * 4096
    starts = [s0]
    if gps - s0 < 80:
        starts.insert(0, s0 - 4096)
    if s0 + 4096 - gps < 80:
        starts.append(s0 + 4096)
    hs, qs, t0, paths = [], [], None, []
    for st in starts:
        p = fetch(det, st)
        if p is None:
            return None
        paths.append(p)
        with h5py.File(p, 'r') as f:
            hs.append(f['strain/Strain'][:])
            dq = f['quality/simple/DQmask'][:]
            names = [n.decode() if isinstance(n, bytes) else n for n in f['quality/simple/DQShortnames'][:]]
            bits = (1 << names.index('DATA')) | (1 << names.index('CBC_CAT1'))
            qs.append((dq & bits) == bits)
            if t0 is None:
                t0 = f['meta/GPSstart'][()]
    return np.concatenate(hs), np.concatenate(qs), float(t0), paths


def band_power(seg, fs, psd_f, psd, bands):
    w = tukey(len(seg), 0.2)
    X = np.fft.rfft(seg * w)
    f = np.fft.rfftfreq(len(seg), 1 / fs)
    P = np.abs(X) ** 2 / np.interp(f, psd_f, psd) / (fs * np.sum(w ** 2) / 2)
    vals = [P[(f >= a) & (f <= b)].mean() for a, b in bands if np.any((f >= a) & (f <= b))]
    return float(np.mean(vals)) if vals else np.nan


def analyse(e, det, rng):
    got = load_strain(det, e['gps'])
    if got is None:
        return dict(status='no_data')
    h, ok, t0, paths = got
    fs = 4096.0
    tc = e['gps'] - t0
    n = int(round((W1 - W0) * fs))
    i_on = int(round((tc + W0) * fs))
    sec_ok = lambda a, b: ok[int(np.floor(a)):int(np.ceil(b))].all() and np.isfinite(h[int(a * fs):int(b * fs)]).all()
    if not sec_ok(tc - 64, tc + 64):
        return dict(status='dq_fail_around_event')
    # Amendment A2: high-pass (8th-order Butterworth, 15 Hz, zero phase) before cutting windows, so
    # sub-15 Hz seismic power cannot leak through the window sidelobes into the analysis bands.
    h = sosfiltfilt(butter(8, 15.0, 'highpass', fs=fs, output='sos'), np.nan_to_num(h))
    pa = welch(h[int((tc - 64) * fs):int((tc - 4) * fs)], fs, nperseg=int(4 * fs))
    pb = welch(h[int((tc + 4) * fs):int((tc + 64) * fs)], fs, nperseg=int(4 * fs))
    psd_f, psd = pa[0], 0.5 * (pa[1] + pb[1])
    def bands(rs, M):
        out = []
        for r in rs:
            fr = r * 2200.0 / M
            if 20 <= fr <= 1000:
                hw = max(0.05 * fr, 4.0)
                out.append((fr - hw, fr + hw))
        return out

    def seg_power(i0, nn, bl):
        seg = h[i0:i0 + nn]
        w = tukey(len(seg), 0.2)
        X = np.fft.rfft(seg * w)
        f = np.fft.rfftfreq(len(seg), 1 / fs)
        P = np.abs(X) ** 2 / np.interp(f, psd_f, psd) / (fs * np.sum(w ** 2) / 2)
        vals = [P[(f >= a) & (f <= b)].mean() for a, b in bl if np.any((f >= a) & (f <= b))]
        return float(np.mean(vals)) if vals else np.nan

    npost = int(round(0.25 * fs))
    i_post = int(round((tc + 0.05) * fs))
    i_pre = int(round((tc - 0.25) * fs))
    tot = len(h) / fs
    starts = []
    tries = 0
    while len(starts) < NOFF and tries < 50000:
        tries += 1
        s0 = rng.uniform(1, tot - 2)
        if abs(s0 - tc) > 8 and sec_ok(s0 - 4, s0 + 4):
            starts.append(s0)
    res = dict(status='ok', n_off=len(starts))
    for name, M in (('src', e['M_src']), ('det', e['M_src'] * (1 + e['z']))):
        bl = bands(RATIOS, M)
        if not bl:
            res[name] = None
            continue
        on = seg_power(i_post, npost, bl)
        off = np.array([seg_power(int(s0 * fs), npost, bl) for s0 in starts])
        ctrl = np.column_stack([np.ones(500), rng.uniform(1.3, 2.9, size=(500, 3))])
        sc = np.array([seg_power(i_post, npost, bands(c, M)) for c in ctrl])
        sc = sc[np.isfinite(sc)]
        res[name] = dict(S_post=on, p_offsource=float((1 + np.sum(off >= on)) / (1 + len(off))),
                         p_random_comb=float((1 + np.sum(sc >= on)) / (1 + len(sc))), bands=bl)
    bb = [(30.0, 300.0)]
    on = seg_power(i_pre, npost, bb)
    off = np.array([seg_power(int(s0 * fs), npost, bb) for s0 in starts])
    res['control_pre_merger_30_300'] = dict(S=on, p=float((1 + np.sum(off >= on)) / (1 + len(off))))
    return res


def main():
    sel = load_events()
    json.dump(sel, open(os.path.join(OUT, 'seleccion.json'), 'w'), indent=1)
    for e in sel:
        out = os.path.join(OUT, e['event'] + '.json')
        if os.path.exists(out):
            continue
        r = dict(e)
        for j, det in enumerate(('H1', 'L1')):
            rng = np.random.default_rng([SEED + 2, int(e['event'][2:8] + e['event'][9:15]), j])
            try:
                r[det] = analyse(e, det, rng)
            except Exception as x:  # record, never hide
                r[det] = dict(status='error', error=repr(x)[:300])
        json.dump(r, open(out, 'w'), indent=1)
        for f in os.listdir(STRAIN):  # keep disk bounded
            os.remove(os.path.join(STRAIN, f))
        print(e['event'], {d: r[d].get('status') for d in ('H1', 'L1')}, flush=True)
    print('FIN', flush=True)


if __name__ == '__main__':
    main()
