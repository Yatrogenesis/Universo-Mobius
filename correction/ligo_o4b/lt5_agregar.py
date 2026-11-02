#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Francisco Molina-Burgos, Avermex Research Division
"""Aggregate LT5 (verifier quality) and amended LT4 (decision) per amendments B1/B2."""
import glob, json, sys
import numpy as np
from scipy.stats import binomtest, combine_pvalues
rows = [json.load(open(f)) for f in sorted(glob.glob(sys.argv[1] + '/GW*.json'))]
units = [(r['event'], d, r[d]) for r in rows for d in ('H1', 'L1') if r[d].get('status') == 'ok']
rho = {}
for _, _, u in units:
    for k, v in u['rho'].items():
        a = rho.setdefault(k, [0, 0, 0])
        a[0] += v['detected']; a[1] += v['n']; a[2] += v['detected_with_comb']
fa = [u['noise_false_alarm'] for _, _, u in units]
fa_pa = sum(x['pa_lt_0p05'] for x in fa); fa_n = sum(x['n'] for x in fa)
ev = []
for r in rows:
    ps = [r[d]['lt4_amended']['p_offsource'] for d in ('H1', 'L1') if r[d].get('status') == 'ok']
    pc = [r[d]['lt4_amended']['p_comb'] for d in ('H1', 'L1') if r[d].get('status') == 'ok']
    if ps:
        p = ps[0] if len(ps) == 1 else combine_pvalues(ps, method='fisher').pvalue
        q = pc[0] if len(pc) == 1 else combine_pvalues(pc, method='fisher').pvalue
        ev.append(dict(event=r['event'], snr=r['snr'], p_offsource=float(p), p_comb=float(q)))
k = sum(e['p_offsource'] < 0.05 for e in ev); n = len(ev)
kc = sum(e['p_comb'] < 0.05 for e in ev)
pb = binomtest(k, n, 0.05, alternative='greater').pvalue
det = {k_: dict(rate=v[0] / v[1], rate_with_comb=v[2] / v[1], n=v[1]) for k_, v in rho.items()}
informative = any(v['rate'] >= 0.90 for k_, v in det.items() if float(k_) <= 12) and fa_pa / fa_n <= 0.05
verdict = ('DESCARTADO' if pb > 0.01 else ('SOSTENIDO' if k / n >= 0.60 else 'EXCESO_SIN_UMBRAL')) if informative else 'NO_INFORMATIVO'
res = dict(units=len(units), verifier=dict(detection_by_rho=det, false_alarm_pa=fa_pa / fa_n, false_alarm_n=fa_n,
                                             informative=bool(informative)),
           lt4_amended=dict(n_events=n, k_p_offsource_lt_0p05=k, frac=k / n, p_binom_vs_5pct=float(pb),
                            fisher_all=float(combine_pvalues([e['p_offsource'] for e in ev], method='fisher').pvalue),
                            k_p_comb_lt_0p05_secondary=kc, verdict=verdict,
                            top_snr=sorted([(e['event'], round(e['snr'], 1), round(e['p_offsource'], 3)) for e in ev], key=lambda x: -x[1])[:8]))
json.dump(dict(res, events=ev), open(sys.argv[2], 'w'), indent=1)
print(json.dumps(res, indent=1))
