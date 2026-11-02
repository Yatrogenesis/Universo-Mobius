#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Francisco Molina-Burgos, Avermex Research Division
"""Aggregate LT6 per its preregistration (quality gate, then decision)."""
import glob, json, sys
from scipy.stats import binomtest, combine_pvalues
rows = [json.load(open(f)) for f in sorted(glob.glob(sys.argv[1] + '/GW*.json'))]
units = [r[d] for r in rows for d in ('H1', 'L1') if r[d].get('status') == 'ok']
fa_k = sum(u['false_alarm']['k'] for u in units); fa_n = sum(u['false_alarm']['n'] for u in units)
det = {k: sum(u['rho'][k]['detected'] for u in units) / sum(u['rho'][k]['n'] for u in units) for k in units[0]['rho']}
fa_p = binomtest(fa_k, fa_n, 0.05, alternative='greater').pvalue
quality = det['12'] >= 0.90 and fa_p > 0.01
ev = []
for r in rows:
    ps = [r[d]['p_on'] for d in ('H1', 'L1') if r[d].get('status') == 'ok']
    if ps:
        ev.append(dict(event=r['event'], snr=r['snr'], p=float(ps[0] if len(ps) == 1 else combine_pvalues(ps, method='fisher').pvalue)))
k = sum(e['p'] < 0.05 for e in ev); n = len(ev)
pb = binomtest(k, n, 0.05, alternative='greater').pvalue
if not quality: verdict = 'NO_INFORMATIVO'
elif pb > 0.01: verdict = 'DESCARTADO'
elif k / n >= 0.60: verdict = 'SOSTENIDO'
else: verdict = 'EXCESO_SIN_UMBRAL'
res = dict(units=len(units), detection=det, false_alarm=dict(k=fa_k, n=fa_n, rate=fa_k / fa_n, p_vs_5pct=float(fa_p)),
           quality_ok=bool(quality), n_events=n, k_p_lt_0p05=k, frac=k / n, p_binom=float(pb),
           fisher_all=float(combine_pvalues([e['p'] for e in ev], method='fisher').pvalue), verdict=verdict,
           status_counts={s: sum(r[d].get('status') == s for r in rows for d in ('H1', 'L1')) for s in {r[d].get('status') for r in rows for d in ('H1', 'L1')}},
           events=sorted([(e['event'], round(e['snr'], 1), round(e['p'], 3)) for e in ev], key=lambda x: -x[1]))
json.dump(res, open(sys.argv[2], 'w'), indent=1)
print(json.dumps(res, indent=1))
