#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Francisco Molina-Burgos, Avermex Research Division
"""Aggregate LT2 per PREREGISTRO_LIGO_O4B.md: per-event Fisher over detectors, fraction p<0.05, global Fisher."""
import glob, json, sys
import numpy as np
from scipy.stats import combine_pvalues, binomtest
rows = [json.load(open(f)) for f in sorted(glob.glob(sys.argv[1] + '/GW*.json'))]
out = {}
for test in ('ocht_src', 'ocht_det', 'gr_f220', 'diag_broadband_30_300'):
    ev = []
    for r in rows:
        ps = [r[d][test]['p'] for d in ('H1', 'L1') if r[d].get('status') == 'ok' and r[d].get(test)]
        if ps:
            p = ps[0] if len(ps) == 1 else combine_pvalues(ps, method='fisher').pvalue
            ev.append(dict(event=r['event'], snr=r['snr'], M=r['M_src'], p=float(p)))
    p = np.array([e['p'] for e in ev])
    k = int(np.sum(p < 0.05)); n = len(p)
    out[test] = dict(n_events=n, n_p_lt_0p05=k, frac=k / n if n else None,
                     p_binom_vs_5pct=float(binomtest(k, n, 0.05, alternative='greater').pvalue) if n else None,
                     fisher_all=float(combine_pvalues(p, method='fisher').pvalue) if n else None,
                     corr_p_mass=float(np.corrcoef(p, [e['M'] for e in ev])[0, 1]) if n > 2 else None,
                     top_snr=[(e['event'], round(e['snr'], 1), round(e['p'], 3)) for e in sorted(ev, key=lambda e: -e['snr'])[:6]])
out['status_counts'] = {s: sum(r[d].get('status') == s for r in rows for d in ('H1', 'L1')) for s in
                        set(r[d].get('status') for r in rows for d in ('H1', 'L1'))}
json.dump(out, open(sys.argv[2], 'w'), indent=1)
print(json.dumps(out, indent=1))
