#!/usr/bin/env python3
"""MR summary: multiple-testing control, sensitivity concordance and axis contrasts.

Primary family = univariable IVW of every consortium / UK Biobank exposure on the eight FinnGen outcomes (mrlib.TO_FINNGEN
x FinnGen); Benjamini-Hochberg FDR over that family. For each pair: weighted-median and MR-Egger estimates, Egger intercept
p, Cochran Q p, and replication in UK Biobank-based outcomes where the design allows (FinnGen exposures and exposures
without UK Biobank samples). Axis contrast per exposure: inverse-variance pooled log-odds over the tenosynovial outcomes
(CTS, TF, DQ) minus that over Dupuytren and frozen shoulder; FinnGen M13 outcomes share controls, so the contrast SE
ignores their covariance and is labelled approximate. Tables: mr_summary, mr_axis_contrast.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import numpy as np
import pandas as pd
from project import Project
from scipy import stats

P = Project(__file__)
TEN, FIB = ['CTS', 'TF', 'DQ'], ['DUP', 'FS']


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1)
    return out


with P.connect(read_only=True) as con:
    r = con.execute('select * from mr_results').df()
with P.run(__file__, params={'fdr': 'BH over primary IVW family (exposures x 8 FinnGen outcomes)'}) as run:
    w = r.pivot_table(index=['exposure', 'outcome', 'cohort', 'axis', 'n_snp'], columns='method',
                      values=['estimate', 'se', 'p'], aggfunc='first')
    w.columns = [f'{a}_{b}'.replace(' ', '_').replace('-', '') for a, b in w.columns]
    w = w.reset_index()
    q = r[r.method == 'IVW'][['exposure', 'outcome', 'cohort', 'heterogeneity_p']]
    w = w.merge(q, on=['exposure', 'outcome', 'cohort'], how='left')
    prim = w.cohort == 'FinnGen'
    w['family'] = np.where(prim, 'primary', 'replication')
    w['fdr'] = np.nan
    w.loc[prim, 'fdr'] = bh(w.loc[prim, 'p_IVW'])
    w['or_ivw'] = np.exp(w.estimate_IVW)
    w['or_lo'] = np.exp(w.estimate_IVW - 1.96 * w.se_IVW); w['or_hi'] = np.exp(w.estimate_IVW + 1.96 * w.se_IVW)
    w['concordant'] = (np.sign(w.estimate_weighted_median) == np.sign(w.estimate_IVW)) & \
                      (np.sign(w.estimate_MREgger) == np.sign(w.estimate_IVW))
    run.table('mr_summary', w, note='Univariable MR per exposure x outcome x cohort: IVW OR (95% CI), BH-FDR over the primary '
                                    'FinnGen family, weighted median / MR-Egger estimates, Egger intercept p, Cochran Q p')
    rows = []
    for e, d in w[prim].groupby('exposure'):
        def pool(names):
            x = d[d.outcome.isin(names)]
            wt = 1 / x.se_IVW ** 2
            return (wt * x.estimate_IVW).sum() / wt.sum(), np.sqrt(1 / wt.sum())
        (bt, st), (bf, sf) = pool(TEN), pool(FIB)
        z = (bt - bf) / np.hypot(st, sf)
        rows.append(dict(exposure=e, tenosynovial_logor=bt, tenosynovial_se=st, fibrotic_logor=bf, fibrotic_se=sf,
                         diff=bt - bf, diff_se_approx=np.hypot(st, sf), z=z, p_approx=2 * stats.norm.sf(abs(z))))
    ax = pd.DataFrame(rows)
    run.table('mr_axis_contrast', ax, note='Per-exposure IVW log-odds pooled over tenosynovial (CTS, TF, DQ) vs fibrotic (DUP, FS) '
                                           'FinnGen outcomes; difference with approximate SE (covariance of shared controls ignored)')
    show = w[prim & (w.fdr < 0.05)].sort_values('fdr')[['exposure', 'outcome', 'n_snp', 'or_ivw', 'or_lo', 'or_hi', 'fdr',
                                                        'estimate_weighted_median', 'p_MREgger_intercept', 'concordant']]
    print(show.round(3).to_string(index=False))
    print(ax.round(3).to_string(index=False))
