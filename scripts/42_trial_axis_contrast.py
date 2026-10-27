#!/usr/bin/env python3
"""Axis-level effects of incretin receptor agonists in the serious-event tables of placebo-controlled trials
(trial_axis_cells, 40_trial_axis.py).

Axes are the core disorders of each genetic factor: tenosynovial = carpal tunnel syndrome + trigger finger,
fibro-capsular = frozen shoulder + Dupuytren's disease. Axis effect = Peto one-step OR pooled over strata
(trial x disorder); axis contrast = ratio of the two ORs with a z test (its trial-level bootstrap interval is in
43_trial_robustness.py). Table: trial_axis_contrast.
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
AXES = {'tenosynovial': ['CTS', 'TF'], 'fibrocapsular': ['FS', 'DUP']}
RANK = {'incretin': 0}


def peto_info(t):
    """Peto O - E and information V summed over strata (rows a1, n1, a0, n0)."""
    n = t.n1 + t.n0; m = t.a1 + t.a0
    e = m * t.n1 / n
    v = m * (n - m) * t.n1 * t.n0 / (n ** 2 * (n - 1))
    return float((t.a1 - e).sum()), float(v.sum())


with P.run(__file__, params={'axes': AXES, 'analysis': 'serious'}) as run:
    with P.connect(read_only=True) as con:
        cells = con.execute("select * from trial_axis_cells where analysis = 'serious'").df()
    run.input('trial_axis_cells', note='per-trial serious events, incretin vs placebo')
    cells = cells[(cells.a1 + cells.a0) > 0]
    rows = []
    for cls in RANK:
        for ax_, ds in AXES.items():
            t = cells[(cells.drug_class == cls) & cells.outcome.isin(ds)]
            oe, v = peto_info(t)
            lor, se = oe / v, 1 / np.sqrt(v)
            rows.append(dict(drug_class=cls, axis=ax_, n_strata=len(t), n_trials=t.nct.nunique(),
                             events_active=int(t.a1.sum()), events_placebo=int(t.a0.sum()), log_or=lor, se=se,
                             or_=np.exp(lor), or_lo=np.exp(lor - 1.96 * se), or_hi=np.exp(lor + 1.96 * se),
                             p=2 * stats.norm.sf(abs(lor / se))))
    a = pd.DataFrame(rows).rename(columns={'or_': 'peto_or'})
    w = a.pivot(index='drug_class', columns='axis')
    con_rows = []
    for cls in RANK:
        d = w.loc[cls, ('log_or', 'fibrocapsular')] - w.loc[cls, ('log_or', 'tenosynovial')]
        s = np.hypot(w.loc[cls, ('se', 'fibrocapsular')], w.loc[cls, ('se', 'tenosynovial')])
        con_rows.append(dict(drug_class=cls, axis='ratio_fibro_vs_teno', log_or=d, se=s, peto_or=np.exp(d),
                             or_lo=np.exp(d - 1.96 * s), or_hi=np.exp(d + 1.96 * s), p=2 * stats.norm.sf(abs(d / s))))
    a = pd.concat([a, pd.DataFrame(con_rows)], ignore_index=True)
    run.table('trial_axis_contrast', a, note='Axis-level Peto OR under incretin receptor agonists (serious events; strata = '
                                             'trial x disorder) and the fibro-capsular / tenosynovial ratio of ORs with z test')
    print(a.round(3).to_string(index=False))
