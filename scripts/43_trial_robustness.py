#!/usr/bin/env python3
"""Serious adverse events under incretin receptor agonists in placebo-controlled trials, with sensitivity analyses.

Population: every randomised trial on ClinicalTrials.gov with a placebo group and a group containing an incretin
receptor agonist (as 40_trial_axis.py) whose serious adverse event table was posted for both (number at risk > 0);
trials with zero events are kept, because serious tables list every event without a frequency threshold.
Per trial and outcome: participants affected per group (max over the set's terms; groups of one role summed), at risk =
serious-table denominators. Estimators per outcome:
  Peto one-step OR (strata with events), Mantel-Haenszel OR without continuity correction (Robins-Breslow-Greenland
  variance) and Mantel-Haenszel risk difference per 1,000 participants (Greenland-Robins variance), all trials as strata.
Outcomes add rotator cuff syndrome, the control events and two axis sets (participants with an event of either disorder
of the axis).
Axis contrast under incretins: ratio of Peto ORs (fibro-capsular FS + DUP vs tenosynovial CTS + TF) with a 95% interval
from 4,000 trial-level bootstrap resamples (both axes share each trial's placebo group).
Tables: trial_robustness, trial_axis_bootstrap, trial_incretin_molecules, trial_incretin_trials.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import duckdb
import numpy as np
import pandas as pd
import config as C
from project import Project
from scipy import stats

P = Project(__file__)
CLASSES = {
    'incretin': r'exenatide|liraglutide|lixisenatide|albiglutide|dulaglutide|semaglutide|efpeglenatide|taspoglutide|'
                r'tirzepatide|retatrutide|survodutide|mazdutide|orforglipron|cagrilintide|danuglipron|cotadutide',
}
OUTCOMES = {
    'CTS': ['carpal tunnel syndrome', 'carpal tunnel decompression'],
    'TF': ['trigger finger', 'tenosynovitis stenosans'],
    'DUP': ["dupuytren's contracture", "dupuytren's contracture operation"],
    'FS': ['periarthritis', 'frozen shoulder', 'adhesive capsulitis'],
    'POS_gallbladder': ['cholelithiasis', 'cholecystitis', 'cholecystitis acute'],
    'RC': ['rotator cuff syndrome'],
    'NEG_cataract': ['cataract'],
    'NEG_bcc': ['basal cell carcinoma'],
    'NEG_bph': ['benign prostatic hyperplasia'],
}
# axis sets for absolute effects: participants with an event of either disorder (maximum over the axis's terms per group)
OUTCOMES['AX_tenosynovial'] = OUTCOMES['CTS'] + OUTCOMES['TF']
OUTCOMES['AX_fibrocapsular'] = OUTCOMES['FS'] + OUTCOMES['DUP']
AXES = {'tenosynovial': ['CTS', 'TF'], 'fibrocapsular': ['FS', 'DUP']}


def peto(a1, n1, a0, n0):
    n = n1 + n0; m = a1 + a0
    k = m > 0
    e = m * n1 / n; v = m * (n - m) * n1 * n0 / (n ** 2 * (n - 1))
    oe, vv = (a1 - e)[k].sum(), v[k].sum()
    if vv <= 0:
        return np.nan, np.nan, np.nan, np.nan
    lor, se = oe / vv, 1 / np.sqrt(vv)
    return np.exp(lor), np.exp(lor - 1.96 * se), np.exp(lor + 1.96 * se), 2 * stats.norm.sf(abs(lor / se))


def mh_or(a1, n1, a0, n0):
    """Mantel-Haenszel OR, no continuity correction; Robins-Breslow-Greenland variance of log OR"""
    a, b, c, d = a1, n1 - a1, a0, n0 - a0
    N = n1 + n0
    R, S = a * d / N, b * c / N
    Pp, Q = (a + d) / N, (b + c) / N
    if R.sum() <= 0 or S.sum() <= 0:
        return np.nan, np.nan, np.nan, np.nan
    orr = R.sum() / S.sum()
    var = ((Pp * R).sum() / (2 * R.sum() ** 2) + (Pp * S + Q * R).sum() / (2 * R.sum() * S.sum())
           + (Q * S).sum() / (2 * S.sum() ** 2))
    se = np.sqrt(var); lor = np.log(orr)
    return orr, np.exp(lor - 1.96 * se), np.exp(lor + 1.96 * se), 2 * stats.norm.sf(abs(lor / se))


def mh_rd(a1, n1, a0, n0):
    """Mantel-Haenszel risk difference (drug minus placebo) per 1,000, Greenland-Robins variance"""
    N = n1 + n0
    w = n1 * n0 / N
    rd = (w * (a1 / n1 - a0 / n0)).sum() / w.sum()
    b1, b0 = n1 - a1, n0 - a0
    var = (a1 * b1 * n0 ** 3 + a0 * b0 * n1 ** 3).sum() / (n1 * n0 * N ** 2).sum() / 1.0
    var = ((a1 * b1 * n0 ** 3 + a0 * b0 * n1 ** 3) / (n1 * n0 * N ** 2)).sum() / (w.sum() ** 2)
    se = np.sqrt(var)
    return 1000 * rd, 1000 * (rd - 1.96 * se), 1000 * (rd + 1.96 * se), 2 * stats.norm.sf(abs(rd / se)) if se > 0 else np.nan


term2sets = {}
for k, ts in OUTCOMES.items():
    for t in ts:
        term2sets.setdefault(t, []).append(k)
with P.run(__file__, params={'classes': CLASSES, 'outcomes': OUTCOMES, 'bootstrap': 4000, 'seed': 20261025}) as run:
    con = duckdb.connect(str(C.WORK / 'ctgov' / 'ctgov.duckdb'), read_only=True)
    eg = con.execute("""
        SELECT e.nct, e.eg_id, e.role, e.ingredients, coalesce(e.n_risk_serious, 0) AS n_risk
        FROM ct_eventgroup e JOIN ct_trial t USING (nct)
        WHERE t.allocation = 'RANDOMIZED' AND e.role IN ('active', 'placebo')""").df()
    terms = sorted(term2sets)
    ae = con.execute(f"""
        SELECT nct, eg_id, lower(term) AS term, num_affected, num_at_risk FROM ct_ae_term
        WHERE kind = 'serious' AND lower(term) IN ({', '.join("'" + t.replace("'", "''") + "'" for t in terms)})""").df()
    con.close()
    run.input('ClinicalTrials.gov', note='results (API v2), downloaded October 2026')
    eg['ing'] = eg.ingredients.map(lambda x: ' '.join(x) if x is not None else '').str.lower()
    ae = ae.assign(set=ae.term.map(term2sets)).explode('set')
    aff = ae.groupby(['nct', 'eg_id', 'set']).agg(aff=('num_affected', 'max'), risk=('num_at_risk', 'max')).reset_index()
    rows, cells_all = [], {}
    for cls, pat in CLASSES.items():
        r = re.compile(pat)
        g = eg.copy(); g['has'] = g.ing.str.contains(r)
        g = g[((g.role == 'active') & g.has) | ((g.role == 'placebo') & ~g.has)]
        ok = g.groupby('nct').apply(lambda d: (d.role == 'active').any() and (d.role == 'placebo').any()
                                    and d[d.role == 'active'].n_risk.sum() > 0 and d[d.role == 'placebo'].n_risk.sum() > 0,
                                    include_groups=False)
        g = g[g.nct.isin(set(ok[ok].index))]
        for s in OUTCOMES:
            x = g.merge(aff[aff.set == s][['nct', 'eg_id', 'aff', 'risk']], on=['nct', 'eg_id'], how='left')
            x['aff'] = x.aff.fillna(0)
            x['risk'] = x[['risk', 'n_risk']].max(axis=1)
            t = x.groupby(['nct', 'role']).agg(aff=('aff', 'sum'), risk=('risk', 'sum')).unstack('role')
            t.columns = [f'{c}_{r_}' for c, r_ in t.columns]
            t = t.rename(columns={'aff_active': 'a1', 'risk_active': 'n1', 'aff_placebo': 'a0', 'risk_placebo': 'n0'}).dropna()
            t = t[(t.n1 > 0) & (t.n0 > 0) & (t.a1 <= t.n1) & (t.a0 <= t.n0)]
            cells_all[(cls, s)] = t
            a1, n1, a0, n0 = (t[c].to_numpy(float) for c in ('a1', 'n1', 'a0', 'n0'))
            pe, mo, rd = peto(a1, n1, a0, n0), mh_or(a1, n1, a0, n0), mh_rd(a1, n1, a0, n0)
            rows.append(dict(drug_class=cls, outcome=s, n_trials=len(t), n_trials_events=int(((a1 + a0) > 0).sum()),
                             events_drug=int(a1.sum()), n_drug=int(n1.sum()), events_placebo=int(a0.sum()), n_placebo=int(n0.sum()),
                             peto_or=pe[0], peto_lo=pe[1], peto_hi=pe[2], peto_p=pe[3],
                             mh_or=mo[0], mh_lo=mo[1], mh_hi=mo[2], mh_p=mo[3],
                             rd_per_1000=rd[0], rd_lo=rd[1], rd_hi=rd[2], rd_p=rd[3]))
    res = pd.DataFrame(rows)
    # incretin trials by molecule (data scale for Figure 1): the class ingredients named in each trial's drug groups
    inc = re.compile(CLASSES['incretin'])
    g = eg[eg.nct.isin(set(cells_all[('incretin', 'CTS')].index))].copy()
    mol = (g[g.role == 'active'].ing.map(lambda x: '+'.join(sorted(set(inc.findall(x))))).replace('', np.nan).dropna()
           .groupby(g.nct).agg(lambda v: ' / '.join(sorted(set(v)))))
    t = cells_all[('incretin', 'CTS')][['n1', 'n0']].join(mol.rename('molecule'))
    t['molecule'] = t.molecule.str.split(' / ').str[0].str.split('+').map(lambda v: v[0] if len(v) == 1 else '+'.join(v))
    mols = t.reset_index().groupby('molecule').agg(trials=('nct', 'size'), n_drug=('n1', 'sum'), n_placebo=('n0', 'sum'))
    run.table('trial_incretin_molecules', mols.reset_index(), note='Placebo-controlled incretin trials with serious adverse '
                                                                   'event tables, by molecule (first named incretin)')
    print(mols.sort_values('trials', ascending=False))
    # the incretin trials pooled in Fig. 5, one row per trial: molecule, participants at risk, participants with each event
    per = t[['molecule', 'n1', 'n0']].copy()
    for s in ['CTS', 'TF', 'FS', 'DUP', 'RC', 'POS_gallbladder', 'NEG_cataract', 'NEG_bcc', 'NEG_bph']:
        c = cells_all[('incretin', s)]
        per[f'{s}_drug'], per[f'{s}_placebo'] = c.a1.reindex(per.index), c.a0.reindex(per.index)
    run.table('trial_incretin_trials', per.reset_index(), note='Placebo-controlled incretin trials pooled in Fig. 5: '
                                                               'participants at risk and with each event per group')
    print('trials listed:', len(per))
    # trial-level bootstrap of the incretin axis ratio
    rng = np.random.default_rng(20261025)
    tri = sorted(set().union(*[set(cells_all[('incretin', s)].index) for s in ['CTS', 'TF', 'FS', 'DUP']]))

    def axis_lor(idx):
        out = {}
        for ax_, ds in AXES.items():
            oe = vv = 0.0
            for s in ds:
                t = cells_all[('incretin', s)].reindex(idx).dropna()
                a1, n1, a0, n0 = (t[c].to_numpy(float) for c in ('a1', 'n1', 'a0', 'n0'))
                n = n1 + n0; m = a1 + a0; k = m > 0
                oe += (a1 - m * n1 / n)[k].sum(); vv += (m * (n - m) * n1 * n0 / (n ** 2 * (n - 1)))[k].sum()
            out[ax_] = oe / vv if vv > 0 else np.nan
        return out
    base = axis_lor(tri)
    boots = []
    for _ in range(4000):
        idx = list(rng.choice(tri, size=len(tri), replace=True))
        b = axis_lor(idx)
        boots.append(b['fibrocapsular'] - b['tenosynovial'])
    boots = np.array(boots); boots = boots[np.isfinite(boots)]
    d0 = base['fibrocapsular'] - base['tenosynovial']
    bt = pd.DataFrame([dict(drug_class='incretin', n_trials=len(tri), ratio=np.exp(d0), lo=np.exp(np.percentile(boots, 2.5)),
                            hi=np.exp(np.percentile(boots, 97.5)), share_le_0=float((boots <= 0).mean()), n_boot=len(boots))])
    run.table('trial_robustness', res, note='Serious events, incretin vs placebo, all trials with posted serious tables '
                                            '(zero-event trials kept): Peto OR, Mantel-Haenszel OR (no continuity correction) '
                                            'and MH risk difference per 1,000')
    run.table('trial_axis_bootstrap', bt, note='Incretin fibro-capsular / tenosynovial Peto OR ratio with trial-level bootstrap '
                                               '95% interval (4,000 resamples)')
    pd.set_option('display.width', 250)
    print(res[res.drug_class == 'incretin'].round(3).to_string(index=False)); print(bt.round(4).to_string(index=False))
