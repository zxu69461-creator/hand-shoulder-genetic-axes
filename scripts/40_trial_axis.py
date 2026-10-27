#!/usr/bin/env python3
"""Soft-tissue serious adverse events under incretin receptor agonists in placebo-controlled randomised trials
(ClinicalTrials.gov results; tables built by 39_ctgov.py).

Trials: randomised, with an event group whose drugs include an incretin receptor agonist (GLP-1 receptor agonists and
GIP / glucagon co-agonists) and a placebo group that does not. Outcomes are exact MedDRA preferred terms (lower-cased)
per disorder; positive control = gallbladder events, negative controls = cataract, basal cell carcinoma and benign
prostatic hyperplasia. Per trial and group: participants affected = maximum over the set's terms in the serious-event
table (reported without a frequency threshold); at risk = maximum reported. Pooled per outcome with the Peto one-step
odds ratio (rare events, trials as strata) and the Mantel-Haenszel risk ratio; leave-one-trial-out range of the Peto OR.
Tables: trial_axis_cells, trial_axis_pooled.
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
    'PF': ['plantar fasciitis'],
    'RC': ['rotator cuff syndrome'],
    'POS_gallbladder': ['cholelithiasis', 'cholecystitis', 'cholecystitis acute'],
    'NEG_cataract': ['cataract'],
    'NEG_bcc': ['basal cell carcinoma'],
    'NEG_bph': ['benign prostatic hyperplasia'],
}


def peto(t):
    """t: rows with a1, n1, a0, n0 per trial -> pooled Peto OR, 95% CI, p."""
    n = t.n1 + t.n0; m = t.a1 + t.a0
    e = m * t.n1 / n
    v = m * (n - m) * t.n1 * t.n0 / (n ** 2 * (n - 1))
    oe, vv = (t.a1 - e).sum(), v.sum()
    if vv <= 0:
        return np.nan, np.nan, np.nan, np.nan
    lor = oe / vv; se = 1 / np.sqrt(vv)
    return np.exp(lor), np.exp(lor - 1.96 * se), np.exp(lor + 1.96 * se), 2 * stats.norm.sf(abs(lor / se))


def mh_rr(t):
    n = t.n1 + t.n0
    num, den = (t.a1 * t.n0 / n).sum(), (t.a0 * t.n1 / n).sum()
    return num / den if den > 0 else np.nan


rx = {k: re.compile(v) for k, v in CLASSES.items()}
term2set = {t: k for k, ts in OUTCOMES.items() for t in ts}
with P.run(__file__, params={'classes': CLASSES, 'outcomes': OUTCOMES, 'pool': 'Peto OR + MH RR, trials as strata'}) as run:
    con = duckdb.connect(str(C.WORK / 'ctgov' / 'ctgov.duckdb'), read_only=True)
    eg = con.execute("""
        SELECT e.nct, e.eg_id, e.role, e.ingredients, greatest(coalesce(e.n_risk_serious,0), coalesce(e.n_risk_other,0)) AS n_risk
        FROM ct_eventgroup e JOIN ct_trial t USING (nct)
        WHERE t.allocation = 'RANDOMIZED' AND e.role IN ('active', 'placebo')""").df()
    ae = con.execute(f"""
        SELECT nct, eg_id, kind, lower(term) AS term, num_affected, num_at_risk FROM ct_ae_term
        WHERE lower(term) IN ({', '.join("'" + t.replace("'", "''") + "'" for t in term2set)})""").df()
    con.close()
    eg['ing'] = eg.ingredients.map(lambda x: ' '.join(x) if x is not None else '').str.lower()
    ae['set'] = ae.term.map(term2set)
    agg = {'serious': ae[ae.kind == 'serious']}
    cells = []
    for analysis, sub in agg.items():
        a = sub.groupby(['nct', 'eg_id', 'set']).agg(aff=('num_affected', 'max'), risk=('num_at_risk', 'max')).reset_index()
        for cls, r in rx.items():
            g = eg.copy()
            g['has'] = g.ing.str.contains(r)
            ok = g.groupby('nct').apply(lambda d: d[(d.role == 'active')].has.any() and (d[d.role == 'placebo'].shape[0] > 0)
                                        and not d[d.role == 'placebo'].has.any(), include_groups=False)
            trials = set(ok[ok].index)
            g = g[g.nct.isin(trials) & (((g.role == 'active') & g.has) | (g.role == 'placebo'))]
            for s in OUTCOMES:
                x = g.merge(a[a.set == s], on=['nct', 'eg_id'], how='inner')
                if x.empty:
                    continue
                reported = set(x.nct)
                gg = g[g.nct.isin(reported)].merge(a[a.set == s][['nct', 'eg_id', 'aff', 'risk']], on=['nct', 'eg_id'], how='left')
                gg['aff'] = gg.aff.fillna(0)
                gg['risk'] = gg[['risk', 'n_risk']].max(axis=1)
                t = gg.groupby(['nct', 'role']).agg(aff=('aff', 'sum'), risk=('risk', 'sum')).unstack('role')
                t.columns = [f'{c}_{r}' for c, r in t.columns]
                t = t.dropna().reset_index().rename(columns={'aff_active': 'a1', 'risk_active': 'n1', 'aff_placebo': 'a0', 'risk_placebo': 'n0'})
                t = t[(t.n1 > 0) & (t.n0 > 0)]
                t.insert(0, 'outcome', s); t.insert(0, 'drug_class', cls); t.insert(0, 'analysis', analysis)
                cells.append(t)
    cells = pd.concat(cells, ignore_index=True)
    rows = []
    for (an, c, s), t in cells.groupby(['analysis', 'drug_class', 'outcome']):
        t = t[(t.a1 + t.a0) > 0]
        orr, lo, hi, p = peto(t) if len(t) else (np.nan,) * 4
        loo = [peto(t.drop(i))[0] for i in t.index] if len(t) > 2 else []
        rows.append(dict(analysis=an, drug_class=c, outcome=s, n_trials=len(t),
                         loo_min=np.nanmin(loo) if loo else np.nan, loo_max=np.nanmax(loo) if loo else np.nan,
                         events_active=int(t.a1.sum()), n_active=int(t.n1.sum()),
                         events_placebo=int(t.a0.sum()), n_placebo=int(t.n0.sum()), peto_or=orr, or_lo=lo, or_hi=hi, p=p,
                         mh_rr=mh_rr(t) if len(t) else np.nan))
    pooled = pd.DataFrame(rows)
    run.table('trial_axis_cells', cells, note='Per placebo-controlled trial: participants affected / at risk in incretin vs '
                                              'placebo groups, per soft-tissue outcome and control event (serious-event tables)')
    run.table('trial_axis_pooled', pooled, note='Peto one-step OR and MH RR per outcome (trials with >=1 event)')
    print(pooled.round(3).to_string(index=False))
