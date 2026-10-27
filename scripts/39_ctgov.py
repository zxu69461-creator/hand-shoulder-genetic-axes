#!/usr/bin/env python3
"""ClinicalTrials.gov results (API v2) -> trial, event-group and adverse-event tables for the trial analyses.

Download (into data/ctgov/):
  default    the trials listed in results/source_data/fig_5_included_trials.csv (the analysed set);
  --search   every study with posted results whose interventions name an incretin receptor agonist, to re-derive the
             selection (the inclusion criteria are applied in 40_trial_axis.py and 43_trial_robustness.py);
  --pages D  parse an existing bulk download instead (API v2 JSON pages, page_*.json.gz).
Event groups are linked to the registered arms: title similarity >= 0.55, otherwise containment of the arm label in the
group title and description >= 0.75. A linked group takes the role of its arm (placebo or sham comparator -> placebo;
experimental, active comparator or other -> active) and the drug names of the arm's drug or biological interventions,
normalised (doses, formulations and salts removed); an unlinked group is classed from its title. Groups that summarise
all participants are marked 'total'.
Output: work/ctgov/ctgov.duckdb with tables ct_trial (nct, allocation), ct_eventgroup (nct, eg_id, title, role,
ingredients, participants at risk in the serious and other tables) and ct_ae_term (nct, eg_id, kind, term,
participants affected and at risk).
"""
import argparse
import difflib
import gzip
import json
import re
import subprocess
import sys
import time
import urllib.parse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import duckdb
import pandas as pd
import config as C

BASE = 'https://clinicaltrials.gov/api/v2/studies'
FIELDS = 'NCTId,DesignModule,ArmsInterventionsModule,AdverseEventsModule'
INCRETINS = ['exenatide', 'liraglutide', 'lixisenatide', 'albiglutide', 'dulaglutide', 'semaglutide', 'efpeglenatide',
             'taspoglutide', 'tirzepatide', 'retatrutide', 'survodutide', 'mazdutide', 'orforglipron', 'cagrilintide',
             'danuglipron', 'cotadutide']
DOSE_RX = re.compile(r'\b\d+(\.\d+)?\s*(mg|mcg|µg|ug|g|ml|iu|units?|%|mg/kg|mg/m2|mg/m²)\b(/\s*(day|kg|m2|dose|week|d|wk))?', re.I)
FORM_RX = re.compile(r'\b(tablets?|capsules?|injection|injectable|solution|oral|iv|i\.v\.|intravenous|subcutaneous|sc|s\.c\.|infusion|cream|gel|patch|once daily|twice daily|bid|qd|tid|q\d+w|weekly|daily|extended.release|xr|er|sr|low dose|high dose|dose|placebo matched)\b', re.I)
SALT_RX = re.compile(r'\b(hydrochloride|hcl|sodium|potassium|calcium|mesylate|maleate|tartrate|citrate|sulfate|sulphate|acetate|succinate|fumarate|besylate|phosphate|bromide|tosylate|dihydrate|monohydrate)\b', re.I)
PLACEBO_RX = re.compile(r'\b(placebo|vehicle|sham|dummy|saline|no treatment|untreated|control)\b', re.I)
CONTROL_RX = re.compile(r'\b(control|standard of care|usual care|observation|no intervention)\b', re.I)
TOTAL_RX = re.compile(r'^(total|all participants|overall|all patients|all subjects|total participants)$', re.I)


def norm(s):
    s = (s or '').lower()
    s = re.sub(r'^(drug|biological|other)\s*:\s*', '', s)
    s = DOSE_RX.sub(' ', s)
    s = FORM_RX.sub(' ', s)
    s = re.sub(r'[\(\)\[\],;:®™]', ' ', s)
    s = SALT_RX.sub(' ', s)
    return re.sub(r'\s+', ' ', s).strip(' -+/')


def drug_name(name):
    """normalised drug name, or None for placebo-like or empty names"""
    n = norm(name)
    return None if not n or PLACEBO_RX.search(n) else n


def sim(a, b):
    return difflib.SequenceMatcher(None, (a or '').lower(), (b or '').lower()).ratio()


def containment(label, text):
    toks = [t for t in re.findall(r'[a-z0-9]+', (label or '').lower()) if t not in ('dose', 'arm', 'group', 'the', 'and', 'with', 'of', 'mg')]
    low = (text or '').lower()
    if not toks:
        return 0.0, 10 ** 9
    hit = [t for t in toks if re.search(r'\b' + re.escape(t) + r'\b', low)]
    distinct = sorted(hit, key=lambda t: (not t.isdigit(), -len(t)))
    pos = low.find(distinct[0]) if distinct else 10 ** 9
    return len(hit) / len(toks), pos


class Trial:
    """Arms, their drug names, and the link of each adverse-event group to an arm."""
    def __init__(self, arms, ivs):
        self.arms, self.cache = arms, {}
        self.arm_ing = {}
        for a in arms:
            names = [iv.get('name') for iv in ivs if a.get('label') in (iv.get('armGroupLabels') or [])
                     and iv.get('type') in ('DRUG', 'BIOLOGICAL')]
            self.arm_ing[a.get('label')] = sorted({n for n in map(drug_name, names) if n})

    def link(self, g):
        title, desc = g.get('title') or '', g.get('description') or ''
        key = (title, desc[:200])
        if key in self.cache:
            return self.cache[key]
        if TOTAL_RX.search(title.strip()):
            r = self.cache[key] = ('total', [])
            return r
        best, bscore = None, 0.0
        for a in self.arms:
            sc = max(sim(title, a.get('label')), sim(title, ' '.join(a.get('interventionNames') or [])))
            if sc > bscore:
                best, bscore = a, sc
        if (best is None or bscore < 0.55) and desc:
            cands = []
            for a in self.arms:
                cf, pos = containment(a.get('label'), title + ' ' + desc)
                cands.append((cf, -pos, a))
            cands.sort(key=lambda x: (x[0], x[1]), reverse=True)
            if cands and cands[0][0] >= 0.75:
                best, bscore = cands[0][2], cands[0][0]
        if best is not None and bscore >= 0.55:
            at = best.get('type') or ''
            role = 'placebo' if at in ('PLACEBO_COMPARATOR', 'SHAM_COMPARATOR') else 'active' if at in ('EXPERIMENTAL', 'ACTIVE_COMPARATOR', 'OTHER') else 'other'
            ing = self.arm_ing.get(best.get('label'), [])
        else:
            role = 'placebo' if PLACEBO_RX.search(title) else 'other' if CONTROL_RX.search(title.strip()) else 'active'
            n = drug_name(title)
            ing = [n] if n else []
        r = self.cache[key] = (role, list(ing))
        return r


def num(x):
    try:
        v = float(x)
        return v if v == v else None
    except (TypeError, ValueError):
        return None


def get(url):
    for attempt in range(8):
        r = subprocess.run(['curl', '-s', '--max-time', '600', url], capture_output=True)
        try:
            d = json.loads(r.stdout)
            if 'studies' in d:
                return d
        except json.JSONDecodeError:
            pass
        time.sleep(15 * (attempt + 1))
    raise RuntimeError(f'request failed: {url[:120]}')


def download(query, out):
    """all pages of an API v2 query -> one gzipped JSON list of studies"""
    studies, token = [], None
    while True:
        url = f'{BASE}?{query}&fields={FIELDS}&pageSize=1000&format=json'
        if token:
            url += f'&pageToken={token}'
        d = get(url)
        studies += d['studies']
        token = d.get('nextPageToken')
        if not token:
            break
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, 'wt') as fh:
        json.dump({'studies': studies}, fh)
    print(len(studies), 'studies ->', out)


def studies_from(files):
    for f in files:
        for st in json.loads(gzip.open(f, 'rt').read())['studies']:
            yield st


INC_RX = re.compile('|'.join(INCRETINS))


def parse(files, keep=None):
    T, EG, AE = [], [], []
    for st in studies_from(files):
        ps, rs = st.get('protocolSection', {}), st.get('resultsSection', {}) or {}
        nct = ps.get('identificationModule', {}).get('nctId')
        if not nct or (keep is not None and nct not in keep):
            continue
        ae = rs.get('adverseEventsModule', {}) or {}
        aim = ps.get('armsInterventionsModule', {}) or {}
        if not ae or not INC_RX.search((json.dumps(aim) + json.dumps(ae.get('eventGroups', []))).lower()):
            continue                                    # no incretin named in the arms, interventions or event groups
        tr = Trial(aim.get('armGroups', []) or [], aim.get('interventions', []) or [])
        T.append(dict(nct=nct, allocation=((ps.get('designModule', {}) or {}).get('designInfo', {}) or {}).get('allocation')))
        for g in ae.get('eventGroups', []) or []:
            role, ing = tr.link(g)
            EG.append(dict(nct=nct, eg_id=g.get('id'), title=(g.get('title') or '')[:300], role=role, ingredients=ing,
                           n_risk_serious=num(g.get('seriousNumAtRisk')), n_risk_other=num(g.get('otherNumAtRisk'))))
        for kind in ('seriousEvents', 'otherEvents'):
            for ev in ae.get(kind, []) or []:
                for s in ev.get('stats', []) or []:
                    AE.append(dict(nct=nct, eg_id=s.get('groupId'), kind='serious' if kind == 'seriousEvents' else 'other',
                                   term=(ev.get('term') or '').strip(), num_affected=num(s.get('numAffected')),
                                   num_at_risk=num(s.get('numAtRisk'))))
    return pd.DataFrame(T), pd.DataFrame(EG), pd.DataFrame(AE)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--search', action='store_true')
    ap.add_argument('--pages', type=Path)
    a = ap.parse_args()
    raw = C.CTGOV / ('studies_search.json.gz' if a.search else 'studies_included.json.gz')
    keep = None
    if a.pages:
        files = sorted(a.pages.glob('page_*.json.gz'))
    else:
        if a.search:
            q = ('filter.advanced=' + urllib.parse.quote('AREA[ResultsFirstPostDate]RANGE[MIN,MAX]')
                 + '&query.intr=' + urllib.parse.quote(' OR '.join(INCRETINS)))
        else:
            ids = pd.read_csv(ROOT / 'results' / 'source_data' / 'fig_5_included_trials.csv').astype(str)
            keep = sorted({v for col in ids for v in ids[col] if re.fullmatch(r'NCT\d{8}', v)})
            q = 'filter.ids=' + ','.join(keep)
        if not raw.exists():
            download(q, raw)
        files = [raw]
    T, EG, AE = parse(files, keep)
    out = C.WORK / 'ctgov'; out.mkdir(parents=True, exist_ok=True)
    with duckdb.connect(str(out / 'ctgov.duckdb')) as con:
        for name, df in (('ct_trial', T), ('ct_eventgroup', EG), ('ct_ae_term', AE)):
            con.register('_d', df)
            con.execute(f'CREATE OR REPLACE TABLE {name} AS SELECT * FROM _d')
            con.unregister('_d')
    print(f'{len(T)} trials, {len(EG)} event groups, {len(AE)} adverse-event rows -> {out / "ctgov.duckdb"}')


if __name__ == '__main__':
    main()
