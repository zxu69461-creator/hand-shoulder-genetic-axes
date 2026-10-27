#!/usr/bin/env python3
"""Two-sample MR of endocrine / metabolic exposures on each hand/shoulder soft-tissue disorder.

Design (no sample overlap between exposure and outcome; mrlib.TO_FINNGEN / TO_UKB):
  - consortium / UK Biobank exposures  -> FinnGen R13 outcomes;
  - exposures without UK Biobank samples -> UK Biobank-based outcomes.
Instruments from `mr_instruments` (30_mr_instruments.py). Harmonised on rsID and alleles; palindromic SNPs dropped.
Outcome effects are log-odds; BOLT-LMM linear-scale effects are converted with beta / (mu (1 - mu)), mu = case fraction.
Estimation in R (MendelianRandomization: IVW random effects, weighted median, MR-Egger) via scripts/31_mr.R.
Tables: mr_harmonised (SNP level), mr_results (pair x method).
"""
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import pandas as pd
import config as C
import mrlib as L
from project import Project

P = Project(__file__)
with P.connect(read_only=True) as con:
    inst = con.execute('select * from mr_instruments').df()
H = []
with P.run(__file__, params={'design': 'no exposure-outcome sample overlap', 'palindromic': 'dropped'}) as run:
    for label, cohort, axis, kind, path, mu in L.OUTCOMES:
        if not Path(path).exists():
            print(f'[skip] outcome {label}/{cohort}: not available yet'); continue
        ex = inst[inst.exposure.isin(L.TO_FINNGEN if cohort == 'FinnGen' else L.TO_UKB)]
        if ex.empty:
            continue
        y = L.outcome_rows(kind, Path(path).resolve(), set(ex.rsid), mu)
        run.input(f'{label}:{cohort}', path=str(Path(path).resolve()), note=f'MR outcome ({axis})')
        for e, x in ex.groupby('exposure'):
            m = L.harmonise(x, y)
            H.append(pd.DataFrame({'exposure': e, 'outcome': label, 'cohort': cohort, 'axis': axis, 'rsid': m.rsid,
                                   'ea': m.ea_x, 'oa': m.oa_x, 'bx': m.beta_x, 'bxse': m.se_x, 'by': m.beta_y,
                                   'byse': m.se_y}))
            print(f'{e:>10} -> {label:<8} {cohort:<16} {len(m)}/{len(x)} SNPs', flush=True)
    if not H:
        sys.exit('no exposure-outcome pair available')
    h = pd.concat(H, ignore_index=True)
    with tempfile.TemporaryDirectory() as tmp:
        fi, fo = Path(tmp) / 'h.tsv', Path(tmp) / 'r.tsv'
        h.to_csv(fi, sep='\t', index=False)
        r = subprocess.run([C.RSCRIPT, str(ROOT / 'scripts' / '31_mr.R'), str(fi), str(fo)], capture_output=True, text=True)
        if r.returncode or not fo.exists():
            sys.exit(f'R MR failed: {r.stdout[-800:]} {r.stderr[-800:]}')
        res = pd.read_csv(fo, sep='\t')
    run.table('mr_harmonised', h, note='SNP-level harmonised exposure-outcome effects (log-odds outcomes; palindromic dropped)')
    run.table('mr_results', res, note='Two-sample MR per exposure x outcome x cohort: IVW (random effects), weighted '
                                      'median, MR-Egger (MendelianRandomization); estimates per SD / log-odds of exposure')
    print(res[res.method == 'IVW'].sort_values('p').head(25).to_string(index=False))
