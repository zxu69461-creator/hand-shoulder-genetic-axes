#!/usr/bin/env python3
"""Adiposity contrast between carpal tunnel syndrome and Dupuytren's disease from non-overlapping outcome samples, plus the
instrument and heritability summaries reported with the genetic analyses.

  mr_bmi_independent  IVW log OR per SD BMI on CTS in one cohort minus that on Dupuytren's disease in the other (FinnGen vs
                      UK Biobank, no shared participants or controls), z test with independent standard errors
  mr_instrument_strength  instruments and F statistics per exposure (clumped genome-wide significant variants)
  structure_h2        observed-scale SNP heritability, standard error, Z and LDSC intercept of the eight FinnGen endpoints
                      (LD score regression, lib/fastldsc.py; 1000 Genomes European LD scores, 200-block jackknife)
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import numpy as np
import pandas as pd
import config as C
import fastldsc as FL
from project import Project
from scipy import stats

P = Project(__file__)
FG = {'CTS': 'G6_CARPTU', 'TF': 'M13_TRIGGERFINGER', 'DQ': 'M13_DEQUERVAIN', 'DUP': 'M13_DUPUTRYEN',
      'FS': 'M13_ADHCAPSULITIS', 'PF': 'M13_PLANTARFASC', 'ENTH': 'M13_ENTESOPATHYOTH', 'RC': 'M13_ROTATORCUFF'}
EXP = ['BMI_GIANT', 'IGF1', 'SHBG_F', 'BIOT_F', 'T2D_XUE']
with P.run(__file__) as run:
    with P.connect(read_only=True) as con:
        s = con.execute("select exposure, outcome, cohort, estimate_IVW, se_IVW from mr_summary "
                        "where exposure = 'BMI_GIANT'").df().set_index(['outcome', 'cohort'])
        ins = con.execute('select exposure, count(*) n, avg(F) mean_f, min(F) min_f from mr_instruments where exposure in ('
                          + ', '.join(f"'{e}'" for e in EXP) + ') group by 1').df()
    run.input('mr_summary', note='IVW estimates per exposure x outcome x cohort')
    rows = []
    for (c1, c2, dd) in (('FinnGen', 'UKB', 'DUP'), ('UKB', 'FinnGen', 'DUP'), ('FinnGen', 'non-Finnish meta', 'DUP_META')):
        a, b = s.loc[('CTS', c1)], s.loc[(dd, c2)]
        d = a.estimate_IVW - b.estimate_IVW; se = np.hypot(a.se_IVW, b.se_IVW)
        rows.append(dict(cts_cohort=c1, dd_cohort=c2, cts_or=np.exp(a.estimate_IVW), dd_or=np.exp(b.estimate_IVW),
                         ratio=np.exp(d), lo=np.exp(d - 1.96 * se), hi=np.exp(d + 1.96 * se), p=2 * stats.norm.sf(abs(d / se))))
    ld, M = FL.load_ld(C.LDSC_LD); hm3 = pd.read_csv(C.HM3, sep='\t')
    eps = list(FG.values())
    Z, N = FL.load_traits([C.MUNGED / f'{e}.sumstats.gz' for e in eps], ld, hm3)
    r = FL.run(Z, N, ld.L2.values, ld.L2.values, M, tile=64)
    h2 = pd.DataFrame({'ep': eps, 'h2_obs': r['h2'], 'h2_se': r['h2_se'], 'intercept': np.diag(r['intercept'])})
    h2['z'] = h2.h2_obs / h2.h2_se
    h2.insert(0, 'disorder', h2.ep.map({v: k for k, v in FG.items()}))
    run.input('FinnGen R13 endpoints', note='HapMap 3 summary statistics (10_munge.py)')
    run.table('mr_bmi_independent', pd.DataFrame(rows), note='BMI effect on CTS vs Dupuytren disease from non-overlapping '
                                                           'cohorts (ratio of ORs, z test)')
    run.table('mr_instrument_strength', ins, note='Instruments and F statistics per exposure')
    run.table('structure_h2', h2, note='Observed-scale SNP heritability, SE, Z and intercept of the eight FinnGen endpoints')
    print(pd.DataFrame(rows).round(4).to_string(index=False)); print(h2.round(3).to_string(index=False)); print(ins.round(1))
