#!/usr/bin/env python3
"""Genomic SEM of the disorder GWAS: genetic correlations, exploratory and confirmatory factor models.

Trait sets: fg_core (six FinnGen disorders), ukb_core (five UK Biobank-based GWAS) and xc (both cohorts together, with
the disease-factor models of 20_gsem_models_xc.R). The two-factor models load the tenosynovial disorders (TEN: carpal
tunnel syndrome, trigger finger, de Quervain tenosynovitis) and the fibro-capsular disorders (FIB: Dupuytren's
disease, plantar fascial disorder, frozen shoulder, enthesopathy) on separate factors.
Runs scripts/20_gsem_structure.R, then writes tables gsem_rg, gsem_efa, gsem_cfa_fit and gsem_cfa_est.
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import pandas as pd
import config as C
from project import Project

P = Project(__file__)
STAGE, OUT = C.WORK / 'gsem', C.WORK / 'gsem' / 'out'
M = C.MUNGED
FG = {'CTS': ('G6_CARPTU', 0), 'TF': ('M13_TRIGGERFINGER', 1), 'DQ': ('M13_DEQUERVAIN', 1), 'DUP': ('M13_DUPUTRYEN', 1),
      'PF': ('M13_PLANTARFASC', 1), 'FS': ('M13_ADHCAPSULITIS', 1)}
UKB = {'CTS': 'UKB_CTS', 'TF': 'UKB_TF', 'DUP': 'UKB_DUP', 'FS': 'UKB_FS', 'ENTH': 'UKB_ENTHESO'}
AXIS = {'CTS': 'TEN', 'TF': 'TEN', 'DQ': 'TEN', 'DUP': 'FIB', 'PF': 'FIB', 'FS': 'FIB', 'ENTH': 'FIB'}
SETS = ['fg_core', 'ukb_core', 'xc']
NOTES = {
    'gsem_rg': 'Genomic SEM LDSC genetic correlations (standardised S, all chromosomes) with jackknife SE, per trait set',
    'gsem_efa': 'Exploratory factor loadings (promax) on the odd-chromosome standardised genetic covariance, 1-3 factors',
    'gsem_cfa_fit': 'Confirmatory model fit on even chromosomes (DWLS): m1 one factor, m2 tenosynovial + fibro-capsular '
                    'factors, m3 = m2 + orthogonal method factor on the FinnGen endpoints sharing one control group',
    'gsem_cfa_est': 'Confirmatory model estimates (standardised loadings, factor correlations) on even chromosomes',
}


def rows(name):
    if name == 'fg_core':
        return [(k, M / f'{FG[k][0]}.sumstats.gz', FG[k][1]) for k in ['CTS', 'TF', 'DQ', 'DUP', 'PF', 'FS']]
    if name == 'ukb_core':
        return [(k, M / f'{UKB[k]}.sumstats.gz', 0) for k in ['CTS', 'TF', 'DUP', 'FS', 'ENTH']]
    if name == 'xc':
        return ([(f'{k}_FG', M / f'{FG[k][0]}.sumstats.gz', FG[k][1]) for k in ['CTS', 'TF', 'DQ', 'DUP', 'PF', 'FS']]
                + [(f'{k}_UKB', M / f'{UKB[k]}.sumstats.gz', 0) for k in ['CTS', 'TF', 'DUP', 'FS', 'ENTH']])
    raise SystemExit(f'unknown set {name}')


STAGE.mkdir(parents=True, exist_ok=True)
lines = ['label\tfile\tset\tm1\tm2\tctrl']
for s in SETS:
    for label, f, ctrl in rows(s):
        if not f.exists():
            raise SystemExit(f'{s}: {f} missing (run 10_munge.py)')
        lines.append('\t'.join([label, str(f), s, 'SOFT', AXIS[label.split('_')[0]], str(ctrl)]))
(STAGE / 'traits.tsv').write_text('\n'.join(lines) + '\n')
subprocess.run([C.RSCRIPT, str(ROOT / 'scripts' / '20_gsem_structure.R'), str(STAGE / 'traits.tsv'), str(C.LDSC_LD),
                str(OUT), str(ROOT / 'scripts' / '20_gsem_models_xc.R')], check=True)

with P.run(__file__, params={'sets': SETS}) as run:
    for tab, pat in [('gsem_rg', 'rg_*.tsv'), ('gsem_efa', 'efa_*.tsv'), ('gsem_cfa_fit', 'cfa_fit_*.tsv'),
                     ('gsem_cfa_est', 'cfa_est_*.tsv')]:
        files = [f for f in sorted(OUT.glob(pat)) if f.stat().st_size > 0]
        if not files:
            continue
        df = pd.concat([pd.read_csv(f, sep='\t') for f in files], ignore_index=True)
        df.columns = [c.replace('.', '_') for c in df.columns]
        run.table(tab, df, note=NOTES[tab])
