#!/usr/bin/env python3
"""Multivariable MR: each endocrine / metabolic exposure together with BMI, on every FinnGen soft-tissue outcome.

Pairs: BMI (GIANT 2015, no UK Biobank) + one of SHBG_F, BIOT_F, TT_F, T2D_XUE, HBA1C, IGF1, ANM (when available).
Instruments: union of both exposures' clumped instruments (mr_instruments), re-clumped on the smaller p (r2<0.001,
10 Mb, 1000 Genomes Europeans); effects of every instrument looked up in both exposure files and in the FinnGen
outcome (no sample overlap with the exposures). Estimation: MendelianRandomization::mr_mvivw (scripts/32_mvmr.R).
Table: mvmr_results (direct effect of each exposure conditional on BMI, per outcome; conditional F statistics).
"""
import io
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
PARTNERS = ['SHBG_F', 'BIOT_F', 'TT_F', 'T2D_XUE', 'HBA1C', 'IGF1', 'ANM']
FG_OUT = [o for o in L.OUTCOMES if o[1] == 'FinnGen']

with P.connect(read_only=True) as con:
    inst = con.execute('select * from mr_instruments').df()
res = []
with P.run(__file__, params={'base': 'BMI_GIANT', 'reclump': 'r2<0.001 10Mb', 'method': 'mr_mvivw random'}) as run, \
        tempfile.TemporaryDirectory() as tmp:
    for partner in PARTNERS:
        if partner not in set(inst.exposure):
            print(f'[skip] {partner}: no instruments yet'); continue
        u = inst[inst.exposure.isin(['BMI_GIANT', partner])].sort_values('p').drop_duplicates('rsid')
        u, _ = L.clump(u.drop(columns=['exposure', 'F']), tmp, f'mv_{partner}')
        snps = set(u.rsid)
        a = L.exposure_lookup('BMI_GIANT', snps)[['rsid', 'ea', 'oa', 'beta', 'se', 'eaf']]
        b = L.harmonise(a, L.exposure_lookup(partner, snps)[['rsid', 'ea', 'oa', 'beta', 'se', 'eaf']])
        x = b.rename(columns={'ea_x': 'ea', 'oa_x': 'oa', 'beta_x': 'b1', 'se_x': 's1', 'beta_y': 'b2', 'se_y': 's2'})
        x = x[['rsid', 'ea', 'oa', 'b1', 's1', 'b2', 's2']]
        for label, cohort, axis, kind, path, mu in FG_OUT:
            y = L.outcome_rows(kind, Path(path).resolve(), set(x.rsid), mu)
            h = L.harmonise(x.rename(columns={'b1': 'beta', 's1': 'se'}), y)
            d = h.rename(columns={'beta_x': 'b1', 'se_x': 's1', 'beta_y': 'by', 'se_y': 'sy'})[['rsid', 'b1', 's1', 'b2', 's2', 'by', 'sy']]
            f = Path(tmp) / f'{partner}_{label}.tsv'
            d.to_csv(f, sep='\t', index=False)
            r = subprocess.run([C.RSCRIPT, str(ROOT / 'scripts' / '32_mvmr.R'), str(f)], capture_output=True, text=True)
            if r.returncode:
                print(f'[fail] {partner} {label}: {r.stderr[-300:]}'); continue
            o = pd.read_csv(io.StringIO(r.stdout), sep='\t')
            o['exposure'] = o.exposure.replace({'partner': partner})
            o.insert(0, 'outcome', label); o.insert(1, 'axis', axis); o.insert(2, 'model', f'BMI + {partner}')
            res.append(o)
        run.input(f'mvmr:{partner}', path=str(Path(L.EXPO[partner][2]).resolve()), note=f'MVMR BMI + {partner}')
        print(f'[done] BMI + {partner}: {len(x)} instruments', flush=True)
    out = pd.concat(res, ignore_index=True)
    run.table('mvmr_results', out, note='Multivariable MR (BMI + one endocrine/metabolic exposure) on FinnGen soft-tissue '
                                        'outcomes: direct effects (log-odds per unit exposure), conditional F')
    print(out.assign(OR=lambda d: d.estimate.map(lambda v: round(2.718281828 ** v, 3))).to_string(index=False))
