#!/usr/bin/env python3
"""LDSC genetic correlations of the two axis factor GWAS with each other, with the seven disorder GWAMAs and with the
endocrine and metabolic exposures (multivariate LD score regression, lib/fastldsc.py; 1000 Genomes European LD
scores, HapMap 3 SNPs, 200-block jackknife). Table: factor_rg.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import numpy as np
import pandas as pd
from scipy import stats
import config as C
import fastldsc as FL
from project import Project

P = Project(__file__)
FAC = ['FAC_TEN', 'FAC_FIB']
DIS = ['CTS', 'TF', 'DQ', 'ENTH', 'PF', 'FS', 'DUP']
EXP = ['BMI_GIANT', 'IGF1', 'SHBG_F', 'BIOT_F', 'TT_F', 'T2D_XUE', 'HBA1C', 'ANM']
files = ([C.FACTOR / f'{x}.sumstats.gz' for x in FAC] + [C.GWAMA / f'{d}.sumstats.gz' for d in DIS]
         + [C.MUNGED / f'{x}.sumstats.gz' for x in EXP])
names = FAC + DIS + EXP
with P.run(__file__, params={'ld': 'eur_w_ld_chr', 'blocks': 200}) as run:
    for n, f in zip(names, files):
        run.input(n, path=str(f.resolve()), note='munged HapMap3 summary statistics')
    ld, M = FL.load_ld(C.LDSC_LD); hm3 = pd.read_csv(C.HM3, sep='\t')
    Z, N = FL.load_traits(files, ld, hm3)
    r = FL.run(Z, N, ld.L2.values, ld.L2.values, M, tile=64)
    rows = [dict(factor=FAC[i], trait=names[j], kind='factor' if j < 2 else 'disorder' if j < 2 + len(DIS) else 'exposure',
                 rg=r['rg'][i, j], se=r['rg_se'][i, j], intercept=r['intercept'][i, j])
            for i in range(2) for j in range(len(names)) if j != i and not (i == 1 and j == 0)]
    out = pd.DataFrame(rows)
    out['p'] = 2 * stats.norm.sf(np.abs(out.rg / out.se))
    run.table('factor_rg', out, note='LDSC rg of the axis factor GWAS with each other, the disorder GWAMAs and the exposures')
    print(out.round(3).to_string(index=False))
