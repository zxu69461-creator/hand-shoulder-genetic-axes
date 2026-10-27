#!/usr/bin/env python3
"""Summarise the axis factor GWAS (25_factor_gwas.sh) and run Mendelian randomization on the latent axes.

1. factor_gwas (Parquet view; input work/factor/factor_<AXIS>.tsv.gz from 25_factor_gwas.sh): per axis (TEN = carpal tunnel syndrome, trigger finger, de Quervain tenosynovitis;
   FIB = Dupuytren's disease, frozen shoulder, plantar fascial disorder) the SNP effect on the common factor (est, SE,
   Z, P) and the heterogeneity statistic Q_SNP (SNP effects not mediated by the factor).
   factor_gwas_summary: SNPs, effective N (GenomicSEM recommendation: mean of (Z/est)^2 / (2 MAF (1 - MAF)) over SNPs
   with MAF 0.1-0.4), genome-wide significant lead SNPs (P < 5e-8, clumped r2 < 0.1 within 1 Mb), lambda_GC.
   work/factor/FAC_<AXIS>.sumstats.gz: SNP A1 A2 Z N for MAGMA (51_magma.py) and LDSC (30_factor_rg.py).
2. mr_factor_results: two-sample MR of the five exposures on each factor (instruments P < 5e-8 among the factor GWAS SNPs,
   r2 < 0.001 within 10 Mb; IVW random effects, weighted median, MR-Egger via scripts/31_mr.R); estimates are SD change
   of the latent axis per SD of exposure (per log-odds of type 2 diabetes liability).
"""
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import numpy as np
import pandas as pd
from scipy import stats
import config as C
import mrlib as L
from project import Project

P = Project(__file__)
OUT = C.FACTOR
OUT.mkdir(parents=True, exist_ok=True)
EXP = ['BMI_GIANT', 'IGF1', 'SHBG_F', 'BIOT_F', 'T2D_XUE']
with P.run(__file__, params={'exposures': EXP}) as run, tempfile.TemporaryDirectory() as tmp:
    allf, summ, harm = [], [], []
    for ax in ('TEN', 'FIB'):
        src = OUT / f'factor_{ax}.tsv.gz'
        if not src.exists():
            print(ax, 'not available yet'); continue
        run.input(f'factor_gwas:{ax}', path=str(src), note='GenomicSEM commonfactorGWAS (DWLS)')
        f = pd.read_csv(src, sep='\t').rename(columns={'se_c': 'SE'})
        f = f[(f.SE > 0) & np.isfinite(f.est)]
        if 'fail' in f:
            f = f[f.fail.astype(str).isin(['0', 'False', 'FALSE'])]
        f = f.rename(columns={'Z_Estimate': 'Z', 'Pval_Estimate': 'P'})
        mid = f[(f.MAF >= 0.1) & (f.MAF <= 0.4)]
        neff = float(np.mean((mid.Z / mid.est) ** 2 / (2 * mid.MAF * (1 - mid.MAF))))
        lam = float(np.median(f.Z ** 2) / stats.chi2.ppf(0.5, 1))
        f[['SNP', 'A1', 'A2', 'Z']].assign(N=neff).to_csv(OUT / f'FAC_{ax}.sumstats.gz', sep='\t', index=False)
        # lead SNPs
        sig = f[f.P < 5e-8][['SNP', 'P']].rename(columns={'P': 'P'})
        n_lead = 0
        if len(sig):
            a = Path(tmp) / f'{ax}.assoc'; sig.to_csv(a, sep='\t', index=False)
            subprocess.run([L.PLINK, '--bfile', L.REF, '--clump', str(a), '--clump-p1', '5e-8', '--clump-p2', '1e-5',
                            '--clump-r2', '0.1', '--clump-kb', '1000', '--out', str(Path(tmp) / ax)], capture_output=True)
            cl = Path(tmp) / f'{ax}.clumped'
            n_lead = len(pd.read_csv(cl, sep=r'\s+')) if cl.exists() else 0
        summ.append(dict(axis=ax, n_snp=len(f), n_eff=neff, lambda_gc=lam, n_gws=int((f.P < 5e-8).sum()), n_lead=n_lead,
                         n_q_sig=int((f.Q_pval < 5e-8).sum())))
        allf.append(f[['SNP', 'CHR', 'BP', 'MAF', 'A1', 'A2', 'est', 'SE', 'Z', 'P', 'Q', 'Q_df', 'Q_pval']].assign(axis=ax))
        # MR on the factor
        y = f[['SNP', 'A1', 'A2', 'est', 'SE']].rename(columns={'SNP': 'rsid', 'A1': 'ea', 'A2': 'oa', 'est': 'beta', 'SE': 'se'})
        snps = set(y.rsid)
        for ex in EXP:
            g = L.exposure_gws(ex)
            g = g[g.rsid.isin(snps)]
            x, _ = L.clump(g, tmp, f'{ex}_{ax}')
            mm = L.harmonise(x[['rsid', 'ea', 'oa', 'beta', 'se']], y)
            harm.append(mm.assign(exposure=ex, outcome=f'FACTOR_{ax}', cohort='GenomicSEM', axis=ax))
    fac = pd.concat(allf, ignore_index=True)
    run.parquet('factor_gwas', fac, note='GenomicSEM common-factor GWAS of the tenosynovial (TEN) and fibro-capsular (FIB) '
                                         'axes: SNP effect on the factor and Q_SNP heterogeneity')
    run.table('factor_gwas_summary', pd.DataFrame(summ), note='Factor GWAS: SNPs, effective N, lambda_GC, genome-wide '
                                                               'significant SNPs and clumped leads, Q_SNP-significant SNPs')
    h = pd.concat(harm, ignore_index=True).rename(columns={'beta_x': 'bx', 'se_x': 'bxse', 'beta_y': 'by', 'se_y': 'byse'})
    fi, fo = Path(tmp) / 'h.tsv', Path(tmp) / 'r.tsv'
    h[['exposure', 'outcome', 'cohort', 'axis', 'rsid', 'bx', 'bxse', 'by', 'byse']].to_csv(fi, sep='\t', index=False)
    r = subprocess.run([C.RSCRIPT, str(ROOT / 'scripts' / '31_mr.R'), str(fi), str(fo)], capture_output=True, text=True)
    if r.returncode:
        raise SystemExit(r.stderr[-800:])
    res = pd.read_csv(fo, sep='\t')
    run.table('mr_factor_results', res, note='MR of five exposures on the latent axes (SD of factor per SD exposure)')
    # TEN minus FIB per exposure (IVW); the two factor GWAS share samples with a positive cross-trait LDSC intercept, so treating
    # the estimates as independent overstates the SE of the difference (conservative)
    iv = res[res.method == 'IVW'].pivot(index='exposure', columns='axis', values=['estimate', 'se'])
    if {'TEN', 'FIB'} <= set(iv.columns.get_level_values(1)):
        c = pd.DataFrame({'diff': iv['estimate', 'TEN'] - iv['estimate', 'FIB'],
                          'se': np.hypot(iv['se', 'TEN'], iv['se', 'FIB'])}).reset_index()
        c['lo'], c['hi'] = c['diff'] - 1.96 * c.se, c['diff'] + 1.96 * c.se
        c['p'] = 2 * stats.norm.sf(np.abs(c['diff'] / c.se))
        run.table('mr_factor_contrast', c, note='IVW effect on the tenosynovial minus the fibro-capsular factor per exposure '
                                                '(SD per SD; estimates treated as independent, conservative)')
        print(c.round(4).to_string(index=False))
    pd.set_option('display.width', 200)
    print(pd.DataFrame(summ).to_string(index=False))
    print(res[res.method.isin(['IVW', 'weighted median'])].round(4).to_string(index=False))
