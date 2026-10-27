#!/usr/bin/env python3
"""Sex-stratified Mendelian randomization without sample overlap between exposure and outcome.

(a) BMI (GIANT, Locke 2015; no UK Biobank participants) -> UK Biobank outcomes in women and in men separately (Neale
    lab round 2, hospital main-diagnosis ICD-10: G56 mononeuropathies of the upper limb, chiefly carpal tunnel syndrome;
    M72 fibroblastic disorders, chiefly Dupuytren's disease; M65 synovitis and tenosynovitis). Linear-probability
    effects are converted to log-odds with the sex-specific case fraction; variants are matched on GRCh37 position and
    alleles (rsID -> position from the 1000 Genomes EUR reference).
(b) Sex hormones measured in men (Ruth 2020, UK Biobank: SHBG GCST90012109, bioavailable testosterone GCST90012103,
    total testosterone GCST90012113) -> FinnGen outcomes, beside the women's instruments used in the main analysis.
Instruments: P < 5e-8, r2 < 0.001 within 10 Mb (mrlib.clump). Estimators: IVW (random effects), weighted median,
MR-Egger (scripts/31_mr.R). Tables: mr_sex_results, mr_sex_harmonised.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import duckdb
import numpy as np
import pandas as pd
import config as C
import mrlib as L
from project import Project

P = Project(__file__)
NEALE = {(c, s): C.NEALE / f'{c}.gwas.imputed_v3.{s}.tsv.bgz' for c in ('G56', 'M72') for s in ('female', 'male')}
CASES = {('G56', 'female'): (5492, 188682), ('G56', 'male'): (2638, 164382), ('M72', 'female'): (801, 193373),
         ('M72', 'male'): (2392, 164628), ('M65', 'female'): (1704, 192470), ('M65', 'male'): (1108, 165912)}
OUTNAME = {'G56': 'CTS', 'M72': 'DUP', 'M65': 'TENOSYN'}
MALE = {'SHBG_M': C.GWAS / '32042192-GCST90012109-EFO_0004696.h.tsv.gz',
        'BIOT_M': C.GWAS / '32042192-GCST90012103-EFO_0004908.h.tsv.gz',
        'TT_M': C.GWAS / '32042192-GCST90012113-EFO_0004908.h.tsv.gz'}


def neale_rows(path, pos37, mu):
    """rows at the wanted GRCh37 positions; Neale beta is for the alt allele (variant = chr:pos:ref:alt)"""
    con = duckdb.connect(); con.execute(f'SET threads={C.THREADS}')
    con.register('want', pos37)
    df = con.execute(f"""
        WITH v AS (SELECT string_split(variant, ':') s, beta, se, minor_AF, low_confidence_variant AS lc
                   FROM read_csv('{path}', delim='\t', header=true, compression='gzip'))
        SELECT w.rsid, upper(s[4]) AS ea, upper(s[3]) AS oa, beta, se, minor_AF AS eaf, lc
        FROM v JOIN want w ON s[1] = w.chr AND CAST(s[2] AS BIGINT) = w.pos""").df()
    df = df[~df.lc.astype(str).str.lower().eq('true')]
    df = df.assign(beta=df.beta / (mu * (1 - mu)), se=df.se / (mu * (1 - mu)))
    return df.dropna(subset=['beta', 'se']).query('se > 0').drop_duplicates('rsid')


with P.run(__file__, params={'neale': [str(p) for p in NEALE.values()], 'male': {k: str(v) for k, v in MALE.items()}}) as run, \
        tempfile.TemporaryDirectory() as tmp:
    bim = pd.read_csv(L.REF + '.bim', sep='\t', header=None, names=['chr', 'rsid', 'cm', 'pos', 'a1', 'a2'])
    bim['chr'] = bim.chr.astype(str)
    harm = []
    # (a) BMI -> sex-specific UK Biobank outcomes
    x, _ = L.clump(L.exposure_gws('BMI_GIANT'), tmp, 'BMI_GIANT')
    pos = bim[bim.rsid.isin(x.rsid)][['rsid', 'chr', 'pos']]
    for (code, sex), path in NEALE.items():
        ca, co = CASES[(code, sex)]
        y = neale_rows(path, pos, ca / (ca + co))
        m = L.harmonise(x[['rsid', 'ea', 'oa', 'beta', 'se']], y[['rsid', 'ea', 'oa', 'beta', 'se']])
        harm.append(m.assign(exposure='BMI_GIANT', outcome=OUTNAME[code], cohort=f'UKB {sex}', axis=code))
        run.input(f'neale:{code}:{sex}', path=str(path), note=f'Neale lab UKB round 2 {code} {sex} ({ca} cases)')
    # (b) male sex-hormone instruments -> FinnGen outcomes
    for ex, path in MALE.items():
        g = L.read_gwas(path, where='p < 5e-8')
        xm, _ = L.clump(g, tmp, ex)
        run.input(ex, path=str(path), note='Ruth 2020, UK Biobank men')
        for out, coh, axis, kind, opath, mu in L.OUTCOMES:
            if coh != 'FinnGen' or out not in ('CTS', 'TF', 'DQ', 'DUP', 'FS', 'PF'):
                continue
            y = L.outcome_rows(kind, opath, set(xm.rsid), mu)
            m = L.harmonise(xm[['rsid', 'ea', 'oa', 'beta', 'se']], y[['rsid', 'ea', 'oa', 'beta', 'se']])
            harm.append(m.assign(exposure=ex, outcome=out, cohort='FinnGen', axis=axis))
    h = pd.concat(harm, ignore_index=True).rename(columns={'beta_x': 'bx', 'se_x': 'bxse', 'beta_y': 'by', 'se_y': 'byse'})
    f_in, f_out = Path(tmp) / 'h.tsv', Path(tmp) / 'r.tsv'
    h[['exposure', 'outcome', 'cohort', 'axis', 'rsid', 'bx', 'bxse', 'by', 'byse']].to_csv(f_in, sep='\t', index=False)
    r = subprocess.run([C.RSCRIPT, str(ROOT / 'scripts' / '31_mr.R'), str(f_in), str(f_out)], capture_output=True, text=True)
    if r.returncode:
        raise SystemExit(r.stderr[-800:])
    res = pd.read_csv(f_out, sep='\t')
    res['odds_ratio'] = np.exp(res.estimate)
    res['or_lo'] = np.exp(res.estimate - 1.96 * res.se); res['or_hi'] = np.exp(res.estimate + 1.96 * res.se)
    run.table('mr_sex_harmonised', h, note='harmonised SNP effects for sex-stratified MR')
    run.table('mr_sex_results', res, note='Sex-stratified MR: BMI -> UK Biobank outcomes by sex (Neale round 2); male '
                                          'sex-hormone instruments -> FinnGen outcomes; IVW / weighted median / MR-Egger')
    pd.set_option('display.width', 220)
    print(res[res.method == 'IVW'][['exposure', 'outcome', 'cohort', 'n_snp', 'odds_ratio', 'or_lo', 'or_hi', 'p']].round(4).to_string(index=False))
