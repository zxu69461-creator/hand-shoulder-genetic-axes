#!/usr/bin/env python3
"""Disease-level cross-cohort GWAS meta-analysis (FinnGen R13 + non-Finnish GWAS) on HapMap3 SNPs.

Per disease: sample-size-weighted Z (Stouffer, weights sqrt(N_eff)), alleles aligned to the first file; the cohorts do not
share samples (FinnGen vs UK Biobank / non-Finnish meta), so no overlap correction is needed. Single-cohort disorders
are passed through. Output: work/gwama/<D>.sumstats.gz (SNP A1 A2 Z N; LDSC/MAGMA ready) and table gwama_inputs.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import duckdb
import pandas as pd
import config as C
from project import Project

P = Project(__file__)
FGM = SFM = C.MUNGED
OUT = C.GWAMA
DISEASES = {
    'CTS': [FGM / 'G6_CARPTU.sumstats.gz', SFM / 'UKB_CTS.sumstats.gz'],
    'TF': [FGM / 'M13_TRIGGERFINGER.sumstats.gz', SFM / 'UKB_TF.sumstats.gz'],
    'DQ': [FGM / 'M13_DEQUERVAIN.sumstats.gz'],
    'DUP': [FGM / 'M13_DUPUTRYEN.sumstats.gz', SFM / 'DD_META.sumstats.gz'],
    'FS': [FGM / 'M13_ADHCAPSULITIS.sumstats.gz', SFM / 'UKB_FS.sumstats.gz'],
    'PF': [FGM / 'M13_PLANTARFASC.sumstats.gz'],
    'ENTH': [FGM / 'M13_ENTESOPATHYOTH.sumstats.gz', SFM / 'UKB_ENTHESO.sumstats.gz'],
}
OUT.mkdir(parents=True, exist_ok=True)
rows = []
with P.run(__file__, params={'method': 'sqrt(N)-weighted Stouffer, allele-aligned, no overlap correction'}) as run:
    con = duckdb.connect(); con.execute(f"SET threads={C.THREADS}; SET enable_progress_bar=false")
    for d, files in DISEASES.items():
        sel = []
        for i, f in enumerate(files):
            run.input(f.name.replace('.sumstats.gz', ''), path=str(f.resolve()), note=f'{d} input {i + 1}')
            sel.append(f"SELECT {i} AS k, SNP, upper(A1) AS A1, upper(A2) AS A2, Z, N FROM read_csv('{f.resolve()}', delim='\t', header=true)")
        con.execute(f"CREATE OR REPLACE TEMP TABLE s AS {' UNION ALL '.join(sel)}")
        con.execute("""CREATE OR REPLACE TEMP TABLE ref AS SELECT SNP, arg_min(A1, k) AS A1, arg_min(A2, k) AS A2 FROM s GROUP BY SNP""")
        df = con.execute("""
            SELECT r.SNP, r.A1, r.A2,
                   sum(sqrt(s.N) * CASE WHEN s.A1 = r.A1 AND s.A2 = r.A2 THEN s.Z ELSE -s.Z END) / sqrt(sum(s.N)) AS Z,
                   sum(s.N) AS N, count(*) AS n_cohorts
            FROM s JOIN ref r USING (SNP)
            WHERE (s.A1 = r.A1 AND s.A2 = r.A2) OR (s.A1 = r.A2 AND s.A2 = r.A1)
            GROUP BY r.SNP, r.A1, r.A2""").df()
        f = OUT / f'{d}.sumstats.gz'
        df[['SNP', 'A1', 'A2', 'Z', 'N']].to_csv(f, sep='\t', index=False, compression='gzip')
        lam = (df.Z ** 2).median() / 0.4549
        rows.append(dict(disease=d, n_cohorts=len(files), n_snp=len(df), n_snp_all_cohorts=int((df.n_cohorts == len(files)).sum()),
                         max_n=float(df.N.max()), lambda_gc=round(lam, 3), n_gws=int((df.Z.abs() > 5.45).sum()), path=str(f)))
        print(rows[-1], flush=True)
    run.table('gwama_inputs', pd.DataFrame(rows), note='Disease-level cross-cohort GWAMA (FinnGen + non-Finnish) on HapMap3 SNPs')
