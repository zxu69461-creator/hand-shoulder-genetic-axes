#!/usr/bin/env python3
"""Colocalisation of disorder signals with their endocrine / metabolic drivers at candidate loci.

Pairs (gene window +-500 kb, GRCh38): disorder = FinnGen R13 endpoint (complete regional statistics);
trait = GWAS Catalog harmonised exposure file (hm_pos GRCh38). Variants matched on chromosome, position and alleles.
coloc.abf (default priors p1 = p2 = 1e-4, p12 = 1e-5) in R (scripts/35_coloc.R). Table: coloc_results.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import duckdb
import pandas as pd
import config as C
import mrlib as L
from project import Project

P = Project(__file__)
GENE38 = {'SHBG': (17, 7_610_000), 'IGF1R': (15, 98_650_000), 'FTO': (16, 53_800_000), 'PAPPA2': (1, 176_500_000),
          'WNT7B': (22, 45_950_000), 'SFRP4': (7, 37_900_000), 'EPDR1': (7, 37_950_000)}
PAIRS = [('SHBG', 'CTS', 'SHBG_F'), ('SHBG', 'TF', 'SHBG_F'), ('IGF1R', 'CTS', 'IGF1'), ('PAPPA2', 'CTS', 'IGF1'),
         ('PAPPA2', 'TF', 'IGF1'), ('FTO', 'CTS', 'BMI_GIANT'), ('FTO', 'TF', 'BMI_GIANT'), ('FTO', 'DUP', 'BMI_GIANT')]
FG = {o[0]: o for o in L.OUTCOMES if o[1] == 'FinnGen'}
CASES = {'CTS': (32951, 431627), 'TF': (6528, 350578), 'DUP': (6972, 350578)}
N_EXP = {'SHBG_F': 189473, 'IGF1': 340567, 'BMI_GIANT': 322154}
W = 500_000
rows = []
with P.run(__file__, params={'window': W, 'priors': 'coloc.abf defaults'}) as run, tempfile.TemporaryDirectory() as tmp:
    for gene, dis, expo in PAIRS:
        ch, pos = GENE38[gene]
        con = duckdb.connect(); con.execute(f"SET threads={C.THREADS}; SET enable_progress_bar=false")
        d = con.execute(f"""SELECT chrom, pos, upper(ref) AS ref, upper(alt) AS alt, beta AS b1, se AS s1, af AS maf1
            FROM read_parquet('{FG[dis][4]}') WHERE chrom = {ch} AND pos BETWEEN {pos - W} AND {pos + W}""").df()
        _, kind, src, acc, _ = L.EXPO[expo]
        x = L.read_gwas(Path(src).resolve(), where=f'chrom = {ch} AND pos BETWEEN {pos - W} AND {pos + W}')
        m = d.merge(x, left_on=['chrom', 'pos'], right_on=['chrom', 'pos'])
        same = (m.alt == m.ea) & (m.ref == m.oa); flip = (m.alt == m.oa) & (m.ref == m.ea)
        m = m[same | flip].copy(); m.loc[flip[same | flip].values, 'beta'] *= -1
        m = m.dropna(subset=['b1', 's1', 'beta', 'se']).query('s1 > 0 and se > 0')
        m['maf'] = m.maf1.where(m.maf1 <= 0.5, 1 - m.maf1)
        m = m[m.maf > 0.001].drop_duplicates(['chrom', 'pos'])
        f = Path(tmp) / f'{gene}_{dis}.tsv'
        m[['chrom', 'pos', 'b1', 's1', 'beta', 'se', 'maf']].to_csv(f, sep='\t', index=False)
        ca, co = CASES[dis]
        r = subprocess.run([C.RSCRIPT, str(ROOT / 'scripts' / '35_coloc.R'), str(f), str(ca / (ca + co)), str(ca + co), str(N_EXP[expo])],
                           capture_output=True, text=True)
        if r.returncode:
            print(f'[fail] {gene} {dis}: {r.stderr[-300:]}'); continue
        pp = dict(kv.split('=') for kv in r.stdout.strip().split() if '=' in kv)
        top = pp.pop('top_snp'); tpos = int(top.split(':')[1])
        rs = con.execute(f"SELECT rsid FROM read_parquet('{FG[dis][4]}') WHERE chrom = {ch} AND pos = {tpos} LIMIT 1").fetchone()
        rows.append(dict(gene=gene, disease=dis, exposure=expo, n_snp=len(m), **{k: float(v) for k, v in pp.items()},
                         top_pos_grch38=tpos, top_rsid=rs[0] if rs else None))
        print(rows[-1], flush=True)
        run.input(f'coloc:{gene}:{dis}:{expo}', note=f'{acc} x FinnGen {dis}')
    run.table('coloc_results', pd.DataFrame(rows), note='coloc.abf of FinnGen disorder signals with exposure GWAS at candidate loci '
                                                       '(+-500 kb, GRCh38, default priors): posterior probabilities H0-H4 and the variant with the largest '
                                                       'SNP-level posterior under H4')
