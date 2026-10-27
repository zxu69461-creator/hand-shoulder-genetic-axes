#!/usr/bin/env python3
"""MAGMA results (51_magma.py) into the results database.

Disorder meta-analyses: magma_genes (disease x gene: Entrez, symbol, chr, n SNPs, Z, P), magma_gsa (disease x atlas x
cell state: beta, SE, P; Benjamini-Hochberg FDR within disease x atlas). Axis factor GWAS (traits FAC_TEN, FAC_FIB):
factor_magma_genes and factor_magma_gsa, same layout.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import numpy as np
import pandas as pd
import config as C
from project import Project

P = Project(__file__)
sym = pd.read_csv(C.GENE_LOC, sep='\t', header=None, names=['entrez', 'chr', 'start', 'end', 'strand', 'symbol']).set_index('entrez').symbol


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1)
    return out


def read(factor):
    genes, gsa = [], []
    keep = (lambda t: t.startswith('FAC_')) if factor else (lambda t: not t.startswith('FAC_'))
    strip = (lambda t: t[4:]) if factor else (lambda t: t)
    for f in sorted(x for x in C.MAGMA_OUT.glob('*.genes.out') if x.name.count('.') == 2):
        t = f.name.split('.')[0]
        if keep(t):
            g = pd.read_csv(f, sep=r'\s+'); g.insert(0, 'disease', strip(t)); genes.append(g)
    for f in sorted(C.MAGMA_OUT.glob('*.gsa.out')):
        t, ds = f.name.split('.')[0], f.name.split('.')[1]
        if not keep(t):
            continue
        g = pd.read_csv(f, sep=r'\s+', comment='#')
        g = g[g.TYPE == 'SET'] if 'TYPE' in g else g
        name = g.FULL_NAME if 'FULL_NAME' in g else g.VARIABLE
        gsa.append(pd.DataFrame({'disease': strip(t), 'dataset': ds, 'cell_type': name.str.split(':', n=1).str[1].str.replace('_', ' '),
                                 'n_genes': g.NGENES, 'beta': g.BETA, 'beta_std': g.BETA_STD, 'se': g.SE, 'p': g.P}))
    return genes, gsa


with P.run(__file__) as run:
    for factor, prefix, what in ((False, 'magma', 'disorder meta-analyses'), (True, 'factor_magma', 'axis factor GWAS')):
        genes, gsa = read(factor)
        if genes:
            g = pd.concat(genes, ignore_index=True)
            g['symbol'] = g.GENE.map(sym)
            g = g.rename(columns={'GENE': 'entrez', 'CHR': 'chr', 'NSNPS': 'n_snps', 'ZSTAT': 'z', 'P': 'p'})
            run.table(f'{prefix}_genes', g[['disease', 'entrez', 'symbol', 'chr', 'n_snps', 'N', 'z', 'p']],
                      note=f'MAGMA gene analysis (SNP-wise mean, 35 kb up / 10 kb down, 1000 Genomes Europeans) of the {what}')
        if gsa:
            s = pd.concat(gsa, ignore_index=True)
            s['fdr'] = s.groupby(['disease', 'dataset']).p.transform(bh)
            run.table(f'{prefix}_gsa', s, note=f'MAGMA competitive gene-set analysis of the {what}: top-10% cell-state-specific '
                                               'genes per atlas; BH-FDR within trait x atlas')
