#!/usr/bin/env python3
"""MAGMA v1.10 gene and competitive gene-set analysis of the disorder meta-analyses and the two axis factor GWAS.

Annotation: 1000 Genomes European reference (GRCh37) and NCBI build 37.3 gene locations, window 35 kb upstream and
10 kb downstream. Gene analysis: SNP-wise mean model on HapMap 3 summary statistics (SNP, Z, N; P = 2 Phi(-|Z|),
floored at 1e-300) with the per-SNP sample size. Gene-set analysis: every trait against every cell-state gene-set file
written by 60_sc_specificity.py (work/magma_sets/<atlas>.sets).
Inputs: work/gwama/<D>.sumstats.gz (50_gwama.py) and work/factor/FAC_<AXIS>.sumstats.gz (26_factor_downstream.py).
Outputs in work/magma: <trait>.genes.out / .genes.raw and <trait>.<atlas>.gsa.out (read by 52_magma_ingest.py).
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import numpy as np
import pandas as pd
from scipy.stats import norm
import config as C

OUT = C.MAGMA_OUT
OUT.mkdir(parents=True, exist_ok=True)
ANNOT = OUT / 'g1000_eur_35_10'


def magma(*args):
    r = subprocess.run([C.MAGMA, *map(str, args)], capture_output=True, text=True)
    if r.returncode:
        raise SystemExit(f'magma {" ".join(map(str, args[:2]))} failed: {r.stdout[-600:]} {r.stderr[-400:]}')


if not Path(f'{ANNOT}.genes.annot').exists():
    magma('--annotate', 'window=35,10', '--snp-loc', f'{C.G1000}.bim', '--gene-loc', C.GENE_LOC, '--out', ANNOT)

traits = {f.name.replace('.sumstats.gz', ''): f for f in sorted(C.GWAMA.glob('*.sumstats.gz'))}
traits.update({f.name.replace('.sumstats.gz', ''): f for f in sorted(C.FACTOR.glob('FAC_*.sumstats.gz'))})
for t, src in traits.items():
    out = OUT / t
    if Path(f'{out}.genes.out').exists():
        continue
    d = pd.read_csv(src, sep='\t', usecols=['SNP', 'Z', 'N']).dropna()
    d = d[np.isfinite(d.Z)]
    d['P'] = np.clip(2 * norm.sf(d.Z.abs()), 1e-300, 1)
    pv = OUT / f'{t}.pval'
    d[['SNP', 'P', 'N']].to_csv(pv, sep='\t', index=False)
    magma('--bfile', C.G1000, '--pval', pv, 'use=SNP,P', 'ncol=N', '--gene-annot', f'{ANNOT}.genes.annot', '--out', out)
    pv.unlink()
    print(t, 'gene analysis done', flush=True)

for t in traits:
    for s in sorted(C.MAGMA_SETS.glob('*.sets')):
        out = OUT / f'{t}.{s.stem}'
        if not Path(f'{out}.gsa.out').exists():
            magma('--gene-results', OUT / f'{t}.genes.raw', '--set-annot', s, '--out', out)
    print(t, 'gene-set analysis done', flush=True)
