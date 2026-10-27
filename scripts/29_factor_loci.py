#!/usr/bin/env python3
"""Lead SNPs and loci of the two axis factor GWAS, and whether each locus is new relative to the disorder GWAS.

Leads: PLINK clumping of factor SNPs with P < 5e-8 (r2 < 0.1 within 1 Mb, 1000 Genomes EUR; as 26_factor_downstream.py).
Loci: leads merged when within 500 kb of each other. A locus is new when no SNP within 500 kb of any of its leads reaches
P < 5e-8 in any of the axis's three indicator GWAS (the cross-cohort meta-analyses that entered the factor model, same
HapMap 3 SNP set), so novelty is judged on identical SNPs. Nearest gene: NCBI build 37.3 gene locations (MAGMA), the
gene whose body is closest to the lead (distance 0 inside the gene).
Tables: factor_leads (axis x lead SNP), factor_loci (axis x locus).
"""
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import numpy as np
import pandas as pd
import config as C
import mrlib as L
from project import Project

P = Project(__file__)
STAGE = C.WORK / 'factor_stage'
AXES = {'TEN': ['CTS', 'TF', 'DQ'], 'FIB': ['DUP', 'FS', 'PF']}
LOC = C.GENE_LOC
genes = pd.read_csv(LOC, sep='\t', header=None, names=['entrez', 'chr', 'start', 'end', 'strand', 'symbol'])
genes['chr'] = genes.chr.astype(str)
WIN = 500_000


def nearest(c, bp):
    g = genes[genes.chr == str(c)]
    d = np.maximum(0, np.maximum(g.start - bp, bp - g.end))
    i = int(np.argmin(d.values))
    return g.symbol.iloc[i], int(d.iloc[i])


with P.run(__file__, params={'win': WIN}) as run, tempfile.TemporaryDirectory() as tmp:
    with P.connect() as con:
        fac = {ax: con.execute(f"select SNP, CHR, BP, A1, A2, est, SE, P, Q_pval from factor_gwas where axis = '{ax}'").df()
               for ax in AXES}
    leads, loci = [], []
    for ax, ind in AXES.items():
        f = fac[ax]
        sig = f[f.P < 5e-8]
        a = Path(tmp) / f'{ax}.assoc'; sig[['SNP', 'P']].to_csv(a, sep='\t', index=False)
        subprocess.run([L.PLINK, '--bfile', L.REF, '--clump', str(a), '--clump-p1', '5e-8', '--clump-p2', '1e-5',
                        '--clump-r2', '0.1', '--clump-kb', '1000', '--out', str(Path(tmp) / ax)], capture_output=True, check=True)
        cl = pd.read_csv(Path(tmp) / f'{ax}.clumped', sep=r'\s+')
        ld = f.set_index('SNP').loc[cl.SNP].reset_index().sort_values(['CHR', 'BP'])
        # merge leads within WIN into loci
        ld['locus'] = (ld.groupby('CHR').BP.diff().fillna(np.inf) > WIN).groupby(ld.CHR).cumsum()
        ld['locus'] = ld.CHR.astype(str) + ':' + ld.locus.astype(int).astype(str)
        # indicator GWAS: genome-wide significant SNPs, placed on the factor SNP positions (same HapMap 3 set)
        pos = f[['SNP', 'CHR', 'BP']]
        ind_hits = []
        for t in ind:
            s = pd.read_csv(STAGE / f'{t}_ss.tsv.gz', sep='\t', usecols=['SNP', 'P'])
            run.input(f'indicator:{t}', path=str(STAGE / f'{t}_ss.tsv.gz'), note='cross-cohort meta-analysis, HapMap 3 (24_factor_stage.py)')
            ind_hits.append(s[s.P < 5e-8].merge(pos, on='SNP').assign(trait=t))
        ih = pd.concat(ind_hits, ignore_index=True)
        rows = []
        for r in ld.itertuples():
            near = ih[(ih.CHR == r.CHR) & (abs(ih.BP - r.BP) <= WIN)]
            g, dist = nearest(r.CHR, r.BP)
            rows.append(dict(axis=ax, snp=r.SNP, chr=int(r.CHR), bp=int(r.BP), a1=r.A1, a2=r.A2, beta=r.est, se=r.SE, p=r.P,
                             q_p=r.Q_pval, locus=r.locus, nearest_gene=g, gene_dist=dist,
                             min_p_indicator=float(near.P.min()) if len(near) else np.nan,
                             indicators_gws=','.join(sorted(set(near.trait)))))
        lt = pd.DataFrame(rows)
        leads.append(lt)
        lo = lt.sort_values('p').groupby('locus', sort=False).agg(
            chr=('chr', 'first'), lead_snp=('snp', 'first'), bp=('bp', 'first'), p=('p', 'first'), beta=('beta', 'first'),
            nearest_gene=('nearest_gene', 'first'), n_leads=('snp', 'size'), start=('bp', 'min'), end=('bp', 'max'),
            indicators_gws=('indicators_gws', lambda x: ','.join(sorted(set(','.join(x).split(',')) - {''})))).reset_index()
        lo['new'] = lo.indicators_gws == ''
        loci.append(lo.assign(axis=ax))
        print(ax, 'leads', len(lt), 'loci', len(lo), 'new loci', int(lo.new.sum()))
    leads = pd.concat(leads, ignore_index=True); loci = pd.concat(loci, ignore_index=True)
    run.table('factor_leads', leads, note='Factor GWAS lead SNPs (clumped r2 < 0.1, 1 Mb), nearest gene, minimum P within 500 kb '
                                          'in the three indicator GWAS and which indicators are genome-wide significant there')
    run.table('factor_loci', loci, note='Factor GWAS loci (leads within 500 kb merged); new = no indicator GWAS SNP with '
                                        'P < 5e-8 within 500 kb of any lead')
    pd.set_option('display.width', 220)
    print(loci[loci.new].sort_values(['axis', 'p'])[['axis', 'lead_snp', 'chr', 'bp', 'p', 'beta', 'nearest_gene']].to_string(index=False))
