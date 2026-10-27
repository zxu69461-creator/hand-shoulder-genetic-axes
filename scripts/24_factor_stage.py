#!/usr/bin/env python3
"""Inputs for the axis factor GWAS (scripts/25_factor_gwas.R).

Indicators: the cross-cohort meta-analyses of each disorder (50_gwama.py; HapMap 3 SNPs, Z and N) for the tenosynovial
(CTS, TF, DQ) and fibro-capsular (DUP, FS, PF) axes. Writes into work/factor_stage: <T>.sumstats.gz (LDSC input),
<T>_ss.tsv.gz (SNP A1 A2 Z P N, for GenomicSEM::sumstats with linprob) and ref_maf.txt (SNP CHR BP MAF A1 A2 from the
1000 Genomes European reference, restricted to the indicator SNPs).
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import pandas as pd
from scipy import stats
import config as C

GW = C.GWAMA
REF = C.G1000
TRAITS = ['CTS', 'TF', 'DQ', 'DUP', 'FS', 'PF']
out = C.WORK / 'factor_stage'; out.mkdir(parents=True, exist_ok=True)
snps = None
for t in TRAITS:
    d = pd.read_csv(GW / f'{t}.sumstats.gz', sep='\t').dropna()
    snps = set(d.SNP) if snps is None else snps | set(d.SNP)
keep = sorted(snps)
keep_s = set(keep)
for t in TRAITS:
    d = pd.read_csv(GW / f'{t}.sumstats.gz', sep='\t').dropna()
    d = d[d.SNP.isin(keep_s)]
    d.to_csv(out / f'{t}.sumstats.gz', sep='\t', index=False)
    d.assign(P=2 * stats.norm.sf(d.Z.abs()))[['SNP', 'A1', 'A2', 'Z', 'P', 'N']].to_csv(out / f'{t}_ss.tsv.gz', sep='\t', index=False)
tmp = out / 'freq'
(out / 'keep.txt').write_text('\n'.join(keep) + '\n')
subprocess.run([C.PLINK, '--bfile', REF, '--extract', str(out / 'keep.txt'), '--freq', '--out', str(tmp)], check=True,
               capture_output=True)
f = pd.read_csv(str(tmp) + '.frq', sep=r'\s+')
bim = pd.read_csv(REF + '.bim', sep='\t', header=None, names=['CHR', 'SNP', 'cm', 'BP', 'a1', 'a2'])
ref = f.merge(bim[['SNP', 'BP']], on='SNP')[['SNP', 'CHR', 'BP', 'MAF', 'A1', 'A2']].dropna()
ref = ref[(ref.MAF > 0) & (ref.MAF < 1)]
ref.to_csv(out / 'ref_maf.txt', sep='\t', index=False)
for p in out.glob('freq*'):
    p.unlink()
(out / 'keep.txt').unlink()
print(len(keep), 'SNPs;', len(ref), 'with reference MAF')
