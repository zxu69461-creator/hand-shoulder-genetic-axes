#!/usr/bin/env python3
"""Disease-level genetic correlations and their factor structure.

Source: gsem_cfa_est, set xc, model 'dis' (each disorder measured in FinnGen and in UK Biobank-based GWAS forms a disease
factor; an orthogonal method factor absorbs the control group shared by the FinnGen M13 endpoints). The disease-factor
correlations are written as `disease_rg`; exploratory factor analysis (1-2 factors, promax; R factanal on the
nearest positive-definite matrix) as `disease_efa`, with eigenvalues in `disease_eigen`.
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
from project import Project

P = Project(__file__)
ORDER = ['dCTS', 'dTF', 'dDQ', 'dENTH', 'dPF', 'dFS', 'dDUP']
R_CODE = r'''
a <- commandArgs(TRUE); R <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE))
R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
out <- do.call(rbind, lapply(1:2, function(k) {
  f <- factanal(covmat = R, factors = k, rotation = if (k > 1) "promax" else "none", n.obs = 1000)
  L <- unclass(f$loadings)
  do.call(rbind, lapply(seq_len(ncol(L)), function(j) data.frame(k = k, factor = paste0("F", j), disease = rownames(L), loading = L[, j])))
}))
write.table(out, a[2], sep = "\t", quote = FALSE, row.names = FALSE)
'''
with P.run(__file__, params={'set': 'xc', 'model': 'dis', 'efa': 'factanal promax, nearPD'}) as run:
    with P.connect(read_only=True) as con:
        e = con.execute("""select lhs, rhs, STD_Genotype as rg, STD_Genotype_SE as se, p_value from gsem_cfa_est
                           where set='xc' and model='dis' and op='~~' and lhs<>rhs and lhs like 'd%' and rhs like 'd%'""").df()
    assert len(e) == 21, f'expected 21 disease pairs, got {len(e)}'
    R = pd.DataFrame(np.eye(len(ORDER)), index=ORDER, columns=ORDER)
    for r in e.itertuples():
        R.loc[r.lhs, r.rhs] = R.loc[r.rhs, r.lhs] = r.rg
    eig = np.sort(np.linalg.eigvalsh(R.values))[::-1]
    with tempfile.TemporaryDirectory() as tmp:
        fi, fo, fr = Path(tmp) / 'R.tsv', Path(tmp) / 'efa.tsv', Path(tmp) / 'efa.R'
        R.to_csv(fi, sep='\t'); fr.write_text(R_CODE)
        r = subprocess.run([C.RSCRIPT, str(fr), str(fi), str(fo)], capture_output=True, text=True)
        if r.returncode:
            sys.exit(r.stderr[-800:])
        efa = pd.read_csv(fo, sep='\t')
    run.table('disease_rg', e.assign(disease1=e.lhs.str[1:], disease2=e.rhs.str[1:])[['disease1', 'disease2', 'rg', 'se', 'p_value']],
              note='Disease-level genetic correlations (FinnGen + UK Biobank disease factors; FinnGen shared-control method factor removed)')
    run.table('disease_efa', efa.assign(disease=efa.disease.str[1:]),
              note='Exploratory factor loadings of the disease-level genetic correlation matrix (1-2 factors, promax)')
    run.table('disease_eigen', pd.DataFrame({'component': range(1, len(eig) + 1), 'eigenvalue': eig}),
              note='Eigenvalues of the disease-level genetic correlation matrix')
    print(efa[efa.k == 2].pivot(index='disease', columns='factor', values='loading').round(2))
    print('eigenvalues', np.round(eig, 2))
