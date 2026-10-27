#!/usr/bin/env python3
"""Two-factor EFA of the disease-level genetic correlations (as 22_disease_level.py: factanal, promax, nearest
positive-definite matrix): uniquenesses, communalities and the factor correlation, plus a sensitivity solution with the
uniquenesses bounded at 0.05 instead of factanal's default 0.005. Table: disease_efa_detail.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import pandas as pd
import config as C
from project import Project

P = Project(__file__)
V = ['CTS', 'TF', 'DQ', 'ENTH', 'PF', 'FS', 'DUP']
R_CODE = r'''
a <- commandArgs(TRUE); d <- read.delim(a[1]); v <- strsplit(a[3], ",")[[1]]
R <- diag(length(v)); dimnames(R) <- list(v, v)
for (i in seq_len(nrow(d))) { R[d$disease1[i], d$disease2[i]] <- d$rg[i]; R[d$disease2[i], d$disease1[i]] <- d$rg[i] }
R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
out <- do.call(rbind, lapply(c(0.005, 0.05), function(lb) {
  f <- factanal(covmat = R, factors = 2, rotation = "promax", n.obs = 1000, lower = lb)
  L <- unclass(f$loadings); Phi <- solve(crossprod(f$rotmat))
  data.frame(lower_bound = lb, disease = rownames(L), loading_F1 = L[, 1], loading_F2 = L[, 2],
             uniqueness = f$uniquenesses, communality = rowSums((L %*% Phi) * L), factor_correlation = Phi[1, 2])
}))
write.table(out, a[2], sep = "\t", quote = FALSE, row.names = FALSE)
'''
with P.run(__file__, params={'lower_bounds': [0.005, 0.05]}) as run, tempfile.TemporaryDirectory() as tmp:
    with P.connect(read_only=True) as con:
        con.execute(f"copy (select * from disease_rg) to '{tmp}/rg.tsv' (delimiter '\t', header)")
    (Path(tmp) / 'efa.R').write_text(R_CODE)
    subprocess.run([C.RSCRIPT, f'{tmp}/efa.R', f'{tmp}/rg.tsv', f'{tmp}/out.tsv', ','.join(V)], check=True)
    out = pd.read_csv(f'{tmp}/out.tsv', sep='\t')
    run.table('disease_efa_detail', out, note='Two-factor promax EFA: loadings, uniquenesses, communalities, factor correlation; '
                                              'default uniqueness bound 0.005 and sensitivity bound 0.05')
    print(out.round(3).to_string(index=False))
