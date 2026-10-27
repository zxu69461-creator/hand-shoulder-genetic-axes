#!/usr/bin/env python3
"""Cell-type expression specificity of public connective-tissue single-cell atlases, for MAGMA cell-type enrichment.

Per dataset and annotated cell type (author annotation, types with >= 50 cells/nuclei): mean of library-size-normalised
counts (counts per 10,000, linear scale) -> specificity = mean in type / sum of means over types (Bryois et al. 2020).
Per type, the 10% most specific genes among genes expressed in that type (mean >= 0.1 CP10k, i.e. detectable) and mapped
to NCBI Entrez IDs (MAGMA NCBI37.3 annotation, by symbol) form its gene set.
Inputs: CELLxGENE h5ad files (raw counts in .raw), or a per-type mean matrix exported from a Seurat object
(scripts/62_dd_means.R). Outputs: table sc_specificity and the MAGMA set file
work/magma_sets/<dataset>.sets. The six runs of the study are listed in scripts/60_sc_sets.sh.

Usage: 60_sc_specificity.py <dataset> <h5ad | means.tsv.gz> [annotation column] [regex of cell types to keep]
  The optional regex restricts the specificity to a subset of types (e.g. fibroblast states only), so that a state's
  gene set reflects what distinguishes it from the other states rather than the shared fibroblast programme.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import config as C
from project import Project

P = Project(__file__)
LOC = C.GENE_LOC
MIN_CELLS, MIN_MEAN, TOP = 50, 0.1, 0.10
name, src = sys.argv[1], Path(sys.argv[2])
col = sys.argv[3] if len(sys.argv) > 3 else 'author_cell_type'
keep = sys.argv[4] if len(sys.argv) > 4 else None


def means_from_h5ad(path, col):
    a = ad.read_h5ad(path)
    X = a.raw.X if a.raw is not None else a.X
    genes = (a.raw.var if a.raw is not None else a.var)
    sym = genes['feature_name'].astype(str).values if 'feature_name' in genes else genes.index.values
    X = sp.csr_matrix(X)
    lib = np.asarray(X.sum(1)).ravel()
    X = sp.diags(1e4 / np.maximum(lib, 1)) @ X
    lab = a.obs[col].astype(str).values
    out, counts = {}, {}
    for t in pd.unique(lab):
        m = lab == t
        if m.sum() < MIN_CELLS:
            continue
        out[t] = np.asarray(X[m].mean(0)).ravel(); counts[t] = int(m.sum())
    M = pd.DataFrame(out, index=sym)
    return M.groupby(level=0).sum(), counts


if src.suffix == '.h5ad':
    M, counts = means_from_h5ad(src, col)
else:
    M = pd.read_csv(src, sep='\t', index_col=0)
    n = pd.read_csv(str(src).replace('.means.tsv.gz', '.counts.tsv'), sep='\t', index_col=0).iloc[:, 0]
    counts = {k: int(v) for k, v in n.items() if v >= MIN_CELLS}
    M = M[[c for c in M.columns if c in counts]]
if keep:
    M = M[[c for c in M.columns if pd.Series([c]).str.contains(keep, case=False, regex=True).iloc[0]]]
loc = pd.read_csv(LOC, sep='\t', header=None, names=['entrez', 'chr', 'start', 'end', 'strand', 'symbol'])
entrez = loc.drop_duplicates('symbol').set_index('symbol').entrez
S = M.div(M.sum(1).replace(0, np.nan), axis=0)
rows, sets = [], []
for t in M.columns:
    d = pd.DataFrame({'gene': M.index, 'mean_cp10k': M[t].values, 'specificity': S[t].values})
    d = d[(d.mean_cp10k >= MIN_MEAN) & d.gene.isin(entrez.index)].dropna()
    cut = d.specificity.quantile(1 - TOP)
    d['top10'] = d.specificity >= cut
    d.insert(0, 'cell_type', t); d.insert(0, 'dataset', name); d['n_cells'] = counts[t]
    d['entrez'] = d.gene.map(entrez).astype(int)
    rows.append(d)
    sets.append(f"{name}:{t.replace(' ', '_')} " + ' '.join(map(str, d[d.top10].entrez)))
res = pd.concat(rows, ignore_index=True)
sd = C.MAGMA_SETS; sd.mkdir(parents=True, exist_ok=True)
(sd / f'{name}.sets').write_text('\n'.join(sets) + '\n')
with P.run(__file__, params={'dataset': name, 'annotation': col, 'subset': keep, 'min_cells': MIN_CELLS, 'min_mean_cp10k': MIN_MEAN, 'top': TOP}) as run:
    run.input(name, path=str(src.resolve()), note=f'single-cell reference ({col})')
    with P.connect(read_only=True) as con:
        have = {r[0] for r in con.execute("select table_name from duckdb_tables()").fetchall()} | \
               {r[0] for r in con.execute("select view_name from duckdb_views()").fetchall()}
        old = con.execute('select * from sc_specificity').df() if 'sc_specificity' in have else pd.DataFrame()
    allr = pd.concat([old[old.dataset != name], res], ignore_index=True) if len(old) else res
    run.parquet('sc_specificity', allr, note='Per dataset x cell type x gene: mean CP10k, specificity (share of summed type '
                                             'means), top-10% flag, Entrez ID; gene sets in work/magma_sets')
    print(name, {t: counts[t] for t in M.columns}, f'{len(res)} rows')
