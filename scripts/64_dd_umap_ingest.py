#!/usr/bin/env python3
"""The authors' UMAP of the Dupuytren's / palmar fascia / dermis atlas (GSE173252; written by 63_dd_umap.R to
work/sc/dd_umap.tsv.gz) with the labels of 62_dd_means.R. Table dd_umap (one row per cell; doublets flagged, not dropped)."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import pandas as pd
import config as C
from project import Project

P = Project(__file__)
with P.run(__file__) as run:
    run.input('GSE173252', path=str(C.SC / 'GSE173252_dd.rds.gz'), note="authors' Seurat object")
    u = pd.read_csv(C.WORK / 'sc' / 'dd_umap.tsv.gz', sep='\t')
    u['doublet'] = u.doublet.astype(str).str.upper() == 'TRUE'
    run.table('dd_umap', u, note="GSE173252 authors' UMAP per cell: broad annotation, tissue source (DD / SF fascia / "
                                 "dermis), sample, Seurat cluster, mesenchymal state names (62_dd_means.R), doublet flag")
    print(u.groupby('source').size(), u.doublet.sum())
