#!/usr/bin/env python3
"""Disorder specificity of a cell state's heritability enrichment: for a state, the MAGMA gene-set beta of one disorder
minus that of each other disorder, z = difference / sqrt(se1^2 + se2^2) (the disorders' GWAMAs share some samples, so the
z is approximate). States: lesional POSTN+ fibroblasts (palmar fascia atlas) with Dupuytren's disease as the reference,
COMP+ THBS4+ fascial fibroblasts with carpal tunnel syndrome. Table: cell_contrast."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import numpy as np
import pandas as pd
from project import Project
from scipy import stats

P = Project(__file__)
REF = [('dupuytren_mes', 'POSTN+ Dupuytren fibroblasts', 'DUP'), ('dupuytren_mes', 'COMP+ THBS4+ fascial fibroblasts', 'CTS')]
with P.run(__file__, params={'states': REF}) as run:
    with P.connect(read_only=True) as con:
        g = con.execute('select * from magma_gsa').df()
    run.input('magma_gsa', note='MAGMA competitive gene-set results per atlas cell state x disorder')
    rows = []
    for ds, ct, ref in REF:
        s = g[(g.dataset == ds) & (g.cell_type == ct)].set_index('disease')
        for d in s.index:
            if d == ref:
                continue
            diff = s.loc[ref, 'beta'] - s.loc[d, 'beta']
            se = np.hypot(s.loc[ref, 'se'], s.loc[d, 'se'])
            rows.append(dict(dataset=ds, cell_type=ct, reference=ref, other=d, diff=diff, se=se, z=diff / se,
                             p=2 * stats.norm.sf(abs(diff / se))))
    out = pd.DataFrame(rows)
    run.table('cell_contrast', out, note='Difference in MAGMA cell-state enrichment beta between a reference disorder and '
                                         'each other disorder (approximate z; GWAMAs partly share samples)')
    print(out.round(4).to_string(index=False))
