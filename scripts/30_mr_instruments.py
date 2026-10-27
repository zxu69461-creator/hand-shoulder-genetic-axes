#!/usr/bin/env python3
"""Genetic instruments for the endocrine and metabolic exposures.

Per exposure (registry in mrlib.EXPOSURES): genome-wide significant variants (p < 5e-8, MAF >= 0.01, outside the MHC
25-34 Mb) -> PLINK 1.9 clumping against 1000 Genomes Europeans (r2 < 0.001 within 10 Mb; rsID-matched) -> table
`mr_instruments` with per-SNP F statistic. Exposures whose file is not downloaded are skipped.
"""
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import pandas as pd
import config as C
import mrlib as L
from project import Project

P = Project(__file__)
rows = []
with P.run(__file__, params={'p': L.P_INST, 'r2': L.R2, 'kb': L.KB, 'ref': '1000 Genomes Europeans', 'mhc': 'excluded'}) as run, \
        tempfile.TemporaryDirectory() as tmp:
    for label, kind, src, acc, note in L.EXPOSURES:
        if not L.available(label):
            print(f'[skip] {label}: not downloaded yet'); continue
        inst, n_gws = L.clump(L.exposure_gws(label), tmp, label)
        inst.insert(0, 'exposure', label)
        rows.append(inst)
        run.input(acc, path=str(Path(src).resolve()), note=note)
        print(f'[done] {label}: {n_gws} genome-wide variants -> {len(inst)} instruments, median F {inst.F.median():.0f}', flush=True)
    if not rows:
        sys.exit('no exposure available')
    run.table('mr_instruments', pd.concat(rows, ignore_index=True),
              note='Clumped genome-wide significant instruments (r2<0.001, 10 Mb, 1000G EUR; MHC excluded) per exposure, with F')
