#!/usr/bin/env python3
"""Tracks for the FTO regional plot (Figure 4B): LD with the shared variant, gene models and recombination rate.

fto_ld: r2 (and D') of every 1000 Genomes European variant within 1 Mb of rs1558902 (PLINK --r2 --ld-snp on the
  reference used for clumping, GRCh37); joined to the plotted GWAS variants by rsID.
fto_gene_models: GENCODE v46 (GRCh38) Ensembl-canonical transcripts of protein-coding genes overlapping the plotted window,
  one row per exon, with gene strand.
fto_recombination: recombination rate (cM/Mb) between consecutive markers of the Beagle GRCh38 PLINK genetic map
  (HapMap II-derived, lifted to GRCh38) across the window.
Window: chr16:53.70-54.00 Mb (GRCh38), as plotted.
"""
import gzip
import io
import re
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import numpy as np
import pandas as pd
import config as C
import mrlib as L
from project import Project

P = Project(__file__)
LEAD, CH, LO, HI = 'rs1558902', 16, 53_700_000, 54_000_000
GTF = C.GENCODE_GTF
MAP = C.GENETIC_MAP

with P.run(__file__, params={'lead': LEAD, 'window_grch38': [CH, LO, HI]}) as run, tempfile.TemporaryDirectory() as tmp:
    # LD with the shared variant (1000 Genomes EUR, the clumping reference)
    out = Path(tmp) / 'ld'
    subprocess.run([L.PLINK, '--bfile', L.REF, '--r2', 'dprime', '--ld-snp', LEAD, '--ld-window-kb', '1000',
                    '--ld-window', '99999', '--ld-window-r2', '0', '--out', str(out)], capture_output=True, check=True)
    ld = pd.read_csv(f'{out}.ld', sep=r'\s+')
    ld = ld.rename(columns={'SNP_B': 'rsid', 'BP_B': 'bp_grch37', 'R2': 'r2', 'DP': 'dprime'})[['rsid', 'bp_grch37', 'r2', 'dprime']]
    run.input('g1000_eur', path=f'{L.REF}.bim', note='1000 Genomes phase 3 European PLINK reference')
    run.table('fto_ld', ld, note=f'LD (r2, D\') of 1000 Genomes European variants within 1 Mb of {LEAD} (GRCh37)')
    # gene models: Ensembl-canonical transcripts of protein-coding genes overlapping the window
    rows = []
    with gzip.open(GTF, 'rt') as fh:
        for line in fh:
            if not line.startswith(f'chr{CH}\t'):
                continue
            f = line.rstrip('\n').split('\t')
            if f[2] != 'exon':
                continue
            s, e_ = int(f[3]), int(f[4])
            if e_ < LO - 200_000 or s > HI + 200_000:
                continue
            a = f[8]
            if 'Ensembl_canonical' not in a or 'gene_type "protein_coding"' not in a:
                continue
            g = re.search(r'gene_name "([^"]+)"', a).group(1)
            t = re.search(r'transcript_id "([^"]+)"', a).group(1)
            n = int(re.search(r'exon_number (\d+)', a).group(1))
            rows.append(dict(gene=g, transcript=t, strand=f[6], exon=n, start=s, end=e_))
    gm = pd.DataFrame(rows)
    span = gm.groupby('gene').agg(g_start=('start', 'min'), g_end=('end', 'max'))
    keep = span[(span.g_end >= LO) & (span.g_start <= HI)].index
    gm = gm[gm.gene.isin(keep)].sort_values(['gene', 'start']).reset_index(drop=True)
    run.input('gencode_v46', path=str(GTF), note='GENCODE v46 basic annotation, GRCh38')
    run.table('fto_gene_models', gm, note='GENCODE v46 Ensembl-canonical protein-coding transcripts overlapping '
                                          f'chr{CH}:{LO}-{HI} (GRCh38), one row per exon')
    # recombination rate from the Beagle GRCh38 genetic map
    with zipfile.ZipFile(MAP) as z:
        mp = pd.read_csv(io.BytesIO(z.read(f'no_chr_in_chrom_field/plink.chr{CH}.GRCh38.map')), sep=r'\s+', header=None,
                         names=['chrom', 'id', 'cM', 'bp'])
    mp = mp[(mp.bp >= LO - 50_000) & (mp.bp <= HI + 50_000)].sort_values('bp').reset_index(drop=True)
    rec = pd.DataFrame({'start': mp.bp.iloc[:-1].to_numpy(), 'end': mp.bp.iloc[1:].to_numpy(),
                        'rate_cM_per_Mb': np.diff(mp.cM.to_numpy()) / np.diff(mp.bp.to_numpy()) * 1e6})
    run.input('beagle_genetic_map_GRCh38', path=str(MAP.resolve()), url='https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip',
              note='Beagle PLINK-format genetic maps (HapMap II-derived), GRCh38')
    run.table('fto_recombination', rec, note=f'Recombination rate between consecutive genetic-map markers, chr{CH}:{LO}-{HI} '
                                             '(+-50 kb), Beagle GRCh38 map')
    print(len(ld), 'LD rows; r2>=0.8:', int((ld.r2 >= 0.8).sum()))
    print(gm.groupby(['gene', 'strand']).agg(n_exon=('exon', 'size'), start=('start', 'min'), end=('end', 'max')))
    print(rec.rate_cM_per_Mb.describe().round(2).to_dict())
