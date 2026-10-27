#!/usr/bin/env python3
"""Harmonise every GWAS used in the study to HapMap 3 variants (LD score regression input: SNP A1 A2 Z N).

Sources and output format (MHC removed; N = effective N 4 / (1/cases + 1/controls), or the sample size of quantitative
traits):
  - FinnGen R13 endpoints (released summary statistics); each is also written as Parquet with every variant of allele
    frequency between 0.001 and 0.999 (outcome effects for Mendelian randomisation and colocalisation);
  - UK Biobank-based soft-tissue GWAS (replication) and the endocrine and metabolic exposure GWAS (GWAS Catalog
    harmonised files, the BOLT-LMM carpal tunnel GWAS and the ReproGen menopause meta-analysis).
Table: munged_inputs.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'lib')]
import duckdb
import pandas as pd
import config as C
import mrlib
from project import Project

P = Project(__file__)
HM3 = C.HM3
GW = C.GWAS
FINNGEN = {   # endpoint: (cases, controls); the musculoskeletal (M13) endpoints share one control group
    'G6_CARPTU': (32951, 431627), 'M13_TRIGGERFINGER': (6528, 350578), 'M13_DEQUERVAIN': (2381, 350578),
    'M13_DUPUTRYEN': (6972, 350578), 'M13_ADHCAPSULITIS': (7586, 350578), 'M13_PLANTARFASC': (4210, 350578),
    'M13_ENTESOPATHYOTH': (15115, 350578), 'M13_ROTATORCUFF': (34158, 350578)}

# label, source kind, file, accession, cases (or N), controls, note
STUDIES = [(ep, 'finngen', C.FINNGEN_RAW / f'finngen_R13_{ep}.gz', f'FinnGen R13 {ep}', ca, co, f'FinnGen R13 {ep}')
           for ep, (ca, co) in FINNGEN.items()] + [
    ('UKB_TF', 'harmonised', GW / 'GCST90104907.h.tsv.gz', 'GCST90104907', 2908, 436579, 'Patel 2022 trigger finger, UK Biobank'),
    ('DD_META', 'harmonised', GW / 'GCST90301252.h.tsv.gz', 'GCST90301252', 11320, 47023,
     "Riesmeijer 2024 Dupuytren's disease meta-analysis (UK, Netherlands, Germany; no FinnGen)"),
    ('UKB_CTS', 'bolt', GW / 'WibergA_2019_UKBB.txt', 'GCST007581', 12312, 389344, 'Wiberg 2019 carpal tunnel syndrome, UK Biobank'),
    ('UKB_FS', 'harmonised', GW / '34111113-GCST90000512-EFO_1000941.h.tsv.gz', 'GCST90000512', 10104, 440995,
     'Green 2021 frozen shoulder, UK Biobank only'),
    ('UKB_TENOSYN', 'harmonised', GW / '34737426-GCST90044562-EFO_1001435.h.tsv.gz', 'GCST90044562', 3692, 452656,
     'Jiang 2021 PheCode 727.1 synovitis and tenosynovitis, UK Biobank'),
    ('UKB_DUP', 'harmonised', GW / '34737426-GCST90044569-EFO_1000438.h.tsv.gz', 'GCST90044569', 3728, 452620,
     "Jiang 2021 PheCode 728.6 Dupuytren's disease, UK Biobank"),
    ('UKB_ENTHESO', 'harmonised', GW / '34737426-GCST90044558-EFO_0009666.h.tsv.gz', 'GCST90044558', 6439, 449909,
     'Jiang 2021 PheCode 726.1 enthesopathy, UK Biobank'),
    ('SHBG_F', 'harmonised', GW / '32042192-GCST90012107-EFO_0004696.h.tsv.gz', 'GCST90012107', 189473, None,
     'Ruth 2020 SHBG, UK Biobank women'),
    ('BIOT_F', 'harmonised', GW / '32042192-GCST90012102-EFO_0004908.h.tsv.gz', 'GCST90012102', 188507, None,
     'Ruth 2020 bioavailable testosterone, UK Biobank women'),
    ('TT_F', 'harmonised', GW / '32042192-GCST90012112-EFO_0004908.h.tsv.gz', 'GCST90012112', 230454, None,
     'Ruth 2020 total testosterone, UK Biobank women'),
    ('IGF1', 'harmonised', GW / 'GCST90019511.h.tsv.gz', 'GCST90019511', 340567, None,
     'Sinnott-Armstrong 2021 IGF-1, UK Biobank European ancestry'),
    ('HBA1C', 'harmonised', GW / '34059833-GCST90002244-EFO_0004541.h.tsv.gz', 'GCST90002244', 146806, None,
     'Chen 2021 MAGIC HbA1c, European ancestry (no UK Biobank)'),
    ('ANM', 'reprogen', GW / 'reprogen_ANM_201K_170621.txt.gz', 'ReproGen 2021', 201323, None,
     'ReproGen 2021 age at natural menopause (201,323 women; per-SNP N)'),
    ('T2D_XUE', 'harmonised', GW / '30054458-GCST006867-EFO_0001360.h.tsv.gz', 'GCST006867', 61714, 593952,
     'Xue 2018 type 2 diabetes (DIAGRAM, GERA, UK Biobank; no FinnGen)'),
    ('BMI_GIANT', 'harmonised', GW / '25673413-GCST002783-EFO_0004340.h.tsv.gz', 'GCST002783', 322154, None,
     'Locke 2015 GIANT body mass index (no UK Biobank)'),
]


def num(x):
    return f"try_cast(nullif(nullif({x},'NA'),'') as double)"


def _con():
    con = duckdb.connect()
    con.execute(f"SET threads={C.THREADS}; SET memory_limit='8GB'")
    return con


def munge_finngen(src, neff, out, slim):
    con = _con()
    con.execute(f"""CREATE TEMP TABLE g0 AS
        SELECT CAST(replace(CAST("#chrom" AS VARCHAR),'X','23') AS INTEGER) AS chrom, pos, ref, alt, rsids AS rsid,
               nearest_genes, beta, sebeta AS se, pval AS p, mlogp, af_alt AS af, af_alt_cases AS af_cases,
               af_alt_controls AS af_controls
        FROM read_csv('{src}', delim='\t', header=true, compression='gzip', nullstr=['NA',''], sample_size=200000)
        WHERE af_alt > 0.001 AND af_alt < 0.999""")
    con.execute(f"COPY g0 TO '{slim}' (FORMAT PARQUET, COMPRESSION ZSTD)")
    con.execute("""CREATE TEMP TABLE g AS SELECT chrom, pos, ref, alt, split_part(rsid, ',', 1) AS rsid, beta, se, af
                   FROM g0 WHERE se > 0""")
    nvar = con.execute("SELECT count(*) FROM g").fetchone()[0]
    con.execute(f"""CREATE TEMP TABLE m AS
        WITH h AS (SELECT * FROM read_csv('{HM3}', delim='\t', header=true))
        SELECT h.SNP, g.alt AS A1, g.ref AS A2, g.beta / g.se AS Z, {neff} AS N
        FROM g JOIN h ON g.rsid = h.SNP
        WHERE g.rsid LIKE 'rs%' AND g.af BETWEEN 0.01 AND 0.99
          AND NOT (g.chrom = 6 AND g.pos BETWEEN 25000000 AND 34000000)
          AND ((g.alt = h.A1 AND g.ref = h.A2) OR (g.alt = h.A2 AND g.ref = h.A1))
        QUALIFY row_number() OVER (PARTITION BY h.SNP ORDER BY abs(g.beta / g.se) DESC) = 1""")
    return _finish(con, nvar, out)


def munge_harmonised(src, neff, out):
    """GWAS Catalog harmonised file (hm_* columns or GWAS-SSF); SE from the 95% CI of the OR when not reported."""
    con = _con()
    c = {r[0] for r in con.execute(f"describe select * from read_csv('{src}', delim='\t', header=true, "
                                   f"all_varchar=true, sample_size=2000)").fetchall()}
    if 'hm_rsid' in c:
        rs, ch, pos, ea, oa = 'hm_rsid', 'hm_chrom', 'hm_pos', 'hm_effect_allele', 'hm_other_allele'
        beta = f"coalesce({num('hm_beta')}, ln({num('hm_odds_ratio')}))"
        eaf = num('hm_effect_allele_frequency')
    else:
        rs, ch, pos, ea, oa = 'rsid', 'chromosome', 'base_pair_location', 'effect_allele', 'other_allele'
        parts = [num('beta')] * ('beta' in c) + [f"ln({num('odds_ratio')})"] * ('odds_ratio' in c)
        beta = f"coalesce({', '.join(parts)})" if len(parts) > 1 else parts[0]
        eaf = num('effect_allele_frequency') if 'effect_allele_frequency' in c else 'null::double'
    se_ci = (f"(ln({num('ci_upper')}) - ln({num('ci_lower')})) / 3.919928" if {'ci_upper', 'ci_lower'} <= c
             else 'null::double')
    se = f"coalesce({num('standard_error')}, {se_ci})" if 'standard_error' in c else se_ci
    con.execute(f"""CREATE TEMP TABLE g AS
        SELECT {rs} AS rsid, try_cast(replace(cast({ch} AS varchar), 'X', '23') AS integer) AS chrom,
               try_cast({pos} AS bigint) AS pos, upper({ea}) AS ea, upper({oa}) AS oa, {beta} AS beta, {se} AS se,
               {eaf} AS eaf, {num('p_value')} AS p
        FROM read_csv('{src}', delim='\t', header=true, all_varchar=true)""")
    nvar = con.execute("SELECT count(*) FROM g WHERE beta IS NOT NULL AND se > 0").fetchone()[0]
    con.execute(f"""CREATE TEMP TABLE m AS
        WITH h AS (SELECT * FROM read_csv('{HM3}', delim='\t', header=true))
        SELECT h.SNP, g.ea AS A1, g.oa AS A2, g.beta / g.se AS Z, {neff} AS N
        FROM g JOIN h ON g.rsid = h.SNP
        WHERE g.beta IS NOT NULL AND g.se > 0 AND isfinite(g.beta / g.se)
          AND (g.eaf IS NULL OR g.eaf BETWEEN 0.01 AND 0.99)
          AND NOT (g.chrom = 6 AND g.pos BETWEEN 25000000 AND 34000000)
          AND ((g.ea = h.A1 AND g.oa = h.A2) OR (g.ea = h.A2 AND g.oa = h.A1))
        QUALIFY row_number() OVER (PARTITION BY h.SNP ORDER BY abs(g.beta / g.se) DESC) = 1""")
    return _finish(con, nvar, out)


def munge_bolt(src, neff, out):
    """BOLT-LMM output (SNP CHR BP ... ALLELE1 ALLELE0 A1FREQ INFO BETA SE P_BOLT_LMM...); GRCh37 positions."""
    con = _con()
    cols = {r[0] for r in con.execute(f"DESCRIBE SELECT * FROM read_csv('{src}', header=true, sample_size=5000)").fetchall()}
    need = {'SNP', 'CHR', 'BP', 'ALLELE1', 'ALLELE0', 'A1FREQ', 'BETA', 'SE'}
    assert need <= cols, f'unexpected BOLT columns: {sorted(cols)}'
    info = 'AND INFO >= 0.8' if 'INFO' in cols else ''
    con.execute(f"""CREATE TEMP TABLE g AS
        SELECT SNP AS rsid, try_cast(CHR AS integer) AS chrom, BP AS pos, upper(ALLELE1) AS ea, upper(ALLELE0) AS oa,
               BETA AS beta, SE AS se, A1FREQ AS af
        FROM read_csv('{src}', header=true, sample_size=5000) WHERE SE > 0 {info}""")
    nvar = con.execute("SELECT count(*) FROM g").fetchone()[0]
    con.execute(f"""CREATE TEMP TABLE m AS
        WITH h AS (SELECT * FROM read_csv('{HM3}', delim='\t', header=true))
        SELECT h.SNP, g.ea AS A1, g.oa AS A2, g.beta / g.se AS Z, {neff} AS N
        FROM g JOIN h ON g.rsid = h.SNP
        WHERE g.af BETWEEN 0.01 AND 0.99 AND NOT (g.chrom = 6 AND g.pos BETWEEN 26000000 AND 34000000)
          AND ((g.ea = h.A1 AND g.oa = h.A2) OR (g.ea = h.A2 AND g.oa = h.A1))
        QUALIFY row_number() OVER (PARTITION BY h.SNP ORDER BY abs(g.beta / g.se) DESC) = 1""")
    return _finish(con, nvar, out)


def munge_reprogen(src, out):
    """ReproGen text (chr:pos identifiers, GRCh37) -> rsIDs via the 1000 Genomes European .bim; per-SNP N."""
    g = mrlib.read_reprogen(src)
    con = _con()
    con.register('g0', g)
    con.execute("CREATE TEMP TABLE g AS SELECT * FROM g0 WHERE se > 0")
    nvar = con.execute("SELECT count(*) FROM g").fetchone()[0]
    con.execute(f"""CREATE TEMP TABLE m AS
        WITH h AS (SELECT * FROM read_csv('{HM3}', delim='\t', header=true))
        SELECT h.SNP, g.ea AS A1, g.oa AS A2, g.beta / g.se AS Z, g.n AS N
        FROM g JOIN h ON g.rsid = h.SNP
        WHERE g.eaf BETWEEN 0.01 AND 0.99 AND NOT (g.chrom = 6 AND g.pos BETWEEN 26000000 AND 34000000)
          AND ((g.ea = h.A1 AND g.oa = h.A2) OR (g.ea = h.A2 AND g.oa = h.A1))""")
    return _finish(con, nvar, out)


def _finish(con, nvar, out):
    nm = con.execute("SELECT count(*) FROM m").fetchone()[0]
    lam = con.execute("SELECT median(Z*Z) / 0.4549 FROM m").fetchone()[0]
    if nvar < 300_000 or nm < 200_000:
        raise RuntimeError(f'implausible counts variants={nvar} hm3={nm}; nothing written')
    con.execute(f"COPY m TO '{out}' (FORMAT CSV, DELIMITER '\t', HEADER, COMPRESSION GZIP)")
    con.close()
    return nvar, nm, lam


C.MUNGED.mkdir(parents=True, exist_ok=True); C.FINNGEN.mkdir(parents=True, exist_ok=True)
rows = []
with P.run(__file__, params={'hm3': str(HM3), 'n': 'effective N = 4/(1/cases+1/controls); quantitative: N'}) as run:
    for label, kind, src, acc, nca, nco, note in STUDIES:
        src = Path(src)
        if not src.exists():
            print(f'[skip] {label}: {src.name} not downloaded'); continue
        neff = 4.0 / (1.0 / nca + 1.0 / nco) if nco else float(nca)
        out = C.MUNGED / f'{label}.sumstats.gz'
        try:
            if kind == 'finngen':
                r = munge_finngen(src, neff, out, C.FINNGEN / f'{label}.parquet')
            elif kind == 'bolt':
                r = munge_bolt(src, neff, out)
            elif kind == 'reprogen':
                r = munge_reprogen(src, out)
            else:
                r = munge_harmonised(src, neff, out)
        except RuntimeError as e:
            print(f'[fail] {label}: {e}', flush=True); continue
        run.input(acc, path=str(src), note=note)
        rows.append(dict(label=label, accession=acc, source=kind, cases=nca, controls=nco, neff=round(neff, 1),
                         n_variants=int(r[0]), n_hm3=int(r[1]), lambda_gc=round(float(r[2]), 3), path=str(out), note=note))
        print(f'[done] {label} hm3={r[1]} lambda={r[2]:.3f}', flush=True)
    if not rows:
        sys.exit('nothing munged: no input present')
    run.table('munged_inputs', pd.DataFrame(rows), note='GWAS harmonised to HapMap 3 variants')
