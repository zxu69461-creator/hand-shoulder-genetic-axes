"""Shared Mendelian randomisation helpers: exposure and outcome registry, GWAS readers, clumping, harmonisation.

Readers return one row per variant with columns rsid, chrom, pos, ea, oa, beta, se, eaf, p (missing ones as null);
pos is GRCh38 for GWAS Catalog harmonised files and FinnGen, GRCh37 for ReproGen (used only for the MHC filter).
GWAS Catalog harmonised files come in two layouts (hm_* columns, or GWAS-SSF); both are handled.
"""
import subprocess
import sys
from pathlib import Path

import duckdb
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C  # noqa: E402

GW = C.GWAS
FGS = C.FINNGEN
PLINK = C.PLINK
REF = C.G1000
P_INST, R2, KB = 5e-8, 0.001, 10000
PALIN = {('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C')}

# label, source kind, file, accession, note
EXPOSURES = [
    ('SHBG_F', 'gwas', GW / '32042192-GCST90012107-EFO_0004696.h.tsv.gz', 'GCST90012107', 'SHBG, UK Biobank women (Ruth 2020)'),
    ('BIOT_F', 'gwas', GW / '32042192-GCST90012102-EFO_0004908.h.tsv.gz', 'GCST90012102', 'bioavailable testosterone, UK Biobank women (Ruth 2020)'),
    ('TT_F', 'gwas', GW / '32042192-GCST90012112-EFO_0004908.h.tsv.gz', 'GCST90012112', 'total testosterone, UK Biobank women (Ruth 2020)'),
    ('IGF1', 'gwas', GW / 'GCST90019511.h.tsv.gz', 'GCST90019511', 'IGF-1, UK Biobank (Sinnott-Armstrong 2021)'),
    ('HBA1C', 'gwas', GW / '34059833-GCST90002244-EFO_0004541.h.tsv.gz', 'GCST90002244', 'HbA1c, MAGIC European ancestry (Chen 2021)'),
    ('ANM', 'reprogen', GW / 'reprogen_ANM_201K_170621.txt.gz', 'ReproGen 2021', 'age at natural menopause, ReproGen 2021 (201,323 women)'),
    ('T2D_XUE', 'gwas', GW / '30054458-GCST006867-EFO_0001360.h.tsv.gz', 'GCST006867', 'type 2 diabetes (Xue 2018)'),
    ('BMI_GIANT', 'gwas', GW / '25673413-GCST002783-EFO_0004340.h.tsv.gz', 'GCST002783', 'body mass index, GIANT (Locke 2015)'),
]
EXPO = {e[0]: e for e in EXPOSURES}

# label, cohort, axis, source kind, path, linear-scale case fraction (None = already log-odds)
OUTCOMES = [
    ('CTS', 'FinnGen', 'tenosynovial', 'fg', FGS / 'G6_CARPTU.parquet', None),
    ('TF', 'FinnGen', 'tenosynovial', 'fg', FGS / 'M13_TRIGGERFINGER.parquet', None),
    ('DQ', 'FinnGen', 'tenosynovial', 'fg', FGS / 'M13_DEQUERVAIN.parquet', None),
    ('DUP', 'FinnGen', 'fibro-capsular', 'fg', FGS / 'M13_DUPUTRYEN.parquet', None),
    ('PF', 'FinnGen', 'fibro-capsular', 'fg', FGS / 'M13_PLANTARFASC.parquet', None),
    ('FS', 'FinnGen', 'fibro-capsular', 'fg', FGS / 'M13_ADHCAPSULITIS.parquet', None),
    ('RC', 'FinnGen', 'fibro-capsular', 'fg', FGS / 'M13_ROTATORCUFF.parquet', None),
    ('ENTH', 'FinnGen', 'fibro-capsular', 'fg', FGS / 'M13_ENTESOPATHYOTH.parquet', None),
    ('CTS', 'UKB', 'tenosynovial', 'bolt', GW / 'WibergA_2019_UKBB.txt', 12312 / (12312 + 389344)),
    ('TF', 'UKB', 'tenosynovial', 'gwas', GW / 'GCST90104907.h.tsv.gz', None),
    ('TENOSYN', 'UKB', 'tenosynovial', 'gwas', GW / '34737426-GCST90044562-EFO_1001435.h.tsv.gz', None),
    ('DUP', 'UKB', 'fibro-capsular', 'gwas', GW / '34737426-GCST90044569-EFO_1000438.h.tsv.gz', None),
    ('DUP_META', 'non-Finnish meta', 'fibro-capsular', 'gwas', GW / 'GCST90301252.h.tsv.gz', None),
    ('FS', 'UKB', 'fibro-capsular', 'gwas', GW / '34111113-GCST90000512-EFO_1000941.h.tsv.gz', None),
    ('ENTH', 'UKB', 'fibro-capsular', 'gwas', GW / '34737426-GCST90044558-EFO_0009666.h.tsv.gz', None),
]
# outcome cohorts each exposure may be paired with (no overlap with the exposure sample); the Benjamini-Hochberg FDR
# of 33_mr_summary.py is computed over every exposure paired with the eight FinnGen outcomes (TO_FINNGEN x 8)
TO_FINNGEN = {'SHBG_F', 'BIOT_F', 'TT_F', 'IGF1', 'HBA1C', 'ANM', 'T2D_XUE', 'BMI_GIANT'}
TO_UKB = {'HBA1C', 'BMI_GIANT'}


def num(x):
    return f"try_cast(nullif(nullif(cast({x} as varchar),'NA'),'') as double)"


def _con():
    con = duckdb.connect()
    con.execute(f"SET threads={C.THREADS}; SET memory_limit='8GB'; SET enable_progress_bar=false")
    return con


def read_gwas(src, where='true', rsids=None):
    """GWAS Catalog harmonised file -> standard columns; `where` filters on the standard columns (e.g. 'p < 5e-8')."""
    con = _con()
    c = {r[0] for r in con.execute(f"DESCRIBE SELECT * FROM read_csv('{src}', delim='\t', header=true, all_varchar=true, "
                                   f"sample_size=2000)").fetchall()}
    if 'hm_rsid' in c:
        rs, ch, ea, oa = 'hm_rsid', 'hm_chrom', 'hm_effect_allele', 'hm_other_allele'
        beta = f"coalesce({num('hm_beta')}, ln({num('hm_odds_ratio')}))"
        eaf = num('hm_effect_allele_frequency')
    else:
        rs, ch, ea, oa = 'rsid', 'chromosome', 'effect_allele', 'other_allele'
        parts = [num('beta')] * ('beta' in c) + [f"ln({num('odds_ratio')})"] * ('odds_ratio' in c)
        beta = f"coalesce({', '.join(parts)})" if len(parts) > 1 else parts[0]
        eaf = num('effect_allele_frequency') if 'effect_allele_frequency' in c else 'null::double'
    pos = 'hm_pos' if 'hm_pos' in c else 'base_pair_location'
    join = ''
    if rsids is not None:
        con.register('want', pd.DataFrame({'w': sorted(rsids)})); join = f'JOIN want ON o.{rs} = want.w'
    df = con.execute(f"""
        SELECT * FROM (
          SELECT o.{rs} AS rsid, try_cast(o.{ch} AS integer) AS chrom, try_cast(o.{pos} AS bigint) AS pos,
                 upper(o.{ea}) AS ea, upper(o.{oa}) AS oa, {beta} AS beta, {num('standard_error')} AS se, {eaf} AS eaf,
                 {num('p_value')} AS p
          FROM read_csv('{src}', delim='\t', header=true, all_varchar=true) o {join}) WHERE {where}""").df()
    con.close()
    return df


def read_reprogen(src, where='true', rsids=None):
    """ReproGen 2021 age at natural menopause (SNP CHR POS Effect_Allele Other_Allele EAF Effect SE Pval N; GRCh37;
    effect in years per allele). Variant IDs are chr:pos, so rsIDs come from the 1000 Genomes European reference .bim
    (GRCh37) matched on chromosome, position and the unordered allele pair."""
    con = _con()
    con.execute(f"""CREATE TEMP VIEW bim AS SELECT column0::integer AS chr, column1 AS rsid, column3::bigint AS bp,
                           upper(column4) AS a1, upper(column5) AS a2
                    FROM read_csv('{REF}.bim', delim='\t', header=false, columns={{'column0':'VARCHAR','column1':'VARCHAR',
                         'column2':'VARCHAR','column3':'VARCHAR','column4':'VARCHAR','column5':'VARCHAR'}})
                    WHERE column0 SIMILAR TO '[0-9]+'""")
    want = ''
    if rsids is not None:
        con.register('want', pd.DataFrame({'w': sorted(rsids)})); want = 'AND b.rsid IN (SELECT w FROM want)'
    df = con.execute(f"""
        SELECT * FROM (
          SELECT b.rsid, o.CHR::integer AS chrom, o.POS::bigint AS pos, upper(o.Effect_Allele) AS ea,
                 upper(o.Other_Allele) AS oa, o.Effect::double AS beta, o.SE::double AS se, o.EAF::double AS eaf,
                 o.Pval::double AS p, o.N::double AS n
          FROM read_csv('{src}', delim='\t', header=true) o
          JOIN bim b ON b.chr = o.CHR AND b.bp = o.POS
           AND ((b.a1 = upper(o.Effect_Allele) AND b.a2 = upper(o.Other_Allele))
             OR (b.a2 = upper(o.Effect_Allele) AND b.a1 = upper(o.Other_Allele)))
          WHERE b.rsid LIKE 'rs%' {want}
          QUALIFY row_number() OVER (PARTITION BY b.rsid ORDER BY o.N DESC) = 1) WHERE {where}""").df()
    con.close()
    return df


def exposure_gws(label):
    _, kind, src, _, _ = EXPO[label]
    return (read_gwas if kind == 'gwas' else read_reprogen)(Path(src).resolve(), f'p < {P_INST}')


def exposure_lookup(label, rsids):
    """Effects of the given rsIDs in an exposure GWAS, any P."""
    _, kind, src, _, _ = EXPO[label]
    return (read_gwas if kind == 'gwas' else read_reprogen)(Path(src).resolve(), rsids=rsids)


def available(label):
    return Path(EXPO[label][2]).exists()


def clump(df, tmp, label):
    """Filter (rsID, se > 0, MAF >= 0.01, no MHC) and PLINK-clump on P against 1000 Genomes Europeans; adds per-SNP F."""
    df = df[df.rsid.str.startswith('rs', na=False) & (df.se > 0) & df.beta.notna()]
    df = df[(df.eaf.isna()) | ((df.eaf >= 0.01) & (df.eaf <= 0.99))]
    df = df[~((df.chrom == 6) & df.pos.between(25_000_000, 34_000_000))]
    df = df.sort_values('p').drop_duplicates('rsid')
    f = Path(tmp) / f'{label}.assoc'
    df[['rsid', 'p']].rename(columns={'rsid': 'SNP', 'p': 'P'}).to_csv(f, sep='\t', index=False)
    r = subprocess.run([PLINK, '--bfile', REF, '--clump', str(f), '--clump-p1', str(P_INST), '--clump-p2', '1',
                        '--clump-r2', str(R2), '--clump-kb', str(KB), '--out', str(Path(tmp) / label)],
                       capture_output=True, text=True)
    cl = Path(tmp) / f'{label}.clumped'
    if r.returncode or not cl.exists():
        raise RuntimeError(f'{label}: plink clump failed: {r.stdout[-400:]} {r.stderr[-400:]}')
    keep = pd.read_csv(cl, sep=r'\s+').SNP
    out = df[df.rsid.isin(keep)].copy()
    out['F'] = (out.beta / out.se) ** 2
    return out, len(df)


def outcome_rows(kind, path, rsids, mu):
    """Outcome effects (log-odds) at the given rsIDs."""
    if kind == 'gwas':
        df = read_gwas(path, rsids=rsids)
    else:
        con = _con()
        con.register('want', pd.DataFrame({'rsid': sorted(rsids)}))
        if kind == 'fg':
            df = con.execute(f"""SELECT w.rsid, upper(alt) AS ea, upper(ref) AS oa, beta, se, af AS eaf
                FROM (SELECT *, unnest(string_split(rsid, ',')) AS rs1 FROM read_parquet('{path}')) o
                JOIN want w ON o.rs1 = w.rsid""").df()
        else:   # BOLT-LMM text
            df = con.execute(f"""SELECT SNP AS rsid, upper(ALLELE1) AS ea, upper(ALLELE0) AS oa, BETA AS beta, SE AS se,
                                        A1FREQ AS eaf
                FROM read_csv('{path}', header=true, sample_size=5000) o JOIN want w ON o.SNP = w.rsid""").df()
        con.close()
    df = df[['rsid', 'ea', 'oa', 'beta', 'se', 'eaf']]
    if mu is not None:
        df = df.assign(beta=df.beta / (mu * (1 - mu)), se=df.se / (mu * (1 - mu)))
    return df.dropna(subset=['beta', 'se']).query('se > 0').drop_duplicates('rsid')


def harmonise(x, y):
    """Align y (outcome or second exposure) to x's effect allele; drop allele mismatches and palindromic SNPs."""
    m = x.merge(y, on='rsid', suffixes=('_x', '_y'))
    same = (m.ea_x == m.ea_y) & (m.oa_x == m.oa_y)
    flip = (m.ea_x == m.oa_y) & (m.oa_x == m.ea_y)
    m = m[same | flip].copy()
    m.loc[flip[same | flip].values, 'beta_y'] *= -1
    return m[[(a, b) not in PALIN for a, b in zip(m.ea_x, m.oa_x)]]
