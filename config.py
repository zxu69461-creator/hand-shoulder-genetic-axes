"""Locations of inputs, intermediate files and external programs.

Every path can be overridden with an environment variable of the same name (e.g. DATA_DIR=/data/hand python3 ...).
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def _path(name, default):
    return Path(os.environ.get(name, default))


# downloaded public inputs (scripts/00_fetch_inputs.sh)
DATA = _path('DATA_DIR', ROOT / 'data')
GWAS = DATA / 'gwas'                      # GWAS Catalog, ReproGen and other summary statistics
FINNGEN_RAW = DATA / 'finngen_raw'        # FinnGen R13 summary statistics as released (finngen_R13_<endpoint>.gz)
NEALE = DATA / 'neale'                    # UK Biobank sex-stratified GWAS (Neale lab round 2)
SC = DATA / 'sc'                          # single-cell atlases
CTGOV = DATA / 'ctgov'                    # ClinicalTrials.gov study records (API v2, JSON)

# reference panels and annotation
REF = _path('REF_DIR', DATA / 'ref')
LDSC_LD = _path('LDSC_LD', REF / 'eur_w_ld_chr')                        # LDSC European LD scores
HM3 = _path('HM3_SNPLIST', REF / 'w_hm3.snplist')                       # HapMap 3 SNP list
G1000 = str(_path('G1000_EUR', REF / 'g1000_eur' / 'g1000_eur'))       # 1000 Genomes European PLINK files (prefix)
GENE_LOC = _path('GENE_LOC', REF / 'NCBI37.3.gene.loc')                 # MAGMA gene locations (NCBI build 37.3)
GENCODE_GTF = _path('GENCODE_GTF', REF / 'gencode.v46.annotation.gtf.gz')
GENETIC_MAP = _path('GENETIC_MAP', REF / 'plink.GRCh38.map.zip')        # PLINK-format HapMap II genetic map, GRCh38

# intermediate files
WORK = _path('WORK_DIR', ROOT / 'work')
FINNGEN = WORK / 'finngen'                # per-endpoint Parquet (all variants with 0.001 < AF < 0.999)
MUNGED = WORK / 'munged'                  # HapMap 3 summary statistics (SNP A1 A2 Z N)
GWAMA = WORK / 'gwama'                    # cross-cohort meta-analyses per disorder
FACTOR = WORK / 'factor'                  # axis factor GWAS
MAGMA_SETS = WORK / 'magma_sets'          # cell-state gene sets
MAGMA_OUT = WORK / 'magma'                # MAGMA gene and gene-set results

# results database (created on first run)
DB = _path('RESULTS_DB', ROOT / 'db' / 'results.duckdb')

# external programs
PLINK = os.environ.get('PLINK', 'plink')          # PLINK 1.9
MAGMA = os.environ.get('MAGMA', 'magma')          # MAGMA 1.10
RSCRIPT = os.environ.get('RSCRIPT', 'Rscript')
THREADS = int(os.environ.get('THREADS', '6'))
