# Genetic axes of hand and shoulder soft-tissue disease

Analysis code and source data for:

> Xu Z, Cheng J, Chen X, Luo J, Zhu X, Guo X. Genetic axes of hand and shoulder soft-tissue disease and their divergent
> responses to adiposity and incretin therapy: a cross-biobank genomic and randomised trial study. *medRxiv* 2026,
> doi:[10.64898/2026.10.27.26366575](https://doi.org/10.64898/2026.10.27.26366575) (preprint: not peer reviewed).

Carpal tunnel syndrome, trigger finger, de Quervain tenosynovitis, Dupuytren's disease, frozen shoulder and plantar
fascial disorder cluster in people with obesity and diabetes. Using genome-wide data from FinnGen release 13 and UK
Biobank-based GWAS, the study shows that these disorders form two heritable axes, a tenosynovial axis and a
fibro-capsular axis, each with its own genome-wide loci (genomic structural equation modelling and common-factor GWAS).
Mendelian randomisation tests five endocrine and metabolic exposures on each disorder and on each axis, gene-level and
cell-state enrichment analyses map the axes onto three single-cell atlases of human tendon and palmar fascia, and
serious adverse events from 210 placebo-controlled trials of incretin receptor agonists test how pharmacological weight
loss affects each axis. All inputs are public.

## Source data

`results/Source_Data.xlsx` holds the data behind every figure and the result tables quoted in the text; its first sheet
lists the contents. The same sheets are in `results/source_data/` as CSV files.

## Public inputs

| Input | Source | Used in |
|---|---|---|
| FinnGen R13 endpoints G6_CARPTU, M13_TRIGGERFINGER, M13_DEQUERVAIN, M13_DUPUTRYEN, M13_ADHCAPSULITIS, M13_PLANTARFASC, M13_ENTESOPATHYOTH, M13_ROTATORCUFF | [FinnGen R13 summary statistics](https://www.finngen.fi/en/access_results) | all genetic analyses |
| UK Biobank-based GWAS: carpal tunnel syndrome (GCST007581), trigger finger (GCST90104907), frozen shoulder (GCST90000512), Dupuytren's disease (GCST90044569), enthesopathy (GCST90044558), synovitis and tenosynovitis (GCST90044562); Dupuytren's disease meta-analysis without FinnGen (GCST90301252) | [GWAS Catalog](https://www.ebi.ac.uk/gwas/) | replication, meta-analysis |
| Exposures: body mass index (GCST002783), IGF-1 (GCST90019511), SHBG, bioavailable and total testosterone in women (GCST90012107, GCST90012102, GCST90012112) and men (GCST90012109, GCST90012103, GCST90012113), type 2 diabetes (GCST006867), HbA1c (GCST90002244), age at natural menopause (ReproGen 2021) | GWAS Catalog; [ReproGen](https://www.reprogen.org) | Mendelian randomisation, colocalisation, genetic correlation |
| UK Biobank GWAS by sex, hospital main diagnosis ICD-10 G56 and M72 | [Neale lab round 2](http://www.nealelab.is/uk-biobank) | sex-stratified Mendelian randomisation |
| Palmar fascia atlas (Dupuytren's disease, non-pathogenic palmar fascia, dermis) | GEO [GSE173252](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173252) | cell-state enrichment, Fig. 1D |
| Achilles tendon and hamstring tendon atlases | [CELLxGENE Discover](https://cellxgene.cziscience.com) | cell-state enrichment |
| Trial results | [ClinicalTrials.gov](https://clinicaltrials.gov) (API v2) | Fig. 5 |
| GENCODE release 46 annotation; HapMap phase II genetic map (GRCh38, PLINK format) | [GENCODE](https://www.gencodegenes.org); [Beagle genetic maps](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/) | Fig. 4B tracks |

`scripts/00_fetch_inputs.sh` downloads all of these into `data/`.

### Reference panels

Place these under `data/ref/` (or point the environment variables in `config.py` to existing copies):

- `eur_w_ld_chr/` and `w_hm3.snplist`: European LD scores and the HapMap 3 SNP list distributed with
  [LDSC](https://github.com/bulik/ldsc);
- `g1000_eur/g1000_eur.{bed,bim,fam}` and `NCBI37.3.gene.loc`: 1000 Genomes European reference and gene locations
  from the [MAGMA website](https://cncr.nl/research/magma/).

## Software

| Analysis | Software |
|---|---|
| Genomic SEM, factor GWAS | R 4.6.1; GenomicSEM 0.0.5, lavaan 0.7.3, data.table 1.18.6.1, Matrix 1.7.5 |
| Mendelian randomisation, colocalisation, palmar fascia atlas | R 4.5.2; MendelianRandomization 0.10.0, coloc 5.2.3, Seurat 5.6.0, SeuratObject 5.4.0, data.table 1.18.6.1, Matrix 1.7.6 |
| Data handling, LD score regression, meta-analysis, trials | Python 3.14; duckdb 1.5.6, pandas 3.0.6, numpy 2.4.3, scipy 1.17.1 |
| Tendon atlases | Python 3.11; anndata 0.12.19, scanpy 1.11.5 |
| Clumping, LD | PLINK v1.9.0-b.7.7 |
| Gene and gene-set analysis | MAGMA v1.10 |

`environment.yml` and `r_packages.txt` list the packages. Paths, programs and thread counts are set in `config.py` and
can be overridden with environment variables (`DATA_DIR`, `REF_DIR`, `WORK_DIR`, `RESULTS_DB`, `PLINK`, `MAGMA`,
`RSCRIPT`, `THREADS`).

## Repository layout

```
config.py    paths to inputs, intermediate files and programs
lib/         project.py (results database), mrlib.py (Mendelian randomisation helpers),
             fastldsc.py (multivariate LD score regression with block-jackknife standard errors)
scripts/     analysis scripts, numbered by analysis layer
results/     Source_Data.xlsx and the same sheets as CSV
```

Scripts write their result tables to one DuckDB database (`db/results.duckdb`, created on the first run) and their
intermediate files to `work/`. Each script states its inputs, method and output tables in its header.

## Run order

```
bash scripts/00_fetch_inputs.sh                 # public inputs
python3 scripts/10_munge.py                     # HapMap 3 summary statistics; FinnGen endpoints as Parquet
python3 scripts/21_gsem.py                      # genetic correlations, factor models (runs 20_gsem_structure.R)
python3 scripts/22_disease_level.py             # disease-level genetic correlations and factor analysis
python3 scripts/23_efa_detail.py
python3 scripts/50_gwama.py                     # cross-cohort meta-analysis per disorder
python3 scripts/24_factor_stage.py              # factor GWAS inputs
bash scripts/25_factor_gwas.sh                  # common-factor GWAS of both axes (25_factor_gwas.R)
python3 scripts/26_factor_downstream.py         # factor GWAS summary, Mendelian randomisation on the axes
python3 scripts/29_factor_loci.py               # lead variants and loci of the axes
python3 scripts/30_factor_rg.py                 # genetic correlations of the axes
python3 scripts/30_mr_instruments.py            # Mendelian randomisation
python3 scripts/31_mr.py
python3 scripts/33_mr_summary.py
python3 scripts/32_mvmr.py
python3 scripts/36_mr_independent.py            # cross-cohort contrast, instrument strength, heritability
python3 scripts/38_sex_mr.py
python3 scripts/35_coloc.py
bash scripts/60_sc_sets.sh                      # cell-state gene sets, palmar fascia UMAP
python3 scripts/51_magma.py                     # gene and cell-state analyses (disorders and axes)
python3 scripts/52_magma_ingest.py
python3 scripts/53_cell_contrast.py
python3 scripts/39_ctgov.py                     # trial records (default: the analysed trials; --search to re-select)
python3 scripts/40_trial_axis.py
python3 scripts/42_trial_axis_contrast.py
python3 scripts/43_trial_robustness.py
python3 scripts/65_fto_locus_tracks.py          # LD, gene and recombination tracks at FTO
```

## Where each result comes from

| Figure or result | Result tables | Scripts |
|---|---|---|
| Fig. 1C (cases per disorder) | `munged_inputs` | `10_munge.py` |
| Fig. 1D (palmar fascia atlas, atlas sizes) | `dd_umap`, `sc_specificity` | `60_sc_sets.sh`, `62_dd_means.R`, `63_dd_umap.R`, `64_dd_umap_ingest.py` |
| Fig. 1E (trials by molecule) | `trial_incretin_molecules` | `43_trial_robustness.py` |
| Fig. 2A, 2B (disease-level structure) | `disease_rg`, `disease_efa`, `disease_eigen`, `disease_efa_detail` | `21_gsem.py`, `22_disease_level.py`, `23_efa_detail.py` |
| Fig. 2C, 2D (cross-cohort genetic correlations); model fit | `gsem_rg`, `gsem_cfa_fit`, `gsem_cfa_est` | `21_gsem.py` |
| Fig. 2E (axis factor GWAS) | `factor_gwas`, `factor_gwas_summary`, `factor_leads`, `factor_loci` | `24`, `25`, `26`, `29` |
| Heritability; genetic correlations of the axes | `structure_h2`, `factor_rg` | `36_mr_independent.py`, `30_factor_rg.py` |
| Fig. 3A, 3B (Mendelian randomisation) | `mr_instruments`, `mr_results`, `mr_summary` | `30_mr_instruments.py`, `31_mr.py`, `33_mr_summary.py` |
| Fig. 3C (cross-cohort contrast); instrument strength | `mr_bmi_independent`, `mr_instrument_strength` | `36_mr_independent.py` |
| Fig. 3D (effects on the axes) | `mr_factor_results`, `mr_factor_contrast` | `26_factor_downstream.py` |
| Fig. 3E (by sex; hormones measured in men) | `mr_sex_results` | `38_sex_mr.py` |
| Fig. 3F (multivariable) | `mvmr_results` | `32_mvmr.py` |
| Fig. 4A (gene-level association) | `magma_genes`, `factor_magma_genes` | `50_gwama.py`, `51_magma.py`, `52_magma_ingest.py` |
| Fig. 4B (FTO) | `coloc_results`, `fto_ld`, `fto_gene_models`, `fto_recombination` | `35_coloc.py`, `65_fto_locus_tracks.py` |
| Fig. 4C, 4D (cell states) | `magma_gsa`, `factor_magma_gsa`, `cell_contrast` | `60_sc_sets.sh`, `51_magma.py`, `52_magma_ingest.py`, `53_cell_contrast.py` |
| Fig. 5 (trials) | `trial_robustness`, `trial_axis_contrast`, `trial_axis_bootstrap`, `trial_incretin_trials` | `39_ctgov.py`, `40_trial_axis.py`, `42_trial_axis_contrast.py`, `43_trial_robustness.py` |

The figures were drawn from these tables; their plotted values are in the Source Data.

## Licence

Code: MIT licence (`LICENSE`). Source data (`results/`): [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
The input data remain under the terms of their providers.
