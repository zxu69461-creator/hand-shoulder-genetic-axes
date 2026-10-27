#!/usr/bin/env bash
# Download the public inputs into $DATA_DIR (default: ./data). Re-running skips files already present.
# Reference panels that are distributed as archives by their authors (LDSC European LD scores and HapMap 3 SNP list;
# MAGMA 1000 Genomes European reference and NCBI 37.3 gene locations) are checked, not downloaded: see README.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DATA=${DATA_DIR:-$ROOT/data}
REF=${REF_DIR:-$DATA/ref}
get() {  # url dest
  [[ -s "$2" ]] && return 0
  mkdir -p "$(dirname "$2")"
  curl -fL --retry 5 -C - -o "$2.part" "$1" && mv "$2.part" "$2"
}
EBI=https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics

# GWAS Catalog: UK Biobank-based soft-tissue GWAS and the exposure GWAS
while read -r path; do
  get "$EBI/$path" "$DATA/gwas/$(basename "$path")"
done <<'LIST'
GCST90104001-GCST90105000/GCST90104907/harmonised/GCST90104907.h.tsv.gz
GCST90301001-GCST90302000/GCST90301252/harmonised/GCST90301252.h.tsv.gz
GCST007001-GCST008000/GCST007581/WibergA_2019_UKBB.txt
GCST90000001-GCST90001000/GCST90000512/harmonised/34111113-GCST90000512-EFO_1000941.h.tsv.gz
GCST90044001-GCST90045000/GCST90044562/harmonised/34737426-GCST90044562-EFO_1001435.h.tsv.gz
GCST90044001-GCST90045000/GCST90044569/harmonised/34737426-GCST90044569-EFO_1000438.h.tsv.gz
GCST90044001-GCST90045000/GCST90044558/harmonised/34737426-GCST90044558-EFO_0009666.h.tsv.gz
GCST90012001-GCST90013000/GCST90012107/harmonised/32042192-GCST90012107-EFO_0004696.h.tsv.gz
GCST90012001-GCST90013000/GCST90012102/harmonised/32042192-GCST90012102-EFO_0004908.h.tsv.gz
GCST90012001-GCST90013000/GCST90012112/harmonised/32042192-GCST90012112-EFO_0004908.h.tsv.gz
GCST90012001-GCST90013000/GCST90012109/harmonised/32042192-GCST90012109-EFO_0004696.h.tsv.gz
GCST90012001-GCST90013000/GCST90012103/harmonised/32042192-GCST90012103-EFO_0004908.h.tsv.gz
GCST90012001-GCST90013000/GCST90012113/harmonised/32042192-GCST90012113-EFO_0004908.h.tsv.gz
GCST90019001-GCST90020000/GCST90019511/harmonised/GCST90019511.h.tsv.gz
GCST90002001-GCST90003000/GCST90002244/harmonised/34059833-GCST90002244-EFO_0004541.h.tsv.gz
GCST006001-GCST007000/GCST006867/harmonised/30054458-GCST006867-EFO_0001360.h.tsv.gz
GCST002001-GCST003000/GCST002783/harmonised/25673413-GCST002783-EFO_0004340.h.tsv.gz
LIST

# ReproGen 2021 age at natural menopause
get https://www.reprogen.org/reprogen_ANM_201K_170621.txt.gz "$DATA/gwas/reprogen_ANM_201K_170621.txt.gz"

# FinnGen release 13 endpoints
for ep in G6_CARPTU M13_TRIGGERFINGER M13_DEQUERVAIN M13_DUPUTRYEN M13_ADHCAPSULITIS M13_PLANTARFASC \
          M13_ENTESOPATHYOTH M13_ROTATORCUFF; do
  get "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_$ep.gz" \
      "$DATA/finngen_raw/finngen_R13_$ep.gz"
done

# UK Biobank, Neale lab round 2, by sex (hospital main diagnosis ICD-10 G56 and M72)
for f in G56.gwas.imputed_v3.female G56.gwas.imputed_v3.male M72.gwas.imputed_v3.female M72.gwas.imputed_v3.male; do
  get "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/$f.tsv.bgz" "$DATA/neale/$f.tsv.bgz"
done

# single-cell atlases: Achilles and hamstring tendon (CELLxGENE Discover), palmar fascia (GEO GSE173252)
get https://datasets.cellxgene.cziscience.com/accd7d2e-d67d-4b91-9ea1-f2d1eee78ebf.h5ad "$DATA/sc/achilles_tendon.h5ad"
get https://datasets.cellxgene.cziscience.com/9fbe9e66-2649-4bc7-b729-f89e37423552.h5ad "$DATA/sc/hamstring_tendon.h5ad"
get https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173252/suppl/GSE173252_dd.rds.gz "$DATA/sc/GSE173252_dd.rds.gz"

# annotation for the FTO regional plot
get https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz \
    "$REF/gencode.v46.annotation.gtf.gz"
get https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip "$REF/plink.GRCh38.map.zip"

# reference panels supplied by their authors
for f in "$REF/eur_w_ld_chr/1.l2.ldscore.gz" "$REF/w_hm3.snplist" "$REF/g1000_eur/g1000_eur.bed" "$REF/NCBI37.3.gene.loc"; do
  [[ -s "$f" ]] || echo "missing reference file: $f (see README, 'Reference panels')"
done
echo "inputs in $DATA"
