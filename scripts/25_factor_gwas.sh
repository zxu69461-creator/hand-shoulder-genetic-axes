#!/usr/bin/env bash
# Axis factor GWAS (scripts/25_factor_gwas.R) for both axes: LDSC, sumstats in chunks, NPAR chunk processes at a time
# (re-entrant: finished chunks are skipped), then the join. About 0.5 GB of memory per chunk process.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
WORK=${WORK_DIR:-$ROOT/work}
IN=$WORK/factor_stage
OUT=$WORK/factor
LD=${LDSC_LD:-${REF_DIR:-${DATA_DIR:-$ROOT/data}/ref}/eur_w_ld_chr}
RS=${RSCRIPT:-Rscript}
NPAR=${NPAR:-8}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1   # one core per chunk process
mkdir -p "$OUT"
R="$RS $ROOT/scripts/25_factor_gwas.R $IN $LD $OUT"
for AX in TEN FIB; do
  $R ldsc "$AX"
  $R prep "$AX"
  n=$(ls "$OUT/chunks/$AX"/ss_*.rds | wc -l | tr -d ' ')
  seq 1 "$n" | xargs -P "$NPAR" -I{} $R run "$AX" {}
  $R join "$AX"
done
echo done
