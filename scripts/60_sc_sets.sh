#!/usr/bin/env bash
# Cell-state gene sets of the three atlases (60_sc_specificity.py) and the palmar fascia UMAP (63_dd_umap.R).
# Achilles and hamstring tendon: all author cell types, and the fibroblast (and chondrocyte) states alone, so that a
# state's gene set reflects what distinguishes it from the other fibroblast states. Palmar fascia (GSE173252): tissue
# compartments (mesenchyme split by tissue source) and the mesenchymal states.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DATA=${DATA_DIR:-$ROOT/data}
WORK=${WORK_DIR:-$ROOT/work}
PY=${PYTHON:-python3}
RS=${RSCRIPT:-Rscript}
S=$ROOT/scripts
mkdir -p "$WORK/sc"
$PY "$S/60_sc_specificity.py" hamstring "$DATA/sc/hamstring_tendon.h5ad" author_cell_type
$PY "$S/60_sc_specificity.py" hamstring_fib "$DATA/sc/hamstring_tendon.h5ad" author_cell_type 'fibroblast'
$PY "$S/60_sc_specificity.py" achilles "$DATA/sc/achilles_tendon.h5ad" author_cell_type
$PY "$S/60_sc_specificity.py" achilles_fib "$DATA/sc/achilles_tendon.h5ad" author_cell_type 'fibroblast|chondrocyte'
$RS "$S/62_dd_means.R" "$DATA/sc/GSE173252_dd.rds.gz" "$WORK/sc/dupuytren"
$PY "$S/60_sc_specificity.py" dupuytren_tissue "$WORK/sc/dupuytren_tissue.means.tsv.gz"
$PY "$S/60_sc_specificity.py" dupuytren_mes "$WORK/sc/dupuytren_mes.means.tsv.gz"
$RS "$S/63_dd_umap.R" "$DATA/sc/GSE173252_dd.rds.gz" "$WORK/sc"
$PY "$S/64_dd_umap_ingest.py"
