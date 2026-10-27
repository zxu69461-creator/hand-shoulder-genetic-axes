# Dupuytren's disease / non-pathogenic palmar fascia / dermis atlas (GSE173252; authors' Seurat object, updated to Seurat v5):
# per-cell-type mean expression (counts per 10,000, linear) for scripts/60_sc_specificity.py, two label sets:
#   tissue — authors' broad annotation, with mesenchyme split by tissue source (Dupuytren's / fascia / dermis);
#   mes    — authors' clusters within mesenchyme (≥ 50 cells), named from their tissue composition and marker means
#            (printed): 3 = Dupuytren-only POSTN+ PDPN+ FAP+ SCX+; 7 = COMP+ THBS4+ (fascia and Dupuytren's);
#            14 = APOD+ (93% fascia); 0 = PI16+ CD34+ universal; 2 = ACTA2+ MYH11+ RGS5+ mural.
# Usage: Rscript 62_dd_means.R <GSE173252_dd.rds.gz> <out prefix>
suppressMessages({library(Seurat); library(Matrix); library(data.table)})
a <- commandArgs(TRUE)
obj <- UpdateSeuratObject(readRDS(if (grepl("\\.gz$", a[1])) gzcon(file(a[1], "rb")) else a[1]))
md <- obj@meta.data
src <- c(DD = "Dupuytren's", SF = "fascia", Dermis = "dermis")[as.character(md$source)]
tissue <- ifelse(md$annotation == "Mesenchyme", paste0("Mesenchyme (", src, ")"), as.character(md$annotation))
mes_name <- c("0" = "PI16+ universal fibroblasts", "2" = "Mural cells", "3" = "POSTN+ Dupuytren fibroblasts",
              "7" = "COMP+ THBS4+ fascial fibroblasts", "14" = "APOD+ fascia fibroblasts")
mes <- ifelse(md$annotation == "Mesenchyme", mes_name[as.character(md$seurat_clusters)], NA)
cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")
lib <- Matrix::colSums(cnt)
cnt <- cnt %*% Diagonal(x = 1e4 / pmax(lib, 1)); dimnames(cnt) <- dimnames(GetAssayData(obj, assay = "RNA", layer = "counts"))
export <- function(lab, tag) {
  keep <- !is.na(lab) & as.character(md$doublet_prediction) != "TRUE"
  types <- names(which(table(lab[keep]) >= 50))
  M <- sapply(types, function(t) Matrix::rowMeans(cnt[, keep & lab == t, drop = FALSE]))
  rownames(M) <- rownames(cnt)
  fwrite(data.table(gene = rownames(M), M), paste0(a[2], "_", tag, ".means.tsv.gz"), sep = "\t")
  fwrite(data.table(cell_type = types, n = as.integer(table(lab[keep])[types])), paste0(a[2], "_", tag, ".counts.tsv"), sep = "\t")
  M
}
invisible(export(tissue, "tissue"))
Mm <- export(mes, "mes")
mk <- c("ACTA2", "PDPN", "FAP", "POSTN", "COL1A1", "COL3A1", "COL11A1", "PI16", "DPP4", "CD34", "APOD", "CXCL12", "COMP",
        "MYH11", "RGS5", "PRG4", "SFRP4", "WNT7B", "EPDR1", "THBS4", "SCX", "TNMD", "MKX", "COL15A1", "PPARG")
mk <- intersect(mk, rownames(Mm))
comp <- table(mes[md$annotation == "Mesenchyme"], md$source[md$annotation == "Mesenchyme"])
print(round(prop.table(comp, 1), 2))
print(round(t(Mm[mk, ]), 2))
cat("done\n")
