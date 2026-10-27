# Dupuytren's / palmar fascia / dermis atlas (GSE173252, authors' Seurat object): the authors' UMAP coordinates with the
# labels used in 62_dd_means.R (broad annotation, tissue source, mesenchymal cluster names), for the Figure 1 atlas map.
# Usage: Rscript 63_dd_umap.R <GSE173252_dd.rds.gz> <output directory>   (writes dd_umap.tsv.gz)
suppressMessages({library(SeuratObject); library(data.table)})
a <- commandArgs(TRUE); OUT <- a[2]; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
obj <- readRDS(gzcon(file(a[1], "rb")))
md <- obj@meta.data
red <- names(obj@reductions); cat("reductions:", red, "\n")
key <- intersect(c("umap", "UMAP", "tsne"), red)[1]
E <- obj@reductions[[key]]@cell.embeddings
mes_name <- c("0" = "PI16+ universal fibroblasts", "2" = "Mural cells", "3" = "POSTN+ Dupuytren fibroblasts",
              "7" = "COMP+ THBS4+ fascial fibroblasts", "14" = "APOD+ fascia fibroblasts")
out <- data.table(cell = rownames(md), u1 = E[rownames(md), 1], u2 = E[rownames(md), 2],
                  annotation = as.character(md$annotation), source = as.character(md$source),
                  cluster = as.character(md$seurat_clusters), doublet = as.character(md$doublet_prediction),
                  sample = as.character(md$orig.ident))
out[, mes := ifelse(annotation == "Mesenchyme", mes_name[cluster], NA_character_)]
fwrite(out, file.path(OUT, "dd_umap.tsv.gz"), sep = "\t")
print(out[, .N, by = .(source)]); print(out[, .N, by = .(annotation)]); print(out[, uniqueN(sample)])
cat("embedding:", key, "\ndone\n")
