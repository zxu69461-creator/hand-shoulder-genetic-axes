# Axis factor GWAS with GenomicSEM. For each axis: LDSC on its three indicators (observed scale), SNP effects from Z
# statistics (sumstats, linprob), then commonfactorGWAS (DWLS) with the LDSC sampling covariance accounting for sample
# overlap between indicators and genomic control by the LDSC intercepts (GC = "standard"). The SNPs are processed in
# chunks of 20,000, one R process per chunk; because genomic control uses the LDSC intercepts in the covariance
# structure, chunking does not change any estimate.
# Usage: Rscript 25_factor_gwas.R <stage dir> <LD score dir> <output dir> <mode> <axis> [chunk]
#   ldsc  AX     covariance structure of the axis indicators -> ldsc_AX.rds
#   prep  AX     sumstats for the axis indicators -> chunks/AX/ss_<k>.rds; prints the number of chunks
#   run   AX k   commonfactorGWAS on chunk k -> chunks/AX/res_<k>.tsv.gz (skipped when present)
#   join  AX     all chunks in order -> factor_AX.tsv.gz
suppressMessages({library(GenomicSEM); library(data.table)})
a <- commandArgs(TRUE)
IN <- a[1]; ld <- paste0(normalizePath(a[2]), "/"); OUT <- a[3]; mode <- a[4]; ax <- a[5]
axes <- list(TEN = c("CTS", "TF", "DQ"), FIB = c("DUP", "FS", "PF"))
tr <- axes[[ax]]
CH <- file.path(OUT, "chunks", ax); dir.create(CH, recursive = TRUE, showWarnings = FALSE)
if (mode == "ldsc") {
  cov <- ldsc(traits = file.path(IN, paste0(tr, ".sumstats.gz")), sample.prev = rep(NA, 3), population.prev = rep(NA, 3),
              ld = ld, wld = ld, trait.names = tr, ldsc.log = file.path(OUT, paste0("ldsc_", ax)))
  saveRDS(cov, file.path(OUT, paste0("ldsc_", ax, ".rds")))
  quit(save = "no")
}
cov <- readRDS(file.path(OUT, paste0("ldsc_", ax, ".rds")))
if (mode == "prep") {
  ss <- sumstats(files = file.path(IN, paste0(tr, "_ss.tsv.gz")), ref = file.path(IN, "ref_maf.txt"), trait.names = tr,
                 se.logit = rep(FALSE, 3), OLS = rep(FALSE, 3), linprob = rep(TRUE, 3), maf.filter = 0.01,
                 parallel = TRUE, cores = 3)
  k <- split(seq_len(nrow(ss)), ceiling(seq_len(nrow(ss)) / 20000))
  for (i in seq_along(k)) saveRDS(ss[k[[i]], ], file.path(CH, sprintf("ss_%03d.rds", i)))
  cat(ax, "SNPs", nrow(ss), "chunks", length(k), "\n")
} else if (mode == "run") {
  f <- file.path(CH, sprintf("res_%03d.tsv.gz", as.integer(a[6])))
  if (!file.exists(f)) {
    res <- suppressWarnings(commonfactorGWAS(covstruc = cov, SNPs = readRDS(file.path(CH, sprintf("ss_%03d.rds", as.integer(a[6])))),
                                             estimation = "DWLS", parallel = FALSE, GC = "standard"))
    fwrite(res, paste0(f, ".tmp.gz"), sep = "\t"); file.rename(paste0(f, ".tmp.gz"), f)
    cat(ax, "chunk", a[6], nrow(res), "SNPs\n")
  }
} else if (mode == "join") {
  fs <- sort(list.files(CH, "^res_[0-9]+\\.tsv\\.gz$", full.names = TRUE))
  stopifnot(length(fs) == length(list.files(CH, "^ss_[0-9]+\\.rds$")))
  res <- rbindlist(lapply(fs, fread)); fwrite(res, file.path(OUT, paste0("factor_", ax, ".tsv.gz")), sep = "\t")
  cat(ax, "joined", length(fs), "chunks,", nrow(res), "SNPs\n")
}
