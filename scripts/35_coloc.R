# coloc.abf for scripts/35_coloc.py: args = region TSV (chrom pos b1 s1 beta se maf), case fraction, N disorder, N trait.
# Prints the posterior probabilities and the variant with the largest SNP-level posterior under H4.
suppressMessages({library(coloc); library(data.table)})
a <- commandArgs(TRUE)
d <- fread(a[1])
d[, snp := paste(chrom, pos, sep = ":")]
r <- coloc.abf(dataset1 = list(beta = d$b1, varbeta = d$s1^2, type = "cc", s = as.numeric(a[2]), N = as.numeric(a[3]),
                               MAF = d$maf, snp = d$snp),
               dataset2 = list(beta = d$beta, varbeta = d$se^2, type = "quant", N = as.numeric(a[4]), MAF = d$maf, snp = d$snp))
s <- r$summary
top <- r$results[which.max(r$results$SNP.PP.H4), ]                     # most likely shared causal variant under H4
cat(sprintf("nsnps=%d PP0=%.4f PP1=%.4f PP2=%.4f PP3=%.4f PP4=%.4f top_snp=%s top_snp_pp4=%.4f\n", as.integer(s["nsnps"]),
            s["PP.H0.abf"], s["PP.H1.abf"], s["PP.H2.abf"], s["PP.H3.abf"], s["PP.H4.abf"], top$snp, top$SNP.PP.H4))
