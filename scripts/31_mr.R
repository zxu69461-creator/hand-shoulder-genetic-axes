# MR estimation for scripts/31_mr.py: args = harmonised TSV in, results TSV out.
suppressMessages({library(MendelianRandomization); library(data.table)})
a <- commandArgs(TRUE)
h <- fread(a[1])
res <- rbindlist(lapply(split(h, by = c("exposure", "outcome", "cohort")), function(d) {
  k <- nrow(d)
  if (k < 3) return(NULL)
  inp <- mr_input(bx = d$bx, bxse = d$bxse, by = d$by, byse = d$byse, snps = d$rsid)
  ivw <- mr_ivw(inp, model = "random")
  wm <- mr_median(inp, weighting = "weighted", iterations = 10000, seed = 1)
  eg <- mr_egger(inp)
  data.table(exposure = d$exposure[1], outcome = d$outcome[1], cohort = d$cohort[1], axis = d$axis[1], n_snp = k,
             method = c("IVW", "weighted median", "MR-Egger", "MR-Egger intercept"),
             estimate = c(ivw@Estimate, wm@Estimate, eg@Estimate, eg@Intercept),
             se = c(ivw@StdError, wm@StdError, eg@StdError.Est, eg@StdError.Int),
             p = c(ivw@Pvalue, wm@Pvalue, eg@Causal.pval, eg@Pleio.pval),
             heterogeneity_q = c(ivw@Heter.Stat[1], NA, NA, NA), heterogeneity_p = c(ivw@Heter.Stat[2], NA, NA, NA))
}))
cat("MendelianRandomization", as.character(packageVersion("MendelianRandomization")), "\n")
fwrite(res, a[2], sep = "\t")
