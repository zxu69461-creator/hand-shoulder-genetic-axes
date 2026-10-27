# MVMR for scripts/32_mvmr.py: arg = TSV (rsid, b1, s1 [BMI], b2, s2 [partner], by, sy); prints results TSV to stdout.
suppressMessages({library(MendelianRandomization); library(data.table)})
d <- fread(commandArgs(TRUE)[1])
d <- d[complete.cases(d)]
inp <- mr_mvinput(bx = as.matrix(d[, .(b1, b2)]), bxse = as.matrix(d[, .(s1, s2)]), by = d$by, byse = d$sy, snps = d$rsid,
                  exposure = c("BMI", "partner"))
r <- mr_mvivw(inp, model = "random")
Fc <- tryCatch(r@CondFstat, error = function(e) c(NA, NA))
out <- data.table(exposure = c("BMI", "partner"), n_snp = nrow(d), estimate = r@Estimate, se = r@StdError, p = r@Pvalue,
                  cond_F = as.numeric(Fc))
fwrite(out, file = "", sep = "\t")
