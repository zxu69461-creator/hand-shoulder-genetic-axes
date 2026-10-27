# Genomic SEM structure of the hand and shoulder soft-tissue disorders (called by 21_gsem.py).
# Usage: Rscript 20_gsem_structure.R <traits.tsv> <LD score directory> <output directory> [extra models .R]
#   traits.tsv: label, file, set[, m1, m2, ctrl]; m1/m2 = factor of each trait in the one- and two-factor models
#   ('' = not loaded); ctrl = 1 for FinnGen musculoskeletal endpoints that share one control group, giving model
#   m3 = m2 + an orthogonal method factor CTRL on those traits.
# Writes per set: LDSC objects (all, odd, even chromosomes), EFA loadings (odd), CFA fits and estimates (even),
# genetic correlations with jackknife SE (all).
suppressMessages({library(GenomicSEM); library(data.table)})
a <- commandArgs(TRUE)
tr <- fread(a[1], colClasses = "character")
ld <- paste0(normalizePath(a[2]), "/")
OUT <- a[3]; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
extra <- if (length(a) > 3) a[4] else NA
cat("GenomicSEM", as.character(packageVersion("GenomicSEM")), "lavaan", as.character(packageVersion("lavaan")), "\n")

cfa_spec <- function(lab, fac, ctrl = NULL) {
  g <- split(lab, fac); g <- g[names(g) != ""]
  if (!is.null(ctrl)) g$CTRL <- lab[ctrl == "1"]
  s <- sapply(names(g), function(f) paste0(f, " =~ NA*", paste(g[[f]], collapse = " + ")))
  s <- c(s, paste0(names(g), " ~~ 1*", names(g)))
  if (!is.null(ctrl)) s <- c(s, paste0("CTRL ~~ 0*", setdiff(names(g), "CTRL")))
  paste(s, collapse = "\n")
}

fit_models <- function(cov, S0, specs, tag) {
  fits <- list(); ests <- list()
  for (m in names(specs)) {
    cat("== CFA", S0, m, "\n", specs[[m]], "\n")
    r <- tryCatch(usermodel(cov, estimation = "DWLS", model = specs[[m]], CFIcalc = TRUE, std.lv = TRUE),
                  error = function(e) {cat("CFA error:", conditionMessage(e), "\n"); NULL})
    if (is.null(r)) next
    fits[[m]] <- data.table(set = S0, model = m, r$modelfit)
    ests[[m]] <- data.table(set = S0, model = m, r$results)
  }
  fwrite(rbindlist(fits, fill = TRUE), file.path(OUT, paste0("cfa_fit_", S0, tag, ".tsv")), sep = "\t")
  fwrite(rbindlist(ests, fill = TRUE), file.path(OUT, paste0("cfa_est_", S0, tag, ".tsv")), sep = "\t")
}

for (S0 in unique(tr$set)) {
  tt <- tr[tr$set == S0]
  res <- list()
  for (sel in c("ALL", "ODD", "EVEN")) {
    cat("== ldsc", S0, sel, "\n")
    res[[sel]] <- ldsc(traits = tt$file, sample.prev = rep(NA, nrow(tt)), population.prev = rep(NA, nrow(tt)),
                       ld = ld, wld = ld, trait.names = tt$label, stand = TRUE,
                       select = if (sel == "ALL") FALSE else sel,
                       ldsc.log = file.path(OUT, paste0("ldsc_", S0, "_", sel)))
  }
  saveRDS(res, file.path(OUT, paste0("ldsc_", S0, ".rds")))

  # exploratory factor analysis on odd chromosomes
  Ss <- as.matrix(Matrix::nearPD(res$ODD$S_Stand, corr = TRUE)$mat)
  dimnames(Ss) <- list(tt$label, tt$label)
  efa <- rbindlist(lapply(1:min(3, nrow(tt) - 3), function(k) {
    f <- tryCatch(factanal(covmat = Ss, factors = k, rotation = if (k > 1) "promax" else "none"), error = function(e) NULL)
    if (is.null(f)) return(NULL)
    L <- unclass(f$loadings); colnames(L) <- paste0("F", seq_len(ncol(L)))
    data.table(set = S0, k = k, trait = rownames(L), L)
  }), fill = TRUE)
  fwrite(efa, file.path(OUT, paste0("efa_", S0, ".tsv")), sep = "\t")

  # confirmatory models on even chromosomes
  if (all(c("m1", "m2") %in% names(tt)) && any(!is.na(tt$m1) & tt$m1 != "")) {
    specs <- list(m1 = cfa_spec(tt$label, tt$m1), m2 = cfa_spec(tt$label, tt$m2))
    if ("ctrl" %in% names(tt) && sum(tt$ctrl == "1") >= 3) specs$m3 <- cfa_spec(tt$label, tt$m2, tt$ctrl)
    fit_models(res$EVEN, S0, specs, "")
  }
  # additional lavaan models: a file defining `extra_models`, a list keyed by set of named model strings
  if (!is.na(extra) && file.exists(extra)) {
    extra_models <- NULL; source(extra)
    if (length(extra_models[[S0]])) fit_models(res$EVEN, S0, extra_models[[S0]], "_extra")
  }

  k <- nrow(tt); SE <- matrix(0, k, k); SE[lower.tri(SE, diag = TRUE)] <- sqrt(diag(res$ALL$V_Stand))
  SE <- SE + base::t(SE) - diag(diag(SE))
  out <- data.table(expand.grid(trait1 = tt$label, trait2 = tt$label, stringsAsFactors = FALSE))
  out[, rg := as.vector(res$ALL$S_Stand)][, se := as.vector(SE)][, set := S0]
  fwrite(out, file.path(OUT, paste0("rg_", S0, ".tsv")), sep = "\t")
}
cat("done\n")
