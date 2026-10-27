# Additional Genomic SEM models for the cross-cohort set (xc): each disorder measured in FinnGen and in UK Biobank-based
# GWAS forms a disease factor (cohort-specific variance left in the residuals); single-cohort disorders enter as
# single-indicator factors. The factor correlations are the disease-level genetic correlations.
# 'dis' adds an orthogonal method factor for the control group shared by the FinnGen musculoskeletal endpoints.
dis_core <- "
dCTS =~ NA*CTS_FG + CTS_UKB
dTF =~ NA*TF_FG + TF_UKB
dDUP =~ NA*DUP_FG + DUP_UKB
dFS =~ NA*FS_FG + FS_UKB
dDQ =~ NA*DQ_FG
dPF =~ NA*PF_FG
dENTH =~ NA*ENTH_UKB
DQ_FG ~~ 0*DQ_FG
PF_FG ~~ 0*PF_FG
ENTH_UKB ~~ 0*ENTH_UKB
"
ctrl <- "
CTRL =~ NA*TF_FG + DQ_FG + DUP_FG + PF_FG + FS_FG
CTRL ~~ 0*dCTS + 0*dTF + 0*dDUP + 0*dFS + 0*dDQ + 0*dPF + 0*dENTH
"
extra_models <- list(xc = list(dis = paste0(dis_core, ctrl), dis_noctrl = dis_core))
