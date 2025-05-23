library(TwoSampleMR)
library(pleio)
library(readxl)
library(dplyr)
library(ieugwasr)
Sys.setenv(OPENGWAS_JWT="")
R <- extract_instruments(
  outcomes = "eqtl-a-", 
  p1 = 5e-06, 
  clump = TRUE, 
  r2 = 0.001,
  kb=100)

outcome_dat1 <- extract_outcome_data(snps = R$SNP, outcomes = 'ebi-a-')
dat1<- harmonise_data(exposure_dat =R,  outcome_dat = outcome_dat1)
res1<- mr(dat1)
pleio1<- mr_pleiotropy_test(dat1)
OR1<-generate_odds_ratios(res1)
het1 <- mr_heterogeneity(dat1)

