library(TwoSampleMR)
library(pleio)
library(readxl)
library(dplyr)
library(ieugwasr)
Sys.setenv(OPENGWAS_JWT="")
ROMO1 <- extract_instruments(
  outcomes = "eqtl-a-ENSG00000125995", 
  p1 = 5e-06, 
  clump = TRUE, 
  r2 = 0.001,
  kb=100)

outcome_dat1 <- extract_outcome_data(snps = ROMO1$SNP, outcomes = 'ebi-a-GCST90001935')
dat1<- harmonise_data(exposure_dat =ROMO1,  outcome_dat = outcome_dat1)
res1<- mr(dat1)
pleio1<- mr_pleiotropy_test(dat1)
OR1<-generate_odds_ratios(res1)
het1 <- mr_heterogeneity(dat1)

IL2RA<- extract_instruments(outcomes = "ebi-a-GCST90001935",p1=5e-06)
outcome_dat2 <- extract_outcome_data(snps = IL2RA$SNP, outcomes = 'prot-c-5356_2_3')
dat2<- harmonise_data(exposure_dat =IL2RA,  outcome_dat = outcome_dat2)
res2<- mr(dat2)
pleio2<- mr_pleiotropy_test(dat2)
OR2<-generate_odds_ratios(res2)
het2 <- mr_heterogeneity(dat2)

MIF <- extract_instruments(outcomes = "prot-c-5356_2_3",p1=5e-06)
outcome_dat3 <- extract_outcome_data(snps = MIF$SNP, outcomes = 'finn-b-RHEUMA_SEROPOS_OTH')
dat3<- harmonise_data(exposure_dat =MIF,  outcome_dat = outcome_dat3)
res3<- mr(dat3)
pleio3<- mr_pleiotropy_test(dat3)
OR3<-generate_odds_ratios(res3)
het3 <- mr_heterogeneity(dat3)


