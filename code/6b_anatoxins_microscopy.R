#### Comparing microscopy data with regard to anatoxin concentrations
### Jordan Zabrecky
## last edited: 01.20.2026

# This script examines how communities as identified by microscopy
# change with increasing anatoxin concentrations by <INSERT>

# also currently has all environmental parameters with an RDA
# go through RDA variable selection and adjust
# https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html

# likely will move dbRDA to supplemental script
# also explore collinearity there

#### (1) Loading libraries & data ####

# set seed for reproducibility
set.seed(2025)

# libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot", "indicspecies"), require, character.only = T)

# read in relative abundance files (data transformed in previous script, "4a_amongrivers_microscopy.R")
data <- lapply(list.files(path = "./data/morphological/transformed/", pattern = ".csv"),
               function(x) read.csv(paste("./data/morphological/transformed/", x, sep = "")))
sample_types <- c("nt", "tac", "tm")
names(data) <- sample_types

# read in environmental covariates & toxin data
env_tog <- read.csv("./data/field_and_lab/env_tox_standardized_together.csv")
env_sep <- read.csv("./data/field_and_lab/env_tox_standardized_byriver.csv")

# lastly, merge with our community data to get rows formatted the same
data_tog <- lapply(data, function(x) left_join(x, env_tog, by = c("site", "site_reach", "field_date")))
data_sep <- lapply(data, function(x) left_join(x, env_sep, by = c("site", "site_reach", "field_date")))

# for microcoleus & anabaena remove rows where there are NAs for anatoxins
which(is.na(data_tog$tm$TM_ATX_all_ug_orgmat_g)) #22
which(is.na(data_tog$tac$TAC_ATX_all_ug_orgmat_g))

# TM 22 was a "not sure if this is microcoleus sample" and there was not enough material
# to analyze for anatoxins
data_tog$tm <- data_tog$tm[-22,]
data_sep$tm <- data_sep$tm[-22,]

# split into lists by river
data_sep_list <- list()
for(i in 1:length(data_sep)) {
  data_sep_list[[i]] <- split(data_sep[[i]], data_sep[[i]]$site)
}
names(data_sep_list) <- names(data_sep)

#### (2) Functions for Analyses ####

# load community analyses functions from other script
source("./code/supplemental_code/S4a_community_analyses_func.R")

# set start column of community data
start_col <- 5

#### (3) PERMANOVA ####

# Do communities significantly change with anatoxin concentrations?

## (a) rivers together
runPERMANOVA(data_tog$tm, start_col, end_col = ncol(data$tm), 
             group = data_tog$tm$`TM_ATX_all_ug_orgmat_g`)
# TM: significant **
runPERMANOVA(data_tog$tac, start_col, end_col = ncol(data$tac), 
             group = data_tog$tac$`TAC_ATX_all_ug_orgmat_g`)
# TAC: significant *
runPERMANOVA(data_tog$nt, start_col, end_col = ncol(data$nt), 
             group = data_tog$nt$`mean_ATX_all_ug_orgmat_g`, na.action = "na.omit")
# NT: significant **

## (b) rivers separately
# (omitting Salmon because only one sample was toxic and very minorly)
runPERMANOVA(data_sep_list$tm$`SFE-M`, start_col, end_col = ncol(data$tm),
             group = data_sep_list$tm$`SFE-M`$TM_ATX_all_ug_orgmat_g)
# SFE TM: significant: **
runPERMANOVA(data_sep_list$tac$`SFE-M`, start_col, end_col = ncol(data$tac),
             group = data_sep_list$tac$`SFE-M`$TAC_ATX_all_ug_orgmat_g)
# SFE TAC: significant ***
runPERMANOVA(data_sep_list$tac$`RUS`, start_col, end_col = ncol(data$tac),
             group = data_sep_list$tac$`RUS`$TAC_ATX_all_ug_orgmat_g)
# RUS TAC: significant *

lapply(data_sep_list$nt, function(x) runPERMANOVA(x, start_col, end_col = ncol(data$nt),
                                                  group = x$`mean_ATX_all_ug_orgmat_g`,
                                                  na.action = "na.omit"))
# NT: SAL NS, RUS NS, SFE-M *

# How about if we have grouping of anatoxin concentrations, such as:
# <insert grouping criteria>

# look at un-transformed anatoxin concentrations
atx <- read.csv("")

#### (4) (db)RDA analyses ####

## (a) variable selection

# check which variables are correlated
correlations <- cor(env_tog[,4:ncol(env_tog)], use = "complete.obs")
view(which(correlations > 0.7 & correlations != 1, arr.ind = TRUE))
# correlated: DO & temp, assumed pH & DO, ammonium & nitrate & DIN, 
# conductivity & Cl Na K & Mg, DOC & Br amm DIN, SO4 & Cl NA
# we are missing some conductivity observations, so choose most correlated instead (Mg)
# will select one from each correlated group:
# included in RDA analyses: <INSERT>

## (b) analyzing rivers together with data standardized together

# run dbRDA for all sample types with all rivers together
dbRDA_tog <- lapply(sample_types, function(x) run_dbRDA(data_tog[[x]], 
                                                       start_col = 5, 
                                                       end_col = ncol(data[[x]]),
                                                       mat_atx = x))
names(dbRDA_tog) <- sample_types

# plot
lapply(dbRDA_tog, function(x) plot(x$object))

# are models significant?
lapply(dbRDA_tog, function(x) print(x$model_sig))
# nt ***, tac ***, tm ***

# what is the rsquared?
lapply(dbRDA_tog, function(x) print(x$rsquared))
# adjusted: nt TBD, tac 0.231, tm 0.457

# significant variables?
lapply(dbRDA_tog, function(x) print(x$sig_variables))
# NT:
# TAC: TAC ATX**, oPhos***
# TM: TM ATX***, oPhos*, TDC**, Mg**

## (c) separated standardization

# run dbRDA for all sample types with all rivers together
# but covariates standardized for each river
dbRDA_sep <- lapply(sample_types, function(x) run_dbRDA(data_sep[[x]], 
                                                        start_col = 5, 
                                                        end_col = ncol(data[[x]]),
                                                        mat_atx = x,
                                                        na.action = "na.omit"))
names(dbRDA_sep) <- sample_types

# plot
lapply(dbRDA_sep, function(x) plot(x$object))

# are models significant?
lapply(dbRDA_sep, function(x) print(x$model_sig))
# nt NS, tac ***, tm NS

# what is the rsquared?
lapply(dbRDA_tog, function(x) print(x$rsquared))
# adjusted: nt TBD, tac 0.265, tm 0.121

# significant variables?
lapply(dbRDA_tog, function(x) print(x$sig_variables))
# NT:
# TAC: TAC ATX***, oPhos**, DOC *
# TM: none

## (b) apply to each river separately with separate standardization

# run dbRDA for all sample types with all rivers separately - NT
# want to figure out what to do with ATX
data_sep$nt

# russian TAC
russian_TAC = run_dbRDA(data_sep$tac$RUS, 
                        start_col = 5, 
                        end_col = ncol(data$tac),
                        mat_atx = "tac")

# sfkeel TAC
sfkeel_TAC = run_dbRDA(data_sep$tac$`SFE-M`,
                       start_col = 5,
                       end_col = ncol(data$tac),
                       mat_atx = "tac")

# sfkeel TM
sfkeel_TM = run_dbRDA(data_sep_list$tm$`SFE-M`,
                      start_col = 5,
                      end_col = ncol(data$tm),
                      mat_atx = "tm")

# add into list
individual_t_rdas <- list(russian_TAC, sfkeel_TAC, sfkeel_TM)
names(individual_t_rdas) <- c("russian_TAC", "sfkeel_TAC", "sfkeel_TM")

# plot
lapply(individual_t_rdas, function(x) plot(x$object))

# are models significant?
lapply(individual_t_rdas, function(x) print(x$model_sig))
# russian tac no, sfkeel tac yes, sfkeel tm no

# what is the rsquared?
lapply(individual_t_rdas, function(x) print(x$rsquared))
# adjusted: nt TBD, tac 0.265, tm 0.121

# significant variables?
lapply(individual_t_rdas, function(x) print(x$sig_variables))
# <insert>

#### (5) Species Indicator Analyses ####







#### OLD CODE BELOW: #####


# lastly, merge with our community data to get rows formatted the same
data_tog <- lapply(data, function(x) left_join(x, env_tog, by = c("site", "site_reach", "field_date")))
data_sep <- lapply(data, function(x) left_join(x, env_sep, by = c("site", "site_reach", "field_date")))

# for microcoleus & anabaena remove rows where there are NAs for anatoxins
which(is.na(data_tog$tm$TM_ATX_all_ug_orgmat_g)) #22
which(is.na(data_tog$tac$TAC_ATX_all_ug_orgmat_g)) #none!

# 22 was a "not sure if this is microcoleus sample" and there was not enough material
# to analyze for anatoxins
data_tog$tm <- data_tog$tm[-22,]
data_sep$tm <- data_sep$tm[-22,]

# split into lists by river
for(i in 1:length(data_sep)) {
  data_sep[[i]] <- split(data_sep[[i]], data_sep[[i]]$site)
}




#### (3) (Distance-Based) Redundancy Analysis ####

## ADD IN VARIABLE SELECTION

## (a) TAC

# all samples together
tac_all_rda <- dbrda(data = data_tog$tac, 
                     data_tog$tac[,5:ncol(data$tac)] ~ 
                       TAC_ATX_all_ug_orgmat_g + nitrate_mg_N_L + oPhos_ug_P_L + ammonium_mg_N_L +
                       temp_C + cond_uS_cm + TDC_mg_L + DOC_mg_L, na.action = "na.exclude")
plot(tac_all_rda)
anova.cca(tac_all_rda, step = 1000, by = "term") # significant for atx*** and orthophosphate*

# Russian River
tac_rus_rda <- dbrda(data = data_sep$tac$RUS, 
                     data_sep$tac$RUS[,5:ncol(data$tac)] ~ 
                       TAC_ATX_all_ug_orgmat_g + nitrate_mg_N_L + oPhos_ug_P_L + ammonium_mg_N_L +
                       temp_C + cond_uS_cm + TDC_mg_L + DOC_mg_L, na.action = "na.exclude")
plot(tac_rus_rda)
anova.cca(tac_rus_rda, step = 1000, by = "term") # significant for atx*

# South Fork Eel River
tac_sfe_rda <- dbrda(data = data_sep$tac$`SFE-M`, 
                     data_sep$tac$`SFE-M`[,5:ncol(data$tac)] ~ 
                       TAC_ATX_all_ug_orgmat_g + nitrate_mg_N_L + oPhos_ug_P_L + ammonium_mg_N_L +
                       temp_C + cond_uS_cm + TDC_mg_L + DOC_mg_L, na.action = "na.exclude")
plot(tac_sfe_rda)
anova.cca(tac_sfe_rda, step = 1000, by = "term") # not significant for anything

## (b) TM

# all samples together
tm_all_rda <- dbrda(data = data_tog$tm, 
                     data_tog$tm[,5:ncol(data$tm)] ~ 
                       TM_ATX_all_ug_orgmat_g + nitrate_mg_N_L + oPhos_ug_P_L + ammonium_mg_N_L +
                       temp_C + cond_uS_cm, na.action = "na.exclude")
plot(tm_all_rda)
anova.cca(tm_all_rda, step = 1000, by = "term") # significant: ATX***, ophos*, cond*

# salmon river only (EXCLUDING AS THERE ARE ONLY TWO TIME POINTS)
tm_sal_rda <- dbrda(data = data_sep$tm$SAL, 
                    data_sep$tm$SAL[,5:ncol(data$tm)] ~ 
                      TM_ATX_all_ug_orgmat_g + nitrate_mg_N_L + oPhos_ug_P_L + 
                      ammonium_mg_N_L, na.action = "na.exclude")
plot(tm_sal_rda)
anova.cca(tm_sal_rda, step = 1000, by = "term") # no signficant

# south fork eel only
tm_sfe_rda <- dbrda(data = data_sep$tm$`SFE-M`, 
                    data_sep$tm$`SFE-M`[,5:ncol(data$tm)] ~ 
                      TM_ATX_all_ug_orgmat_g + nitrate_mg_N_L + oPhos_ug_P_L + ammonium_mg_N_L +
                      temp_C + cond_uS_cm, na.action = "na.exclude")
plot(tm_sfe_rda)
anova.cca(tm_sfe_rda, step = 1000, by = "term") # anatoxins significant *, cond significant *

## (c) NT

# will complete after TITAN analysis

#### (4) TITAN ####

## (a) TM

# first establish that anatoxin is associated with significant changes in composition
adonis2(data = data_tog$tm, data_tog$tm[,5:ncol(data$tm)] ~ 
          TM_ATX_all_ug_orgmat_g) # yes**

# rivers separately
adonis2(data = data_sep$tm$`SFE-M`, data_sep$tm$`SFE-M`[,5:ncol(data$tm)] ~ 
          TM_ATX_all_ug_orgmat_g) # yes**
adonis2(data = data_sep$tm$`SAL`, data_sep$tm$`SAL`[,5:ncol(data$tm)] ~ 
          TM_ATX_all_ug_orgmat_g)
# no, but only one sample had small amounts of atx so dropping this site

## (b) TAC

# first establish that anatoxin is associated with significant changes in composition
adonis2(data = data_tog$tac, data_tog$tac[,5:ncol(data$tac)] ~ 
          TAC_ATX_all_ug_orgmat_g) # yes*

# rivers separately
adonis2(data = data_sep$tac$`SFE-M`, data_sep$tac$`SFE-M`[,5:ncol(data$tac)] ~ 
          TAC_ATX_all_ug_orgmat_g) # yes**
adonis2(data = data_sep$tac$`RUS`, data_sep$tac$`RUS`[,5:ncol(data$tac)] ~ 
          TAC_ATX_all_ug_orgmat_g) # yes*
# no, but only one sample had small amounts of atx so dropping this site


# still need to look at dispersion

# Our data is too sparse and imperfect for TITAN :(
# will move attempts to external script as well as RDA
titan.PDIR <- titan(env = test$TAC_ATX_all_ug_orgmat_g,
                    txa =  test[,5:14])

plot_taxa_ridges(titan.PDIR,
                 xlabel = "PDIR")

?txa.screen()
TITAN2::txa.screen(data_tog$tac[,5:ncol(data$tac)])
view(data_tog$tac)

titan.PDIR$sppmax %>%
  as.data.frame() %>%
  select(zenv.cp, freq, maxgrp, IndVal, purity, reliability, filter) %>%
  filter(filter != 0) %>%
  arrange(maxgrp, zenv.cp)

# remove taxa that occur in less than 3 samples
test <- data_tog$tac[,-which(colSums(data_tog$tac[,5:ncol(data$tac)]) < 3)]

test <- data_tog$tac[,-(4 + which(apply(data_tog$tac[,5:ncol(data$tac)], 2, function(x) length(which(x != 0))) < 3))]
which(colSums())

glades.titan <- titan(glades.env, glades.taxa,
                      minSplt = 5, numPerm = 250, boot = TRUE, nBoot = 500, imax = FALSE,
                      ivTot = FALSE, pur.cut = 0.95, rel.cut = 0.95, ncpus = 8, memory = FALSE
)

## (c) NT

# to be completed for both TM atx and TAC atx

# maybe TITAN will work for bacterial data

##### (5) Species Indicator Analysis for Groupings ####

# determine categories for ATX:

# based on river separately
ggplot(data = data_sep$tac$RUS, aes(x = TAC_ATX_all_ug_orgmat_g)) +
  geom_histogram() +
  facet_wrap(~site)

# maybe don't want to use log-transformed data here but will do a prelim
# <-1 no anatoxins detected
# 0-1 is low anatoxins
# >1 is high anatoxins

test <- lapply(data_tog, function(x) x <- x %>% 
                                    mutate(atx_group = case_when(TAC_ATX_all_ug_orgmat_g < -0.5 ~ "non-detect",
                                                                 TAC_ATX_all_ug_orgmat_g < 0 ~ "low",
                                                                 TAC_ATX_all_ug_orgmat_g < 1 ~ "medium",
                                                                 TAC_ATX_all_ug_orgmat_g > 1 ~ "high")))
test2 <- lapply(data_tog, function(x) x <- x %>% 
                 mutate(atx_group = case_when(TM_ATX_all_ug_orgmat_g < -0.5 ~ "non-detect",
                                              TM_ATX_all_ug_orgmat_g < 0 ~ "low",
                                              TM_ATX_all_ug_orgmat_g < 1 ~ "medium",
                                              TM_ATX_all_ug_orgmat_g > 1 ~ "high")))

summary(multipatt(test$tac[,5:ncol(data$tac)], test$tac$atx_group, func = "r.g", control = how(nperm = 999)))
summary(multipatt(test2$tm[,5:ncol(data$tm)], test2$tm$atx_group, func = "r.g", control = how(nperm = 999)))
# epithemia associated with high and medium
# non epithemia diatoms associated with low and non-detect

# how about within a river
test3 <- lapply(data_sep$tac, function(x) x <- x %>% 
                  mutate(atx_group = case_when(TAC_ATX_all_ug_orgmat_g < -0.5 ~ "non-detect",
                                               TAC_ATX_all_ug_orgmat_g < 0 ~ "low",
                                               TAC_ATX_all_ug_orgmat_g < 1 ~ "medium",
                                               TAC_ATX_all_ug_orgmat_g > 1 ~ "high")))
summary(multipatt(test3$`SFE-M`[,5:ncol(data$tac)], test3$`SFE-M`$atx_group, func = "r.g", control = how(nperm = 999)))
# high is microcoleus while low is nodularia, nostoc, low-medium is epithemia

summary(multipatt(test3$`RUS`[,5:ncol(data$tac)], test3$`RUS`$atx_group, func = "r.g", control = how(nperm = 999)))
# medium + non-detect
# aka unclear

test4 <- lapply(data_sep$tm, function(x) x <- x %>% 
                  mutate(atx_group = case_when(TM_ATX_all_ug_orgmat_g < -0.5 ~ "non-detect",
                                               TM_ATX_all_ug_orgmat_g < 0 ~ "low",
                                               TM_ATX_all_ug_orgmat_g < 1 ~ "medium",
                                               TM_ATX_all_ug_orgmat_g > 1 ~ "high")))
summary(multipatt(test4$`SFE-M`[,5:ncol(data$tac)], test4$`SFE-M`$atx_group, func = "r.g", control = how(nperm = 999)))
# maybe try differential abundance analysis
# epithemia and leptolyngbya for high
