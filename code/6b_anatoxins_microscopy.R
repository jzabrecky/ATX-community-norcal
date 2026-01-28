#### Comparing microscopy data with regard to anatoxin concentrations
### Jordan Zabrecky
## last edited: 01.27.2026

# This script examines how communities as identified by microscopy
# change with increasing anatoxin concentrations by <INSERT>

# also currently has all environmental parameters with an RDA
# go through RDA variable selection and adjust
# https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html

# likely will move dbRDA to supplemental script
# also explore collinearity there

#### (1) Loading libraries & data ####

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

#### (3) PERMANOVA (& BETADISP) ####

## Do communities significantly change with anatoxin concentrations?
## Using log-transformed and cenetered anatoxin concentrations here...

## (a) rivers together
set.seed(1)
runPERMANOVA(data_tog$tm, start_col, end_col = ncol(data$tm), 
             group = data_tog$tm$`TM_ATX_all_ug_orgmat_g`, 
             strata = data_tog$tm$site)
anova(betadisper(vegdist(data_tog$tm[,start_col:ncol(data$tm)], method = "bray"), 
                 data_tog$tm$`TM_ATX_all_ug_orgmat_g`))
# TM: PERMANOVA (barely) significant*, but dispersion ***
# which makes me think that visually, these would not be different
set.seed(1)
runPERMANOVA(data_tog$tac, start_col, end_col = ncol(data$tac), 
             group = data_tog$tac$`TAC_ATX_all_ug_orgmat_g`,
             strata = data_tog$tac$site)
anova(betadisper(vegdist(data_tog$tac[,start_col:ncol(data$tac)], method = "bray"), 
                 data_tog$tac$`TAC_ATX_all_ug_orgmat_g`))
# TAC: not significant, also no significant dispersion differences
set.seed(1)
runPERMANOVA(data_tog$nt, start_col, end_col = ncol(data$nt), 
<<<<<<< HEAD
             group = data_tog$nt$`mean_ATX_all_ug_orgmat_g`, na.action = "na.omit")
# NT: significant **
# also, check dispersion
anova(betadisper(vegdist(data_tog$tm[,start_col:ncol(data$tm)], method = "bray"), 
                 data_tog$tm$TM_ATX_all_ug_orgmat_g))
# TM: significant dispersion differences ***
anova(betadisper(vegdist(data_tog$tac[,start_col:ncol(data$tac)], method = "bray"), 
                 data_tog$tac$TAC_ATX_all_ug_orgmat_g))
# TAC: no significant dispersion differences
anova(betadisper(vegdist(data_tog$nt[,start_col:ncol(data$nt)], method = "bray"), 
                 data_tog$nt$mean_ATX_all_ug_orgmat_g))
# NT: significant dispersion differences *

# view NMDS to visualize differences
NMDS_tog <- lapply(names(data_tog), function(x) getNMDSdata(data_tog[[x]], start_col, 
                                                            end_col = ncol(data[[x]])))
names(NMDS_tog) <- sample_types
makeNMDSplot(NMDS_tog$tm, FALSE, FALSE, color = "TM_atx_category", shape = "site")
=======
             group = data_tog$nt$`mean_ATX_all_ug_orgmat_g`,
             strata = data_tog$nt$`site`[-which(is.na(data_tog$nt$mean_ATX_all_ug_orgmat_g))],
             na.action = "na.omit")
anova(betadisper(vegdist(data_tog$nt[,start_col:ncol(data$nt)], method = "bray"), 
                 data_tog$nt$`mean_ATX_all_ug_orgmat_g`))
# NT: significant * but also significant dispersion *
>>>>>>> aac9ab7456d7548273e268593b9919b453cf0ef2

## (b) rivers separately
# (omitting Salmon because only one sample was toxic and very minorly)
set.seed(1)
runPERMANOVA(data_sep_list$tm$`SFE-M`, start_col, end_col = ncol(data$tm),
             group = data_sep_list$tm$`SFE-M`$TM_ATX_all_ug_orgmat_g)
anova(betadisper(vegdist(data_sep_list$tm$`SFE-M`[,start_col:ncol(data$tm)], method = "bray"), 
                 data_sep_list$tm$`SFE-M`$`TM_ATX_all_ug_orgmat_g`))
# SFE TM: significant: **, so is dispersion **
set.seed(1)
runPERMANOVA(data_sep_list$tac$`SFE-M`, start_col, end_col = ncol(data$tac),
             group = data_sep_list$tac$`SFE-M`$TAC_ATX_all_ug_orgmat_g)
anova(betadisper(vegdist(data_sep_list$tac$`SFE-M`[,start_col:ncol(data$tac)], method = "bray"), 
                 data_sep_list$tac$`SFE-M`$`TAC_ATX_all_ug_orgmat_g`))
# SFE TAC: significant ***, dispersion ***
set.seed(1)
runPERMANOVA(data_sep_list$tac$`RUS`, start_col, end_col = ncol(data$tac),
             group = data_sep_list$tac$`RUS`$TAC_ATX_all_ug_orgmat_g)
<<<<<<< HEAD
# RUS TAC: significant *
=======
anova(betadisper(vegdist(data_sep_list$tac$`RUS`[,start_col:ncol(data$tac)], method = "bray"), 
                 data_sep_list$tac$`RUS`$`TAC_ATX_all_ug_orgmat_g`))
# RUS TAC: significant **, but dispersion is not

set.seed(1)
>>>>>>> aac9ab7456d7548273e268593b9919b453cf0ef2
lapply(data_sep_list$nt, function(x) runPERMANOVA(x, start_col, end_col = ncol(data$nt),
                                                  group = x$`mean_ATX_all_ug_orgmat_g`,
                                                  na.action = "na.omit"))
lapply(data_sep_list$nt, function(x) 
  anova(betadisper(vegdist(x[,start_col:ncol(data$nt)], method = "bray"), 
                  x$`mean_ATX_all_ug_orgmat_g`)))
# NT: SAL NS, RUS NS, SFE-M * with dispersion only significantly different in Russian

<<<<<<< HEAD
# To-do: PERMDISP, view NMDS to visualize differences

## How about just comparing samples with and without detectable toxin?
=======

## How about if we do it by grouping of anatoxin concentrations?
>>>>>>> aac9ab7456d7548273e268593b9919b453cf0ef2

## (a) rivers together
set.seed(1)
runPERMANOVA(data_tog$tm, start_col, end_col = ncol(data$tm), 
             group = data_tog$tm$`TM_atx_category`,
             strata = data_tog$tm$site)
# TM: not significant (but very close- results may be influenced by Salmon?)
set.seed(1)
runPERMANOVA(data_tog$tac, start_col, end_col = ncol(data$tac), 
             group = data_tog$tac$`TAC_atx_category`,
             strata = data_tog$tac$site)
# TAC: not significant, could be do to separate scaling equalizing both SFE-M and RUS
set.seed(1)
runPERMANOVA(data_tog$nt, start_col, end_col = ncol(data$nt), 
             group = data_tog$nt$`NT_atx_category`, 
             strata = data_tog$nt$site[-which(is.na(data_tog$nt$mean_ATX_all_ug_orgmat_g))],
             na.action = "na.omit")
# NT: not significant (same note as above)

## (b) rivers separately
# (omitting Salmon because only one sample was toxic and very minorly)
set.seed(1)
runPERMANOVA(data_sep_list$tm$`SFE-M`, start_col, end_col = ncol(data$tm),
             group = data_sep_list$tm$`SFE-M`$TM_atx_category)
anova(betadisper(vegdist(data_sep_list$tm$`SFE-M`[,start_col:ncol(data$tm)], method = "bray"), 
                 data_sep_list$tm$`SFE-M`$TM_atx_category))
# SFE TM: significant: **, dispersion is not!
set.seed(1)
runPERMANOVA(data_sep_list$tac$`SFE-M`, start_col, end_col = ncol(data$tac),
             group = data_sep_list$tac$`SFE-M`$TAC_atx_category)
anova(betadisper(vegdist(data_sep_list$tac$`SFE-M`[,start_col:ncol(data$tac)], method = "bray"), 
                 data_sep_list$tac$`SFE-M`$TAC_atx_category))
# SFE TAC: significant *, dispersion is not!
set.seed(1)
runPERMANOVA(data_sep_list$tac$`RUS`, start_col, end_col = ncol(data$tac),
             group = data_sep_list$tac$`RUS`$TAC_atx_category)
anova(betadisper(vegdist(data_sep_list$tac$`RUS`[,start_col:ncol(data$tac)], method = "bray"), 
                 data_sep_list$tac$`RUS`$TAC_atx_category))
# RUS TAC: significant *, dispersion is not!

lapply(data_sep_list$nt, function(x) runPERMANOVA(x, start_col, end_col = ncol(data$nt),
                                                  group = x$`NT_atx_category`,
                                                  na.action = "na.omit"))
# not for any river, especially Russian and Salmon

## Lastly, what about with a binary?

# (a) just doing single within rivers
# (omitting Salmon because only one sample was toxic and very minorly)
set.seed(1)
runPERMANOVA(data_sep_list$tm$`SFE-M`, start_col, end_col = ncol(data$tm),
             group = data_sep_list$tm$`SFE-M`$TM_atx_detected)
anova(betadisper(vegdist(data_sep_list$tm$`SFE-M`[,start_col:ncol(data$tm)], method = "bray"), 
                 data_sep_list$tm$`SFE-M`$TM_atx_detected))
# SFE TM: significant: ***, dispersion is not!
set.seed(1)
runPERMANOVA(data_sep_list$tac$`SFE-M`, start_col, end_col = ncol(data$tac),
             group = data_sep_list$tac$`SFE-M`$TAC_atx_detected)
anova(betadisper(vegdist(data_sep_list$tac$`SFE-M`[,start_col:ncol(data$tac)], method = "bray"), 
                 data_sep_list$tac$`SFE-M`$TAC_atx_detected))
# SFE TAC: significant *, dispersion is not!
set.seed(1)
runPERMANOVA(data_sep_list$tac$`RUS`, start_col, end_col = ncol(data$tac),
             group = data_sep_list$tac$`RUS`$TAC_atx_detected)
anova(betadisper(vegdist(data_sep_list$tac$`RUS`[,start_col:ncol(data$tac)], method = "bray"), 
                 data_sep_list$tac$`RUS`$TAC_atx_detected))
# RUS TAC: significant *, dispersion is not!

#### (4) NMDS Visualization ####

# South Fork Eel Microcoleus
makeNMDSplot(getNMDSdata(data_sep_list$tm$`SFE-M`, start_col, end_col = ncol(data$tm)),
             color = "TM_atx_detected", shape = "TM_atx_category", loading = FALSE)

# Russian Anabaena
makeNMDSplot(getNMDSdata(data_sep_list$tac$`RUS`, start_col, end_col = ncol(data$tac)),
             color = "TAC_atx_detected", shape = "TAC_atx_category", loading = FALSE)

# Russian Anabaena
makeNMDSplot(getNMDSdata(data_sep_list$tac$`SFE-M`, start_col, end_col = ncol(data$tac)),
             color = "TAC_atx_detected", shape = "TAC_atx_category", loading = FALSE)
# difference between those with toxins detected and those without?

#### (5) (db)RDA analyses ####

## (a) variable selection

# check which variables are correlated
correlations <- cor(env_tog[,10:ncol(env_tog)], use = "complete.obs")
view(which(correlations > 0.7 & correlations != 1, arr.ind = TRUE))
# correlated: DO & temp, assumed pH & DO, ammonium & nitrate & DIN, 
# conductivity & Cl Na K & Mg, DOC & Br amm DIN, SO4 & Cl NA
# we are missing some conductivity observations, so choose most correlated instead (Mg)
# will select one from each correlated group:
# included in RDA analyses: <INSERT>

## (b) analyzing rivers together with data standardized together

# run dbRDA for all sample types with all rivers together
dbRDA_tog <- lapply(sample_types, function(x) run_dbRDA(data_tog[[x]], 
                                                       start_col = start_col, 
                                                       end_col = ncol(data[[x]]),
                                                       mat_atx = x,
                                                       na.action = "na.omit"))
names(dbRDA_tog) <- sample_types

# plot
lapply(dbRDA_tog, function(x) plot(x$object))

plot(dbRDA_tog$tm$object, type = "n", scaling = 3)
points(dbRDA_tog$tm$object, display = "sites", pch=20, cex=0.7, col="gray32", scaling=3)

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
russian_TAC = run_dbRDA(data_sep_list$tac$`SFE-M`, 
                        start_col = 5, 
                        end_col = ncol(data$tac),
                        mat_atx = "tac")

# sfkeel TAC
sfkeel_TAC = run_dbRDA(data_sep_list$tac$`SFE-M`,
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
# russian tac yes, sfkeel tac yes, sfkeel tm no

# what is the rsquared?
lapply(individual_t_rdas, function(x) print(x$rsquared))
# adjusted: nt TBD, tac 0.265, tm 0.121

# significant variables?
lapply(individual_t_rdas, function(x) print(x$sig_variables))

#### (5) Species Indicator Analyses ####

<<<<<<< HEAD
# with low, medium, high categories
=======
## (a) within river, ATX groupings

# South Fork Eel Microcoleus
summary(multipatt(data_sep_list$tm$`SFE-M`[,start_col:ncol(data$tm)], 
                  data_sep_list$tm$`SFE-M`$TM_atx_category,
                  func = "r.g", control = how(nperm = 999)))
# uniquely assigned to "high"- Epithemia and Leptolyngbya
>>>>>>> aac9ab7456d7548273e268593b9919b453cf0ef2

# South Fork Eel Anabaena
summary(multipatt(data_sep_list$tac$`SFE-M`[,start_col:ncol(data$tac)], 
                  data_sep_list$tac$`SFE-M`$TAC_atx_category,
                  func = "r.g", control = how(nperm = 999)))
# nothing uniquely identified here

# Russian Anabaena
summary(multipatt(data_sep_list$tac$`RUS`[,start_col:ncol(data$tac)], 
                  data_sep_list$tac$`RUS`$TAC_atx_category,
                  func = "r.g", control = how(nperm = 999)))
# medium (which is highest for the Russian)- oscillatoria

# Non-Target within each River
lapply(data_sep_list$nt, function(x) {
       
  # remove rows with NA
  y <- x %>% 
    filter(!is.na(NT_atx_category))
  
  print(x$site[1])
  summary(multipatt(y[,start_col:ncol(data$nt)], y$NT_atx_category,
                    func = "r.g", control = how(nperm = 999)))})
# RUS: medium: phomidium*, none: scenedesmus*, low+none: leptolynbya*
# SFE: high: rhopalodia, epithemia, anabaena
# SAL: nothing (which to be expected, only one sample has low anatoxins)

## (b) within river ATX binary

# South Fork Eel Microcoleus
summary(multipatt(data_sep_list$tm$`SFE-M`[,start_col:ncol(data$tm)], 
                  data_sep_list$tm$`SFE-M`$TM_atx_detected,
                  func = "r.g", control = how(nperm = 999)))
# coccoids** and geilerinema* for detected samples, nostoc** for undetected

# South Fork Eel Anabaena
summary(multipatt(data_sep_list$tac$`SFE-M`[,start_col:ncol(data$tac)], 
                  data_sep_list$tac$`SFE-M`$TAC_atx_detected,
                  func = "r.g", control = how(nperm = 999)))
# nostoc* with none that's it

# Russian Anabaena
summary(multipatt(data_sep_list$tac$`RUS`[,start_col:ncol(data$tac)], 
                  data_sep_list$tac$`RUS`$TAC_atx_detected,
                  func = "r.g", control = how(nperm = 999)))
# nothing here either