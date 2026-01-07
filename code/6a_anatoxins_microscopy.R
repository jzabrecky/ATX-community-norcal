#### Comparing microscopy data with regard to anatoxin concentrations
### Jordan Zabrecky
## last edited: 01.06.2026

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
lapply(c("tidyverse", "plyr", "vegan", "cowplot", "TITAN2"), require, character.only = T)

# read in relative abundance files (data transformed in previous script, "4a_amongrivers_microscopy.R")
data <- lapply(list.files(path = "./data/morphological/transformed/", pattern = ".csv"),
               function(x) read.csv(paste("./data/morphological/transformed/", x, sep = "")))
names(data) <- c("nt", "tac", "tm")

# read in environmental covariates & toxin data
env <- read.csv("data/field_and_lab/environmental_covariates_and_toxins.csv") %>% 
  # analyze congeners all together, so removed separted congeners
  select(!c("pH", "TM_ATXa_ug_g", "TAC_ATXa_ug_g", "TM_dhATXa_ug_g", 
            "TAC_dhATXa_ug_g", "TM_HTXa_ug_g", "TAC_HTXa_ug_g", "TM_Chla_ug_g",
            "TAC_Chla_ug_g", "TM_Pheo_ug_g", "TAC_Pheo_ug_g", "TM_percent_organic_matter",
            "TAC_percent_organic_matter", "TM_ATX_all_ug_chla_ug", "TAC_ATX_all_ug_chla_ug",
            "TM_ATX_all_ug_g", "TAC_ATX_all_ug_g")) %>% 
  filter(year(ymd(field_date)) == 2022)

#### (2) Predictor Data Transformation ####

# replace 0's with small value for log-transformation (to avoid outweighing outliers)
# DO NOT replace NAs because in those circumstances, no mat was present for us to sample
# first look at histogram
hist(env$TM_ATX_all_ug_orgmat_g, breaks = 30)
hist(env$TAC_ATX_all_ug_orgmat_g, breaks = 30)

# lowest non-zero value is 0.03, how about 0.005 as zero replacement?
log(0.03) # -3.506558
log(0.01) # -4.60517

# replace zeros
env$TM_ATX_all_ug_orgmat_g <- replace(env$TM_ATX_all_ug_orgmat_g , 
                                      env$TM_ATX_all_ug_orgmat_g == 0, 0.01)
env$TAC_ATX_all_ug_orgmat_g <- replace(env$TAC_ATX_all_ug_orgmat_g , 
                                       env$TAC_ATX_all_ug_orgmat_g == 0, 0.01)

# log-transform ATX variables 
env[,4:ncol(env)] <- log(env[,4:ncol(env)])

# center & scale
# data all-together (standardized across rivers)
env_tog <- env
env_tog[,4:ncol(env_tog)] <- apply(env_tog[,4:ncol(env_tog)], 2, scale)
# rivers separately (group-level variable)
env_sep <- env
env_sep[,4:ncol(env_sep)] <- apply(env_sep[,4:ncol(env_sep)], 2,
                                   function(x) ave(x, env_sep$site, FUN = scale))

# compare histograms of the two
hist(env_tog$TM_ATX_all_ug_orgmat_g)
hist(env_sep$TM_ATX_all_ug_orgmat_g)
hist(env_tog$TAC_ATX_all_ug_orgmat_g)
hist(env_sep$TAC_ATX_all_ug_orgmat_g)

# save to avoid doing this in future scripts (run once)
#write.csv(env_tog, "./data/field_and_lab/env_tox_standardized_together.csv", 
#          row.names = FALSE)
#write.csv(env_sep, "./data/field_and_lab/env_tox_standardized_byriver.csv", 
#          row.names = FALSE)

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
