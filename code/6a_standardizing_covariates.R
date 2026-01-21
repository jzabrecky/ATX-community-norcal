#### Standardizing anatoxins & environmental covariates for analyses
### Jordan Zabrecky
## last edited: 01.20.2026

# This script standardizes anatoxin concentrations and environmental covariates
# to use in analyses such as an RDA

#### (1) Loading libraries & data ####

# libraries
lapply(c("tidyverse"), require, character.only = T)
       
# read in environmental covariates & toxin data
env <- read.csv("data/field_and_lab/environmental_covariates_and_toxins.csv") %>% 
  # will analyze congeners all together standardized by OM, 
  # so removing individual congeners and standardization by chl-a
  # also removing pH as it has a lot of NAs for 2022
  select(!c("pH", "TM_ATXa_ug_g", "TAC_ATXa_ug_g", "TM_dhATXa_ug_g", 
            "TAC_dhATXa_ug_g", "TM_HTXa_ug_g", "TAC_HTXa_ug_g", "TM_Chla_ug_g",
            "TAC_Chla_ug_g", "TM_Pheo_ug_g", "TAC_Pheo_ug_g", "TM_percent_organic_matter",
            "TAC_percent_organic_matter", "TM_ATX_all_ug_chla_ug", "TAC_ATX_all_ug_chla_ug",
            "TM_ATX_all_ug_g", "TAC_ATX_all_ug_g")) %>% 
  filter(year(ymd(field_date)) == 2022) # focusing on 2022 data for this analysis

#### (4) Come up with ATX groupings ####

# to use species indicator analysis, we need categorical groupings
# >50 ug is high
# 5-10 ug is medium
# <5 is low
# 0 is none

# need to make a column for NT samples to compare to toxin concentrations
# (in case we have both taxa or only one taxa) will take the average of two or the one value
env <- env %>% 
  mutate(mean_ATX_all_ug_orgmat_g = case_when(is.na(TM_ATX_all_ug_orgmat_g) & 
                                                is.na(TAC_ATX_all_ug_orgmat_g) ~ NA,
                                              is.na(TM_ATX_all_ug_orgmat_g) ~ TAC_ATX_all_ug_orgmat_g,
                                              is.na(TAC_ATX_all_ug_orgmat_g) ~ TM_ATX_all_ug_orgmat_g,
                                              TRUE ~ (TAC_ATX_all_ug_orgmat_g + TM_ATX_all_ug_orgmat_g) / 2))

# add in categorical grouping
env <- env %>% 
  mutate(TM_atx_category = case_when(TM_ATX_all_ug_orgmat_g >= 10 ~ "high",
                                     TM_ATX_all_ug_orgmat_g < 10 & TM_ATX_all_ug_orgmat_g >= 1 ~ "medium",
                                     TM_ATX_all_ug_orgmat_g <= 1 & TM_ATX_all_ug_orgmat_g > 0 ~ "low",
                                     TM_ATX_all_ug_orgmat_g == 0 ~ "none",
                                     TRUE ~ NA),
         TAC_atx_category = case_when(TAC_ATX_all_ug_orgmat_g >= 10 ~ "high",
                                      TAC_ATX_all_ug_orgmat_g < 10 & TAC_ATX_all_ug_orgmat_g >= 1 ~ "medium",
                                      TAC_ATX_all_ug_orgmat_g <= 1 & TAC_ATX_all_ug_orgmat_g > 0 ~ "low",
                                      TAC_ATX_all_ug_orgmat_g == 0 ~ "none",
                                      TRUE ~ NA),
         NT_atx_category = case_when(mean_ATX_all_ug_orgmat_g >= 10 ~ "high",
                                     mean_ATX_all_ug_orgmat_g < 10 & mean_ATX_all_ug_orgmat_g >= 1 ~ "medium",
                                     mean_ATX_all_ug_orgmat_g <= 1 & mean_ATX_all_ug_orgmat_g > 0 ~ "low",
                                     mean_ATX_all_ug_orgmat_g == 0 ~ "none",
                                     TRUE ~ NA)) %>% 
  relocate(TM_atx_category, .before = temp_C) %>% 
  relocate(TAC_atx_category, .before = temp_C) %>% 
  relocate(NT_atx_category, .before = temp_C)

#### (3) Log-Transforming ATX data ####

# anatoxins have really high outliers (as below!), so will log transform before scaling
hist(env$TM_ATX_all_ug_orgmat_g, breaks = 100)
hist(env$TAC_ATX_all_ug_orgmat_g, breaks = 100)

# DO NOT replace NAs because in those circumstances, no mat was present for us to sample
# first look at histogram

# need to replace zeros with a value 
# lowest non-zero value is 0.03, how about 0.02 as zero replacement?
log(0.03) # -3.506558
log(0.02) # -3.912

# replace zeros
ATX_columns <- c("TM_ATX_all_ug_orgmat_g", "TAC_ATX_all_ug_orgmat_g", 
                 "mean_ATX_all_ug_orgmat_g")
for(i in 1:length(ATX_columns)) {
  env[ATX_columns[i]] <- replace(env[ATX_columns[i]], 
                                 env[ATX_columns[i]] == 0, 0.02)
}

# log-transform ATX variables 
for(i in 1:length(ATX_columns)) {
  env[ATX_columns[i]] <- log(env[ATX_columns[i]])
}

# it's okay looking- probably little room for improvement...?
hist(env$TM_ATX_all_ug_orgmat_g, breaks = 100)
hist(env$TAC_ATX_all_ug_orgmat_g, breaks = 100)

#### (3) Center & Scale Covariates ####

# standardize data all-together
env_tog <- env
env_tog[,8:ncol(env_tog)] <- apply(env_tog[,8:ncol(env_tog)], 2, scale)

# standardize data within each river separately
env_sep <- env
env_sep[,8:ncol(env_sep)] <- apply(env_sep[,8:ncol(env_sep)], 2,
                                   function(x) ave(x, env_sep$site, FUN = scale))

# save csvs
write.csv(env_tog, "./data/field_and_lab/env_tox_standardized_together.csv", 
          row.names = FALSE)
write.csv(env_sep, "./data/field_and_lab/env_tox_standardized_byriver.csv", 
          row.names = FALSE)
