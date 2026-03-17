#### Log-transforming anatoxins and making categories for analyses
### Jordan Zabrecky
## last edited: 03.12.2026

# This script standardizes anatoxin concentrations to use in Q3 analyses

#### (1) Loading libraries & data ####

# libraries
lapply(c("tidyverse"), require, character.only = T)

# read in environmental covariates & toxin data
atx <- read.csv("data/field_and_lab/environmental_covariates_and_toxins.csv") %>% 
  # will analyze congeners all together standardized by OM, 
  # so removing individual congeners and standardization by chl-a
  select(field_date, site_reach, site, TM_ATX_all_ug_orgmat_g, TAC_ATX_all_ug_orgmat_g)

#### (2) Log-transforming and deciding ATX categories ####

## (a) looking at non-logged data

# see medians, etc. for each
summary(atx$TM_ATX_all_ug_orgmat_g)
summary(atx$TAC_ATX_all_ug_orgmat_g)

# maybe should just use a cut-off of greater than >10 ug?
ggplot(data = atx, aes(y = TM_ATX_all_ug_orgmat_g)) +
  geom_boxplot()
ggplot(data = atx, aes(y = TAC_ATX_all_ug_orgmat_g)) +
  geom_boxplot()

## (b) log-transform values

# need to replace zeros with a value 
# lowest non-zero value is 0.03, how about 0.02 as zero replacement?
log(0.03) # -3.506558
log(0.02) # -3.912

# log values (replace with 0.02 when 0)
atx <- atx %>% 
  mutate(log_TM_ATX_all_ug_orgmat_g = case_when(TM_ATX_all_ug_orgmat_g == 0 ~ log(0.02),
                                                TRUE ~ log(TM_ATX_all_ug_orgmat_g)),
         log_TAC_ATX_all_ug_orgmat_g = case_when(TAC_ATX_all_ug_orgmat_g == 0 ~ log(0.02),
                                                 TRUE ~ log(TAC_ATX_all_ug_orgmat_g)))

# see medians, etc. for each for log values
summary(atx$log_TM_ATX_all_ug_orgmat_g)
summary(atx$log_TAC_ATX_all_ug_orgmat_g)

# maybe should just use a cut-off of greater than >10 ug?
ggplot(data = atx, aes(y = log_TM_ATX_all_ug_orgmat_g)) +
  geom_boxplot()
ggplot(data = atx, aes(y = log_TAC_ATX_all_ug_orgmat_g)) +
  geom_boxplot()

# need to make a column for NT samples to compare to toxin concentrations
# (in case we have both taxa or only one taxa) will take the average of two or the one value
atx <- atx %>% 
  mutate(mean_ATX_all_ug_orgmat_g = case_when(is.na(TM_ATX_all_ug_orgmat_g) & 
                                                is.na(TAC_ATX_all_ug_orgmat_g) ~ NA,
                                              is.na(TM_ATX_all_ug_orgmat_g) ~ TAC_ATX_all_ug_orgmat_g,
                                              is.na(TAC_ATX_all_ug_orgmat_g) ~ TM_ATX_all_ug_orgmat_g,
                                              TRUE ~ (TAC_ATX_all_ug_orgmat_g + TM_ATX_all_ug_orgmat_g) / 2)) %>% 
  # lastly, log that value as above
  mutate(log_NT_ATX_all_ug_orgmat_g = case_when(mean_ATX_all_ug_orgmat_g == 0 ~ log(0.02),
                                                TRUE ~ log(mean_ATX_all_ug_orgmat_g)))

# add in categorical grouping
atx <- atx %>% 
  mutate(TM_atx_category = case_when(TM_ATX_all_ug_orgmat_g >= 10 ~ "high",
                                     TM_ATX_all_ug_orgmat_g > 0 ~ "low",
                                     TM_ATX_all_ug_orgmat_g == 0 ~ "none",
                                     TRUE ~ NA),
         TAC_atx_category = case_when(TAC_ATX_all_ug_orgmat_g >= 10 ~ "high",
                                      TAC_ATX_all_ug_orgmat_g > 0 ~ "low",
                                      TAC_ATX_all_ug_orgmat_g == 0 ~ "none",
                                      TRUE ~ NA),
         NT_atx_category = case_when(mean_ATX_all_ug_orgmat_g >= 10 ~ "high",
                                     mean_ATX_all_ug_orgmat_g > 0 ~ "low",
                                     mean_ATX_all_ug_orgmat_g == 0 ~ "none",
                                     TRUE ~ NA),
         TM_atx_detected = case_when(TM_atx_category == "none" ~ "none",
                                     is.na(TM_atx_category) ~ NA,
                                     TRUE ~ "detected"),
         TAC_atx_detected = case_when(TAC_atx_category == "none" ~ "none",
                                      is.na(TAC_atx_category) ~ NA,
                                      TRUE ~ "detected"),
         NT_atx_detected = case_when(NT_atx_category == "none" ~ "none",
                                     is.na(NT_atx_category) ~ NA,
                                     TRUE ~ "detected"))

#### (3) Saving CSV ####

# save csv
write.csv(atx, "./data/field_and_lab/atx_w_categorical_groupings.csv", 
          row.names = FALSE)
