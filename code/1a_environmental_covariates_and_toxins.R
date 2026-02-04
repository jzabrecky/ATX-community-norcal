#### Gathering environmental covariates for each field date at each reach
### Jordan Zabrecky
## last edited: 02.03.2026

# This script creates a dataframe with environmental covariates (water chemistry
# and anatoxin data) to use with both microscopy and molecular data

# For code processing of water chemistry and target sample anatoxin 
# data from original EDI csv's, see
# https://github.com/jzabrecky/ATX-synchrony-norcal/blob/main/code/2b_water_chemistry.R
# https://github.com/jzabrecky/ATX-synchrony-norcal/blob/main/code/2e_target_sample_anatoxins.R

#### (1) Loading libraries & data ####

# load tidyverse library
library("tidyverse")

# read in processed anatoxin & water chemistry data
atx_target <- read.csv("./data/field_and_lab/cyano_atx.csv") # processed version
water_chemistry <- read.csv("./data/field_and_lab/water_chemistry.csv") %>% 
  select(!c(time, reach))

#### (2) Calculate DIN ####

# calculate DIN (dissolved inorganic nitrogen) as the sum of nitrate and ammonium
water_chemistry <- water_chemistry %>% 
  mutate(DIN_mg_N_L = ammonium_mg_N_L + nitrate_mg_N_L) %>% 
  relocate(DIN_mg_N_L, .after = ammonium_mg_N_L)

#### (2) Joining data and saving csv ####

# TM sample that was taken on 9/8 was to make up for sample on 9/6 that was taken incorrectly
atx_target$field_date[which(atx_target$field_date == "2022-09-08")] <- "2022-09-06"

# add prefix of sample_type to column names to make dataframe wider 
# (i.e., one row per reach visit on a day)
atx_target_wider <- atx_target %>% 
  pivot_wider(names_from = sample_type, values_from = c(4:ncol(atx_target)), 
              names_glue = "{sample_type}_{.value}")

# left join to water chemistry
data <- left_join(water_chemistry, atx_target_wider, by = c("field_date", "site_reach"))

# save csv
write.csv(data, "./data/field_and_lab/environmental_covariates_and_toxins.csv",
          row.names = FALSE)
