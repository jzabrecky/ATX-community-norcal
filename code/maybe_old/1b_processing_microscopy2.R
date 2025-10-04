#### Processing microscopy data further (removing target taxa)
### Jordan Zabrecky
## last edited 05.09.2025

## This code processes target microscopy data from processed from the previous
## code by removing the target taxa (i.e. Microcoleus or Anabaena/Cylindrospermum)
## and recalculates the relative abundance of other taxa

#### (1) Loading libraries ####

# loading libraries
lapply(c("tidyverse", "lubridate"), require, character.only = T)

# read in processed microscopy data
tm <- read.csv("./data/morphological/tm_algalonly_with_covar.csv")
tac <- read.csv("./data/morphological/tac_algalonly_with_covar.csv")

# saving covariate data to add in later
covariates <- colnames(tm[,27:ncol(tm)])
tm_covariates <- tm %>% 
  select(site_reach, field_date, covariates)
tac_covariates <- tac %>% 
  select(site_reach, field_date, covariates)

#### (2) Removing target taxa from each sample ####

## (a) TM // target microcoleus samples

# remove environmental covariate data and microcoleus abundance
tm_temp <- tm %>% 
  select(!c(microcoleus, covariates))

# calculate total abundances w/o microcoleus
tm_temp$total_community_percent <- rowSums(tm_temp[c(5:ncol(tm_temp))])

# recalculate relative abundances w/o microcoleus
tm_processed <- tm_temp %>% 
  relocate(total_community_percent, .after = sample_type) %>%  # quick reorganization
  pivot_longer(cols = c(6:ncol(tm_temp)), values_to = "old_percent", names_to = "taxa") %>% 
  mutate(new_percent = (old_percent / total_community_percent) * 100) %>% 
  select(!c(old_percent, total_community_percent)) %>% 
  pivot_wider(names_from = "taxa", values_from = "new_percent")

# double-check that recalculated
check_tm <- which(rowSums(tm_processed[5:ncol(tm_processed)]) != 100)

# check these values- likely minor decimal issues 
# as computers are bad at dealing with small numbers
rowSums(tm_processed[5:ncol(tm_processed)])[check_tm] # yup, all 100!

## (b) anabaena/cylindrospermum

# remove environmental covariate data and microcoleus abundance
tac_temp <- tac %>% 
  select(!c(anabaena_and_cylindrospermum, covariates))

# calculate total abundances w/o microcoleus
tac_temp$total_community_percent <- rowSums(tac_temp[c(5:ncol(tac_temp))])

# recalculate relative abundances w/o microcoleus
tac_processed <- tac_temp %>% 
  relocate(total_community_percent, .after = sample_type) %>%  # quick reorganization
  pivot_longer(cols = c(6:ncol(tac_temp)), values_to = "old_percent", names_to = "taxa") %>% 
  mutate(new_percent = (old_percent / total_community_percent) * 100) %>% 
  select(!c(old_percent, total_community_percent)) %>% 
  pivot_wider(names_from = "taxa", values_from = "new_percent")

# double-check that recalculated
check_tac <- which(rowSums(tac_processed[5:ncol(tac_processed)]) != 100)

# check these values- likely minor decimal issues 
# as computers are bad at dealing with small numbers
rowSums(tac_processed[5:ncol(tac_processed)])[check_tac] # yup, all 100!

#### (3) Add back in environmental covariate data and save ####

# join back in environmental covariate data
tm_final <- left_join(tm_processed, tm_covariates, by = c("site_reach", "field_date"))
tac_final <- left_join(tac_processed, tac_covariates, by = c("site_reach", "field_date"))

# save!
write.csv(tm_final, "./data/morphological/tm_algalonly_nomicro_with_covar.csv", 
          row.names = FALSE)
write.csv(tac_final, "./data/morphological/tac_algalonly_noanacyl_with_covar.csv", 
          row.names = FALSE)
