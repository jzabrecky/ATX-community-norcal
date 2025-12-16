#### Processing microscopy data further (removing target taxa)
### Jordan Zabrecky
## last edited 12.15.2025

## This code processes target microscopy data from processed from the previous
## code by removing the target taxa (i.e. Microcoleus or Anabaena/Cylindrospermum)
## and recalculates the relative abundance of other taxa. This is also done
## again for TAC samples removing green algae (the substrate)

#### (1) Loading libraries ####

# loading libraries
lapply(c("tidyverse", "lubridate"), require, character.only = T)

# read in processed microscopy data
tm <- read.csv("./data/morphological/tm_algalonly.csv")
tac <- read.csv("./data/morphological/tac_algalonly.csv")

#### (2) Removing target taxa from each sample ####

## (a) TM // target microcoleus samples

# remove  microcoleus abundance
tm_temp <- tm %>% 
  select(!microcoleus)

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

# remove anabaena/cylindrospermum abundance
tac_temp <- tac %>% 
  select(!anabaena_and_cylindrospermum)

# calculate total abundances w/o anabaena/cylindrospermum
tac_temp$total_community_percent <- rowSums(tac_temp[c(5:ncol(tac_temp))])

# recalculate relative abundances w/o anabaena/cylindrospermum
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

## (c) green algae in tac samples

# remove green algae in addition to the above 
tac_temp2 <- tac_processed %>% 
  select(!green_algae)

# calculate total abundances w/o green algae
tac_temp2$total_community_percent <- rowSums(tac_temp2[c(5:ncol(tac_temp2))])

# recalculate relative abundances w/o green algae
tac_processed2 <- tac_temp2 %>% 
  relocate(total_community_percent, .after = sample_type) %>%  # quick reorganization
  pivot_longer(cols = c(6:ncol(tac_temp2)), values_to = "old_percent", names_to = "taxa") %>% 
  mutate(new_percent = (old_percent / total_community_percent) * 100) %>% 
  select(!c(old_percent, total_community_percent)) %>% 
  pivot_wider(names_from = "taxa", values_from = "new_percent")

# double-check that recalculated
check_tac2 <- which(rowSums(tac_processed2[5:ncol(tac_processed2)]) != 100)

# check these values- likely minor decimal issues 
# as computers are bad at dealing with small numbers
rowSums(tac_processed2[5:ncol(tac_processed2)])[check_tac2] # yup, all 100!

#### (3) Save CSV's ####

# save!
write.csv(tm_processed, "./data/morphological/tm_algalonly_nomicro.csv", 
          row.names = FALSE)
write.csv(tac_processed, "./data/morphological/tac_algalonly_noanacyl.csv", 
          row.names = FALSE)
write.csv(tac_processed2, "./data/morphological/tac_algalonly_noanacylgreenalgae.csv",
          row.names = FALSE)
