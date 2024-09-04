#### Assigning sample names and types to sequence reads
### 9.4.2024

# NOTE: our "unknown" IDs from the plates are 108*, 801, 552, 554, and 931 
# *there were 2 108s divided into 1000 and 1001
# our "unknown" IDs from the sample metadata are 104, 115, 208, ??, ??
# ?? are probably numbers > 500 and are unafiliated with this project

# characteristics of unknown sample IDs:
# 104 - Russian River, non-targeted/community composition (8/17/2022)
# 115 - Russian River, targeted-anabaena/cylindrospermum (9/1/2022)
# 208 - Russian River, non-targeted/community composition (9/15/2022)
# ??? - Microcoleus rock from South Fork Eel (for another project) (SFE-M-1S, 9/8/22)
# ??? - Control rock (green algae) from South Fork Eel (for another project) (SFE-M-1S, 9/8/22)

# the known 108 is a RUS-3 fake TM from 8/17/2022

# also have sample 27 which we were unsure if it was Microcoleus or not
# (though I don't think it has a corresponding anatoxin concentration measurement bc sample was so small anyways)

library(tidyverse)
library(lubridate)

## Read in data

# taxonomy data
taxonomy <- read.csv("./data/molecular/Bac_counts_taxonomy_cleanedup.csv")
# taxonomy <- read.csv("Bac_counts_taxonomy_cleanedup.csv")

# mutation from other code...
Bac_counts_taxonomy<- taxonomy %>%
  mutate(Phylum = gsub('_.*','',Phylum))

Bac_counts_taxonomy<- Bac_counts_taxonomy %>%
  mutate(Genus = gsub('_.*','',Genus))

Bac_counts_taxonomy<- Bac_counts_taxonomy %>%
  mutate(Genus = gsub('-.*','',Genus))

# metadata
metadata <- read.csv("./data/molecular/16s_sample_metadata.csv")
# metadata <- read.csv("metadata.csv")

## merge with metadata...
metadata <- metadata %>% 
  rename(Sample_name = vial_ID) %>% 
  mutate(site_reach_date = paste(site_reach, field_date, sep = " "))
all <- left_join(Bac_counts_taxonomy, metadata, by = "Sample_name")
all$field_date <- mdy(all$field_date)
all <- arrange(all, field_date, increasing = TRUE) # order by date for later plots

# separating into groups based on river
sfkeel <- all %>% 
  filter(site == "SFE-M") %>% 
  filter(sample_type != "blank") # do not want to include blanks for now

russian <- all %>% 
  filter(site == "RUS") %>% 
  filter(sample_type != "blank")

salmon <- all %>% 
  filter(site == "SAL") %>% 
  filter(sample_type != "blank")

# thoughts on what to do with blanks

## russian river
rus_NT <- russian %>% 
  filter(sample_type == "NT")

ggplot(rus_NT, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = Phylum)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# these should be reordered w/ time but too lazy to do that at present

# see if any of our "fake" target samples actually contained microcoleus
rus_fakes <- russian %>% 
  filter(sample_type == "TM") %>% 
  filter(fake_target == "y") %>% 
  filter(Phylum == "Cyanobacteria")

ggplot(rus_fakes, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = Genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# mostly synechococcus or something pink
rus_fakes_micro <- rus_fakes %>% 
  filter(Genus == "Microcoleus")
unique(rus_fakes_micro$site_reach_date)

rus_fakes_metdata <- metadata %>% 
  filter(site == "RUS" & fake_target == "y")
# so 6/15 of our "fake" samples actually did contain microcoleus
# "fake" meaning when we could not find microcoleus (which never happened on the Russian)
# we took anything that looked like it could possibly be it (dark, in a riffle, etc.)
# to confirm presence/absence
# small presence has also been confirmed via microscopy

# target anabaena & cylindrospermum
rus_anacylin <- russian %>% 
  filter(sample_type == "TAC") %>% 
  filter(Phylum == "Cyanobacteria")

ggplot(rus_anacylin, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = Genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# feel like I am missing samples here???
# see note about Trichormus below... (is Anabaena essentially)
# interesting the one sample with a lot of blue!
# also lots of NA counts so maybe CyanoSeq is a valuable option
### INVESTIGATE THIS ####

## south fork eel 

# NT/ reach community composition
sfk_NT <- sfkeel %>% 
  filter(sample_type == "NT")

ggplot(sfk_NT, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = Phylum)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# seeing a lot more cyanobacteria here than at the Russian!
# again, feel like I'm missing samples

# TAC, anabaena & cylindrospermum
sfk_TAC <- sfkeel %>% 
  filter(sample_type == "TAC") %>% 
  filter(fake_target == "n") %>% # need to specify fake here bc they took "fake" a week I was gone
  filter(Phylum == "Cyanobacteria")

ggplot(sfk_TAC, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = Genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# lots of synnechoccus and trichormus which I don't think we saw a ton of
# in our microscopy data (saw some synnechococcus but not a ton)
# Trichormus looks very similar to anabaena & cylindrospermum under the scope
# OKAY- looked at my emails- Anabaena was recently separated into 
# Dolichiospermum (planktonic) and Trichormus (benthic) so for our purposes
# those will be grouped in with Anabaena

# TM, microcoleus
sfk_TM <- sfkeel %>% 
  filter(sample_type == "TM") %>% 
  filter(Phylum == "Cyanobacteria")

ggplot(sfk_TM, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = Genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# mostly microcoleus!
# early season sample is more mixed which makes sense
# and then there is a sample that is NA
# honestly site 4 barely had any microcoleus growth that year

# let's look at the one sample that was a maybe????
sfk_TM_maybe <- sfk_TM %>% 
  filter(fake_target == "maybe")
# okay this is one of the mysteriously missing samples

## Lastly, the Salmon
salmon_NT <- salmon %>% 
  filter(sample_type == "NT")

ggplot(salmon_NT, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = Phylum)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

salmon_TM <- salmon %>% 
  filter(sample_type == "TM") %>% 
  filter(Phylum == "Cyanobacteria") 

ggplot(salmon_TM, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = Genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# mostly microcoleus! that last sample is a fake- I think it was mostly diatoms under the scope
# when we looked at it back at the Angelo Reserve
# purple SIO unknown????

# we had one TAC sample
salmon_TAC <- salmon %>% 
  filter(sample_type == "TAC") %>% 
  filter(Phylum == "Cyanobacteria")

ggplot(salmon_TAC, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = Genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# mostly trichormus like other TACs

#### INVESTIGATING UNKNOWN ID SAMPLES

# either 1000 or 1001 is a fake RUS TM- seems like more likely 1001

sample_1000_108 <- all %>% # sample ID is 108 from plate but renamed to 1000 since there were two 108s
  filter(Sample_name == 1000)

ggplot(sample_1000_108, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = Phylum)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# lots of cyanobacteria and proteobacteria; >60000 counts

sample_1000_108_cyano <- sample_1000_108 %>% 
  filter(Phylum == "Cyanobacteria")

ggplot(sample_1000_108_cyano, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = Genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# mostly synechoccus or trichormus???

sample_1001_108 <- all %>% 
  filter(Sample_name == 1001)

ggplot(sample_1001_108, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = Phylum)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# more cyanobacteria and proteobacteria; >100,000 counts

sample_1000_108_cyano <- sample_1001_108 %>% 
  filter(Phylum == "Cyanobacteria")

ggplot(sample_1000_108_cyano, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = Genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# mostly SIO2C1

sample_801 <- all %>% 
  filter(Sample_name == 801)

ggplot(sample_801, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = Phylum)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# more cyanobacteria and proteobacteria

sample_801_cyano <- sample_801 %>% 
  filter(Phylum == "Cyanobacteria")

ggplot(sample_801_cyano, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = Genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# mostly trichormus??

sample_552 <- all %>% 
  filter(Sample_name == 552) # sample not found

sample_554 <- all %>% 
  filter(Sample_name == 554) # sample not found

sample_931 <- all %>% 
  filter(Sample_name == 931)

ggplot(sample_931, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = Phylum)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# mostly cyanobacteria; not as diverse (i.e. less total phyla)

sample_931_cyano <- sample_931 %>% 
  filter(Phylum == "Cyanobacteria")

ggplot(sample_931_cyano, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = Genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# okay this is almost certainly the microcoleus rock from the other study!!

# need to look at ther samples (552 and 554) that are missing to confirm

#### SEEMINGLY MISSING SAMPLES THAT SHOULD NOT BE MISSING
all_sample_list <- as.data.frame(unique(all$Sample_name))
colnames(all_sample_list) <- c("vial_ID")
all_sample_list$test <- "yes"

# testing each plate to see if we are missing a plate
plate_1 <- read.table("./data/molecular/plate1_metadata.txt")
colnames(plate_1) <- c("plate_ID", "vial_ID")
test1 <- left_join(plate_1, all_sample_list, by = "vial_ID")
# only NA for the two IDs that both have 108 so not an issue!
# need to resolve above 1000 and 1001! (which is hopefully possible)

plate_2 <- read.table("./data/molecular/plate2_metadata.txt")
colnames(plate_2) <- c("plate_ID", "vial_ID")
test2 <- left_join(plate_2, all_sample_list, by = "vial_ID")
# this plate appears to be missing!!!

plate_3 <- read.table("./data/molecular/plate3_metadata.txt")
colnames(plate_3) <- c("plate_ID", "vial_ID")
test3 <- left_join(plate_3, all_sample_list, by = "vial_ID")
# vial ID 34 and 18 missing???

#### NEED TO DEAL WITH TRIPLICATES AND BLANKS