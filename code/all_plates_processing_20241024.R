## greengenes full sample analyses

# 10.24.2022

### loading libraries & data ###

library(tidyverse)
library(stringr)
library(lubridate)

# read in sample metadata
metadata <- read.csv("./data/molecular/16s_sample_metadata.csv")

# in this data, the sample number has already been reassigned to the plate
data <- read.csv("data/molecular/mat_counts_unfiltered_all.csv")

# separating out taxonomic ranks
data$phylum <- str_match(data$Taxon, "p__\\s*(.*?)\\s*; c__")[,2]
data$class <- str_match(data$Taxon, "c__\\s*(.*?)\\s*; o__")[,2]
data$order <- str_match(data$Taxon, "o__\\s*(.*?)\\s*; f__")[,2]
data$family <- str_match(data$Taxon, "f__\\s*(.*?)\\s*; g__")[,2]
data$genus <- str_match(data$Taxon, "g__\\s*(.*?)\\s*; s__")[,2]
data$species <- str_extract(data$Taxon, "(?<=s__).*")

# get rid of columns we no longer care about
data <- data %>% 
  dplyr::select(!X & !feature_id & !Taxon)

# get rid of counts = NA
data <- data %>% 
  filter(!is.na(counts))

## merge with metadata...
metadata <- metadata %>% 
  rename(Sample_name = vial_ID) %>% 
  mutate(site_reach_date = paste(site_reach, field_date, sep = " "))
all <- left_join(data, metadata, by = "Sample_name")
all$field_date <- mdy(all$field_date)
all <- arrange(all, field_date, increasing = TRUE) # order by date for later plots

#### looking at data!!

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
  geom_col(aes(x = site_reach_date, y = counts, fill = phylum)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# these should be reordered w/ time but too lazy to do that at present

RUS_NT_cyano <- rus_NT %>% 
  filter(phylum == "Cyanobacteria")

ggplot(RUS_NT_cyano, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# see if any of our "fake" target samples actually contained microcoleus
rus_fakes <- russian %>% 
  filter(sample_type == "TM") %>% 
  filter(fake_target == "y")

ggplot(rus_fakes, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = phylum)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

rus_fake_cyano <- rus_fakes %>% 
  filter(phylum == "Cyanobacteria")

ggplot(rus_fake_cyano, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# mostly synechococcus or something pink; two have lots of microcoleus actually
rus_fakes_micro <- rus_fakes %>% 
  filter(genus == "Microcoleus" & !is.na(counts))
unique(rus_fakes_micro$site_reach_date)

rus_fakes_metdata <- metadata %>% 
  filter(site == "RUS" & fake_target == "y")
# so 10/15 of our "fake" samples actually did contain microcoleus
# "fake" meaning when we could not find microcoleus (which never happened on the Russian)
# we took anything that looked like it could possibly be it (dark, in a riffle, etc.)
# to confirm presence/absence
# small presence has also been confirmed via microscopy

# target anabaena & cylindrospermum
rus_anacylin <- russian %>% 
  filter(sample_type == "TAC") %>% 
  filter(phylum == "Cyanobacteria")

ggplot(rus_anacylin, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = genus)) + 
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
  geom_col(aes(x = site_reach_date, y = counts, fill = phylum)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# seeing a lot more cyanobacteria here than at the Russian!
# again, feel like I'm missing samples

# TAC, anabaena & cylindrospermum
sfk_TAC <- sfkeel %>% 
  filter(sample_type == "TAC") %>% 
  filter(fake_target == "n") %>% # need to specify fake here bc they took "fake" a week I was gone
  filter(phylum == "Cyanobacteria")

# might actually be worth it to look at family level instead, since that distinguishes
# heterocysts

ggplot(sfk_TAC, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = genus)) + 
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
  filter(phylum == "Cyanobacteria")

ggplot(sfk_TM, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# mostly microcoleus!
# early season sample is more mixed which makes sense
# and then there is a sample that is NA
# honestly site 4 barely had any microcoleus growth that year

# let's look at the one sample that was a maybe????
sfk_TM_maybe <- sfk_TM %>% 
  filter(fake_target == "maybe")

ggplot(sfk_TM_maybe, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# mostly looks like anabaena samples so let's go with a no

## Lastly, the Salmon
salmon_NT <- salmon %>% 
  filter(sample_type == "NT")

ggplot(salmon_NT, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = phylum)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

salmon_TM <- salmon %>% 
  filter(sample_type == "TM") %>% 
  filter(phylum == "Cyanobacteria") 

ggplot(salmon_TM, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# mostly microcoleus! that last sample is a fake- I think it was mostly diatoms under the scope
# when we looked at it back at the Angelo Reserve
# purple SIO unknown????

# we had one TAC sample
salmon_TAC <- salmon %>% 
  filter(sample_type == "TAC") %>% 
  filter(phylum == "Cyanobacteria")

ggplot(salmon_TAC, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# mostly trichormus like other TACs

#### INVESTIGATING UNKNOWN ID SAMPLES

# either 1000 or 1001 is a fake RUS TM- seems like more likely 1001

sample_1000_108 <- all %>% # sample ID is 108 from plate but renamed to 1000 since there were two 108s
  filter(Sample_name == 1000)

ggplot(sample_1000_108, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = phylum)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# lots of cyanobacteria and proteobacteria; >60000 counts

sample_1000_108_cyano <- sample_1000_108 %>% 
  filter(phylum == "Cyanobacteria")

ggplot(sample_1000_108_cyano, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# mostly SIO2C1

sample_1001_108 <- all %>% 
  filter(Sample_name == 1001)

ggplot(sample_1001_108, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = phylum)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# more cyanobacteria and proteobacteria; >100,000 counts

sample_1001_108_cyano <- sample_1001_108 %>% 
  filter(phylum == "Cyanobacteria")

ggplot(sample_1001_108_cyano, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# mostly SIO2C1 (identified as Zarconia in cyanoseq- a member of oscillatoria)
# more than 60000 counts

sample_801 <- all %>% 
  filter(Sample_name == 801)

ggplot(sample_801, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = phylum)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# more cyanobacteria and proteobacteria; almost 30000 reads

sample_801_cyano <- sample_801 %>% 
  filter(phylum == "Cyanobacteria")

ggplot(sample_801_cyano, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# mostly SIO2C1

sample_552 <- all %>% 
  filter(Sample_name == 552) # sample not found

sample_554 <- all %>% 
  filter(Sample_name == 554)

ggplot(sample_554, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = phylum)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# mostly cyanobacteria; 15000 reads

sample_554_cyano <- sample_554 %>% 
  filter(phylum == "Cyanobacteria")

ggplot(sample_554_cyano, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# mostly tychonema which makes me think it is a TM or microcoleus rock
# (probably meaghan's microcoleus rock)

sample_931 <- all %>% 
  filter(Sample_name == 931)

ggplot(sample_931, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = phylum)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# mostly cyanobacteria; not as diverse (i.e. less total phyla)

sample_931_cyano <- sample_931 %>% 
  filter(phylum == "Cyanobacteria")

ggplot(sample_931_cyano, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# okay this is almost certainly the microcoleus rock from the other study!!

#### NEED TO DEAL WITH TRIPLICATES AND BLANKS