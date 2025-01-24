#### cyanoseq results


library(tidyverse)

# get rid of any potential dplyr masking
filter <- dplyr::filter

### looking at cyanoseq results // silva
cyanoseq <- read.csv("./data/molecular/Mat_counts_taxonomy_cyanoseq_clean.csv")

# metadata
metadata <- read.csv("./data/molecular/16s_sample_metadata.csv")
# metadata <- read.csv("metadata.csv")

## merge with metadata...
metadata <- metadata %>% 
  dplyr::rename(Sample_name = vial_ID) %>% 
  mutate(site_reach_date = paste(site_reach, field_date, sep = " "))
all <- left_join(cyanoseq, metadata, by = "Sample_name")
all$field_date <- mdy(all$field_date)
all <- arrange(all, field_date, increasing = TRUE) # order by date for later plots

### looking at results ####
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
# lots of ... nostoc?

# TM, microcoleus
sfk_TM <- sfkeel %>% 
  filter(sample_type == "TM") %>% 
  filter(Phylum == "Cyanobacteria")

ggplot(sfk_TM, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = Genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# defintely seems to be reading microcoleus as tychonema
# keith said he had an issue and would look at the other hits

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
  filter(Phylum == "Cyanobacteriota") 

ggplot(salmon_TM, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = Genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# again tychnema issue

# we had one TAC sample
salmon_TAC <- salmon %>% 
  filter(sample_type == "TAC") %>% 
  filter(Phylum == "Cyanobacteriota")

ggplot(salmon_TAC, aes(x = site_reach_date, y = counts)) + 
  geom_col(aes(x = site_reach_date, y = counts, fill = Genus)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# mostly trichormus like other TACs

