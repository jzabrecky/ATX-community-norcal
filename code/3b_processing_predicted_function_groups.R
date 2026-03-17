#### Processing predicted functional groups for samples
### Jordan Zabrecky
## last edited: 03.11.2026

# This script processes the predicted functional groups from PICRUSt2-SC
# standard pipeline, focusing on KEGG orthologs
# Note: This will only be done for rarefied (95% threshold) & filtered .biom data

# TO DO: rerun with 95 rarefied and filtered data! ran with merged FASTA file!!

#### (1) Loading libraries & data ####

# loading libraries
lapply(c("tidyverse", "ggpicrust2"), require, character.only = T)

## (a) getting metadata

# get list of file names
plate_files <- list.files(path = "./data/molecular/metadata/", pattern = "plate")

# create empty list (different plates have same IDs so want to keep them separate!)
plate_data <- list()

# fill in list with dataframes
for(i in 1:length(plate_files)) {
  plate_data[[i]] <- read.table(paste("data/molecular/metadata/", plate_files[i], sep = ""))
  colnames(plate_data[[i]]) <- c("plate_ID", "vial_ID")
}

# preface plate ID's with p2 or p3 if 2nd or 3rd plate respectively
plate_data[[2]]$plate_ID <- paste("p2", plate_data[[2]]$plate_ID, sep = "")
plate_data[[3]]$plate_ID <- paste("p3", plate_data[[3]]$plate_ID, sep = "")

# put into one dataframe
plate_data_final <- plate_data[[1]]
for(i in 2:length(plate_data)) {
  plate_data_final <- rbind(plate_data_final, plate_data[[i]])
}


# read in sample metadata & add into plate data
metadata <- read.csv("./data/molecular/metadata/16s_sample_metadata.csv")

# join in plate_data_final with metadata
metadata <- left_join(metadata, plate_data_final, by = c("vial_ID"))

## (b) reading in PICRUSt2-SC estimated 16s rRNA copy number values
picrust_predfunc <- read_tsv("./data/molecular/picrust2_outputs/nochimera_rarefied95/pred_metagenome_unstrat_descrip2.tsv") %>% 
  mutate(`function` = str_remove(`function`, "ko:")) %>% 
  dplyr::rename(ko_id = `function`)

## (c) reading in PICRUST2-SC weighted NSTI values
picrust_nsti <- read_tsv("./data/molecular/picrust2_outputs/nochimera_rarefied95/weighted_nsti2.tsv") %>% 
                         dplyr::rename(plate_ID = sample)

## (d) reading in KO_2_Kegg reference from "ggpicrust2"
ref <- ko_to_kegg_reference %>% 
  #select(pathway_number, pathway_name, ko_id, ko_description, level1) %>% 
  filter(ko_id %in% picrust_predfunc$ko_id) %>% 
  # select only pathway groups we'd realistically have pathways in
  filter(level1%in% c("09100 Metabolism", "09120 Genetic Information Processing",
                      "09130 Environmental Information Processing","09140 Cellular Processes"))

#### (2) Joining Data Together ####

# pivot longer and join with plate metadata
data = left_join(picrust_predfunc %>% pivot_longer(cols = c(3:ncol(picrust_predfunc)) ,  names_to = "plate_ID",
                                                              values_to = "abundance") %>% filter(abundance != 0),
                 metadata, by = "plate_ID")

# what samples are we missing?
intersect(metadata$vial_ID)
which(is)

#### (3) Adding Sample Annotations ####

# deciding to only care about pathways we care about
important_refs <- refs %>% 
  # nitrogen metabolism, lipid metabolism, thiamine synthesis, 
  filter(pathway_number %in% c(00910, 99983, 00730, ))

# in cases where multiple ko_ids match multiple pathways, ggpicrust2 was making some weird choices,
# (e.g., assigning stuff to worm longetivity pathway instead of cellular processes so going to do this manually ... :)
# and mostly pay attention to those with the largest number of each
















# view most abundant genes
summary = data %>% 
  dplyr::group_by(ko_id) %>% 
  dplyr::summarize(median = median(abundance))

# get references that are NOT duplicated this is not working IMO
single_instance <- ref %>% 
  group_by(ko_id) %>% 
  dplyr::summarize(count = length(ko_id)) %>% 
  ungroup() %>% 
  filter(count == 1) %>% 
  select(ko_id)
not_duplicates <- ref %>% 
  filter(ko_id %in% single_instance$ko_id)

# join in with summary 
summary <- left_join(summary, ref_not_duplicates, by = c("ko_id"))

#
duplicated_refs = ref[which(duplicated(ref$ko_id)),]

# for the remainder, need to determine which 

#### Misc. Prelim Data Exploration ####
ggplot(data = data, aes(x = sample_type, y = abundance, fill = `function`)) +
  geom_bar(stat = "identity")

# just look at nitrogen fixers
nfixers <- c("ko:K02588", "ko:K02586", "ko:K02591")

# filter out
nfixersonly <- data %>% filter(`function` %in% nfixers)


#### LEFT OFF HERE: WANT TO DO CHECKING IN 3a FIRST!

# remove "fake" targets 
# (taken when target taxa was absent- won't be analyzing these but sometimes found 
# target taxa even if it was not present! see script "2c_processing_QIIME_outputs.R")

#### (3) Relativizing Abundances ####

#### (4) Processing Triplicates ####

#### (5) Adding KEGG Information and Making Broader Groupings ####

# this should be added with python script?










#### draft script for processing PICRUSt2-SC data

# things to think of: we removed low quality reads and mitochondria, chloroplasts from others
# does this mean we need stratified samples?

# ideal overview: read in data, match to 95% rarefied data, 
# figure out broader functional grouping and add in

# see how abundance matches on each dataframe?
# consider saving outputs from script 2d elsewhere
# and having a script 3c where predicted functional groups are verified


### OLD below

library(ggpicrust2)
library(tidyverse)

# read in files (following ggpicrust2 tutorial)
abundance_file <- "./data/molecular/picrust2_outputs/pred_metagenome_unstrat1.tsv"

# get list of file names
plate_files <- list.files(path = "./data/molecular/metadata/", pattern = "plate")

# create empty list (different plates have same IDs so want to keep them separate!)
plate_data <- list()

# fill in list with dataframes
for(i in 1:length(plate_files)) {
  plate_data[[i]] <- read.table(paste("data/molecular/metadata/", plate_files[i], sep = ""))
  colnames(plate_data[[i]]) <- c("plate_ID", "vial_ID")
}

test = read_tsv("./data/molecular/picrust2_outputs/pred_metagenome_unstrat1.tsv") %>% 
  mutate(`function` = str_remove(`function`, "ko:"))
test2 = read_tsv("./data/molecular/picrust2_outputs/pred_metagenome_unstrat2.tsv") %>% 
  mutate(`function` = str_remove(`function`, "ko:"))
#test = ko2kegg_abundance("./data/molecular/picrust2_outputs/pred_metagenome_unstrat1.tsv")
#test = cbind(rownames(test), test)
#test2 = ko2kegg_abundance("./data/molecular/picrust2_outputs/pred_metagenome_unstrat2.tsv")
#test2 = cbind(rownames(test2), test2)
#colnames(test)[1] = "pathway_id"
#colnames(test2)[1] = "pathway_id"

view(test)

# pivot longer
test_longer <- test %>% 
  pivot_longer(c(2:ncol(test)), values_to = "abundance", names_to = "plate_ID")
test_longer2 <- test2 %>% 
  pivot_longer(c(2:ncol(test2)), values_to = "abundance", names_to = "plate_ID")

together1 <- left_join(test_longer, plate_data[[1]],
                       by = "plate_ID")
together2 <- left_join(test_longer2, plate_data[[2]],
                       by = "plate_ID")

together_all <- rbind(together1, together2)
together_all <- left_join(together_all, metadata,
                          by = "vial_ID")
colnames(together_all)[1] <- "ko_id"

# need to figure out what each kegg pathway is; use ggpicrust2
view(ko_to_kegg_reference)

together_all <- left_join(together_all, ko_to_kegg_reference, by = "ko_id")

# struggling to match the ko and kegg abundances??

# things to do:
# look at NMDS plots

test = pathway_annotation(file = "./data/molecular/picrust2_outputs/nochimera_rarefied95/pred_metagenome_unstrat2.tsv",
                          pathway = "KO",
                          ko_to_kegg = FALSE)
test = ko2kegg_abundance(data = picrust_predfunc %>% select(!description))
test = read_tsv

kegg_data = read_csv("./data/molecular/picrust2_outputs/kegg_pathways.csv", col_names = FALSE)
colnames(kegg_data) = c("kegg_pathway", "pathway_name")


#### Playing with KEGG data ####