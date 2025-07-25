#### Reading in and putting together QIIME2 output that were rarefied
### Jordan Zabrecky
## last edited 06.20.2025

## NOTE CURRENTLY SEEMS TO BE AN ISSUE WITH FEATURE TABLE 1 FOR RAREFIED PROCESSING

## This code reads in QIIME2 outputs (sequence abundances and SILVA taxonomy assignment;
## for processing that included rarefaction) and matches them with metadata. 
## The script detects what samples were missing (either due to label misreadings or 
## low quality reads that were filtered in QIIME2) and summarizes that. It then saves 
## the final dataframe to be further processed in another script.

#### (1) Loading in libraries and data ####

# loading libraries
lapply(c("tidyverse", "plyr", "biomformat"), require, character.only = T)

## (a) read in .biom files

# get file path names 
# (instead of plates separate here, all plates are together; two for INSERT)
biom_files <- list.files(path = "./data/molecular/raw_files/rarefied", pattern = ".biom")

# create empty list 
# (need to keep them separate because IDs repeat across different plates but aren't same sample)
biom_data <- list()

# fill in list with dataframes
for(i in 1:length(biom_files)) {
  # (some slight lexical error in read_biom but seems to be reading in fine)
  raw <- read_biom(paste("./data/molecular/raw_files/rarefied/", biom_files[i], sep = ""))
  temp <-  as(biom_data(raw), 'TsparseMatrix')
  df <- data.frame(feature_ID = temp @Dimnames[[1]][ temp @i + 1],
                   plate_ID = temp @Dimnames[[2]][ temp @j + 1],
                   abundance = temp @x)
  biom_data[[i]] <- df
}

## (b) read in sequencing plate and sample metadata files

##  read in plate_IDs

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
plate_data_metadata <- left_join(plate_data_final, metadata, by = "vial_ID")

## (c) read in Silva taxnomy files

# read in taxonomy file and filter 
taxonomy_data <- read_tsv("./data/molecular/raw_files/rarefied/taxonomy_cy.tsv")
colnames(taxonomy_data) <- c("feature_ID", "taxon_full", "confidence")

# break down full taxonomy assignment
taxonomy_data <- taxonomy_data %>% 
  # if phylum is given, take phrase between d and p, else take entire phrase after d
  mutate(domain = case_when(grepl("c__", taxon_full) ~ 
                              str_match(taxon_full, "d__\\s*(.*?)\\s*; p__")[,2],
                            TRUE ~ str_extract(taxon_full, "(?<=d__).*")),
         # if class is given, take phrase between p and c, else take entire phrase after c
         phylum = case_when(grepl("c__", taxon_full) ~ 
                              str_match(taxon_full, "p__\\s*(.*?)\\s*; c__")[,2],
                            TRUE ~ str_extract(taxon_full, "(?<=p__).*")),
         # if order is given, take phrase between c and o, else take entire phrase after c
         class = case_when(grepl("o__", taxon_full) ~ 
                             str_match(taxon_full, "c__\\s*(.*?)\\s*; o__")[,2],
                           TRUE ~ str_extract(taxon_full, "(?<=c__).*")),
         # if family is given, take phrase between  and f, else take entire phrase after o
         order = case_when(grepl("f__", taxon_full) ~
                             str_match(taxon_full, "o__\\s*(.*?)\\s*; f__")[,2],
                           TRUE ~ str_extract(taxon_full, "(?<=o__).*")),
         # if genus is given, take phrase between f and g, else take entire phrase after f
         family = case_when(grepl("g__", taxon_full) ~
                              str_match(taxon_full, "f__\\s*(.*?)\\s*; g__")[,2],
                            TRUE ~ str_extract(taxon_full, "(?<=f__).*")),
         # if species is given, take phrase between g and s, else take entire phrase after g
         genus = case_when(grepl("s__", taxon_full) ~
                             str_match(taxon_full, "g__\\s*(.*?)\\s*; s__")[,2],
                           TRUE ~ str_extract(df$taxon_full, "(?<=g__).*")),
         # for species there is no further classification so anything after s
         species = str_extract(df$taxon_full, "(?<=s__).*"))

#### (2) Merging dataframes together ####

# empty data frame for merged files
merged_data <- list()

## (a) merging sequence abundances and plate IDs w/ metadata

# have separate no of plates and biom data

# left join in plate metadata
for(i in 1:length(biom_data)) {
  merged_data[[i]] <- left_join(biom_data[[i]], plate_data_metadata, by = "plate_ID")
}

# checking to see if any are missing a vial
merged_data[[2]]$plate_ID[which(is.na(merged_data[[2]]$vial_ID))]
# 9G from plate 2, C5 from plate 3
view(plate_data[[1]])
# 9G is vial 554, C5 is 931
# these were either from another project or labels that got lost in translation

# remove those with missing vial IDs
merged_data <- lapply(merged_data, function(x) {
  x <- x[-which(is.na(x$vial_ID)),]
})

# let's also check that there are no repeats in vial_IDs on plates
# in case things got lost in translation...
plate_test <- append(plate_data[[1]]$vial_ID, plate_data[[2]]$vial_ID)
plate_test <- append(plate_test, plate_data[[3]]$vial_ID)

# check to see that they are same length
eval(length(plate_test) == length(unique(plate_test)))

# we have one repeating sample which is 108 so let's get rid of 
# since there was clearly an error in translation which is a bummer but what can you do
merged_data <- lapply(merged_data, function(x) {
  x <- x %>% filter(vial_ID != 108)
})

# checking to see if any are missing sample names
merged_data[[2]]$vial_ID[which(is.na(merged_data[[2]]$site_reach))]

## (b) adding in taxonomy

# add in taxonomy
for(i in 1:length(merged_data)) {
  merged_data[[i]] <-left_join(merged_data[[i]], taxonomy_data, by = "feature_ID")
}

# check to see each sequence ID got an assignment!
merged_data[[2]]$vial_ID[which(is.na(merged_data[[2]]$taxon_full))]
# yay everything transferred!

## (c) putting all into one dataframe

# NEED TO SEE WHAT THE DIFFERENCE IS BETWEEN THESE TWO FRAMES!!
# believe 1 is _90 and 2 is _95
# temporary
final <- merged_data[[2]]

# final dataframe with all plates
final <- rbind(merged_data[[1]], merged_data[[2]])
final <- rbind(final, merged_data[[3]])

# real quick want to see if we are missing any vials from metadata
eval(length(unique(final$vial_ID)) == length(metadata$vial_ID))
# no, seems like we are missing seven...
setdiff(metadata$vial_ID, unique(final$vial_ID)) # 18, 34, 104, 108, 115, 135, 208

# we know what happened for vial 108 from early (two plates had that ID)
# 208 could be 801 or the second 108 but unsure which would be the correct 108
# so will just have to accept that those samples are gone

# found vial IDs on plates for the following:
# vial 135 is plate 2 ID 8B, vial 18 is plate 3 ID B8, and vial 34 is plate 3 ID A11
# confirmed these sequences were read into the pipeline
# they were samples had very few reads and were ultimately filtered out by QIIME

# seems like 104 must have also just got lost in translation

## so in summary:
## vials lost in translation: 104, 108, 208
# 104 RUS-1S NT 8-17-2022
# 108 RUS-3 TM 8-17-2022 (fake TM meaning we did not see Microcoleus but tried to take a sample anyways)
# ^ (likely throwing those out so that sample is inconsequential)
# 208 RUS-1S NT 9-15-2022
## vials with low quality reads: 18, 34, 135
# 18 RUS-2 NT 7-6-2022 (luckily is a triplicate)
# 34 SAL-3 NT 7-12-2022
# 135 SFE-M-4 9-6-2022 (another fake TM sample)

# NOTE: will remove chloroplasts, unassigned taxa, "fake targets", etc. in future script

#### (3) Re-organizing and saving data ####

# using select on dataframe to quickly reorganize columns
final_tosave <- final %>% 
  select(site_reach, site, field_date, sample_type, material, triplicate, fake_target,
         container, plate_ID, vial_ID, abundance, feature_ID, taxon_full, domain, 
         phylum, class, order, family, genus, species, confidence)

# saving outputs
write.csv(final_tosave, "./data/molecular/16s_nochimera_rarefied.csv", row.names = FALSE)
