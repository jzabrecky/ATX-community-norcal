#### Reading in and putting together QIIME2 output that were rarefied
### Jordan Zabrecky
## last edited 05.16.2025

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
biom_files <- list.files(path = "./data/molecular/raw_files", pattern = "_rarefied")

# create empty list 
# (need to keep them separate because IDs repeat across different plates but aren't same sample)
biom_data <- list()

# fill in list with dataframes
for(i in 1:length(biom_files)) {
  # (some slight lexical error in read_biom but seems to be reading in fine)
  raw <- read_biom(paste("data/molecular/raw_files/", biom_files[i], sep = ""))
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

# read in sample metadata
metadata <- read.csv("./data/molecular/metadata/16s_sample_metadata.csv")

## (c) read in Silva taxnomy files

# get list of file names
taxonomy_files <- list.files(path = "./data/molecular/raw_files/", pattern = "taxonomy")

# create empty list (taxonomy for each plate)
taxonomy_data <- list()

# fill in list with dataframes
for(i in 1:length(taxonomy_files)) {
  df <- read_tsv(paste("./data/molecular/raw_files/", taxonomy_files[i], sep = ""))
  colnames(df) <- c("feature_ID", "taxon_full", "confidence")
  # break down full taxonomy assignment
  df <- df %>% 
    # if class is given, take phrase between p and c, else take entire phrase after c
    mutate(phylum = case_when(grepl("c__", taxon_full) ~ 
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
           species = str_extract(df$taxon_full, "(?<=s__).*")
    )
  
  # add in dataframe to taxonomy data list
  taxonomy_data[[i]] <- df
}

#### (2) Merging dataframes together ####

# empty data frame for merged files
merged_data <- list()

## (a) merging plate metadata and sequence abundances

# left join in plate metadata
for(i in 1:length(plate_files)) {
  merged_data[[i]] <- left_join(biom_data[[i]], plate_data[[i]], by = "plate_ID")
}

# checking to see if any are missing a vial
merged_data[[3]]$plate_ID[which(is.na(merged_data[[3]]$vial_ID))]
# 2E from plate 1, 9G from plate 2, C5 from plate 3
view(plate_data[[1]])
# 2E is vial 801, 9G is vial 554, C5 is 931
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

# join in sample metadata
for(i in 1:length(plate_files)) {
  merged_data[[i]] <- left_join(merged_data[[i]], metadata, by = "vial_ID")
}

# checking to see if any are missing sample names
merged_data[[3]]$vial_ID[which(is.na(merged_data[[3]]$site_reach))]

## (b) adding in taxonomy

# add in taxonomy
for(i in 1:length(plate_files)) {
  merged_data[[i]] <-left_join(merged_data[[i]], taxonomy_data[[i]], by = "feature_ID")
}

# check to see each sequence ID got an assignment!
merged_data[[1]]$vial_ID[which(is.na(merged_data[[1]]$taxon_full))]
# yay everything transferred!

## (c) putting all into one dataframe

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
         container, plate_ID, vial_ID, abundance, feature_ID, taxon_full, phylum, class,
         order, family, genus, species, confidence)

# saving outputs
write.csv(final_tosave, "./data/molecular/16s_nochimera.csv", row.names = FALSE)
