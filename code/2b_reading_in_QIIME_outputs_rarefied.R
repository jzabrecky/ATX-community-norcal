#### Reading in and putting together QIIME2 output that were rarefied
### Jordan Zabrecky
## last edited 10.07.2025

# This code reads in QIIME2 outputs (sequence abundances and SILVA taxonomy assignment;
# for processing that included rarefaction) and matches them with metadata. 
# The script detects what samples were missing (either due to label misreadings or 
# low quality reads that were filtered in QIIME2) and summarizes that. It then saves 
# the final dataframe to be further processed in another script.

# Note: two rarefied dataframes are read in here. "_90" refers to rarefaction sampling
# depth where 90% of samples are retained and "_95" refers to rarefaction sampling 
# depth where 95% of samples are retained

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

# add names to biom_data list to know which is which
names(biom_data) <- biom_files

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
  mutate(domain = case_when(grepl("p__", taxon_full) ~ 
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
                           TRUE ~ str_extract(taxon_full, "(?<=g__).*")),
         # for species there is no further classification so anything after s
         species = str_extract(taxon_full, "(?<=s__).*"))

#### (2) Merging dataframes together ####

## (a) merging sequence abundances and plate IDs w/ metadata

# create list of merged data
merged_data <- list()

# left join in plate metadata
for(i in 1:length(biom_data)) {
  merged_data[[i]] <- left_join(biom_data[[i]], plate_data_metadata, by = "plate_ID")
}

# add in names to merged_data list (distinguish 90 vs 95 rarefied)
names(merged_data) <- names(biom_data)

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
merged_data[[1]]$vial_ID[which(is.na(merged_data[[1]]$site_reach))]
# none missing!

## (b) adding in taxonomy

# add in taxonomy
for(i in 1:length(merged_data)) {
  merged_data[[i]] <-left_join(merged_data[[i]], taxonomy_data, by = "feature_ID")
}

# check to see each sequence ID got an assignment! (change 1 & 2 manually :) )
merged_data[[1]]$vial_ID[which(is.na(merged_data[[1]]$taxon_full))]
# yay everything transferred!

## (c) some final checks

## LEFT OFF HERE 7/24/2025

# want to see if we are missing any vials from metadata (change 1 & 2 manually :) )
eval(length(unique(merged_data[[2]]$vial_ID)) == length(metadata$vial_ID))
# no, seems like we are missing some
setdiff(metadata$vial_ID, unique(merged_data[[1]]$vial_ID)) # missing IDs:
# 4  18  34  40  43  58  86  90  94  98  99 102 104 107 108 115 124 135 137 142 208 215 216
setdiff(metadata$vial_ID, unique(merged_data[[2]]$vial_ID)) # missing IDs:
# 4  18  34  40  43  58  99 102 104 108 115 135 137 142 208 215

# we know from looking at nonrarefied data that:
# 104, 108, and 208 were lost in translation
# 18, 34, and 135 had poor sequencing quality

## lost from 95% sequencing depth rarefaction:
# 4: SAL-3 blank 6-27-2022 (inconsequential!)
# 40: SAL-2 TM 7-26-2022
# 43: SFE-M-3 TAC 7-14-2022
# 58: SFE-M-3 TAC 7-28-2022 (triplicate- inconsequential!)
# 99: RUS-1S TM 8-17-2022 (fake target- inconsequential!)
# 102: RUS-1S blank 8-17-2022 (blank- inconsequential!)
# 115: RUS-2 TAC 9-1-2022
# 137: SFE-M-1S blank 9-6-2022 (blank- inconsequential!)
# 142: RUS-1S TM 9-1-2022 (fake target- inconsequential!)
# 215: RUS-1S TM 9-15-2022 (fake target- inconsequential!)

## additional lost from 90% sequencing depth rarefaction:
# 86: SFE-M-3 TM 8-10-2022
# 90: SFE-M-4 TAC 8-10-2022
# 94: RUS-1S TAC 8-17-2022
# 98: SFE-M-3 TM 8-23-2022 (triplicate- inconsequential!)
# 107: RUS-2 TAC 8-17-2022 
# 124: RUS-3 TAC 9-1-2022
# 216: SAL-1S NT 9-22-2022

# to compare, how many samples did we have total 
# (excluding blanks, fake targets, and counting triplicates as one sample)
nrow(unique(metadata %>% filter(sample_type != "blank" & fake_target == "n") %>% 
  select(site_reach, sample_type, field_date))) # 105
nrow(unique(merged_data[[1]] %>% filter(sample_type != "blank" & fake_target == "n") %>% 
              select(site_reach, sample_type, field_date))) # 93 for 90% rarefaction
nrow(unique(merged_data[[2]] %>% filter(sample_type != "blank" & fake_target == "n") %>% 
              select(site_reach, sample_type, field_date))) # 99 for 95% rarefaction

# NOTE: will remove chloroplasts, unassigned taxa, "fake targets", etc. in future script

#### (3) Re-organizing and saving data ####

# trim dataframe and save!
for(i in 1:length(merged_data)) {
  # keep columns we care about
  final <- merged_data[[i]] %>% select(site_reach, site, field_date, sample_type, material, triplicate, fake_target,
                              container, plate_ID, vial_ID, abundance, feature_ID, taxon_full, domain, 
                              phylum, class, order, family, genus, species, confidence)
  # save!
  write.csv(final, paste("./data/molecular/intermediate_csvs/16s_nochimera_", 
                         str_match(names(merged_data)[i], "table_\\s*(.*?)\\s*.biom")[,2], 
                         "_unfiltered.csv", sep = ""), row.names = FALSE)
}
