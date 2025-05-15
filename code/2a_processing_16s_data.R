#### Processing .biom files
### Jordan Zabrecky
## last edited 05.15.2025

## TBD!?!?
# next script: reading in rarefied data
# next script: a supplemental comparing the two

# code outline- read in no chimera, read in rarefied, 
# analyze triplicate 16s rRNA, average across samples??!?!?!? look at blanks; relativize in this one? 
# remove "maybe" or fake targets

#### (1) Loading in libraries and data ####

# loading libraries
lapply(c("tidyverse", "plyr", "biomformat"), require, character.only = T)

## (a) read in .biom files

# get file path names
biom_files <- list.files(path = "./data/molecular/raw_files", pattern = "_nochim")

# create empty list 
# (need to keep them separate because IDs repeat across different plates but aren't same sample)
biom_data <- list()

# fill in list with dataframes
for(i in 1:length(biom_files)) {
  raw = read_biom(paste("data/molecular/raw_files/", biom_files[i], sep = ""))
  temp <-  as(biom_data(raw), 'TsparseMatrix')
  df <- data.frame(feature_ID = temp @Dimnames[[1]][ temp @i + 1],
                   plate_ID = temp @Dimnames[[2]][ temp @j + 1],
                   abundance = temp @x)
  biom_data[[i]] <- df
}

## (b) read in sequencing plate and sample metadatafiles

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
  df$phylum <- str_match(df$taxon_full, "p__\\s*(.*?)\\s*; c__")[,2]
  df$class <- str_match(df$taxon_full, "c__\\s*(.*?)\\s*; o__")[,2]
  df$order <- str_match(df$taxon_full, "o__\\s*(.*?)\\s*; f__")[,2]
  df$family <- str_match(df$taxon_full, "f__\\s*(.*?)\\s*; g__")[,2]
  df$genus <- str_match(df$taxon_full, "g__\\s*(.*?)\\s*; s__")[,2]
  df$species <- str_extract(df$taxon_full, "(?<=s__).*")
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
# 2E is vial 801, 9G is vial 554, C5 is 931, all for another project

# remove those with missing vial IDs
merged_data <- lapply(merged_data, function(x) {
  x <- x[-which(is.na(x$vial_ID)),]
})

# join in sample metadata
for(i in 1:length(plate_files)) {
  merged_data[[i]] <- left_join(merged_data[[i]], metadata, by = "vial_ID")
}

# checking to see if any are missing sample names
merged_data[[1]]$vial_ID[which(is.na(merged_data[[1]]$site_reach))]

## (b) adding in taxonomy

# add in taxonomy
for(i in 1:length(plate_files)) {
  merged_data[[i]] <-left_join(merged_data[[i]], taxonomy_data[[i]], by = "feature_ID")
}

# check to see each sequence ID got an assignment!

## (c) putting all into one dataframe


# will remove chloroplasts, unassigned taxa, etc. in a future list

#### (3) Re-organizing 

## (c) read in metadata files
# metadata & plate files
# check unique for metadata

#### add in metadata, assign in taxonomy, relative abundance!?!

### repeat this script with rarefied sequences in a second script

# rarefied version (i.e. where low abundance samples get 
# filtered out and existing features get subsampled.)