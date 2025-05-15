#### Processing .biom files
### Jordan Zabrecky
## last edited 03.15.2025

## TBD!?!?

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
  df <- data.frame(taxon_ID = temp @Dimnames[[1]][ temp @i + 1],
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
}

# read in sample metadata
metadata <- read.csv("./data/molecular/metadata/16s_sample_metadata.csv")

## (c) read in Silva taxnomy files

# all feature IDs are unique so can read into one larger dataframe!

taxonomy1 <- read_tsv("./data/molecular/raw_files/taxonomy1_silva.tsv")
taxonomy2 <- read_tsv("./data/molecular/raw_files/taxonomy2_silva.tsv")
merged <- rbind(taxonomy1, taxonomy2)
# can put them all together??
# will there be repeats of IDs?

#### (3) Merging together

## (c) read in metadata files
# metadata & plate files
# check unique for metadata

#### add in metadata, assign in taxonomy, relative abundance!?!

### repeat this script with rarefied sequences in a second script

# rarefied version (i.e. where low abundance samples get 
# filtered out and existing features get subsampled.)