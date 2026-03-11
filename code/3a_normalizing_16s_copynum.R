#### Normalizing relative abundance with predicted 16s gene copy numbers
### Jordan Zabrecky
## last edited: 03.10.2026

# This script pulls in relative abundance based on the predicted 16s gene 
# copy numbers for each ASV obtained via PICRUSt2-SC and compares it to 
# the relative abundance without doing so (QIIME2 outputs obtained and
# processed in step 2 of code)

# Note: This will only be done for rarefied (95% threshold) data
# as that balances samples and does not remove many

# TO DO: rerun QIIME2 scripts and save 95 rarefied instead and change name
# remove anacyl and microcoleus at end

#### (1) Loading libraries & data ####

# loading libraries
lapply(c("tidyverse"), require, character.only = T)

## (a) getting metadata

# get list of plates in order (will read in 1, 2, then 3)
plate_data <- lapply(list.files("./data/molecular/metadata/",
                                  pattern = "plate"), function(x) { 
                                      y = read.table(paste("./data/molecular/metadata/", x, sep = ""))
                                      colnames(y) <- c("plate_ID", "vial_ID")
                                      return(y)
                                  })

# read in sample metadata & add into plate data
metadata <- read.csv("./data/molecular/metadata/16s_sample_metadata.csv")

## (b) reading in PICRUSt2-SC NSTI values (whereby lower values closer to 0 indicate better match)
picrust_nsti <- lapply(list.files("./data/molecular/picrust2_outputs/nochimera_nonrarefied",
                                  pattern = "nsti"), function(x) 
                                  read_tsv(paste("./data/molecular/picrust2_outputs/nochimera_nonrarefied/", x, sep = "")) %>% 
                                             dplyr::rename(plate_ID = sample))

## (c) reading in PICRUSt2-SC estimated 16s rRNA copy number values
picrust_predcopynum <- lapply(list.files("./data/molecular/picrust2_outputs/nochimera_nonrarefied",
                                        pattern = "seqtab"), function(x) 
                                          read_tsv(paste("./data/molecular/picrust2_outputs/nochimera_nonrarefied/", x, sep = "")))

## (d) reading in QIIME2 outputs
qiime2 <- read.csv("./data/molecular/16s_nochimera_rarefied_95.csv")

#### (2) Matching PICRUSt2 outputs to metadata ####

# merge in plate data with nsti values
merged_data <- list()
for(i in 1:length(plate_data)) {
  merged_data[[i]] <- left_join(plate_data[[i]], picrust_nsti[[i]], by = "plate_ID")
}

# use decostand in vegan to relativize abundances from PICRUSt2 that are normalized by 16s copy number
picrust_copynumnormalized = lapply(picrust_predcopynum, function(x) {
  # first transpose to use vegan::decostand which needs taxa as columns and samples as rows
  y = as.data.frame(t(x)) %>% 
    # get plateIDs into a column and move to front
    mutate(plate_ID = rownames(.)) %>% 
    relocate(plate_ID, .before = 1)
  
  # move ASV names to column names
  colnames(y) = y["normalized",]
  
  # remove old row of ASVs and change name of first column to plate_ID
  y = y[-1,] 
  
  # convert abundance values from character to numeric
  y[,2:ncol(y)] = sapply(y[,2:ncol(y)], as.numeric)
  
  # normalize with vegan::decostand
  z = y
  z[,2:ncol(z)] = decostand(y[,2:ncol(y)], method = "total")
  
  # confirm that abundances all add to 1
  print(rowSums(z[,2:ncol(z)]))
  
  # return normalized data
  return(z)
})

# join in with plate data

# join in sample metadata

# will probably need to move relativized to the bottom

#### (4) Dealing with Triplicates ####

#### (3) Comparing ASV normalized and not normalized ####

#### (4) Saving outputs ####

