#### putting together microscopy data for EDI data package
### Jordan Zabrecky
## last edited: 12.09.2024

# This code reads in microscopy data for each sample type, makes sure the data
# is complete, and then saves as combined csvs

#### (1) Loading packages and reading in data #### 

# loading libraries
lapply(c("tidyverse", "plyr"), require, character.only = T)

## loading raw microscopy data by sample type

# target-microcoleus (TM) samples
TM_samples <- ldply(list.files(path = "./data/morphological/raw/", pattern = "TM"), function(filename) {
  d <- read.csv(paste("data/morphological/raw/", filename, sep = ""))
  return(d)
})

# target-anabaena and cylindrospermum (TAC) samples
TAC_samples <- ldply(list.files(path = "./data/morphological/raw/", pattern = "TAC"), function(filename) {
  d <- read.csv(paste("data/morphological/raw/", filename, sep = ""))
  return(d)
})

# non-target (NT) samples
NT_samples <- ldply(list.files(path = "./data/morphological/raw/", pattern = "NT"), function(filename) {
  d <- read.csv(paste("data/morphological/raw/", filename, sep = ""))
  return(d)
})

#### (2) Misc. data checks and save csvs ####

# replace NA with 0's
TM_samples <- replace(TM_samples, is.na(TM_samples), 0)
TAC_samples <- replace(TAC_samples, is.na(TAC_samples), 0)
NT_samples <- replace(NT_samples, is.na(NT_samples), 0)

# checking no multiples of column with mispellings, etc.
colnames(TM_samples[8:29])
colnames(TAC_samples[8:31])
colnames(NT_samples[8:55])

# checking the each sample (row) adds up to 100%
which(rowSums(TM_samples[8:29]) != 100)
which(rowSums(TAC_samples[8:31]) != 100)
which(rowSums(NT_samples[8:55]) != 100)
# all add up to 100!

# function to create site_reach column from site and reach information
site_reach <- function(df) {
  df$site_reach <- paste(df$site, "-", df$reach, sep = "")
  new_df <- df %>% 
    dplyr::select(!site & !reach) %>% # get rid of separate site_reach columns
    relocate(site_reach, .before = sample_type)
  return(new_df)
}

# apply function to final data frames
TM_samples_final <- site_reach(TM_samples)
TAC_samples_final <- site_reach(TAC_samples)
NT_samples_final <- site_reach(NT_samples)

# save each as separate csv
write.csv(TM_samples_final, "./data/EDI_data_package/microscopy_target_microcoleus.csv", row.names = FALSE)
write.csv(TAC_samples_final, "./data/EDI_data_package/microscopy_target_anabaena_cylindrospermum.csv", row.names = FALSE)
write.csv(NT_samples_final, "./data/EDI_data_package/microscopy_non_target.csv", row.names = FALSE)
