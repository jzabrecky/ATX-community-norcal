#### putting together microscopy data for EDI data package
### Jordan Zabrecky
## last edited: 01.27.2025

# This code reads in microscopy data for each sample type, makes sure the data
# is complete, and then saves as combined csvs

#### (1) Loading packages and reading in data #### 

# loading libraries
lapply(c("tidyverse", "plyr", "gtools"), require, character.only = T)

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

# checking no multiples of column with mispellings, repeats, etc.
colnames(TM_samples[8:29])
colnames(TAC_samples[8:30])
colnames(NT_samples[8:54])

# checking the each sample (row) adds up to 100%
which(rowSums(TM_samples[8:29]) != 100)
which(rowSums(TAC_samples[8:30]) != 100)
which(rowSums(NT_samples[8:54]) != 100)
# all add up to 100!

# joining TM and TAC to one target sample data frame
T_samples <- smartbind(TM_samples, TAC_samples)

# replace NA with 0's
T_samples <- replace(T_samples, is.na(T_samples), 0)
NT_samples <- replace(NT_samples, is.na(NT_samples), 0)

# function to create site_reach column from site and reach information
site_reach <- function(df) {
  df$site_reach <- paste(df$site, "-", df$reach, sep = "")
  new_df <- df %>% 
    dplyr::select(!site & !reach) %>% # get rid of separate site_reach columns
    relocate(site_reach, .before = sample_type) %>% 
    mutate(field_date = mdy(field_date), # change date format to yyyy-mm-dd for EDI
           date_analyzed = mdy(date_analyzed))
  return(new_df)
}

# apply function to final data frames
T_samples_final <- site_reach(T_samples)
NT_samples_final <- site_reach(NT_samples)

# save each as separate csv
write.csv(T_samples_final, "./data/EDI_data_package/microscopy_target_samples.csv", row.names = FALSE)
write.csv(NT_samples_final, "./data/EDI_data_package/microscopy_non_target_samples.csv", row.names = FALSE)
