#### Comparing microscopy data among rivers
### Jordan Zabrecky
## last edited: 11.14.2025

## This code compares 16s rRNA (microbial community) data from NT, TM, and TAC samples
## across rivers to answer Q1

#### (1) Loading libraries & data ####

# set seed for reproducibility
set.seed(2025)

# libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot"), require, character.only = T)

# read in files (note: two target taxa csvs- one with target taxa included, other with it excluded)
# doing as a list rather than one csv as each has a different # of columns
files <- list.files(path = "./data/morphological/", pattern = ".csv")
data_wide <- lapply(files, function(x) read.csv(paste("./data/microscopy/", x, sep = "")))
names(data_wide) <- files
