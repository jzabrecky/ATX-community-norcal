#### Comparing microscopy data across time
### Jordan Zabrecky
## last edited: 12.15.2025

## This code compares microscopy data from NT, TM, and TAC samples
## across rivers to answer Q2

## to-do: cluster analysis
## ANOSIM w/ time stamp
## indicator species analsis?

## CONSIDER WHAT TO DO WITH 2023 DATA

# MAYBE MOVE NMDS FUNCTIONS TO SUPPLEMENTAL CODE FILES

#### (1) Loading libraries & data ####

# set seed for reproducibility
set.seed(2025)

# libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot"), require, character.only = T)

# read in files (note: two target taxa csvs- one with target taxa included, other with it excluded)
# doing as a list rather than one csv as each has a different # of columns
files <- list.files(path = "./data/morphological/transformed/", pattern = ".csv")
data_wide <- lapply(files, function(x) read.csv(paste("./data/morphological/transformed/", x, sep = "")))
names(data_wide) <- files

# only doing square-root transformation, as established in previous script,
# "4a_amongrivers_microscopy.R" not transforming or additionally removing rare 
# taxa did not change results of analysis in a meaningful way 
# (and that is already complete as we read in that data)

#### (2) Add  Tag for Sampling Time ####

# we have field dates for sampling distinguishing samples, however, let's change it
# to be the the number of sampling event (where x is data)
add_event_no <- function(data) {
  data %>% 
    mutate(field_date = ymd(field_date)) %>% 
    mutate(event_no = case_when((field_date >= ymd("2022-06-24") & field_date <= ymd("2022-06-29")) ~ 1,
                                (field_date >= ymd("2022-07-06") & field_date <= ymd("2022-07-14")) ~ 2,
                                (field_date >= ymd("2022-07-20") & field_date <= ymd("2022-07-28")) ~ 3,
                                (field_date >= ymd("2022-08-02") & field_date <= ymd("2022-08-10")) ~ 4,
                                (field_date >= ymd("2022-08-17") & field_date <= ymd("2022-08-23")) ~ 5,
                                (field_date >= ymd("2022-09-01") & field_date <= ymd("2022-09-06")) ~ 6,
                                (field_date >= ymd("2022-09-15") & field_date <= ymd("2022-09-22")) ~ 7)) %>% 
    relocate(event_no, .before = "field_date")
}

# apply to all dataframes
data_wide <- lapply(data_wide, function(x) add_event_no(x))

#### (3) test cluster analysis

# test with microcoleus
tm_dist <- vegdist(data_wide$tm_algalonly_nomicro_sqrttransformed.csv[,7:ncol(data_wide$tm_algalonly_nomicro_sqrttransformed.csv)],
                   method = "bray")
tm_dist_eu <- vegdist(data_wide$tm_algalonly_nomicro_sqrttransformed.csv[,7:ncol(data_wide$tm_algalonly_nomicro_sqrttransformed.csv)],
                      method = "euc")

# need to square root bray curtis dissimilarity for use with ward's mininum variance
tm_dist <- sqrt(tm_dist)
apply(tm_dist, 1:2, sqrt)
library("cluster")

# add labels?
attr(tm_dist, "labels") <- data_wide$tm_algalonly_nomicro_sqrttransformed.csv$site

# cluster
ward <- hclust(tm_dist, method = "ward.D2")
plot(ward, main = "Algal Assemblages Microcoleus", 
     labels = paste(data_wide$tm_algalonly_nomicro_sqrttransformed.csv$site_reach, data_wide$tm_algalonly_sqrttransformed.csv$field_date))

cluster_analyses <- function(data, plot_title) {
  dist = vegdist(data[,7:ncol(data)], method = "bray")
  dist_trans = sqrt(dist)
  
  ward = hclust(dist, method = "ward.D2")
  plot(ward, main = plot_title, 
       labels = paste(data$site_reach, " no: ", data$event_no))
}

titles <- c("Non-Target Samples", "Anabaena/Cylindrospermum Samples (including)", 
            "Anabaena/Cylindrospermum Samples (excluding)", 
            "Anabaena/Cylindrospermum Samples (also excluding GA)",
            "Microcoleus Samples (including)",
            "Microcoleus Samples (excluding)")

for(i in 1:length(data_wide)) {
  cluster_analyses(data_wide[[i]], titles[i])
}    
