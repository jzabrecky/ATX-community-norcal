#### Comparing microscopy data across time
### Jordan Zabrecky
## last edited: 12.15.2025

## This code compares microscopy data from NT, TM, and TAC samples
## across rivers to answer Q2

## to-do:
## ANOSIM w/ time stamp
## indicator species analsis? vs. SIMPER?
## look into how other papers are determining which is the species most responsible

## CONSIDER WHAT TO DO WITH 2023 DATA

# MAYBE MOVE NMDS FUNCTIONS TO SUPPLEMENTAL CODE FILES

# ALSO ADD IN MONTH TO THIS

#### (1) Loading libraries & data ####

# set seed for reproducibility
set.seed(2025)

# libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot"), require, character.only = T)

# MAYBE ONLY READ IN SELECT FILES? 

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

#### (3) PERMANOVA ####

# Is there a significant difference in communities across time?
dist_matrix = vegdist(data_wide$tm_algalonly_nomicro_sqrttransformed.csv
                      [,7:ncol(data_wide$tm_algalonly_nomicro_sqrttransformed.csv)], method = "bray")
adonis2(dist_matrix ~ event_no, data = data_wide$tm_algalonly_nomicro_sqrttransformed.csv,
        strata = data_wide$tm_algalonly_nomicro_sqrttransformed.csv$site)

# how about we split data
sfkeel_only <- data_wide$tm_algalonly_nomicro_sqrttransformed.csv %>% 
  filter(site == "SFE-M")
salmon_only <- data_wide$tm_algalonly_nomicro_sqrttransformed.csv %>% 
  filter(site == "SAL")

sfkeeldist <- vegdist(sfkeel_only[,7:ncol(sfkeel_only)], method = "bray")
adonis2(sfkeeldist ~ event_no, data = sfkeel_only) # significant

saldist <- vegdist(salmon_only[,7:ncol(salmon_only)], method = "bray")
adonis2(saldist ~ event_no, data = salmon_only) # not significant

# 
sfkeel_only <- data_wide$tm_algalonly_nomicro_sqrttransformed.csv %>% 
  filter(site == "SFE-M")
salmon_only <- data_wide$tm_algalonly_nomicro_sqrttransformed.csv %>% 
  filter(site == "SAL")

sfkeeldist <- vegdist(sfkeel_only[,7:ncol(sfkeel_only)], method = "bray")
adonis2(sfkeeldist ~ event_no, data = sfkeel_only) # significant

saldist <- vegdist(salmon_only[,7:ncol(salmon_only)], method = "bray")
adonis2(saldist ~ event_no, data = salmon_only) # not significant

# maybe want to look at rivers separately?

#### (3) test indicator species analysis

# Probably want to add this to script 5a WILL GO BACK AND DO THIS NEXT
# maybe try out both???
# can you use indicator species analyses on relative abundances?

# more curious how this will change over time

# does this need to be hellinger transformed and why?
library(indicspecies)

# work with NT data to start
NT <- data_wide$nt_algalonly_sqrttransformed.csv

# func = "r.g" is for relative abundances 
groups <- NT$
test <- multipatt(NT[,7:ncol(NT)], groups, func = "r.g", control = how(nperm = 999))
summary(test)
