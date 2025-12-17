#### Side fun with clustering
### Jordan Zabrecky
## last edited: 12.16.2025

## This code compares clusters communities using the "cluster" package
## Instead, for analyses with are just using PERMANOVA and NMDS

#### (1) Loading libraries & data ####

# set seed for reproducibility
set.seed(2025)

# load libraries
lapply(c("tidyverse", "vegan", "cluster"), require, character.only = T)

# loading files individually because I only want certain ones
nt_micro <- read.csv("./data/morphological/transformed/nt_algalonly_sqrttransformed.csv")
tm_micro <- read.csv("./data/morphological/transformed/tm_algalonly_nomicro_sqrttransformed.csv")
tac_micro <- read.csv("./data/morphological/transformed/tac_algalonly_noanacylgreenalgae_sqrttransformed.csv")
nt_molec <- read.csv("./data/molecular/16s_nochimera_rarefied_95_FINAL.csv") %>% 
  filter(sample_type == "NT")
tm_molec <- read.csv("./data/molecular/16s_nochimera_rarefied_95_TAC_noanacyl.csv")
tac_molec <- read.csv("./data/molecular/16s_nochimera_rarefied_95_TM_nomicro.csv")
# note, the above molecular files have not been altered; wait till work through script 4b

# put all data into a list (keep grouping by category)
microscopy <- list(nt_micro, tm_micro, tac_micro) # nt_molec, tm_molec, tac_molec
names(microscopy) <- c("NT Algal Assemblages", "TM Algal Assemblages", "TAC Algal Assemblages")
# "NT Microbial Community", "TM Microbial Community", "TAC Microbial Community")

# add event no for date
#we have field dates for sampling distinguishing samples, however, let's change it
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
microscopy <- lapply(microscopy, function(x) add_event_no(x))

#### (2) Functions for Clustering ####

# do cluster analyses (returns plot)
# (arguments: data = wide data frame, data_start = col # of abundance data beginning,
# plot_title = title for plot, string)
cluster_analyses <- function(data, data_start, plot_title) {
  dist = vegdist(data[,data_start:ncol(data)], method = "bray")
  dist_trans = sqrt(dist)
  
  ward = hclust(dist, method = "ward.D2")
  plot(ward, main = plot_title, 
       labels = paste(data$site_reach, " date: ", data$event_no))
}

#### (3) Run analyses ####

# microscopy/algal assemblages data
for(i in 1:length(microscopy)) {
  cluster_analyses(microscopy[[i]], 7, names(microscopy)[i])
}

# notes:
# TAC - lots of river intermingling; single Salmon more related to Russian; 
# generally July and August on left and September on the right
# TM - Salmon and South Fork Eel distinct; however time periods for South Fork Eel are intermingled
# NT - distinct early year Salmon branch; some branches are distinctly Russian while
# other South Fork Eel but some intermingled; Salmon September closer to Russian

# subsetting by river for NT...
nt_micro_list <- split(nt_micro, nt_micro$site)
for(i in 1:length(nt_micro_list)) {
  cluster_analyses(nt_micro_list[[i]], 7, names(nt_micro_list)[i])
}
# Meh, I think NMDS is better for this stuff, definitely distinct September in Salmon
# distinct through time seems more apparent with Russian maybe then South Fork Eel?