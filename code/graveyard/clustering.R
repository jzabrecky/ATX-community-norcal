#### Side fun with clustering
### Jordan Zabrecky
## last edited: 12.16.2025

# This code compares clusters communities using the "cluster" package
# Instead, for analyses with are just using PERMANOVA and NMDS

#### (1) Loading libraries & data ####

# set seed for reproducibility
set.seed(2025)

# load libraries
lapply(c("tidyverse", "vegan", "cluster"), require, character.only = T)

# loading files individually because I only want certain ones
nt_micro <- read.csv("./data/morphological/transformed/nt_algalonly_sqrttransformed.csv")
tm_micro <- read.csv("./data/morphological/transformed/tm_algalonly_nomicro_sqrttransformed.csv")
tac_micro <- read.csv("./data/morphological/transformed/tac_algalonly_noanacylgreenalgae_sqrttransformed.csv")
nt_molec <- read.csv("./data/molecular/transformed/16s_nochimera_rarefied_95_NT_sqrttransformed.csv")
tac_molec <- read.csv("./data/molecular/transformed/16s_nochimera_rarefied_95_TAC_noanacyl_sqrttransformed.csv")
tm_molec <- read.csv("./data/molecular/transformed/16s_nochimera_rarefied_95_TM_nomicro_sqrttransformed.csv")
# note, the above molecular files have not been altered; wait till work through script 4b

# put all data into a list (keep grouping by category)
microscopy <- list(nt_micro, tm_micro, tac_micro) # nt_molec, tm_molec, tac_molec
names(microscopy) <- c("NT Algal Assemblages", "TM Algal Assemblages", "TAC Algal Assemblages")
molecular <- list(nt_molec, tm_molec, tac_molec)
names(molecular) <- c("NT Microbial Assemblages", "TM Microbial Assemblages", "TAC Microbial Assemblages")

# change format of molecular data
molecular <- lapply(molecular, function(x) x <- x %>% mutate(field_date = mdy(field_date)))

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
molecular <- lapply(molecular, function(x) add_event_no(x))

#### (2) Functions for Clustering ####

# do cluster analyses (returns plot) using "cluster"
# (arguments: data = wide data frame, data_start = col # of abundance data beginning,
# plot_title = title for plot, string)
cluster_analyses <- function(data, data_start, plot_title) {
  dist = vegdist(data[,data_start:ncol(data)], method = "bray")
  dist_trans = sqrt(dist)
  
  ward = hclust(dist, method = "ward.D2")
  print(plot(ward, main = plot_title, 
             labels = paste(data$site_reach, " sampling no. ", data$event_no)))
  return(ward)
}

# function to identify optimal number of clusters
opt_clusters <- function(data, data_start, cluster_object) {
  tests = numeric(nrow(data)) # empty vector for number of groups
  
  # get our sqrt bray-curtis transformed distances (again)
  dist = vegdist(data[,data_start:ncol(data)], method = "bray")
  dist_trans = sqrt(dist)
  
  for(i in 2:(nrow(data) - 1)) { # start with two groups and end with all minus 1
    silhouette = silhouette(cutree(cluster_object, k = i), dist_trans)
    tests[i] <- summary(silhouette)$avg.width
  }
  
  return(which.max(tests))
}

#### (3) Microscopy Clusters ####

## (a)  broad microscopy/algal assemblages data
for(i in 1:length(microscopy)) {
  cluster_analyses(microscopy[[i]], 7, names(microscopy)[i])
}
# notes:
# TAC - lots of river intermingling; single Salmon more related to Russian; 
# generally July and August on left and September on the right
# TM - Salmon and South Fork Eel distinct; however time periods for South Fork Eel are intermingled
# NT - distinct early year Salmon branch; some branches are distinctly Russian while
# other South Fork Eel but some intermingled; Salmon September closer to Russian

## (b) subsetting by river for NT...
nt_micro_list <- split(microscopy$`NT Algal Assemblages`, microscopy$`NT Algal Assemblages`$site)
for(i in 1:length(nt_micro_list)) {
  nt_cluster = cluster_analyses(nt_micro_list[[i]], 7, names(nt_micro_list)[i])
  print(paste(nt_micro_list[[i]]$site[1], ": ", opt_clusters(nt_micro_list[[i]], nt_cluster),
              sep = ""))
}

# subsetting by river for TM...
tm_micro_list <- split(microscopy$`TM Algal Assemblages`, microscopy$`TM Algal Assemblages`$site)
for(i in 1:length(tm_micro_list)) {
  cluster_analyses(tm_micro_list[[i]], 7, names(tm_micro_list)[i])
}

# subsetting by river for TAC...
tac_micro_list <- split(microscopy$`TAC Algal Assemblages`, microscopy$`TAC Algal Assemblages`$site)
tac_micro_list <- tac_micro_list[-2]
for(i in 1:length(tac_micro_list)) {
  cluster_analyses(tac_micro_list[[i]], 7, names(tac_micro_list)[i])
}

#### (4) Microbial 16s Clusters ####

## (a)  broad microscopy/algal assemblages data
for(i in 1:length(molecular)) {
  cluster_analyses(molecular[[i]], 7, names(molecular)[i])
}
# TAC- SFE and RUS all separate, one SAL TAC matched to RUS
# TM- sites distinct
# NT- sites also more distinct
# much more distinct than microscopy data!

## (b) subsetting by river for NT...
nt_molec_list <- split(molecular$`NT Microbial Assemblages`, molecular$`NT Microbial Assemblages`$site)
for(i in 1:length(nt_molec_list)) {
  nt_cluster = cluster_analyses(nt_molec_list[[i]], 7, names(nt_molec_list)[i])
  print(paste(nt_molec_list[[i]]$site[1], ": ", opt_clusters(nt_molec_list[[i]], 7, nt_cluster),
              sep = ""))
}

# subsetting by river for TM...
tm_molec_list <- split(molecular$`TM Microbial Assemblages`, molecular$`TM Microbial Assemblages`$site)
for(i in 1:length(tm_molec_list)) {
  tm_cluster = cluster_analyses(tm_molec_list[[i]], 7, names(tm_molec_list)[i])
  print(paste(tm_molec_list[[i]]$site[1], ": ", opt_clusters(tm_molec_list[[i]], 7, tm_cluster),
              sep = ""))
}

# subsetting by river for TAC...
tac_molec_list <- split(molecular$`TAC Microbial Assemblages`, molecular$`TAC Microbial Assemblages`$site)
tac_molec_list <- tac_molec_list[-2]
for(i in 1:length(tac_molec_list)) {
  tac_cluster = cluster_analyses(tac_molec_list[[i]], 7, names(tac_molec_list)[i])
  print(paste(tac_molec_list[[i]]$site[1], ": ", opt_clusters(tac_molec_list[[i]], 7, tac_cluster),
              sep = ""))
}
