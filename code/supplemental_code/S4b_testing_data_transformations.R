#### Playing with different data transformations used in community ecology
### Jordan Zabrecky
## last edited: 12.17.2025

## Do different data transformations changes our results?
## This script explores: (1) using un-transformed relative abundances
## (2) squareroot transformed relative abundances (Hellinger-transformation)
## (3) and Hellinger-transformed data with rare taxa removed

## In doing so, we also compare differences in results of targeted samples 
## with the targeted taxa remove, and, in the case of Anabaena samples (TAC), 
## how removing green algae--typically the substrate--alters results

## We look predominantly at PERMANOVAs and Species Indicator Analysis Results
## and then visually assess potential differences in NMDS

#### (1) Loading libraries & data ####

# set seed for reproducibility
set.seed(2025)

# libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot", "indicspecies"), require, character.only = T)

# microscopy data
files <- list.files(path = "./data/morphological/", pattern = ".csv")
data_wide <- lapply(files, function(x) read.csv(paste("./data/morphological/", x, sep = "")) %>% 
                      mutate(field_date = ymd(field_date),
                             month = month(field_date)) %>% # add month tag
                      relocate(month, .before = "field_date") %>% 
                      filter(year(field_date) == 2022)) # filter for 2022 only data
names(data_wide) <- files

# microbial data

# get functions from supplemental code script
source("./code/supplemental_code/S4a_community_analyses_func.R")

#### (2) Data Transformations ####

## (a) subset original relative abundance

# create new dataframe
algal_untransformed <- data_wide

# set number for start column of relative abundance data
start_col <- 6

# remove any taxa where there are no observations at any river for each sample type
# note, need to add 4 to match the subset taken for colSums
algal_notaxa <- lapply(algal_untransformed, function(x) colnames(x)[c((start_col - 1) 
                                                                      + which(colSums(x[,start_col:ncol(x)]) == 0))]) # no aphanothece ever recorded in TM samples

# remove these taxa (that are not present in any samples) from the dataframes
for(i in 1:length(data_wide)) {
  algal_untransformed[[i]] <- algal_untransformed[[i]] %>% 
    dplyr::select(!c(algal_notaxa[[i]]))
}

# histogram of raw relative abundances (highly right-skew)
lapply(algal_untransformed, function(x) hist(pivot_longer(data = x, cols = c(start_col:ncol(x)),
                                           names_to = "taxa", values_to = "percent")$percent))

## (b) square-root transform 
# (Hellinger transformation as it is relative abundance being transformed)

# create new dataframe
algal_hellinger <- algal_untransformed

# sqrt-transform
for(i in 1:length(algal_hellinger)) {
  algal_hellinger[[i]][,start_col:ncol(algal_hellinger[[i]])] <- 
    sqrt(algal_hellinger[[i]][,start_col:ncol(algal_hellinger[[i]])])
}

# histogram of transformed relative abundances (still right-skewed but better!)
lapply(algal_hellinger, function(x) hist(pivot_longer(data = x, cols = c(start_col:ncol(x)),
                                                          names_to = "taxa", values_to = "percent")$percent))

## (c) remove rare taxa

# function to figure out rare-taxa on untransformed data where max percent is cut-off for rare taxa
rare_taxa <- function(data, max_percent) {
  taxa <- data %>% pivot_longer(cols = c(start_col:ncol(data)), names_to = "taxa", values_to = "percent") %>% 
    dplyr::group_by(taxa) %>% 
    dplyr::summarize(max = max(percent),
                     min = min(percent),
                     mean = mean(percent)) %>% 
      filter(max < max_percent)
  # return string vector
  return(taxa$taxa)
}

# apply function
algal_raretaxa <- lapply(algal_hellinger, function(x) rare_taxa(x, max_percent = 1))

# create new list starting from hellinger-transformed data
algal_raretaxaremoved <- algal_hellinger

# remove rare taxa from dataframe
for(i in 1:length(algal_raretaxaremoved)) {
  algal_raretaxaremoved[[i]] <- algal_raretaxaremoved[[i]] %>% 
    select(!c(algal_raretaxa[[i]]))
}


# histogram of transformed relative abundances w/ rare taxa removed 
# (less skew but definitely still skewed!)
lapply(algal_raretaxaremoved, function(x) hist(pivot_longer(data = x, cols = c(start_col:ncol(x)),
                                                      names_to = "taxa", values_to = "percent")$percent))

#### (3) Q1: PERMANOVAs ####

# Is there a significant difference in assemblages among rivers?

# create summary table
Q1_algal_permanovas <- data.frame(data = NA,
                                  transformation = NA,
                                  significant = NA)

# run PERMANOVAs for Q1
for(i in 1:length(algal_untransformed)) {
  # make a temporary dataframe for all PERMANOVAs at index i
  temp = data.frame(data = c(names(algal_untransformed)[i], names(algal_hellinger)[i],
                             names(algal_raretaxaremoved)[i]),
                    transformation = c("none", "hellinger", "hellinger_w_raretaxaremoved"),
                    significant = c(runPERMANOVA(algal_untransformed[[i]], start_col, algal_untransformed[[i]]$site)$`Pr(>F)`[1], 
                                    runPERMANOVA(algal_hellinger[[i]], start_col, algal_hellinger[[i]]$site)$`Pr(>F)`[1],
                                    runPERMANOVA(algal_raretaxaremoved[[i]], start_col, algal_raretaxaremoved[[i]]$site)$`Pr(>F)`[1]))
  
  # add to existing dataframe
  Q1_algal_permanovas <- rbind(Q1_algal_permanovas, temp)
}
view(Q1_algal_permanovas)
# RESULTS: what is significant stays significant and is not stays not regardless of 
# transformation or what is removed (target taxa)

## TO-DO for microbial

#### (4) Q1: Species Indicator Analysis ####

# & species indicator test

#### (4) PERMANOVA & Species Indicator Q2 test ####

#### (5) Visualize NMDS ####

# run function to get data to make NMDS plots
algal_NMDS_list <- list()
algal_NMDS_list$`untransformed` <- lapply(algal_untransformed, function(x) getNMDSdata(x, start_col))
algal_NMDS_list$`hellinger` <- lapply(algal_hellinger, function(x) getNMDSdata(x, start_col))
algal_NMDS_list$`hellinger_raretaxaremoved` <- lapply(algal_raretaxaremoved, function(x) getNMDSdata(x, start_col))

# viewing plots against each other
for(i in 1:length(algal_untransformed)) {
  plots <- list()
  for(j in 1:length(algal_NMDS_list)) {
    plots[[j]] = makeNMDSplot(algal_NMDS_list[[j]][[i]], TRUE, TRUE, shape = "month", color = "site") +
      labs(title = paste(names(algal_untransformed)[i], names(algal_NMDS_list)[j], sep = " "))
  }
  print(plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol = 1))
}
# RESULTS: generally the same across transformations

#### (6) RDA for Q3? ####