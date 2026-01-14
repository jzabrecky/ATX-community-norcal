#### Main figure to show differences in assemblages associated with target taxa among rivers
### Jordan Zabrecky
## last edited: 01.13.2026

# This script creates a main figure to show differences in assemblages associated 
# with Microcoleus and Anabaena among rivers sampled

#### (1) Loading libraries & data ####

# load from analysis script
source("./code/4a_amongrivers_microscopy.R")

# rename data to not get overrun
data_morphological <- data
data_longer_morphological <- data_longer
NMDS_list_morphological <- NMDS_list
NMDS_plots_morphological <- NMDS_plots

# load from second analysis script
source("./code/4b_amongrivers_16s.R")

#### (2) Putting Figures Together ####

# list of desired figures:
# a- NMDS algal
# b- algal taxa groups (bar plot)
# c- NMDS bacterial
# d- alpha diversity (scatterplot)
# e- venn diagram
# f- functional groups (bar plot)

# creating string vector to iterate through sample types
sample_types = c("nt", "tac", "tm")


## (a) NMDS algal

# need to remake plot without envfit taxa
fig_a <- list()
for(i in sample_types) {
  fig_a[[i]] = makeNMDSplot(NMDS_list_morphological[[i]], FALSE, FALSE,
                            color = "site", shape = "site")
}

## (b) algal taxa bar plot

fig_b <- list()
for(i in sample_types) {
 fig_b[[i]] <- barplot_broader_plots[[i]]
}
lapply(fig_b, print)
# need to regroup NT
# check why there are some NA
# color palette?
# remove title

# maybe supplemental figure of each taxa individually?

## (c) NMDS bacterial

# create string vector to iterate through & empty list to hold figures
fig_c <- list()

for(i in sample_types) {
  fig_c[[i]] <- NMDS_plots_molecular[[i]]
}