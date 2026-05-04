### assessing evenness and diversity over time

# calculate evenness and diversity over # of species
# add this into main analyses???

# feel unsure about algal/morphological as taxa not on same level for that

#### algal stuff (as microbial is on current analyses)

# load in algal data
nt <- read.csv("./data/morphological/nt_algalonly.csv") %>% 
  filter(year(ymd(field_date)) == 2022)
tm <- read.csv("./data/morphological/tm_algalonly_nomicro.csv") %>% 
  filter(year(ymd(field_date)) == 2022)
tac <- read.csv("./data/morphological/tac_algalonly_noanacylgreenalgae.csv") %>% 
  filter(year(ymd(field_date)) == 2022)
  
# load from supplemental script
source("./code/supplemental_code/S4b_community_analyses_func.R")


microscopy <- list(nt, tm, tac)
names(microscopy) <- c("nt", "tm", "tac")

microscopy <- lapply(microscopy, add_event_no)

library(vegan)

# function to calculate diversity of ASVs
# @param data is dataframe of wide abundances
# @param start_col is column where abundance data starts
calc_diversity <- function(data, start_col) {
  vector = diversity(data[,start_col:ncol(data)])
}

# another for species number
# @param data is dataframe of wide abundances
# @param start_col is column where abundance data starts
calc_speciesnum <- function(data, start_col) {
  vector = specnumber(data[,start_col:ncol(data)])
}

diversity <- lapply(microscopy, function(x) {
  x = x %>% 
    mutate(shannon_diversity = calc_diversity(x, 7),
           species_num = calc_speciesnum(x, 7),
           evenness = shannon_diversity / log(species_num)) %>% 
    select(field_date, site_reach, site, event_no, shannon_diversity, species_num, evenness)})

# Q1 evenness and diversity across rivers

# plot diversity as boxplots
for(i in 1:length(diversity)) {
  boxplot = ggplot(data = diversity[[i]], aes(x = site, y = shannon_diversity, fill = site)) +
    geom_boxplot()
  print(boxplot)
}

# plot evenness as boxplots
for(i in 1:length(diversity)) {
  boxplot = ggplot(data = diversity[[i]], aes(x = site, y = evenness, fill = site)) +
    geom_boxplot()
  print(boxplot)
}

#### Q2 through time ####
for(i in 1:length(diversity)) {
  plot = ggplot(data = diversity[[i]], aes(x = as.factor(event_no), y = shannon_diversity,
                                           fill = site)) +
    geom_boxplot() +
    facet_wrap(~site, ncol = 1)
  print(plot)
}

for(i in 1:length(diversity)) {
  plot = ggplot(data = diversity[[i]], aes(x = as.factor(event_no), y = evenness,
                                           fill = site)) +
    geom_boxplot() +
    facet_wrap(~site, ncol = 1)
  print(plot)
}

#### Q3 w/ varying ATX concentrations ####

# load in ATX
