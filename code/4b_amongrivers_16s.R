#### Comparing microscopy data among rivers
### Jordan Zabrecky
## last edited: 11.16.2025

## This code compares 16s rRNA (microbial community) data from NT, TM, and TAC samples
## across rivers to answer Q1

#### (1) Loading libraries & data ####

# set seed for reproducibility
set.seed(2025)

# libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot"), require, character.only = T)

# read in files (note: two target taxa csvs- one with target taxa included, other with it excluded)
# doing as a list rather than one csv as each has a different # of columns
files <- list.files(path = "./data/molecular/", pattern = "16s")
data_original <- lapply(files, function(x) read.csv(paste("./data/molecular/", x, sep = "")))
names(data_original) <- files

# separate out sample types from first csv
data <- split(data_original$`16s_nochimera_rarefied_90_FINAL.csv`,
                   data_original$`16s_nochimera_rarefied_90_FINAL.csv`$sample_type)

# add two without target taxa
data[[4]] <- data_original$`16s_nochimera_rarefied_90_TAC_noanacyl.csv`
data[[5]] <- data_original$`16s_nochimera_rarefied_90_TM_nomicro.csv`
names(data)[4:5] <- c("TAC_no_AC", "TM_no_M")

#### (2) Data Transformations ####

## (a) sqrt-transform relative abundances

# histogram of raw percent (VERY right-skew)
lapply(data, function(x) hist(x$relative_abundance))

# square root transform for multivariate analyses since data is highly right-skewed
# also keeps values positive for Bray-Curtis

# sqrt-transform relative abundances
data <- lapply(data, function(x) x <- x %>%  mutate(sqrt_rel_abundance = sqrt(relative_abundance))
               %>% relocate(sqrt_rel_abundance, .after = "relative_abundance"))

# a bit better!
lapply(data, function(x) hist(x$sqrt_rel_abundance, breaks=seq(0,10,l=20)))

### CONSIDER MAKING BROADER GROUPS BUT THAT SHOULD BE SOLVED BY PHYLUM TREE


## (c) consider removing "rare" taxa

# consider removing taxa that compose less than 0.05% across all samples
rare_taxa <- lapply(data, function(x) x %>% 
                      dplyr::group_by(feature_ID) %>% 
                      dplyr::summarize(max = max(relative_abundance),
                                       min = min(relative_abundance),
                                       mean = mean(relative_abundance)) %>% 
                      filter(max < 0.05))

# separate list of dataframes with rare taxa removed
data_filtered <- list()
for(i in 1:length(data)) {
  data_filtered[[i]] <- data[[i]] %>% 
    filter(! feature_ID %in% rare_taxa[[i]]$feature_ID)
}
names(data_filtered) <- names(data)

# histogram of data with "rare" taxa removed (much better but still zero-inflation!)
lapply(data_filtered, function(x) hist(x$sqrt_rel_abundance))

## (c) make wide versions (for NMDS)!

data_wide_sqrt <- lapply(data, function(x) x %>% 
                           select(site_reach, site, field_date, sample_type, 
                                  sqrt_rel_abundance, feature_ID) %>% 
                           pivot_wider(names_from = feature_ID, values_from = sqrt_rel_abundance,
                                       values_fill = 0))
data_wide_unaltered <-lapply(data, function(x) x %>% 
                               select(site_reach, site, field_date, sample_type, 
                                      relative_abundance, feature_ID) %>% 
                               pivot_wider(names_from = feature_ID, values_from = relative_abundance,
                                           values_fill = 0))
data_wide_filtered <- lapply(data_filtered, function(x) x %>% 
                               select(site_reach, site, field_date, sample_type, 
                                      sqrt_rel_abundance, feature_ID) %>% 
                               pivot_wider(names_from = feature_ID, values_from = sqrt_rel_abundance,
                                           values_fill = 0))

#### (3) Functions for Analysis ####

# set universal plot theme
theme_set(theme_bw() + theme(panel.grid = element_blank(),
                             panel.border = element_rect(fill = NA, color = "black"),
                             legend.position = "right"))

#### TO-DO ####
# (1) NEED TO ADD IN A FILTER WHERE THOSE < THAN A CERTAIN AMOUNT GET FILTERED INTO OTHER
# (2) NEED TO CLEAN OUT CLASS
# (3) NEED TO FIGURE OUT Leptolyngbyales vs. Leptolyngbyaceae order within cyanobacteria
# (4) NEED TO, IN GENERAL, FINISH THE NAMING CLEANING
# (5) NITROSPINOTA vs. NITROSPIROTA in PHYLUMS

# bar plot at phylum level (input: data in long format)
barplot_phylum <- function(data) {
  plot = ggplot(data = data, aes(x = site, y = relative_abundance, fill = phylum)) +
    geom_bar(position = "fill", stat = "identity")
  
  return(plot)
}

# bar plot at order level (input: data in long format)
barplot_class <- function(data) {
  plot = ggplot(data = data, aes(x = site, y = relative_abundance, fill = class)) +
    geom_bar(position = "fill", stat = "identity")
  
  return(plot)
}

# bar plot orders within cyanobacteria phylum (input: data in long format)
barplot_cyano_order <- function(data) {
  plot = ggplot(data = data %>% filter(class == "Cyanobacteriia"),
                aes(x = site, y = relative_abundance, fill = order)) + 
    geom_bar(position = "fill", stat = "identity")
  
  return(plot)
}

# bar plot genus within cyanobacteria phylum (input: data in long format)
barplot_cyano_genus <- function(data) {
  plot = ggplot(data = data %>% filter(class == "Cyanobacteriia"),
                aes(x = site, y = relative_abundance, fill = genus)) + 
    geom_bar(position = "fill", stat = "identity")
  
  return(plot)
}

# summarize phylums (input: data in long format)
summarize_phylums <- function(data) {
  summary = data %>% 
    # have to do totals because each phylum has multiple ASVs
    # some of which are 0 which screws with mean
    dplyr::group_by(phylum, site) %>% 
    dplyr::summarize(total = sum(relative_abundance))
  
  # get total abundance for each site
  site_totals = summary %>% 
    dplyr::group_by(site) %>% 
    dplyr::summarize(site_total = sum(total))
  
  # calculate relative abundance of each for site
  final = left_join(summary, site_totals, by = c("site")) %>% 
    dplyr::mutate(relative_abun = (total / site_total) * 100) %>% 
    select(site, relative_abun)
  
  return(final)
}

# summarize classes
summarize_classes <- function(data) {
  summary = data %>% 
    # have to do totals because each class has multiple ASVs
    # some of which are 0 which screws with mean
    dplyr::group_by(class, site) %>% 
    dplyr::summarize(total = sum(relative_abundance))
  
  # get total abundance for each site
  site_totals = summary %>% 
    dplyr::group_by(site) %>% 
    dplyr::summarize(site_total = sum(total))
  
  # calculate relative abundance of each for site
  final = left_join(summary, site_totals, by = c("site")) %>% 
    dplyr::mutate(relative_abun = (total / site_total) * 100) %>% 
    select(site, relative_abun)
  
  return(final)
}

# summarizes orders within cyanobacteria
summarize_cyanoorders <- function(data) {
  summary = data %>% 
    filter(phylum == "Cyanobacteriia") %>% 
    # have to do totals because each order has multiple ASVs
    # some of which are 0 which screws with mean
    dplyr::group_by(order, site) %>% 
    dplyr::summarize(total = sum(relative_abundance))
  
  # get total abundance for each site
  site_totals = summary %>% 
    dplyr::group_by(site) %>% 
    dplyr::summarize(site_total = sum(total))
  
  # calculate relative abundance of each for site
  final = left_join(summary, site_totals, by = c("site")) %>% 
    dplyr::mutate(relative_abun = (total / site_total) * 100) %>% 
    select(site, relative_abun)
  
  return(final)
}

#### (4) 

# put bar plots into lists
barplot_phylum_plots <- lapply(data, function(x) barplot_phylum(x))
barplot_class_plots <- lapply(data, function(x) barplot_class(x))
barplot_cyanoorder_plots <- lapply(data, function(x) barplot_cyano_order(x))
barplot_cyanogenus_plots <- lapply(data, function(x) barplot_cyano_genus(x))
barplot_phylum_plots$TM_no_M
barplot_class_plots$TM
barplot_cyanoorder_plots$TM_no_M
barplot_cyanogenus_plots$TM_no_M
# duplicates on Leptolyngbyales?


#### DIVERSITY METRICS ####

#### (6) Q: What is dominant taxa across samples? ####

# get summaries for each taxa
summaries_phylum <- lapply(data, function(x) summarize_phylums(x))
summaries_class <- lapply(data, function(x) summarize_classes(x))
summaries_cyanoorder <- lapply(data, function(x) summarize_cyanoorders(x))


# why is this not looking the same as plots, is that what is happening on micro communities?

view(summaries_phylum$TM_no_M)
view(summaries_class$TM)
