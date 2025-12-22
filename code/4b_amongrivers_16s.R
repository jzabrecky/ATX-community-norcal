#### Comparing microscopy data among rivers
### Jordan Zabrecky
## last edited: 12.22.2025

# This code compares 16s rRNA (microbial community) data from NT, TM, and TAC samples
# across rivers to answer Q1

# FOLLOW script 4a but also want to look at diversity metrics

#### (1) Loading libraries & data ####

# set seed for reproducibility
set.seed(2025)

## libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot", "indicspecies"), require, character.only = T)

# we will use the NT data, and TM excluding Microcoleus, and TAC excluding Anabaena and Green Algae
nt <- read.csv("./data/molecular/16s_nochimera_rarefied_95_FINAL.csv") %>% 
  filter(sample_type == "NT")
tm <- read.csv("./data/molecular/16s_nochimera_rarefied_95_TM_nomicro.csv")
tac <- read.csv("./data/molecular/16s_nochimera_rarefied_95_TAC_noanacyl.csv")

# add into list
data_long <- list(nt, tm, tac)
names(data_long) <- c("nt", "tm", "tac")

# pivot wider using ASVs for community matrix
data <- lapply(data_long, function(x) {
  # subset columns we care about and pivot wider
  y = x %>% 
    select(site_reach, site, field_date, sample_type, triplicate, feature_ID, relative_abundance) %>% 
    pivot_wider(names_from = feature_ID, values_from = relative_abundance)
  # replace NA (indicating that ASV was not present in sample) with 0
  y[,6:ncol(y)][is.na(y[,6:ncol(y)])] = 0
  return(y)})

# set start col for community matrix
start_col <- 6  

# see our exploration with data transformation in another script, "S4a_testing_data_transformations.R"
# decided on square-root transformation on the relative abundances (Hellinger transformation)
for(i in 1:length(data)) {
  data[[i]][,start_col:ncol(data[[i]])] <- sqrt(data[[i]][,start_col:ncol(data[[i]])])
}


# save this data to use in future scripts (RUN ONCE)
#write.csv(data$nt, "./data/molecular/transformed/16s_nochimera_rarefied_95_NT_sqrttransformed.csv",
#          row.names = FALSE)
#write.csv(data$tm, "./data/molecular/transformed/16s_nochimera_rarefied_95_TM_nomicro_sqrttransformed.csv",
#          row.names = FALSE)
#write.csv(data$tac, "./data/molecular/transformed/16s_nochimera_rarefied_95_TAC_noanacyl_sqrttransformed.csv",
#          row.names = FALSE)

# lastly, add in broader categories to data_long so we don't have several with barplots
# could do this quickly with forcats package, but since we are lumping, may make sense
# to have more customized categories
phylums <- lapply(data_long, function(x) x %>%
                    # calculate average abundance for each phylum
                    group_by(phylum) %>% 
                    dplyr::summarize(mean = mean(relative_abundance)) %>% 
                    # put into "Other" category if NA or % is less than #%
                    mutate(phylum_cat = case_when(is.na(phylum) ~ "Other",
                                                  mean < 0.03  ~ "Other",
                                                  TRUE ~ phylum)) %>% 
                    select(phylum, phylum_cat))
classes <- lapply(data_long, function(x) x %>%
                    # calculate average abundance for each class
                    group_by(class) %>% 
                    dplyr::summarize(mean = mean(relative_abundance)) %>% 
                    # put into "Other" category if NA or % is less than #%
                    mutate(classes_cat = case_when(is.na(class) ~ "Other",
                                                  mean < 0.05  ~ "Other",
                                                  TRUE ~ class)) %>% 
                    select(class, classes_cat))
cyano_order <- lapply(data_long, function(x) x %>%
                        # filter for cyanobacteria phylum
                        filter(phylum == "Cyanobacteria") %>% 
                        # calculate average abundance for each order
                        group_by(order) %>% 
                        dplyr::summarize(mean = mean(relative_abundance)) %>% 
                        # put into "Other" category if NA or % is less than #%
                        mutate(cyano_order_cat = case_when(is.na(order) ~ "Other",
                                                      mean < 0.035  ~ "Other",
                                                      TRUE ~ order)) %>% 
                        select(order, cyano_order_cat))
cyano_genus <- lapply(data_long, function(x) x %>%
                        # filter for cyanobacteria phylum
                        filter(phylum == "Cyanobacteria") %>% 
                        # calculate average abundance for each genus
                        group_by(genus) %>% 
                        dplyr::summarize(mean = mean(relative_abundance)) %>% 
                        # put into "Other" category if NA or % is less than #%
                        mutate(cyano_genus_cat = case_when(is.na(genus) ~ "Other",
                                                      mean < 0.11  ~ "Other",
                                                      TRUE ~ genus)) %>% 
                        select(genus, cyano_genus_cat))
# add these categories into our long dataframes
data_long_broader <- data_long
for(i in 1:length(data_long)) {
  data_long_broader[[i]] <- left_join(data_long_broader[[i]], phylums[[i]], by = "phylum")
  data_long_broader[[i]] <- left_join(data_long_broader[[i]], classes[[i]], by = "class")
  data_long_broader[[i]] <- left_join(data_long_broader[[i]], cyano_order[[i]], by = "order")
  data_long_broader[[i]] <- left_join(data_long_broader[[i]], cyano_genus[[i]], by = "genus")
}


#### (2) Functions for Analyses ####

# load from supplemental script
source("./code/supplemental_code/S4a_community_analyses_func.R")
source("./code/supplemental_code/S4c_barplot_func.R")

# function to calculate diversity of ASVs
# @param data is dataframe of wide abundances
# @param start_col is column where abundance data starts
calc_diversity <- function(data, start_col) {
  vector = diversity(data[,start_col:ncol(data)])
}

#### (3) Relative Abundance Bar Plots ####

# put bar plots into lists
# (note if issues- something is masked; restart R)
barplot_phylum_plots <- lapply(data_long_broader, function(x) barplot(x, x = "site", y = "relative_abundance",
                                                              fill = "phylum_cat"))
barplot_class_plots <- lapply(data_long_broader, function(x) barplot(x, x = "site", y = "relative_abundance",
                                                              fill = "classes_cat"))
barplot_cyanoorder_plots <- lapply(data_long_broader, function(x) barplot(x %>% filter(phylum == "Cyanobacteria"),
                                                             x = "site", y = "relative_abundance",
                                                              fill = "cyano_order_cat"))
barplot_cyanogenus_plots <- lapply(data_long_broader, function(x) barplot(x %>% filter(phylum == "Cyanobacteria"),
                                                                  x = "site", y = "relative_abundance",
                                                                  fill = "cyano_genus_cat"))

# titles for plots
titles <- c("Non-Target Samples", 
            "Microcoleus Samples (excluding M)",
            "Anabaena/Cylindrospermum Samples (excluding AC & GA)")

# view plots (all bacteria categories)
for(i in 1:length(barplot_phylum_plots)) {
  print(barplot_phylum_plots[[i]] + labs(title = titles[i]))
  print(barplot_class_plots[[i]] + labs(title = titles[i]))
}
# could combine TAC categories more!

# view plots (cyanobacteria categories)
for(i in 1:length(barplot_phylum_plots)) {
  print(barplot_cyanoorder_plots[[i]] + labs(title = titles[i]))
  print(barplot_cyanogenus_plots[[i]] + labs(title = titles[i]))
}
# could also probably customize more here

#### (4) Alpha Diversity Metrics ####

# calculate diversity for each dataframe
diversity <- lapply(data, function(x) {
  x = x %>% 
    mutate(shannon_diversity = calc_diversity(x, start_col)) %>% 
    select(field_date, site_reach, site, shannon_diversity)})

# plot diversity as boxplots
for(i in 1:length(diversity)) {
  boxplot = ggplot(data = diversity[[i]], aes(x = site, y = shannon_diversity, fill = site)) +
    geom_boxplot()
  print(boxplot)
}

# Does diversity differ across rivers?
lapply(diversity, function(x) kruskal.test(shannon_diversity~site, data = x))
# not significantly different for any group but close for TM (p = 0.06)

#### (5) NMDS Plots ####

# NEED TO TRY THIS ON REMOTE DESKTOP

# get NMDS for each dataframe (sqrt-transformed!)
NMDS_list <- lapply(data, function(x) getNMDSdata(x, start_col))

# making plots
NMDS_plots <- lapply(NMDS_list, function(x) makeNMDSplot(x, TRUE, TRUE, 
                                                         color = "site", shape = "month"))

lapply(NMDS_plots, print)
# RESULT: TM and NT groups are visually distinct among rivers, but not for TAC

?vegan::richness()

diversity(data$nt[,6:ncol(data$nt)])







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
