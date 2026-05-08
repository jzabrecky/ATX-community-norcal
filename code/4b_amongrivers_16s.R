#### Comparing molecular 16s data among rivers
### Jordan Zabrecky
## last edited: 05.07.2026

# This code compares normalized 16s relative data from NT, TM, and TAC samples
# across rivers to answer Q1. First data is transformed (sqrt).
# Data is analyzed using NMDS and PERMANOVA. We also averaged across all samples
# from a river and created bar plots to visually compare average samples at each river

# IDENTIFY OUTLIER SAMPLE - I THINK IT CREATES ISSUES IN ATX?!??

# okay need to confirm that the relative abundances sum to 100

#### (1) Loading libraries & data ####

## (a) load in data and libraries)

# load libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot", "indicspecies"), require, character.only = T)

# we will use the NT data, and TM excluding Microcoleus, and TAC excluding Anabaena
nt <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_FINAL.csv") %>% 
  filter(sample_type == "NT")
tm <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_TM_nomicro.csv")
tac <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_TAC_noanacyl.csv")

## (b) making a long format dataframe for bar graphs

# add into list
data_long <- list(nt, tm, tac)
names(data_long) <- c("nt", "tm", "tac")

# change format of field_data & add column for both phylum - class
data_long <- lapply(data_long, function(x) x <- x %>% 
                      mutate(field_date = mdy(field_date),
                             phylum_class = paste(phylum, " - ", class)))

## (c) making wider dataframe with ASVs & Hellinger-transforming

# pivot wider using ASVs for community matrix
data <- lapply(data_long, function(x) {
  # subset columns we care about and pivot wider
  y = x %>% 
    select(site_reach, site, field_date, sample_type, triplicate, feature_ID, picrust2_relative_abundance) %>% 
    dplyr::rename(relative_abundance = picrust2_relative_abundance) %>% 
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

## (d) adding broader categories for bar plots to long data format

# lastly, add in broader categories to data_long so we don't have several with barplots
# could do this quickly with forcats package, but since we are lumping, may make sense
# to have more customized categories
data_long <- lapply(data_long, function(x) x %>% dplyr::rename(relative_abundance = picrust2_relative_abundance))
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
                    # calculate average abundance for each (phylum - ) class
                    group_by(phylum_class) %>% 
                    dplyr::summarize(mean = mean(relative_abundance)) %>% 
                    # put into "Other" category if NA or % is less than #%
                    mutate(classes_cat = case_when(is.na(phylum_class) ~ "Other",
                                                  mean < 0.05  ~ "Other",
                                                  TRUE ~ phylum_class)) %>% 
                    select(phylum_class, classes_cat))
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
  data_long_broader[[i]] <- left_join(data_long_broader[[i]], classes[[i]], by = "phylum_class")
  data_long_broader[[i]] <- left_join(data_long_broader[[i]], cyano_order[[i]], by = "order")
  data_long_broader[[i]] <- left_join(data_long_broader[[i]], cyano_genus[[i]], by = "genus")
}

## (e) lastly, make a wide community matrix based on classes
# this is to see if differences are also significant with a broader category rather than just ASVs
data_wide_class <- lapply(data_long, function(x) {
  y <- x %>% 
    select(site, site_reach, field_date, sample_type, triplicate, phylum_class, qiime2_relative_abundance) %>%
    # sum for different ASVs in same class
    dplyr::group_by(site, site_reach, field_date, sample_type, triplicate, phylum_class) %>% 
    dplyr::summarize(qiime2_relative_abundance = sum(qiime2_relative_abundance)) %>% 
    ungroup() %>% 
    # now, can pivot_wider
    pivot_wider(names_from = phylum_class, values_from = qiime2_relative_abundance, values_fill = 0)
})

#### (2) Functions for Analyses ####

# load from supplemental script
source("./code/supplemental_code/S4b_community_analyses_func.R")

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
# could combine TAC categories more! will do so in figure script

# view plots (cyanobacteria categories)
for(i in 1:length(barplot_phylum_plots)) {
  print(barplot_cyanoorder_plots[[i]] + labs(title = titles[i]))
  print(barplot_cyanogenus_plots[[i]] + labs(title = titles[i]))
}

#### (4) Alpha Diversity Metrics ####

# calculate diversity for each dataframe
diversity <- lapply(data, function(x) {
  x = x %>% 
    mutate(shannon_diversity = calc_diversity(x, start_col),
           species_num = calc_speciesnum(x, start_col),
           evenness = shannon_diversity / log(species_num)) %>% 
    select(field_date, site_reach, site, shannon_diversity, species_num, evenness)})

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

# save test results
p_table <- data.frame(sample_type = NA,
                      diversity_kruskal = NA,
                      evenness_kruskal = NA)

for(i in names(diversity)) {
  
  # run tests
  set.seed(1)
  diversity_kruskal = (kruskal.test(shannon_diversity ~ site, data = diversity[[i]]))$p.value
  evenness_kruskal = (kruskal.test(evenness ~ site, data = diversity[[i]]))$p.value
  
  # add to table
  p_table <- rbind(p_table, data.frame(sample_type = i,
                                       diversity_kruskal = diversity_kruskal,
                                       evenness_kruskal = evenness_kruskal))
}

# What is the mean of each?
means_medians <- lapply(diversity, function(x) x %>% 
                          dplyr::group_by(site) %>% 
                          dplyr::summarize(mean = mean(shannon_diversity),
                                           median = median(shannon_diversity),
                                           sd = sd(shannon_diversity)))
view(means_medians$tac)

# save diversity calculations (RUN ONCE)
lapply(names(diversity), function(x) write.csv(diversity[[x]], 
                                               paste("./data/molecular/shannon_diversity/", x, "_diversity.csv", sep = ""),
                                               row.names = FALSE))

# also save test results
write.csv(p_table[-1,], "./data/kruskal_wallis_results/Q1_diversity.csv",
          row.names = FALSE)

#### (5) NMDS Plots ####

## (a) ASV-based (main focus)

# get NMDS for each dataframe (sqrt-transformed!)
set.seed(1) # set seed for reproducibility
NMDS_list <- lapply(data, function(x) getNMDSdata(x, start_col, ASV = TRUE))

# making plots
NMDS_plots <- lapply(NMDS_list, function(x) makeNMDSplot(x, FALSE, FALSE, 
                                                         color = "site", shape = "month"))

lapply(NMDS_plots, print)
# RESULT: TM and NT groups are visually distinct among rivers, but not for TAC

## (b) class-based
# get NMDS for each dataframe (sqrt-transformed!)
set.seed(1) # set seed for reproducibility
NMDS_list_class <- lapply(data_wide_class, function(x) getNMDSdata(x, start_col, ASV = TRUE))

# making plots
NMDS_plots_class <- lapply(NMDS_list_class, function(x) makeNMDSplot(x, FALSE, FALSE, 
                                                         color = "site", shape = "month"))

lapply(NMDS_plots_class, print)
# RESULT: NT definitely plots closer, TM still different seeming, maybe TAC slightly further which is weird

#### (6) Q: Are communities from each river significantly different? (PERMANOVA) ####

# empty table for permanova outputs
p_table <- data.frame(test = NA,
                      subtype = NA,
                      sample_type = NA,
                      p_value = NA,
                      F_stat = NA)

## (a) ASV based (main analysis)

# run PERMANOVAs
set.seed(1)
permanovas <- lapply(data, function(x) runPERMANOVA(x, start_col, group = x$`site`))

# print and add test results to table
for(i in 1:length(permanovas)) {
  
  # print test results to console
  print(names(permanovas[i]))
  print(permanovas[[i]])
  
  # save stats to table
  p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                       subtype = "ASVs",
                                       sample_type = names(permanovas[i]),
                                       p_value = permanovas[[i]]$`Pr(>F)`[1],
                                       F_stat = permanovas[[i]]$`F`[1]))
}
# RESULTS: significant difference for all but not convinced about the TAC
# as visually they plotted on top of each other but had different dispersion

# check dispersion to see if that influences results
for(i in 1:length(data)) {
  set.seed(1)
  anova = anova(betadisper(vegdist(data[[i]][,start_col:ncol(data[[i]])], method = "bray"), 
                           data[[i]]$site))
  
  
  # print results
  print(names(data)[i])
  print(anova)
  
  # add results table
  p_table <- rbind(p_table, data.frame(test = "PERMDISP",
                                       subtype = "ASVs",
                                       sample_type = names(data)[i],
                                       p_value = anova$`Pr(>F)`[1],
                                       F_stat = anova$`F value`[1]))
}
# dispersion not significantly different for NT, a little for TM (*) and very for
# TAC (***)

# due to these results, I am not convinced TAC is necessarily different....

## what if we remove the single salmon sample?
test <- data$tac %>% filter(site != "SAL")
set.seed(1)
print(anova(betadisper(vegdist(test[start_col:ncol(test)], method = "bray"), 
                       test$site)))
# still significantly different, but less so (**)
set.seed(1)
runPERMANOVA(test, start_col, group = test$`site`) # still very significant here (***)

## (b) class based (secondary analysis)

# run PERMANOVAs
set.seed(1)
permanovas_class <- lapply(data_wide_class, function(x) runPERMANOVA(x, start_col, group = x$`site`))

# print and add test results to table
for(i in 1:length(permanovas_class)) {
  
  # print test results to console
  print(names(permanovas_class[i]))
  print(permanovas_class[[i]])
  
  # save stats to table
  p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                       subtype = "classes",
                                       sample_type = names(permanovas_class[i]),
                                       p_value = permanovas_class[[i]]$`Pr(>F)`[1],
                                       F_stat = permanovas_class[[i]]$`F`[1]))
}
# RESULTS: NT and TAC different **, but not TAC
# this confirms my suspicions

# check dispersion to see if that influences results
for(i in 1:length(data_wide_class)) {
  set.seed(1)
  anova_classes = anova(betadisper(vegdist(data_wide_class[[i]][,start_col:ncol(data_wide_class[[i]])], method = "bray"), 
                           data_wide_class[[i]]$site))
  
  
  # print results
  print(names(data_wide_class)[i])
  print(anova_classes)
  
  # add results table
  p_table <- rbind(p_table, data.frame(test = "PERMDISP",
                                       subtype = "classes",
                                       sample_type = names(data)[i],
                                       p_value = anova_classes$`Pr(>F)`[1],
                                       F_stat = anova_classes$`F value`[1]))
}
# NT not significant, TM is *, and TAC is *

# save tests
write.csv(p_table[-1,], "./data/PERMANOVA_results/Q1_molecular.csv", row.names = FALSE)

#### (8) Misc. Questions ####

## (1) What are the top five most abundant phylum-class for each river?
summaries_phyla <- lapply(data_long, function(x) {
  y <- x %>% 
    # first group different ASVs together by phylum class
    dplyr::group_by(site_reach, site, field_date, sample_type, phylum) %>% 
    dplyr::summarize(total_abun = sum(relative_abundance)) %>% 
    ungroup()
  
  z <- y %>% 
    dplyr::group_by(site, sample_type, phylum) %>% 
    dplyr::summarize(mean = mean(total_abun),
                     sd = sd(total_abun))

  return(z)
})
view(summaries_phyla$nt)
summaries_class <- lapply(data_long, function(x) {
  y <- x %>% 
    # first group different ASVs together by phylum class
    dplyr::group_by(site_reach, site, field_date, sample_type, phylum_class) %>% 
    dplyr::summarize(total_abun = sum(relative_abundance)) %>% 
    ungroup()
  
  z <- y %>% 
    dplyr::group_by(site, sample_type, phylum_class) %>% 
    dplyr::summarize(mean = mean(total_abun),
                     sd = sd(total_abun))
  
  return(z)
})
view(summaries_class$nt)

## (2) How about for the sample type together (specifically for TAC)?
sample_summaries <- lapply(data_long, function(x) {
  y <- x %>% 
    # first group different ASVs together by phylum class
    dplyr::group_by(site_reach, site, field_date, sample_type, phylum) %>% 
    dplyr::summarize(total_abun = sum(relative_abundance)) %>% 
    ungroup()
  
  z <- y %>% 
    dplyr::group_by(sample_type, phylum) %>% 
    dplyr::summarize(mean = mean(total_abun),
                     sd = sd(total_abun))
  
  return(z)
})
view(sample_summaries$tac)

## (2) How much of the sample is the target taxa for target samples?
tm_w_target <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_FINAL.csv") %>% 
  filter(sample_type == "TM") %>% 
  filter(genus == "Microcoleus") %>% 
  group_by(site, site_reach, sample_type, field_date) %>% 
  dplyr::summarize(total = sum(picrust2_relative_abundance))

# need to figure out how much
min(tm_w_target$total) # 0.12% (quite low! is SFE-M-3 9/6/2022 - will look at it more later)
max(tm_w_target$total) # 92.02%
mean(tm_w_target$total) # 53.6%
sd(tm_w_target$total) # 30.6%

curious_sample <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_FINAL.csv") %>% 
  filter(site_reach == "SFE-M-3" & field_date == "9/6/2022" & sample_type == "TM")
# seems like mostly proteobacteria

tac_w_target <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_FINAL.csv") %>%
  filter(sample_type == "TAC") %>% 
  filter(genus %in% c("Anabaena","Cylindrospermum","Trichormus","Cylindrospermopsis")) %>% 
  group_by(site, site_reach, sample_type, field_date) %>% 
  dplyr::summarize(total = sum(picrust2_relative_abundance))
min(tac_w_target$total) # 0.04
max(tac_w_target$total) # 29.9
mean(tac_w_target$total) # 8.7%
sd(tac_w_target$total) # 8.6%
