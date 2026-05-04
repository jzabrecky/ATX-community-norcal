#### Comparing molecular 16s data among rivers
### Jordan Zabrecky
## last edited: 05.03.2026

# This code compares normalized 16s relative data from NT, TM, and TAC samples
# across rivers to answer Q1. First data is transformed (sqrt).
# Data is analyzed using NMDS and PERMANOVA. We also averaged across all samples
# from a river and created bar plots to visually compare average samples at each river

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

# Does diversity differ across rivers?
set.seed(1)
lapply(diversity, function(x) kruskal.test(shannon_diversity~site, data = x))
# significantly different for TM (p = 0.02), close but not for NT (0.05), and not for TAC (p = 0.45)

# Does evenness differ across rivers?
set.seed(1)
lapply(diversity, function(x) kruskal.test(evenness~site, data = x))
# significantly different for TM (p = 0.007), close but not for NT (0.007), and... close for TAC (p = 0.0501)


# What is the mean of each?
means_medians <- lapply(diversity, function(x) x %>% 
                          dplyr::group_by(site) %>% 
                          dplyr::summarize(mean = mean(shannon_diversity),
                                           median = median(shannon_diversity)))
view(means_medians$tac)

# save diversity calculations (RUN ONCE)
#lapply(names(diversity), function(x) write.csv(diversity[[x]], 
#                                               paste("./data/molecular/shannon_diversity/", x, "_diversity.csv", sep = ""),
#                                               row.names = FALSE))

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

#### (7) Q: What explains these differences? Species Indicator Analyses ####

# Previously looked at various groupings (orders & genuses w/in cyanobacteria, etc.)
# but have decided phylums and phylum-classes made most sense

## (a) phylums

phylums <- lapply(data_long, function(x) {
    # subset columns we care about and pivot wider
      y = x %>% 
        select(site_reach, site, field_date, sample_type, phylum, relative_abundance) %>% 
        dplyr::group_by(site_reach, site, field_date, sample_type, phylum) %>% 
        # remove multiples of phylums due to differing ASVs
        dplyr::summarize(relative_abundance = sum(relative_abundance)) %>% 
        pivot_wider(names_from = phylum, values_from = relative_abundance)
      # if there is a column NA, remove it
      if("NA" %in% colnames(y)) {
        y = y %>% 
          select(!c("NA"))
      }
      # replace NA (indicating that ASV was not present in sample) with 0
      y[,5:ncol(y)][is.na(y[,5:ncol(y)])] = 0
      # lastly, square-root transform the data
      y[,5:ncol(y)] <- sqrt(y[,5:ncol(y)])
      return(y)})

# (i) NT
set.seed(1)
nt_phylum_test <- summary(multipatt(phylums$nt[,5:ncol(phylums$nt)], phylums$nt$site, func = "r.g", control = how(nperm = 999)))
# lots of unique identified with RUS including ***: Actinobacteriota, Elusimicrobiota,
# Desulfobacteria, Planctomycetota, Verrucomicrobiota
# SAL: Deincoccota (***), Armatimonadota (**)
# SFE-M: Sumerlaeota
# RUS + SAL: Proteobacteria *
# SAL + SFE-M: Cyanobacteria ***
write.csv(nt_phylum_test$sign, "./data/ISA_results/Q1_nt_molecular_phylum.csv")

# (ii) TM
set.seed(1)
tm_phylum_test <- summary(multipatt(phylums$tm[,5:ncol(phylums$tm)], phylums$tm$site, func = "r.g", control = how(nperm = 999)))
# only identified for SFE-M including **: Desulfobacteria, Cyanobacteria, Verrucomicrobiota,
# Chloroflexiota
write.csv(tm_phylum_test$sign, "./data/ISA_results/Q1_tm_molecular_phylum.csv")

# (iii) TAC
set.seed(1)
tac_phylum_test <- summary(multipatt(phylums$tac[,5:ncol(phylums$tac)], phylums$tac$site, func = "r.g", control = how(nperm = 999)))
# only identified for SAL: Fibrobacterota *, Spirochaetota ** 
write.csv(tac_phylum_test$sign, "./data/ISA_results/Q1_tac_molecular_phylum.csv")

## (b) classes

classes <- lapply(data_long, function(x) {
  # subset columns we care about and pivot wider
  y = x %>% 
    select(site_reach, site, field_date, sample_type, phylum_class, relative_abundance) %>% 
    dplyr::group_by(site_reach, site, field_date, sample_type, phylum_class) %>% 
    # remove multiples of phylums due to differing ASVs
    dplyr::summarize(relative_abundance = sum(relative_abundance)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = phylum_class, values_from = relative_abundance)
  # if there is a column NA, remove it
  if("NA" %in% colnames(y)) {
    y = y %>% 
      select(!c("NA"))
  }
  # replace NA (indicating that ASV was not present in sample) with 0
  y[,5:ncol(y)][is.na(y[,5:ncol(y)])] = 0
  # lastly, square-root transform the data
  y[,5:ncol(y)] <- sqrt(y[,5:ncol(y)])
  return(y)})

# (i) NT
set.seed(1)
nt_class_test <- summary(multipatt(classes$nt[,5:ncol(classes$nt)], classes$nt$site, func = "r.g", control = how(nperm = 999)))
write.csv(nt_class_test$sign, "./data/ISA_results/Q1_molecular_nt_molecular_class.csv")
# lots for RUS (not listing)
# SAL: Fimbriimonadia ***, Deinococci ***
# SFE-M: Rhodothermia *
# RUS + SAL: Gammaproteobacteria **, Acidimicrobiia *
# RUS + SFE-M: including Desulfuromonadia **, Desulfovibrionia **
# SAL + SFE-M: Cyanobacteria ***, Gracilibacteria *

# (ii) TM
set.seed(1)
tm_class_test <- summary(multipatt(classes$tm[,5:ncol(classes$tm)], classes$tm$site, func = "r.g", control = how(nperm = 999)))
write.csv(tm_class_test$sign, "./data/ISA_results/Q1_molecular_tm_molecular_class.csv")
# SAL: Bacteroidia *
# SFE-M: many but no with ***-

# (iii) TAC
set.seed(1)
tac_class_test <- summary(multipatt(classes$tac[,5:ncol(classes$tac)], classes$tac$site, func = "r.g", control = how(nperm = 999)))
write.csv(tac_class_test$sign, "./data/ISA_results/Q1_molecular_tac_molecular_class.csv")
# only for SAL including Vampirivibrionia ***

#### (8) Misc. Questions ####

## (1) What are the top five most abundant phylum-class for each river?
summaries <- lapply(data_long, function(x) {
  y <- x %>% 
    dplyr::group_by(site, phylum_class) %>% 
    dplyr::summarize(mean = mean(relative_abundance),
                     sd = sd(relative_abundance))
})
view(summaries$nt)
## (NT)
# South Fork Eel: Cyanobacteria > Alphaproteobacteria > Verrucomicrobiae > Gammaproteobacteria
# Salmon: Cyanobacteria > Deinococci > Alphaproteobacteria > Gammaproteobacteria > Fimbriimonadia
# Russian: Cyanobacteria > Verrucomicrobiae > Gammaproteobacteria > Abditibacteria
