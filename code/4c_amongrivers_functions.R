#### Comparing molecular predicted functional profiles data among rivers
### Jordan Zabrecky
## last edited: 05.07.2026

# This code compares normalized select orthologs/functions predicted via PICRUSt2-SC,
# from NT, TM, and TAC samples across rivers to answer Q1.
# Data is analyzed using Kruskal-Wallis Tests and visualizations

# MAKE SURE FUNCTIONAL GROUPINGS ARE GROUPED FIRST :)

#### (1) Loading libraries & data ####

# load libraries
lapply(c("tidyverse", "plyr", "vegan", "indicspecies"), require, character.only = T)

# load data (for selected functional groups)
nt <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select.csv") %>% 
  filter(sample_type == "NT")
tm <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tm_nomicro.csv")
tac <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tac_noanacyl.csv")

# put all dataframes into a list
data_select <- list(nt, tm, tac)
names(data_select) <- c("nt", "tm", "tac")

# load data for all orthologs
nt <- read.csv("./data/molecular/PICRUSt2_predicted_KO_all_all.csv") %>% 
  filter(sample_type == "NT")
tm <- read.csv("./data/molecular/PICRUSt2_predicted_KO_all_tm_nomicro.csv")
tac <- read.csv("./data/molecular/PICRUSt2_predicted_KO_all_tac_noanacyl.csv")

# put all dataframes into a list
data_all <- list(nt, tm, tac)
names(data_all) <- c("nt", "tm", "tac")

# souce functions for analyses
source("./code/supplemental_code/S4b_community_analyses_func.R")

#### (2) Plotting Select Genes #####

# make box plots!!!
boxplots <- lapply(data_select, function(x) {
  plot = ggplot(x, aes(x = site, y = predicted_gene_abundance, fill = site)) +
    geom_boxplot() +
    facet_wrap(~functional_grouping, scales = "free")
  print(plot)
  return(plot)
})

#### (3) Kruskal Wallis Tests ####

# run tests
kruskal_test_results <- lapply(data_select, function(x) {
  function_groups = unique(x$functional_grouping)
  results = data.frame(function_groups = NA,
                       kruskal_test = NA)
  
  for(i in function_groups) {
    results = rbind(results, data.frame(function_groups = i,
                                        kruskal_test = (kruskal.test(site~predicted_gene_abundance, data = (x %>% filter(functional_grouping == i))))$p.value))
  }
  
  return(results[-1,])
})

lapply(kruskal_test_results, function(x) x[which(x$kruskal_test < 0.05),])
# nitrification significantly different for non-target but nothing else

# save results
lapply(names(kruskal_test_results), function(x) {
  write.csv(kruskal_test_results[[x]], paste("./data/kruskal_wallis_results/Q1_", x, "_selectorthologs.csv"),
            row.names = FALSE)
})

#### (4) NMDS Plots ####

# will try two versions of data: (1) hellinger transform and (2) total predicted abundances
# for both all orthologs and our selected functional groups

## (a) entire functional profile (hellinger transformed)

# pivot wider & Hellinger-transform
data_all_wide_hellinger <- lapply(data_all, function(x){
  
  # pivot wider and then use decostand to relative
  y <- x %>% 
    # need to group same ko id's in same sample
    group_by(site, site_reach, ko_id, field_date, sample_type) %>% 
    dplyr::summarize(total_abundance = sum(predicted_gene_abundance)) %>% 
    pivot_wider(names_from = "ko_id", values_from = "total_abundance", values_fill = 0) %>% 
    mutate(field_date = mdy(field_date))
  y[,5:ncol(y)] <- decostand(y[,5:ncol(y)], method = "hellinger")
  
  return(y)
})

# get NMDS for each dataframe (sqrt-transformed!)
set.seed(1)
NMDS_list_all_hel <- lapply(data_all_wide_hellinger, function(x) getNMDSdata(x, 5, ASV = TRUE))

# making plots
NMDS_plots_all_hel <- lapply(NMDS_list_all_hel, function(x) makeNMDSplot(x, FALSE, FALSE, 
                                                         color = "site", shape = "month"))

lapply(NMDS_plots_all_hel, print)

## (b) entire functional profile (total abundances)

# pivot wider
data_all_wide <- lapply(data_all, function(x){
  
  # pivot wider and then use decostand to relative
  y <- x %>% 
    # need to group same ko id's in same sample
    group_by(site, site_reach, ko_id, field_date, sample_type) %>% 
    dplyr::summarize(total_abundance = sum(predicted_gene_abundance)) %>% 
    pivot_wider(names_from = "ko_id", values_from = "total_abundance", values_fill = 0) %>% 
    mutate(field_date = mdy(field_date))
  
  return(y)
})

# get NMDS for each dataframe
set.seed(1)
NMDS_list_all_total <- lapply(data_all_wide, function(x) getNMDSdata(x, 5, ASV = TRUE))

# making plots
NMDS_plots_all_total <- lapply(NMDS_list_all_total, function(x) makeNMDSplot(x, FALSE, FALSE, 
                                                                         color = "site", shape = "month"))
lapply(NMDS_plots_all_total, print)

## (c) select functions (hellinger transformed)

# pivot wider & Hellinger-transform
data_select_wide_hellinger <- lapply(data_select, function(x){
  
  # pivot wider and then use decostand to relative
  y <- x %>% 
    # need to group same ko id's in same sample
    group_by(site, site_reach, functional_grouping, field_date, sample_type) %>% 
    dplyr::summarize(total_abundance = sum(predicted_gene_abundance)) %>% 
    pivot_wider(names_from = "functional_grouping", values_from = "total_abundance", values_fill = 0) %>% 
    mutate(field_date = mdy(field_date))
  y[,5:ncol(y)] <- decostand(y[,5:ncol(y)], method = "hellinger")
  
  return(y)
})

# get NMDS for each dataframe (hellinger transformed!)
set.seed(1)
NMDS_list_select_rel <- lapply(data_select_wide_hellinger, function(x) getNMDSdata(x, 5))

# making plots
NMDS_plots_all_rel <- lapply(NMDS_list_select_rel, function(x) makeNMDSplot(x, TRUE, TRUE, 
                                                                 color = "site", shape = "month"))
lapply(NMDS_plots_all_rel, print)
# relativism may especially mask differences here? may not investigate this because of that!

## (d) select functions (total abundances!)

# pivot wider
data_select_wide <- lapply(data_select, function(x){
  
  # pivot wider and then use decostand to relative
  y <- x %>% 
    # need to group same ko id's in same sample
    group_by(site, site_reach, functional_grouping, field_date, sample_type) %>% 
    dplyr::summarize(total_abundance = sum(predicted_gene_abundance)) %>% 
    pivot_wider(names_from = "functional_grouping", values_from = "total_abundance", values_fill = 0) %>% 
    mutate(field_date = mdy(field_date))
  
  return(y)
})

# get NMDS for each dataframe (total abundance!)
set.seed(1)
NMDS_list_select_total <- lapply(data_select_wide, function(x) getNMDSdata(x, 5))

# making plots
NMDS_plots_select_total <- lapply(NMDS_list_select_total, function(x) makeNMDSplot(x, TRUE, TRUE, 
                                                                            color = "site", shape = "month"))

lapply(NMDS_plots_select_total, print)
# relativism may especially mask differences here? may not investigate this because of that!

#### (5) PERMANOVA ####

# empty table for permanova outputs
p_table <- data.frame(test = NA,
                      sample_type = NA,
                      all_or_select = NA,
                      data_transformation = NA,
                      p_value = NA,
                      F_stat = NA)

# run PERMANOVAs
set.seed(1)
permanovas_all_hell <- lapply(data_all_wide_hellinger, function(x) runPERMANOVA(data = x, 
                                                    start_col = 5, 
                                                    group = x$`site`))
permanovas_all_total <- lapply(data_all_wide, function(x) runPERMANOVA(data = x, 
                                                                                 start_col = 5, 
                                                                                 group = x$`site`))
permanovas_select_hell <- lapply(data_select_wide_hellinger, function(x) runPERMANOVA(data = x, 
                                                                           start_col = 5, 
                                                                           group = x$`site`))
permanovas_select_total <- lapply(data_select_wide, function(x) runPERMANOVA(data = x, 
                                                                            start_col = 5, 
                                                                            group = x$`site`))
# print and add test results to table
for(i in 1:length(permanovas_all_hell)) {
  
  # save stats to table
  # permanovs for all hellinger transformed
  p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                       sample_type = names(permanovas_all_hell[i]),
                                       all_or_select = "all",
                                       data_transformation = "hellinger",
                                       p_value = permanovas_all_hell[[i]]$`Pr(>F)`[1],
                                       F_stat = permanovas_all_hell[[i]]$`F`[1]))
  
  # permanovas for square-rooted total abundances
  p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                       sample_type = names(permanovas_all_total[i]),
                                       all_or_select = "all",
                                       data_transformation = "total_abundances",
                                       p_value = permanovas_all_total[[i]]$`Pr(>F)`[1],
                                       F_stat = permanovas_all_total[[i]]$`F`[1]))
  
  # permanovas for select hellinger transformed
  p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                       sample_type = names(permanovas_select_hell[i]),
                                       all_or_select = "select",
                                       data_transformation = "hellinger",
                                       p_value = permanovas_select_hell[[i]]$`Pr(>F)`[1],
                                       F_stat = permanovas_select_hell[[i]]$`F`[1]))
  
  # permanovas for square-rooted total abundances
  p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                       sample_type = names(permanovas_select_total[i]),
                                       all_or_select = "select",
                                       data_transformation = "total_abundaces",
                                       p_value = permanovas_select_total[[i]]$`Pr(>F)`[1],
                                       F_stat = permanovas_select_total[[i]]$`F`[1]))
}
view(p_table)
# Okay in order: 
# sample type- all total abun, all hellinger transformed, select total abun, select hellinger transformed
# yes for yes it is significant
# NT: yes, yes, yes, no
# TM: no, yes, no, yes
# TAC: no, no, no, yes
# so dependent on transformation especially for TM

# not sure what to think here ATM

#### (6) Indicator Species Analysis ####

# only doing for select functional groups because there are too many orthologs
# if we look at the full functional profile!

# (i) NT
set.seed(1)
nt_test <- multipatt(data_select_wide$nt[,5:ncol(data_select_wide$nt)], 
                     data_select_wide$nt$site, func = "r.g", control = how(nperm = 999))
summary(nt_test)
#write.csv(nt_test$sign, "./data/ISA_results/Q1_nt_microscopy.csv")
# RUS & SAL: higher nitrification
# SAL & SFE-M: higher pyridoxal

## (ii) TM
set.seed(1)
tm_test <- multipatt(data_select_wide$tm[,5:ncol(data_select_wide$tm)], 
                     data_select_wide$tm$site, func = "r.g", control = how(nperm = 999))
summary(tm_test)
# nothing

## (iii) TAC
set.seed(1)
tac_test <- multipatt(data_select_wide$tac[,5:ncol(data_select_wide$tac)], 
                     data_select_wide$tac$site, func = "r.g", control = how(nperm = 999))
summary(tac_test)
# SAL higher quite a few things