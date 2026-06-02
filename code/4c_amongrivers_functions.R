#### Comparing predicted functional profiles data among rivers
### Jordan Zabrecky
## last edited: 06.01.2026

# This code compares both the full predicted function and select orthologs/functions 
# obtained from PICRUSt2 for NT, TM, and TAC samples across rivers to answer Q1.
# Data is analyzed using PERMANOVA and PCoA (for full predicted function) 
# and Kruskal-Wallis Tests and boxplot visualizations (for select genes)

#### (1) Loading libraries & data ####

# load libraries
lapply(c("tidyverse", "plyr", "vegan", "indicspecies", "dunn.test"), require, character.only = T)

# load data for all orthologs
nt <- read.csv("./data/molecular/PICRUSt2_predicted_KO_all_all.csv") %>% 
  filter(sample_type == "NT")
tm <- read.csv("./data/molecular/PICRUSt2_predicted_KO_all_tm_nomicro.csv")
tac <- read.csv("./data/molecular/PICRUSt2_predicted_KO_all_tac_noanacyl.csv")

# put all dataframes into a list
data_all <- list(nt, tm, tac)
names(data_all) <- c("nt", "tm", "tac")

# load data (for selected functional groups)
nt <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select.csv") %>% 
  filter(sample_type == "NT")
tm <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tm_nomicro.csv")
tac <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tac_noanacyl.csv")

# put all dataframes into a list
data_select <- list(nt, tm, tac)
names(data_select) <- c("nt", "tm", "tac")

# souce functions for analyses
source("./code/supplemental_code/S4b_community_analyses_func.R")

# need to group different orthologs in same functional group for selected genes
data_select <- lapply(data_select, function(x) {
  y = x %>% 
    dplyr::group_by(site, site_reach, field_date, sample_type, functional_grouping) %>% 
    dplyr::summarize(predicted_gene_abundance = sum(predicted_gene_abundance))
})

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

# confirm distributions are not normal to use with ANOVA
shapiro.test((data_select$nt %>% filter(functional_grouping == "nitrogen_fixation" & site == "SFE-M"))$predicted_gene_abundance)
# not normal

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
view(kruskal_test_results) # nothing signficiantly differed

# want to double-check numbers: comparing nitrogen fixation for NT
kruskal.test(site~predicted_gene_abundance, data = data_select$nt %>% filter(functional_grouping == "nitrogen_fixation"))
# 0.4726 which is same as in table
# is phosphatase also the same???
kruskal.test(site~predicted_gene_abundance, data = data_select$nt %>% filter(functional_grouping == "phosphatase_transporters"))
# yes, weird

# save results
lapply(names(kruskal_test_results), function(x) {
  write.csv(kruskal_test_results[[x]], paste("./data/kruskal_wallis_results/Q1_", x, "_selectorthologs.csv"),
            row.names = FALSE)
})

#### (4) PCoA Plots ####

# as decided in the supplemental script, "S4e_testing_data_sqrtformations_predgenes.R",
# we will square-root the predicted gene abundances to minimize impact of high gene counts

# pivot wider & sqrt-transform
data_sqrt <- lapply(data_all, function(x){
  
  # pivot wider and then use decostand to relative
  y <- x %>% 
    # need to group same ko id's in same sample
    group_by(site, site_reach, ko_id, field_date, sample_type) %>% 
    dplyr::summarize(total_abundance = sum(predicted_gene_abundance)) %>% 
    pivot_wider(names_from = "ko_id", values_from = "total_abundance", values_fill = 0) %>% 
    mutate(field_date = mdy(field_date)) %>% 
    ungroup()
  y[,5:ncol(y)] <- sqrt(y[,5:ncol(y)])
  
  return(y)
})

# RUN ONCE: save these in the transformed folder to use in future script
#write.csv(data_sqrt$nt, "./data/molecular/transformed/PICRUSt2_predicted_KO_all_nt_sqrttransformed.csv",
#          row.names = FALSE)
#write.csv(data_sqrt$tm, "./data/molecular/transformed/PICRUSt2_predicted_KO_all_tm_nomicro_sqrttransformed.csv",
#          row.names = FALSE)
#write.csv(data_sqrt$tac, "./data/molecular/transformed/PICRUSt2_predicted_KO_all_tac_noanacyl_sqrttransformed.csv",
#          row.names = FALSE)


# get PCoA for each dataframe 
set.seed(1)
PCoA_list <- lapply(data_sqrt, function(x) getPCoAdata(x, 5))

# making plots
PCoA_plots <- lapply(PCoA_list, function(x) makePCoAplot(x, color = "site", shape = "month"))
lapply(PCoA_plots, print)
# TAC Russian Outlier in July

# tac outlier remover
tac_no_rus_outlier <- data_sqrt$tac %>% filter(!c(site_reach == "RUS-2" & field_date == "2022-07-20"))
tac_no_rus_out_list <- getPCoAdata(tac_no_rus_outlier, 5)
tac_no_rus_out_plot <- makePCoAplot(tac_no_rus_out_list, color = "site", shape = "month")
tac_no_rus_out_plot

#### (5) PERMANOVA ####

# empty table for permanova outputs
p_table <- data.frame(test = NA,
                      sample_type = NA,
                      p_value = NA,
                      F_stat = NA)

# run PERMANOVAs
set.seed(1)
permanovas <- lapply(data_sqrt, function(x) runPERMANOVA(data = x, start_col = 5, group = x$`site`))

# print and add test results to table
for(i in 1:length(permanovas)) {
  
  # print test results to console
  print(names(permanovas[i]))
  print(permanovas[[i]])
  
  # save stats to table
  p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                       sample_type = names(permanovas[i]),
                                       p_value = permanovas[[i]]$`Pr(>F)`[1],
                                       F_stat = permanovas[[i]]$`F`[1]))
}
# significant difference for NT but not TAC nor TM

# check results with outlier removed
runPERMANOVA(data = tac_no_rus_outlier, start_col = 5, group = tac_no_rus_outlier$`site`)
# still not significant!

# check dispersion to see if that influences results
for(i in 1:length(data_sqrt)) {
  set.seed(1)
  anova = anova(betadisper(vegdist(data_sqrt[[i]][,5:ncol(data_sqrt[[i]])], method = "bray"), 
                           data_sqrt[[i]]$site))
  
  
  # print results
  print(names(data_sqrt)[i])
  print(anova)
  
  # add results table
  p_table <- rbind(p_table, data.frame(test = "PERMDISP",
                                       sample_type = names(data_sqrt)[i],
                                       p_value = anova$`Pr(>F)`[1],
                                       F_stat = anova$`F value`[1]))
}
# only TAC had significant dispersion but its PERMANOVA was not significant, so no issues!

# compare PERMANOVA with outlier for TAC Russian Removed

# save results
write.csv(p_table[-1,], "./data/PERMANOVA_results/Q1_predfunction.csv", row.names = FALSE)
