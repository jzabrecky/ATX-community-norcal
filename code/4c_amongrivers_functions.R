#### Comparing predicted functional profiles data among rivers
### Jordan Zabrecky
## last edited: 05.24.2026

# This code compares both the full predicted function and select orthologs/functions 
# obtained from PICRUSt2 for NT, TM, and TAC samples across rivers to answer Q1.
# Data is analyzed using PERMANOVA and NMDS (for full predicted function) 
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
# nitrification significantly different for non-target but nothing else
view(kruskal_test_results)

# perform Dunn's Test as a follow-up to the Kruskal-Wallis Test
dunn.test(x = (data_select$nt %>% filter(functional_grouping == "nitrification"))$predicted_gene_abundance,
          g = (data_select$nt %>% filter(functional_grouping == "nitrification"))$site,
          method = "bonferroni")
# difference exists between SFE-M and RUS and SAL, but not between RUS and SAL

# save results
lapply(names(kruskal_test_results), function(x) {
  write.csv(kruskal_test_results[[x]], paste("./data/kruskal_wallis_results/Q1_", x, "_selectorthologs.csv"),
            row.names = FALSE)
})

#### (4) NMDS Plots ####

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
write.csv(data_sqrt$nt, "./data/molecular/transformed/PICRUSt2_predicted_KO_all_nt_sqrttransformed.csv",
          row.names = FALSE)
write.csv(data_sqrt$tm, "./data/molecular/transformed/PICRUSt2_predicted_KO_all_tm_nomicro_sqrttransformed.csv",
          row.names = FALSE)
write.csv(data_sqrt$tac, "./data/molecular/transformed/PICRUSt2_predicted_KO_all_tac_noanacyl_sqrttransformed.csv",
          row.names = FALSE)


# get NMDS for each dataframe 
set.seed(1)
NMDS_list <- lapply(data_sqrt, function(x) getNMDSdata(x, 5, ASV = TRUE))

# making plots
NMDS_plots <- lapply(NMDS_list, function(x) makeNMDSplot(x, FALSE, FALSE, 
                                                         color = "site", shape = "month"))

lapply(NMDS_plots, print)
# Russian outlier issue
test = data_sqrt$nt %>% filter(!(site_reach == "RUS-1S" & field_date == "2022-07-20"))
getNMDSdata(test, 5, ASV = TRUE)

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

# save results
write.csv(p_table[-1,], "./data/PERMANOVA_results/Q1_predfunction.csv", row.names = FALSE)
