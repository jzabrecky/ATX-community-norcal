

# read in KO data as grouping pivot wider and then just unique KOs
nt_grouping <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select.csv") %>% 
  filter(sample_type == "NT") %>% 
  mutate(month = month(mdy(field_date))) %>% 
  mutate(field_date = mdy(field_date)) %>% 
  dplyr::select(!ko_id) %>% 
  # group and sum for each KO id (because there may be multiple entries of the same) 
  dplyr::group_by(site, site_reach, field_date, month, my_grouping, sample_type) %>% 
  dplyr::summarize(total = sum(predicted_gene_abundance)) %>% 
  pivot_wider(names_from = "my_grouping", values_from = "total", values_fill = 0)
tm_grouping <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tm_nomicro.csv") %>% 
  mutate(month = month(mdy(field_date))) %>% 
  mutate(field_date = mdy(field_date)) %>% 
  # group and sum for each KO id (because there may be multiple entries of the same) 
  dplyr::group_by(site, site_reach, field_date, month, my_grouping, sample_type) %>% 
  dplyr::summarize(total = sum(predicted_gene_abundance)) %>% 
  pivot_wider(names_from = "my_grouping", values_from = "total", values_fill = 0)
tac_grouping <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tac_noanacyl.csv") %>% 
  mutate(month = month(mdy(field_date))) %>% 
  mutate(field_date = mdy(field_date)) %>% 
  dplyr::select(!ko_id) %>% 
  # group and sum for each KO id (because there may be multiple entries of the same) 
  dplyr::group_by(site, site_reach, field_date, month, my_grouping, sample_type) %>% 
  dplyr::summarize(total = sum(predicted_gene_abundance)) %>% 
  pivot_wider(names_from = "my_grouping", values_from = "total", values_fill = 0)
  
nt_ko <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select.csv") %>% 
  filter(sample_type == "NT") %>% 
  mutate(month = month(mdy(field_date))) %>% 
  mutate(field_date = mdy(field_date)) %>% 
  select(!my_grouping) %>% 
  pivot_wider(names_from = "ko_id", values_from = "predicted_gene_abundance", values_fill = 0)
tm_ko <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tm_nomicro.csv") %>% 
  mutate(month = month(mdy(field_date))) %>% 
  mutate(field_date = mdy(field_date)) %>% 
  select(!my_grouping) %>% 
  pivot_wider(names_from = "ko_id", values_from = "predicted_gene_abundance", values_fill = 0)
tac_ko <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tac_noanacyl.csv") %>% 
  mutate(month = month(mdy(field_date))) %>% 
  mutate(field_date = mdy(field_date)) %>% 
  select(!my_grouping) %>% 
  pivot_wider(names_from = "ko_id", values_from = "predicted_gene_abundance", values_fill = 0)




# load in functions for NMDS plot
source("./code/supplemental_code/S4b_community_analyses_func.R")

# make NMDS plot for TM
tm_NMDS_list <- getNMDSdata(tm_ko, start_col = 6, ASV = TRUE)
makeNMDSplot(tm_NMDS_list, FALSE, FALSE, 
             color = "site", shape = "month")
# functional plot
tm_NMDS_list <- getNMDSdata(tm_grouping, start_col = 6, ASV = TRUE)
makeNMDSplot(tm_NMDS_list, FALSE, FALSE, 
             color = "site", shape = "month")


tac_NMDS_list <- getNMDSdata(tac_ko, start_col = 6, ASV = TRUE)
makeNMDSplot(tac_NMDS_list, FALSE, FALSE, 
             color = "site", shape = "month")
tac_NMDS_list <- getNMDSdata(tac_grouping, start_col = 6, ASV = TRUE)
makeNMDSplot(tac_NMDS_list, FALSE, FALSE, 
             color = "site", shape = "month")

nt_NMDS_list <- getNMDSdata(nt_ko, start_col = 6, ASV = TRUE)
makeNMDSplot(nt_NMDS_list, FALSE, FALSE, 
             color = "site", shape = "month")
nt_NMDS_list <- getNMDSdata(nt_grouping, start_col = 6, ASV = TRUE)
makeNMDSplot(nt_NMDS_list, FALSE, FALSE, 
             color = "site", shape = "month")

#### PERMANOVA test ####

# run PERMANOVAs
set.seed(1)
runPERMANOVA(tm_ko, 6, group = tm_ko$`site`) # not significantly different p = 0.321
runPERMANOVA(tm_grouping, 6, group = tm_grouping$`site`) # p = 0.329
runPERMANOVA(tac_ko, 6, group = tac_ko$`site`) # not significantly different p = 0.053 close??
runPERMANOVA(tac_grouping, 6, group = tac_grouping$`site`) # p = 0.059
runPERMANOVA(nt_ko, 6, group = nt_ko$`site`) # not significantly different p = 0.004
runPERMANOVA(nt_grouping, 6, group = nt_grouping$`site`) # p = 0.006
anova(betadisper(vegdist(nt_ko[,6:ncol(nt_ko)], method = "bray"), 
                 nt_ko$site)) # not significant
anova(betadisper(vegdist(nt_grouping[,6:ncol(nt_grouping)], method = "bray"), 
                 nt_grouping$site)) # not significant

anova(betadisper(vegdist(tac_ko[,6:ncol(tac_ko)], method = "bray"), 
                 tac_ko$site)) # not significant
anova(betadisper(vegdist(tac_grouping[,6:ncol(tac_grouping)], method = "bray"), 
                 tac_grouping$site)) # not significant
