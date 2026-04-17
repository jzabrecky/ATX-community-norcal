# to-do: make NMDS plots for each sample type based on KO with shape as month

# read in KO data
nt_ko <- read.csv("./data/molecular/PICRUSt2_predicted_KO_all_all.csv") %>% 
  filter(sample_type == "NT") %>% 
  mutate(month = month(mdy(field_date))) %>% 
  mutate(field_date = mdy(field_date)) %>% 
  # group and sum for each KO id (because there may be multiple entries of the same) 
  dplyr::group_by(site, site_reach, field_date, month, ko_id, sample_type) %>% 
  dplyr::summarize(total = sum(predicted_gene_abundance)) %>% 
  pivot_wider(names_from = "ko_id", values_from = "total", values_fill = 0)
tm_ko <- read.csv("./data/molecular/PICRUSt2_predicted_KO_all_tm_nomicro.csv") %>% 
  mutate(month = month(mdy(field_date))) %>% 
  mutate(field_date = mdy(field_date)) %>% 
  # group and sum for each KO id (because there may be multiple entries of the same) 
  dplyr::group_by(site, site_reach, field_date, month, ko_id, sample_type) %>% 
  dplyr::summarize(total = sum(predicted_gene_abundance)) %>% 
  pivot_wider(names_from = "ko_id", values_from = "total", values_fill = 0)
tac_ko <- read.csv("./data/molecular/PICRUSt2_predicted_KO_all_tac_noanacyl.csv") %>% 
  mutate(month = month(mdy(field_date))) %>% 
  mutate(field_date = mdy(field_date)) %>% 
  # group and sum for each KO id (because there may be multiple entries of the same) 
  dplyr::group_by(site, site_reach, field_date, month, ko_id, sample_type) %>% 
  dplyr::summarize(total = sum(predicted_gene_abundance)) %>% 
  pivot_wider(names_from = "ko_id", values_from = "total", values_fill = 0)
 
# need to pivot wider for NMDS and add in month


# load in functions for NMDS plot
source("./code/supplemental_code/S4b_community_analyses_func.R")

# make NMDS plot for TM
tm_NMDS_list <- getNMDSdata(tm_ko, start_col = 6, ASV = TRUE)
makeNMDSplot(tm_NMDS_list, FALSE, FALSE, 
             color = "site", shape = "month")
# functional plot


tac_NMDS_list <- getNMDSdata(tac_ko, start_col = 6, ASV = TRUE)
makeNMDSplot(tac_NMDS_list, FALSE, FALSE, 
             color = "site", shape = "month")

nt_NMDS_list <- getNMDSdata(nt_ko, start_col = 6, ASV = TRUE)
makeNMDSplot(nt_NMDS_list, FALSE, FALSE, 
             color = "site", shape = "month")

#### removing NT outlier ####
nt_nooutlier <- nt_ko %>% filter(!(field_date == as.Date("2022-08-17") & site_reach == "RUS-2"))


nt_NMDS_list_noout <- getNMDSdata(nt_nooutlier, start_col = 6, ASV = TRUE)
makeNMDSplot(nt_NMDS_list_noout, FALSE, FALSE, 
             color = "site", shape = "month")
