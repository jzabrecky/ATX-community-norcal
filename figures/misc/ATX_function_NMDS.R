
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

# load in atx groupings 
atx <- read.csv("./data/field_and_lab/atx_w_categorical_groupings.csv") %>% 
  mutate(field_date = ymd(field_date))

# join into function
tm_ko <- left_join(tm_ko, atx %>% select(!ATX_all_ug_org_mat), by = c("field_date", "site", "site_reach", "sample_type")) %>% 
  relocate(atx_detected, .before = "field_date") %>% 
  relocate(atx_group, .before = "field_date")
# join into function
tac_ko <- left_join(tac_ko, atx %>% select(!ATX_all_ug_org_mat), by = c("field_date", "site", "site_reach", "sample_type")) %>% 
  relocate(atx_detected, .before = "field_date") %>% 
  relocate(atx_group, .before = "field_date")

# load in functions for NMDS plot
source("./code/supplemental_code/S4b_community_analyses_func.R")

# make NMDS plot for TM
tm_NMDS_list <- getNMDSdata(tm_ko, start_col = 8, ASV = TRUE)
makeNMDSplot(tm_NMDS_list, FALSE, FALSE, 
             color = "atx_group", shape = "site")
# functional plot

# South Fork Eel Only
tm_NMDS_list <- getNMDSdata(tm_ko %>% filter(site == "SFE-M"), start_col = 8, ASV = TRUE)
makeNMDSplot(tm_NMDS_list, FALSE, FALSE, 
             color = "atx_group", shape = "atx_group")

# make NMDS plot for TAC
tac_NMDS_list <- getNMDSdata(tac_ko, start_col = 8, ASV = TRUE)
makeNMDSplot(tac_NMDS_list, FALSE, FALSE, 
             color = "atx_group", shape = "site")
# functional plot

# South Fork Eel Only
tac_NMDS_list <- getNMDSdata(tac_ko %>% filter(site == "SFE-M"), start_col = 8, ASV = TRUE)
makeNMDSplot(tac_NMDS_list, FALSE, FALSE, 
             color = "atx_group", shape = "atx_group")

# Russian River only
tac_NMDS_list <- getNMDSdata(tac_ko %>% filter(site == "RUS"), start_col = 8, ASV = TRUE)
makeNMDSplot(tac_NMDS_list, FALSE, FALSE, 
             color = "atx_group", shape = "atx_group")


#### Species indicator test ####
library(indicspecies)

multipatt(tm_ko[,8:ncol(tm_ko)], tm_ko$atx_group, func = "r.g", control = how(nperm = 999))
# ok this gives like 1000 orthologs lol