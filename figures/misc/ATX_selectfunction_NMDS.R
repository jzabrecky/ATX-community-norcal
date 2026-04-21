
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
tm_grouping <- left_join(tm_grouping, atx %>% select(!ATX_all_ug_org_mat), by = c("field_date", "site", "site_reach", "sample_type")) %>% 
  relocate(atx_detected, .before = "field_date") %>% 
  relocate(atx_group, .before = "field_date")
# join into function
tac_grouping <- left_join(tac_grouping, atx %>% select(!ATX_all_ug_org_mat), by = c("field_date", "site", "site_reach", "sample_type")) %>% 
  relocate(atx_detected, .before = "field_date") %>% 
  relocate(atx_group, .before = "field_date")

# load in functions for NMDS plot
source("./code/supplemental_code/S4b_community_analyses_func.R")

# make NMDS plot for TM all, KO
tm_NMDS_list <- getNMDSdata(tm_ko, start_col = 8, ASV = TRUE)
makeNMDSplot(tm_NMDS_list, FALSE, FALSE, 
             color = "atx_group", shape = "site")

tm_NMDS_list <- getNMDSdata(tm_grouping, start_col = 8, ASV = TRUE)
makeNMDSplot(tm_NMDS_list, FALSE, FALSE, 
             color = "atx_group", shape = "site")

# south fork eel only
tm_NMDS_list <- getNMDSdata(tm_ko %>% filter(site == "SFE-M"), start_col = 8, ASV = TRUE)
makeNMDSplot(tm_NMDS_list, FALSE, FALSE, 
             color = "atx_group", shape = "site")

tm_NMDS_list <- getNMDSdata(tm_grouping %>% filter(site == "SFE-M"), start_col = 8, ASV = TRUE)
makeNMDSplot(tm_NMDS_list, FALSE, FALSE, 
             color = "atx_group", shape = "site")


# make NMDS plot for TAC all, KO
tac_NMDS_list <- getNMDSdata(tac_ko, start_col = 8, ASV = TRUE)
makeNMDSplot(tac_NMDS_list, FALSE, FALSE, 
             color = "atx_group", shape = "atx_group")

tac_NMDS_list <- getNMDSdata(tac_grouping, start_col = 8, ASV = TRUE)
makeNMDSplot(tac_NMDS_list, FALSE, FALSE, 
             color = "atx_group", shape = "atx_group")

# russian river only
tac_NMDS_list <- getNMDSdata(tac_ko %>% filter(site == "RUS"), start_col = 8, ASV = TRUE)
makeNMDSplot(tac_NMDS_list, FALSE, FALSE, 
             color = "atx_group", shape = "atx_group")

tac_NMDS_list <- getNMDSdata(tac_grouping %>% filter(site == "RUS"), start_col = 8, ASV = TRUE)
makeNMDSplot(tac_NMDS_list, FALSE, FALSE, 
             color = "atx_group", shape = "atx_group")

# south fork eel river only
tac_NMDS_list <- getNMDSdata(tac_ko %>% filter(site == "SFE-M"), start_col = 8, ASV = TRUE)
makeNMDSplot(tac_NMDS_list, FALSE, FALSE, 
             color = "atx_group", shape = "atx_group")

tac_NMDS_list <- getNMDSdata(tac_grouping %>% filter(site == "SFE-M"), start_col = 8, ASV = TRUE)
makeNMDSplot(tac_NMDS_list, FALSE, FALSE, 
             color = "atx_group", shape = "atx_group")
