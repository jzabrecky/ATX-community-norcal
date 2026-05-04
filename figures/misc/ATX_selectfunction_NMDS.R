
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

tm_NMDS_list <- getNMDSdata(tm_grouping, start_col = 8, ASV = FALSE)
makeNMDSplot(tm_NMDS_list, TRUE, TRUE, 
             color = "atx_group", shape = "site")
runPERMANOVA(tm_grouping, start_col = 8, group = tm_grouping$atx_group)

# south fork eel only
tm_NMDS_list <- getNMDSdata(tm_ko %>% filter(site == "SFE-M"), start_col = 8, ASV = TRUE)
makeNMDSplot(tm_NMDS_list, FALSE, FALSE, 
             color = "atx_group", shape = "site")

tm_NMDS_list <- getNMDSdata(tm_grouping %>% filter(site == "SFE-M"), start_col = 8, ASV = FALSE)
makeNMDSplot(tm_NMDS_list, TRUE, TRUE, 
             color = "atx_group", shape = "site")
runPERMANOVA(tm_grouping %>% filter(site == "SFE-M"), start_col = 8, group = (tm_grouping %>% filter(site == "SFE-M"))$atx_group)

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

tac_NMDS_list <- getNMDSdata(tac_grouping %>% filter(site == "RUS"), start_col = 8, ASV = FALSE)
makeNMDSplot(tac_NMDS_list, TRUE, FALSE, 
             color = "atx_group", shape = "atx_group")
runPERMANOVA(tac_grouping %>% filter(site == "RUS"), start_col = 8, group = (tac_grouping %>% filter(site == "RUS"))$atx_group)


# south fork eel river only
tac_NMDS_list <- getNMDSdata(tac_ko %>% filter(site == "SFE-M"), start_col = 8, ASV = TRUE)
makeNMDSplot(tac_NMDS_list, FALSE, FALSE, 
             color = "atx_group", shape = "atx_group")

tac_NMDS_list <- getNMDSdata(tac_grouping %>% filter(site == "SFE-M"), start_col = 8, ASV = FALSE)
makeNMDSplot(tac_NMDS_list, TRUE, TRUE, 
             color = "atx_group", shape = "atx_group")
runPERMANOVA(tac_grouping %>% filter(site == "SFE-M"), start_col = 8, group = (tac_grouping %>% filter(site == "SFE-M"))$atx_group)

#### NT addition 

### matching anatoxin values to NT assemblages
atx_nt <- full_join(atx %>% filter(sample_type == "TM") %>% 
                      na.omit() %>% 
                      dplyr::rename(TM_atx_group = atx_group) %>% 
                      select(!c(sample_type, ATX_all_ug_org_mat, atx_detected)),
                    atx %>% filter(sample_type == "TAC") %>% 
                      na.omit() %>% 
                      dplyr::rename(TAC_atx_group = atx_group) %>% 
                      select(!c(sample_type, ATX_all_ug_org_mat, atx_detected)), by = c("field_date", "site", "site_reach")) %>% 
  mutate(nt_atx_group = case_when(is.na(TAC_atx_group) ~ TM_atx_group,
                                  is.na(TM_atx_group) ~ TAC_atx_group,
                                  TAC_atx_group == "high" | TM_atx_group == "high" ~ "high",
                                  TM_atx_group == "medium" | TAC_atx_group == "medium" ~ "medium",
                                  TAC_atx_group == "low" | TM_atx_group == "low" ~ "low",
                                  TRUE ~ "TBD"),
         atx_detected = case_when(nt_atx_group == "none" ~ "n",
                                  TRUE ~ "y")) %>% 
  dplyr::rename(atx_group = nt_atx_group)

# join w/ NT
nt_ko_final <- left_join(nt_ko, atx_nt, by = c("field_date", "site", "site_reach")) %>% 
  relocate(atx_detected, .before = "field_date") %>% 
  relocate(atx_group, .before = "field_date") %>% 
  mutate(atx_detected = case_when(is.na(atx_detected) ~ "not_present",
                                  TRUE ~ atx_detected),
         atx_group = case_when(is.na(atx_group) ~ "not_present",
                               TRUE ~ atx_group))
nt_grouping_final <- left_join(atx_nt, nt_grouping, by = c("field_date", "site", "site_reach")) %>% 
  relocate(atx_detected, .before = "field_date") %>% 
  relocate(atx_group, .before = "field_date") %>% 
  mutate(atx_detected = case_when(is.na(atx_detected) ~ "not_present",
                                  TRUE ~ atx_detected),
         atx_group = case_when(is.na(atx_group) ~ "not_present",
                               TRUE ~ atx_group))

# remove rows w/ missing data
nt_ko_final <- nt_ko_final %>% 
  select(!c(TM_atx_group, TAC_atx_group)) %>% 
  na.omit()
nt_grouping_final <- nt_grouping_final %>% 
  select(!c(TM_atx_group, TAC_atx_group)) %>% 
  na.omit()

# make NMDS plot for NT
nt_NMDS_list <- getNMDSdata(nt_ko_final, start_col = 8, ASV = TRUE)
makeNMDSplot(nt_NMDS_list, FALSE, FALSE, 
             color = "atx_group", shape = "site")
runPERMANOVA(nt_ko_final, start_col = 8, group = nt_ko_final$atx_group)
# functional plot


# south fork eel only
nt_NMDS_list <- getNMDSdata(nt_ko_final %>% filter(site == "SFE-M"), start_col = 8, ASV = TRUE)
makeNMDSplot(nt_NMDS_list, FALSE, FALSE, 
             color = "atx_group", shape = "site")
runPERMANOVA(nt_ko_final %>% filter(site == "SFE-M"), start_col = 8, group = (nt_ko_final
                                                                              %>% filter(site == "SFE-M"))$atx_group)
makeNMDSplot(nt_NMDS_list, TRUE, FALSE, 
             color = "atx_detected", shape = "site")
runPERMANOVA(nt_ko_final %>% filter(site == "SFE-M"), start_col = 8, group = (nt_ko_final
                                                                              %>% filter(site == "SFE-M"))$atx_detected)

nt_NMDS_list <- getNMDSdata(nt_grouping_final %>% filter(site == "SFE-M"), start_col = 8, ASV = FALSE)
makeNMDSplot(nt_NMDS_list, TRUE, FALSE, 
             color = "atx_detected", shape = "site")

runPERMANOVA(nt_grouping_final %>% filter(site == "SFE-M"), start_col = 8, group = (nt_grouping_final
                                                                              %>% filter(site == "SFE-M"))$atx_detected)
# p = 0.109
runPERMANOVA(nt_ko_final %>% filter(site == "SFE-M"), start_col = 8, group = (nt_ko_final
                                                                                    %>% filter(site == "SFE-M"))$atx_detected)

# not quite 0.242


# Russian River only -- isues with this
nt_NMDS_list <- getNMDSdata((nt_ko_final %>% filter(site == "RUS"))[-3,], start_col = 8, ASV = TRUE)
makeNMDSplot(nt_NMDS_list, FALSE, FALSE, 
             color = "atx_group", shape = "site")
runPERMANOVA(nt_ko_final[-30,] %>% filter(site == "RUS"), start_col = 8, group = nt_ko_final[-30,]$atx_group
             %>% filter(site == "RUS"))