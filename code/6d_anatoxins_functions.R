# prelim ATM



# load data 
nt <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select.csv") %>% 
  filter(sample_type == "NT") %>% 
  mutate(field_date = mdy(field_date),
         month = month(field_date)) %>% 
  group_by(my_grouping, sample_type, field_date, site, site_reach) %>% 
  dplyr::summarize(predicted_gene_abundance = sum(predicted_gene_abundance)) %>% 
  ungroup()
tm <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tm_nomicro.csv") %>% 
  mutate(field_date = mdy(field_date),
         month = month(field_date)) %>% 
  group_by(my_grouping, sample_type, field_date, site, site_reach) %>% 
  dplyr::summarize(predicted_gene_abundance = sum(predicted_gene_abundance)) %>% 
  ungroup()
tac <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tac_noanacyl.csv") %>% 
  mutate(field_date = mdy(field_date),
         month = month(field_date)) %>% 
  group_by(my_grouping, sample_type, field_date, site, site_reach) %>% 
  dplyr::summarize(predicted_gene_abundance = sum(predicted_gene_abundance))  %>% 
  ungroup()

atx <- read.csv("./data/field_and_lab/atx_w_categorical_groupings.csv") %>% 
  mutate(field_date = ymd(field_date))

# left join in data
all_tm <- left_join(tm, atx, by = c("field_date", "site", "site_reach", "sample_type")) %>% 
  na.omit() # get rid of 8/23 that did not have enough data for analyses
all_tac <- left_join(tac, atx, by = c("field_date", "site", "site_reach", "sample_type"))
# why am i missing ATX info for 8/23 sfe-

ggplot(all_tm %>% filter(site == "SFE-M"), aes(x = predicted_gene_abundance)) +
  geom_histogram() +
  facet_wrap(~my_grouping, scales = "free")

# look at gene abundances
ggplot(all_tm, aes(x = predicted_gene_abundance)) +
  geom_histogram() +
  facet_wrap(~my_grouping, scales = "free")
ggplot(all_tm, aes(x = log(predicted_gene_abundance))) +
  geom_histogram() +
  facet_wrap(~my_grouping, scales = "free")

# normality test
for(i in 1:length(unique(all_tm$my_grouping))) {
  group = unique(all_tm$my_grouping)[i]
  # compare log and non-log
  
  print(shapiro.test((all_tm %>% filter(my_grouping == group))$predicted_gene_abundance))
}

# compare
ggplot(data = all_tm %>% filter(site == "SFE-M"), aes(x = log(predicted_gene_abundance), y = log_ATX_all_ug_org_mat)) +
  facet_wrap(~my_grouping, scales = "free") +
  geom_smooth(method = "lm") +
  geom_point(aes(shape = site_reach, color = field_date))

ggplot(data = all_tac %>% filter(site == "RUS"), aes(x = log(predicted_gene_abundance), y = log_ATX_all_ug_org_mat)) +
  facet_wrap(~my_grouping, scales = "free") +
  geom_smooth(method = "lm") +
  geom_point(aes(shape = site_reach, color = field_date))


ggplot(data = all_tac %>% filter(site == "SFE-M"), aes(x = log(predicted_gene_abundance), y = log_ATX_all_ug_org_mat)) +
  facet_wrap(~my_grouping, scales = "free") +
  geom_smooth(method = "lm") +
  geom_point(aes(shape = site_reach, color = field_date))

test = lm(log_ATX_all_ug_org_mat ~ log(predicted_gene_abundance), data = all_tm %>% filter(site == "SFE-M") %>% 
     filter(my_grouping == "nitrification"))
summary(test)

test = lm(log_ATX_all_ug_org_mat ~ log(predicted_gene_abundance), data = all_tac %>% filter(site == "SFE-M") %>% 
            filter(my_grouping == "nitrification"))
summary(test)

test = lm(log_ATX_all_ug_org_mat ~ log(predicted_gene_abundance), data = all_tm %>% filter(site == "SFE-M") %>% 
            filter(my_grouping == "thiamine"))
summary(test)

test = lm(log_ATX_all_ug_org_mat ~ log(predicted_gene_abundance), data = all_tac %>% filter(site == "RUS") %>% 
            filter(my_grouping == "phosphatase_transporters"))
summary(test)

# comparing mats w/ toxins detected versus those without
# just looking at plots, probably not worthwhile :/

#### Also want to look at NT

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
                                  TRUE ~ "y"),
  dplyr::rename(atx_group = nt_atx_group))
  
atx_nt$log_ATX_all_ug_org_mat <- pmax(atx_nt$log_ATX_all_ug_org_mat.x, atx_nt$log_ATX_all_ug_org_mat.y, na.rm = TRUE)

# join w/ NT
all_nt <- left_join(nt, atx_nt, by = c("field_date", "site", "site_reach")) %>% 

  select(field_date, site, site_reach, my_grouping, predicted_gene_abundance, log_ATX_all_ug_org_mat) %>% 
  na.omit()

# NT - SFE
ggplot(data = all_nt %>% filter(site == "SFE-M"), aes(x = log(predicted_gene_abundance), y = log_ATX_all_ug_org_mat)) +
  facet_wrap(~my_grouping, scales = "free") +
  geom_smooth(method = "lm") +
  geom_point(aes(shape = site_reach, color = field_date))

# NT - RUS
ggplot(data = all_nt %>% filter(site == "RUS"), aes(x = log(predicted_gene_abundance), y = log_ATX_all_ug_org_mat)) +
  facet_wrap(~my_grouping, scales = "free") +
  geom_smooth(method = "lm") +
  geom_point(aes(shape = site_reach, color = field_date))

test = lm(log_ATX_all_ug_org_mat ~ log(predicted_gene_abundance), data = all_nt %>% filter(site == "SFE-M") %>% 
            filter(my_grouping == "nitrification"))
summary(test)

test = lm(log_ATX_all_ug_org_mat ~ log(predicted_gene_abundance), data = all_nt %>% filter(site == "SFE-M") %>% 
            filter(my_grouping == "nitrogen_fixation"))
summary(test)

test = lm(log_ATX_all_ug_org_mat ~ log(predicted_gene_abundance), data = all_nt %>% filter(site == "SFE-M") %>% 
            filter(my_grouping == "cobalamin_B12"))
summary(test)

test = lm(log_ATX_all_ug_org_mat ~ log(predicted_gene_abundance), data = all_nt %>% filter(site == "SFE-M") %>% 
            filter(my_grouping == "thiamine"))
summary(test)

test = lm(log_ATX_all_ug_org_mat ~ log(predicted_gene_abundance), data = all_nt %>% filter(site == "RUS") %>% 
            filter(my_grouping == "phosphatase_transporters"))
summary(test)
