# prelim ATM



# load data 
data <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select.csv") %>% 
  mutate(field_date = mdy(field_date))

atx <- read.csv("./data/field_and_lab/atx_w_categorical_groupings.csv") %>% 
  mutate(field_date = ymd(field_date))
tm_atx <- atx %>% 
  select(field_date, site_reach, site, log_TM_ATX_all_ug_orgmat_g, TM_atx_detected)
tac_atx <- atx %>% 
  select(field_date, site_reach, site, log_TAC_ATX_all_ug_orgmat_g, TAC_atx_detected)

# left join in data
all_tm <- left_join(data %>% filter(sample_type == "TM"), tm_atx, by = c("field_date", "site", "site_reach")) %>% 
  na.omit() # get rid of 8/23 that did not have enough data for analyses
all_tac <- left_join(data %>% filter(sample_type == "TAC"), tac_atx, by = c("field_date", "site", "site_reach"))
# why am i missing ATX info for 8/23 sfe-

ggplot(all_tm %>% filter(site == "SFE-M"), aes(x = predicted_gene_abundance)) +
  geom_histogram() +
  facet_wrap(~my_grouping, scales = "free")

# compare
ggplot(data = all_tm %>% filter(site == "SFE-M"), aes(x = log(predicted_gene_abundance), y = log_TM_ATX_all_ug_orgmat_g)) +
  facet_wrap(~my_grouping, scales = "free") +
  geom_smooth(method = "lm") +
  geom_point()


ggplot(data = all_tac %>% filter(site == "RUS"), aes(x = log(predicted_gene_abundance), y = log_TAC_ATX_all_ug_orgmat_g)) +
  facet_wrap(~my_grouping, scales = "free") +
  geom_smooth(method = "lm") +
  geom_point()


ggplot(data = all_tac %>% filter(site == "SFE-M"), aes(x = log(predicted_gene_abundance), y = log_TAC_ATX_all_ug_orgmat_g)) +
  facet_wrap(~my_grouping, scales = "free") +
  geom_smooth(method = "lm") +
  geom_point()

test = lm(log_TM_ATX_all_ug_orgmat_g ~ log(predicted_gene_abundance), data = all_tm %>% filter(site == "SFE-M") %>% 
     filter(my_grouping == "nitrification"))
summary(test)
