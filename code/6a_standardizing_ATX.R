#### Making anatoxin categories for analyses
### Jordan Zabrecky
## last edited: 04.20.2026

# This script makes anatoxin concentrations groupings (e.g., high versus low) 
# to use in Q3 analyses

#### (1) Loading libraries & data ####

# libraries
lapply(c("tidyverse"), require, character.only = T)

# read in environmental covariates & toxin data
atx <- read.csv("data/field_and_lab/environmental_covariates_and_toxins.csv") %>% 
  # will analyze congeners all together standardized by OM, 
  # so removing individual congeners and standardization by chl-a
  select(field_date, site_reach, site, TM_ATX_all_ug_orgmat_g, TAC_ATX_all_ug_orgmat_g) %>% 
  filter(year(ymd(field_date)) == 2022)

#### (2) Exploring ATX data ####

# ATX in long format
atx_long <- atx %>% 
  pivot_longer(c("TM_ATX_all_ug_orgmat_g", "TAC_ATX_all_ug_orgmat_g"), values_to = "ATX_all_ug_org_mat",
                    names_to = "taxa_ATX")

# see distribution of all
ggplot(data = atx_long %>% na.omit(), aes(y = ATX_all_ug_org_mat)) +
  geom_boxplot() +
  scale_y_continuous(trans="pseudo_log")

# see distribution of all (with zeros removed)
ggplot(data = atx_long %>% na.omit() %>% filter(ATX_all_ug_org_mat > 0), aes(y = ATX_all_ug_org_mat)) +
  geom_boxplot() +
  scale_y_continuous(trans="pseudo_log")

# get summary (with zeros removed)
summary(atx_long %>% na.omit() %>% filter(ATX_all_ug_org_mat > 0))
# makes sense to do <50% quantile, =< 1.98
# 50-75% quantile =< 5.5 and > 1.98
# 75-100% quantile > 5.5

#### (3) Adding ATX categories to data

# save median and 3rd quantile
med <- median((atx_long %>% na.omit() %>% filter(ATX_all_ug_org_mat > 0))$ATX_all_ug_org_mat)
third_q <- quantile((atx_long %>% na.omit() %>% filter(ATX_all_ug_org_mat > 0))$ATX_all_ug_org_mat)[4]

# add in categorical grouping
atx_long <- atx_long %>% 
  mutate(atx_detected = case_when(ATX_all_ug_org_mat > 0 ~ "y",
                                  TRUE ~ "n"),
         atx_group = case_when(ATX_all_ug_org_mat <= med & ATX_all_ug_org_mat > 0 ~ "low",
                               ATX_all_ug_org_mat <= third_q & ATX_all_ug_org_mat > med ~ "medium",
                               ATX_all_ug_org_mat > third_q ~ "high",
                               TRUE ~ "none")) %>% 
  mutate(sample_type = case_when(taxa_ATX == "TM_ATX_all_ug_orgmat_g" ~ "TM",
                                 taxa_ATX == "TAC_ATX_all_ug_orgmat_g"  ~ "TAC")) %>% 
  select(field_date, site, site_reach, sample_type, ATX_all_ug_org_mat, atx_detected, atx_group)

#### (3) Saving CSV ####

# save csv
write.csv(atx_long, "./data/field_and_lab/atx_w_categorical_groupings.csv", 
          row.names = FALSE)
