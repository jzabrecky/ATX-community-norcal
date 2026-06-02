#### Making anatoxin categories for analyses
### Jordan Zabrecky
## last edited: 06.01.2026

# This script makes anatoxin concentrations groupings (e.g., high versus low) 
# to use in Q2 analyses based on each sample type & river

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

# want to add in category for NT atx (which is mean of TM & TAC when possible)
atx <- atx %>% 
  mutate(NT_ATX_all_ug_orgmat_g = case_when(is.na(TM_ATX_all_ug_orgmat_g) ~ TAC_ATX_all_ug_orgmat_g,
                                           is.na(TAC_ATX_all_ug_orgmat_g) ~ TM_ATX_all_ug_orgmat_g,
                                           TRUE ~ (TM_ATX_all_ug_orgmat_g + TAC_ATX_all_ug_orgmat_g) / 2))

# pivot to long format
atx_long <- atx %>% 
  pivot_longer(c("TM_ATX_all_ug_orgmat_g", "TAC_ATX_all_ug_orgmat_g", "NT_ATX_all_ug_orgmat_g"), 
               values_to = "ATX_all_ug_org_mat",
               names_to = "taxa_ATX") %>% 
  mutate(sample_type = case_when(taxa_ATX == "TM_ATX_all_ug_orgmat_g" ~ "TM",
                                 taxa_ATX == "TAC_ATX_all_ug_orgmat_g" ~ "TAC", 
                                 taxa_ATX == "NT_ATX_all_ug_orgmat_g" ~ "NT")) %>% 
  na.omit() %>% 
  select(!taxa_ATX)

# see distribution of all
ggplot(data = atx_long %>% na.omit() %>% filter(!c(sample_type == "NT")),
       aes(y = ATX_all_ug_org_mat, x = site)) +
  geom_boxplot() +
  scale_y_continuous(trans="pseudo_log")

# see distribution of all (with zeros removed)
ggplot(data = atx_long %>% na.omit() %>% filter(ATX_all_ug_org_mat > 0) %>% filter(!c(sample_type == "NT")), 
       aes(y = ATX_all_ug_org_mat, x = site)) +
  geom_boxplot() +
  scale_y_continuous(trans="pseudo_log")

#### (3) Log-Transforming ATX ####

# need to replace zeros with a value 
# lowest non-zero value is 0.03, how about 0.02 as zero replacement?
log(0.03) # -3.506558
log(0.02) # -3.912

# log values (replace with 0.02 when 0)
atx_long <- atx_long %>% 
  mutate(log_ATX_all_ug_org_mat = case_when(ATX_all_ug_org_mat == 0 ~ log(0.02),
                                                TRUE ~ log(ATX_all_ug_org_mat)))

# view plot
hist(atx_long$log_ATX_all_ug_org_mat) # yes, this is better
hist(atx_long$ATX_all_ug_org_mat)

#### (4) Adding ATX categories to data ####

# makes sense to just decide into upper ("high") and lower values ("low") with our sample sizes
# get medians for each sample type and river
medians <- atx_long %>% 
  # remove samples with no detections
  filter(ATX_all_ug_org_mat > 0) %>% 
  # calculate median for each sample type on each river
  group_by(site, sample_type) %>% 
  dplyr::summarize(med = median(ATX_all_ug_org_mat),
                   total_samples = length(ATX_all_ug_org_mat))

# add in categorical grouping
atx_long <- left_join(atx_long, medians, by = c("site", "sample_type")) %>% 
  mutate(atx_detected = case_when(ATX_all_ug_org_mat > 0 ~ "y",
                                  TRUE ~ "n"),
         atx_group = case_when(ATX_all_ug_org_mat <= med & ATX_all_ug_org_mat > 0 ~ "low",
                               ATX_all_ug_org_mat > med ~ "high",
                               TRUE ~ "none")) %>% 
  select(field_date, site, site_reach, sample_type, ATX_all_ug_org_mat, log_ATX_all_ug_org_mat, atx_group, atx_detected)

#### (3) Saving CSV ####

# save csv
write.csv(atx_long, "./data/field_and_lab/atx_w_categorical_groupings.csv", 
          row.names = FALSE)
