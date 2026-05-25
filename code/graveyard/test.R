# gene abundances normalized for % organic matter of sample ???

# read in OM
OM <- read.csv("./data/field_and_lab/environmental_covariates_and_toxins.csv") %>% 
  select(site_reach, field_date, TM_percent_organic_matter, TAC_percent_organic_matter) %>% 
  mutate(field_date = ymd(field_date))
  filter(year(field_date) == 2022)

# nt OM?
OM_nt <- read.csv("./data/EDI_data_package/non_target_sample_miscellaneous.csv")

# read in selected genes
# load data (for selected functional groups)
nt <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select.csv") %>% 
  filter(sample_type == "NT")
tm <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tm_nomicro.csv")
tac <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tac_noanacyl.csv")

# put all dataframes into a list
data_select <- list(nt, tm, tac)
names(data_select) <- c("nt", "tm", "tac")
data_select <- lapply(data_select, function(x) x %>% mutate(field_date = mdy(field_date)))

# left_join
data_select$tac <- left_join(data_select$tac, OM %>% select(!TM_percent_organic_matter),
                             by = c("site_reach", "field_date"))
data_select$tm <-  left_join(data_select$tm, OM %>% select(!TAC_percent_organic_matter),
                             by = c("site_reach", "field_date"))
