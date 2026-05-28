#### Testing different groupings for ATX groups
### Jordan Zabrecky
## 05.25.2026

# This code tests different grouping for ATX groups (e.g., "high" versus "low")
# and sees how the results of PERMANOVA and ISA vary with these groups
# as well as how NMDS visualizations differ

#### (1) Loading libraries & data ####

# load libraries
lapply(c("tidyverse", "vegan", "indicspecies"), require, character.only = T)

# only loading microscopy & molecular data as functional profiles
# in the initially analysis are not significantly different between
# detected and non-detected for target samples

# microscopy data
microscopy <- lapply(list.files(path = "./data/morphological/transformed/", pattern = ".csv"),
                    function(x) read.csv(paste("./data/morphological/transformed/", x, sep = "")))
names(microscopy) <- c("NT", "TAC", "TM")

# molecular data 
molecular <- lapply(list.files(path = "./data/molecular/transformed/", pattern = "16s_nochimera"),
                     function(x) read.csv(paste("./data/molecular/transformed/", x, sep = "")))
names(molecular) <- c("NT", "TAC", "TM")

# atx data
atx <- read.csv("data/field_and_lab/environmental_covariates_and_toxins.csv") %>% 
  # will analyze congeners all together standardized by OM, 
  # so removing individual congeners and standardization by chl-a
  select(field_date, site_reach, site, TM_ATX_all_ug_orgmat_g, TAC_ATX_all_ug_orgmat_g) %>% 
  filter(year(ymd(field_date)) == 2022)

# add in NT category
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

#### (2) Making ATX groups ####

## (a) option 1
# original groupings from SFS presentation (low, medium, and high determined from all)
med <- median((atx_long %>% na.omit() %>% filter(ATX_all_ug_org_mat > 0 & sample_type != "NT"))$ATX_all_ug_org_mat)
third_q <- quantile((atx_long %>% na.omit() %>% filter(ATX_all_ug_org_mat > 0 & sample_type != "NT"))$ATX_all_ug_org_mat)[4]

# add in categorical grouping
atx_long <- atx_long %>% 
  mutate(atx_group_opt1 = case_when(ATX_all_ug_org_mat <= med & ATX_all_ug_org_mat > 0 ~ "low",
                               ATX_all_ug_org_mat <= third_q & ATX_all_ug_org_mat > med ~ "medium",
                               ATX_all_ug_org_mat > third_q ~ "high",
                               TRUE ~ "none")) %>% 
  select(field_date, site, site_reach, sample_type, ATX_all_ug_org_mat, atx_group_opt1)

## (b) option 2
# high & low determined by ATX from all
atx_long <- atx_long %>% 
  mutate(atx_group_opt2 = case_when(ATX_all_ug_org_mat <= med & ATX_all_ug_org_mat > 0 ~ "low",
                                    ATX_all_ug_org_mat >= med ~ "high",
                                    TRUE ~ "none")) %>% 
  select(field_date, site, site_reach, sample_type, ATX_all_ug_org_mat, atx_group_opt1, atx_group_opt2)

## (c) option 3
# high & low determined for each sample type in each river
medians <- atx_long %>% 
  # remove samples with no detections
  filter(ATX_all_ug_org_mat > 0) %>% 
  # calculate median for each sample type on each river
  group_by(site, sample_type) %>% 
  dplyr::summarize(med = median(ATX_all_ug_org_mat),
                   total_samples = length(ATX_all_ug_org_mat))

# add in option to dataframe
atx_long <- left_join(atx_long, medians, by = c("site", "sample_type")) %>% 
  mutate(atx_group_opt3 = case_when(ATX_all_ug_org_mat <= med & ATX_all_ug_org_mat > 0 ~ "low",
                                    ATX_all_ug_org_mat >= med ~ "high",
                                    TRUE ~ "none")) %>% 
  select(field_date, site, site_reach, sample_type, ATX_all_ug_org_mat, atx_group_opt1, atx_group_opt2, atx_group_opt3)

## (d) option 4
# high & low & medium for each sample type in each river
thirds <- atx_long %>% 
  # remove samples with no detections
  filter(ATX_all_ug_org_mat > 0) %>% 
  # calculate median for each sample type on each river
  group_by(site, sample_type) %>% 
  dplyr::summarize(quantile_1 = quantile(ATX_all_ug_org_mat, .333),
                   quantile_2 = quantile(ATX_all_ug_org_mat, .667),
                   total_samples = length(ATX_all_ug_org_mat))

# add in option to dataframe
atx_long <- left_join(atx_long, thirds, by = c("site", "sample_type")) %>% 
  mutate(atx_group_opt4 = case_when(ATX_all_ug_org_mat <= quantile_1 & ATX_all_ug_org_mat > 0 ~ "low",
                                    ATX_all_ug_org_mat <= quantile_2 & ATX_all_ug_org_mat > quantile_1 ~ "medium",
                                    ATX_all_ug_org_mat > quantile_2 ~ "high",
                                    TRUE ~ "none")) %>%  
  select(field_date, site, site_reach, sample_type, ATX_all_ug_org_mat, atx_group_opt1, atx_group_opt2, atx_group_opt3,
         atx_group_opt4)

#### (3) Join ATX and Microscopy Data ####

# remove "none" as we are only interested in comparing samples with detected atx
atx_long <- atx_long %>% filter(ATX_all_ug_org_mat != 0)

# split atx by sample type
atx_list <- split(atx_long, atx_long$sample_type)

# join microscopy
microscopy_atx <- lapply(microscopy, function(x) left_join(x, atx_long, by = c("site_reach", "site", "field_date", "sample_type")))
molecular_atx <- lapply(molecular, function(x) left_join(x, atx_long, by = c("site_reach", "site", "field_date", "sample_type")))

#### (4) Testing PERMANOVA ####

# load community analyses
source("./code/supplemental_code/S4b_community_analyses_func.R")

# set start_col
micro_start <- 5
molec_start <- 6

# list of different ATX groups
atx_groups <- c("atx_group_opt1", "atx_group_opt2", "atx_group_opt3", "atx_group_opt4")

## (a) SFE-M TM 

# microscopy PERMANOVAs
sfe_tm_micro <- microscopy_atx$TM %>% filter(site == "SFE-M") %>% na.omit()
set.seed(1)
lapply(atx_groups, function(x) {
  test = runPERMANOVA(sfe_tm_micro, micro_start, end_col = ncol(microscopy$TM), group = sfe_tm_micro[[x]])
  print(paste(x, ": ", test$`Pr(>F)`[1], sep = ""))
})
# all not significant: 0.095, 0.099, 0.47, 0.294

# molecular PERAMNOVAS
sfe_tm_molec <- molecular_atx$TM %>% filter(site == "SFE-M") %>% na.omit()
set.seed(1)
lapply(atx_groups, function(x) {
  test = runPERMANOVA(sfe_tm_molec, molec_start, end_col = ncol(molecular$TM), group = sfe_tm_molec[[x]])
  print(paste(x, ": ", test$`Pr(>F)`[1], sep = ""))
})
# all not significant here as well

## (b) SFE-M TAC

# microscopy PERMANOVAs
sfe_tac_micro <- microscopy_atx$TAC %>% filter(site == "SFE-M") %>% na.omit()
set.seed(1)
lapply(atx_groups, function(x) {
  test = runPERMANOVA(sfe_tac_micro, micro_start, end_col = ncol(microscopy$TAC), group = sfe_tac_micro[[x]])
  print(paste(x, ": ", test$`Pr(>F)`[1], sep = ""))
})
# all not significant though opt 3 & 4 get closer to being significant
# will assess via NMDS

# molecular PERAMNOVAS
sfe_tac_molec <- molecular_atx$TAC %>% filter(site == "SFE-M") %>% na.omit()
set.seed(1)
lapply(atx_groups, function(x) {
  test = runPERMANOVA(sfe_tac_molec, molec_start, end_col = ncol(molecular$TAC), group = sfe_tac_molec[[x]])
  print(paste(x, ": ", test$`Pr(>F)`[1], sep = ""))
})
# opt 3 is significant but not the rest

## (c) SFE-M NT
# microscopy PERMANOVAs
sfe_nt_micro <- microscopy_atx$NT %>% filter(site == "SFE-M") %>% na.omit()
set.seed(1)
lapply(atx_groups, function(x) {
  test = runPERMANOVA(sfe_nt_micro, micro_start, end_col = ncol(microscopy$NT), group = sfe_nt_micro[[x]])
  print(paste(x, ": ", test$`Pr(>F)`[1], sep = ""))
})
#  opt 4 is almost significant

# molecular PERAMNOVAS
sfe_nt_molec <- molecular_atx$NT %>% filter(site == "SFE-M") %>% na.omit()
set.seed(1)
lapply(atx_groups, function(x) {
  test = runPERMANOVA(sfe_nt_molec, molec_start, end_col = ncol(molecular$NT), group = sfe_nt_molec[[x]])
  print(paste(x, ": ", test$`Pr(>F)`[1], sep = ""))
})
# opt 2 & 3 is significant but not the rest......

#### (5) Testing NMDS ####

## (a) SFE-M TM

# microscopy
lapply(atx_groups, function(x) {
  # temporary change atx option to atx group
  temp = sfe_tm_micro
  colnames(temp)[which(colnames(temp) == x)] = "atx_group"
  NMDS_data = getNMDSdata(temp, micro_start, end_col = ncol(microscopy$TM), FALSE)
  makeNMDSplot(NMDS_data, loading = TRUE, significant = TRUE, shape = "atx_group", color = "atx_group") +
    ggtitle(paste(x, "TM microscopy", sep = " "))
})

# molecular
lapply(atx_groups, function(x) {
  # temporary change atx option to atx group
  temp = sfe_tm_molec
  colnames(temp)[which(colnames(temp) == x)] = "atx_group"
  NMDS_data = getNMDSdata(temp, molec_start, end_col = ncol(molecular$TM), TRUE)
  makeNMDSplot(NMDS_data, loading = FALSE, significant = FALSE, shape = "atx_group", color = "atx_group") +
    ggtitle(paste(x, "TM molecular", sep = " "))
})

## (b) SFE-M TAC

# microscopy
lapply(atx_groups, function(x) {
  # temporary change atx option to atx group
  temp = sfe_tac_micro
  colnames(temp)[which(colnames(temp) == x)] = "atx_group"
  NMDS_data = getNMDSdata(temp, micro_start, end_col = ncol(microscopy$TAC), FALSE)
  makeNMDSplot(NMDS_data, loading = TRUE, significant = TRUE, shape = "atx_group", color = "atx_group") +
    ggtitle(paste(x, "TAC microscopy", sep = " "))
})

# molecular
lapply(atx_groups, function(x) {
  # temporary change atx option to atx group
  temp = sfe_tac_molec
  colnames(temp)[which(colnames(temp) == x)] = "atx_group"
  NMDS_data = getNMDSdata(temp, molec_start, end_col = ncol(molecular$TAC), TRUE)
  makeNMDSplot(NMDS_data, loading = FALSE, significant = FALSE, shape = "atx_group", color = "atx_group") +
    ggtitle(paste(x, "TAC molecular", sep = " "))
})
# stress is 0 issue here

## (c) SFE-M NT
# microscopy
lapply(atx_groups, function(x) {
  # temporary change atx option to atx group
  temp = sfe_nt_micro
  colnames(temp)[which(colnames(temp) == x)] = "atx_group"
  NMDS_data = getNMDSdata(temp, micro_start, end_col = ncol(microscopy$NT), FALSE)
  makeNMDSplot(NMDS_data, loading = TRUE, significant = TRUE, shape = "atx_group", color = "atx_group") +
    ggtitle(paste(x, "NT microscopy", sep = " "))
})

# molecular
lapply(atx_groups, function(x) {
  # temporary change atx option to atx group
  temp = sfe_nt_molec
  colnames(temp)[which(colnames(temp) == x)] = "atx_group"
  NMDS_data = getNMDSdata(temp, molec_start, end_col = ncol(molecular$NT), TRUE)
  makeNMDSplot(NMDS_data, loading = FALSE, significant = FALSE, shape = "atx_group", color = "atx_group") +
    ggtitle(paste(x, "NT molecular", sep = " "))
})

## (d) RUS TAC
rus_tac_micro <- microscopy_atx$TAC %>% filter(site == "RUS") %>% na.omit()
rus_tac_molec <- molecular_atx$TAC %>% filter(site == "RUS") %>% na.omit()

# microscopy
lapply(atx_groups, function(x) {
  # temporary change atx option to atx group
  temp = rus_tac_micro
  colnames(temp)[which(colnames(temp) == x)] = "atx_group"
  NMDS_data = getNMDSdata(temp, micro_start, end_col = ncol(microscopy$TAC), FALSE)
  makeNMDSplot(NMDS_data, loading = TRUE, significant = TRUE, shape = "atx_group", color = "atx_group") +
    ggtitle(paste(x, "TAC microscopy (Russian)", sep = " "))
})
# zero stress issue here

# molecular
lapply(atx_groups, function(x) {
  # temporary change atx option to atx group
  temp = rus_tac_molec
  colnames(temp)[which(colnames(temp) == x)] = "atx_group"
  NMDS_data = getNMDSdata(temp, molec_start, end_col = ncol(molecular$TAC), TRUE)
  makeNMDSplot(NMDS_data, loading = FALSE, significant = FALSE, shape = "atx_group", color = "atx_group") +
    ggtitle(paste(x, "TAC molecular (Russian)", sep = " "))
})
# zero stress issue here

#### (6) Testing ISA ####

## SFE-M TM
set.seed(1)
lapply(atx_groups, function(x) {
  test = multipatt(sfe_tm_micro[,micro_start:ncol(microscopy$TM)],
            sfe_tm_micro[[x]], func = "r.g", control = how(nperm = 999))
  print(summary(test))
})
# Opt 1: (low & medium) green algae
# Opt 2: (high) leptolyngbya & geitlerinema)
# Opt 3: none
# Opt 4: none

## SFE-M TAC
set.seed(1)
lapply(atx_groups, function(x) {
  test = multipatt(sfe_tac_micro[,micro_start:ncol(microscopy$TAC)],
                   sfe_tac_micro[[x]], func = "r.g", control = how(nperm = 999))
  print(summary(test))
})
# Opt 1: (high & low) geitlerinema
# Opt 2: (low) nodularia
# Opt 3: (high) microcoleus, (low) nodularia and non_e_diatoms
# Opt 4: (high) microcoleus, (low & med) non_e_diatoms, (low) nodularia

## SFE-M NT
set.seed(1)
lapply(atx_groups, function(x) {
  test = multipatt(sfe_nt_micro[,micro_start:ncol(microscopy$NT)],
                   sfe_nt_micro[[x]], func = "r.g", control = how(nperm = 999))
  print(summary(test))
})
# more things here so only writing notable ones
# Opt 1: (high & medium) epithemia
# Opt 2: (high) epithemia, anabaena
# Opt 3: (high) anabaena, epithemia, rhopalodia
# Opt 4: (high & medium) rhopalodia & epithemia

## RUS TAC
set.seed(1)
lapply(atx_groups, function(x) {
  test = multipatt(rus_tac_micro[,micro_start:ncol(microscopy$TAC)],
                   rus_tac_micro[[x]], func = "r.g", control = how(nperm = 999))
  print(summary(test))
})
# nothing!

## RUS NT
rus_nt_micro <- microscopy_atx$NT %>% filter(site == "RUS") %>% na.omit()
set.seed(1)
lapply(atx_groups, function(x) {
  test = multipatt(rus_nt_micro[,micro_start:ncol(microscopy$NT)],
                   rus_nt_micro[[x]], func = "r.g", control = how(nperm = 999))
  print(summary(test))
})
# Opt 1: (low) non_e_r_diatoms
# Opt 2: (low) non_e_r_diatoms
# Opt 3: (low) non_e_r_diatoms, chroococcus
# Opt 4: (high) other coccoids
