#### Comparing molecular predicted functional profiles with varying ATX concentrations
### Jordan Zabrecky
## last edited: 05.11.2026

# This code compares normalized select orthologs/functions predicted via PICRUSt2-SC,
# from NT, TM, and TAC samples with varying anatoxin concentrations for Q3.
# Data is analyzed using Kruskal-Wallis Tests and visualizations

#### (1) Load libraries & data ####

# libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot"), require, character.only = T)

# load data (for selected functional groups)
nt <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select.csv") %>% 
  filter(sample_type == "NT")
tm <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tm_nomicro.csv")
tac <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tac_noanacyl.csv")

# put all dataframes into a list
data_select <- list(nt, tm, tac)
names(data_select) <- c("NT", "TM", "TAC")

# need to log predicted genes!
data_select <- lapply(data_select, function(x) x %>% 
                        group_by(field_date, site, site_reach, sample_type, functional_grouping) %>% 
                        # need to also group as multiple KOs are listed for each functional group
                        dplyr::summarize(predicted_gene_abundance = sum(predicted_gene_abundance)) %>% 
                        ungroup() %>% 
                        mutate(log_predicted_gene_abundance = log(predicted_gene_abundance)))

# load data for all orthologs
nt <- read.csv("./data/molecular/PICRUSt2_predicted_KO_all_all.csv") %>% 
  filter(sample_type == "NT")
tm <- read.csv("./data/molecular/PICRUSt2_predicted_KO_all_tm_nomicro.csv")
tac <- read.csv("./data/molecular/PICRUSt2_predicted_KO_all_tac_noanacyl.csv")

# put all dataframes into a list
data_all <- list(nt, tm, tac)
names(data_all) <- c("NT", "TM", "TAC")

# read in toxin data
atx <- read.csv("./data/field_and_lab/atx_w_categorical_groupings.csv")

# split
atx <- split(atx, atx$sample_type)

# join in w/ functional data
data_all <- lapply(names(data_all), function(x) left_join(data_all[[x]] %>% 
                                                            mutate(field_date = mdy(field_date)),
                                                          atx[[x]] %>% 
                                                            mutate(field_date = ymd(field_date)),
                                                          by = c("sample_type", "site_reach", "site",
                                                                 "field_date")))
data_select <- lapply(names(data_select), function(x) left_join(data_select[[x]] %>% 
                                                                  mutate(field_date = mdy(field_date)),
                                                                atx[[x]] %>% 
                                                            mutate(field_date = ymd(field_date)),
                                                          by = c("sample_type", "site_reach", "site",
                                                                 "field_date")))
names(data_all) <- c("NT", "TM", "TAC")
names(data_select) <- c("NT", "TM", "TAC")

# lastly, source other scripts for functions
source("./code/supplemental_code/S4b_community_analyses_func.R")
source("./code/supplemental_code/S4d_linear_analyses.R")

#### (2) Kruskal Wallis Tests ####

# make empty results table
kruskal_test_results = data.frame(function_groups = NA,
                                  sample_type = NA,
                                  site = NA,
                                  test = NA,
                                  kruskal_test = NA)

# run tests
for(i in c("NT", "TAC", "TM")) {
  for(s in c("SFE-M", "RUS")) {
    if(s == "RUS" & i == "TM") {
    } else {
      for(f in unique(data_select$TM$functional_grouping)) {
        # detected versus not detected
        kruskal_test_results =  rbind(kruskal_test_results,
                                      data.frame(function_groups = f,
                                           sample_type = i,
                                           site = s,
                                           test = "detected",
                                           kruskal_test = kruskal.test(atx_detected~predicted_gene_abundance, 
                                                                        data = (data_select[[i]] %>% filter(functional_grouping == f)))$p.value))
        # atx groupings
        kruskal_test_results =  rbind(kruskal_test_results,
                                      data.frame(function_groups = f,
                                           sample_type = i,
                                           site = s,
                                           test = "atx_groups",
                                           kruskal_test = kruskal.test(atx_group~predicted_gene_abundance, 
                                                                       data = (data_select[[i]] %>% filter(functional_grouping == f)))$p.value))
      }
    }
  }
}

view(kruskal_test_results)
# nitrification significantly different for detected versus non-detected in SFE-M TM but nothing else

# save results
write.csv(kruskal_test_results,
          "./data/kruskal_wallis_results/Q3_selectfunctions.csv", row.names = FALSE)

#### (3) Linear Models ####

## (a) SFE-M NT
predfunctions_sfe_nt <- lapply(unique(data_select$TM$functional_grouping), function(f) {
  return(linear_model(data = data_select$NT %>% filter(site == "SFE-M") %>% 
                        filter(functional_grouping == f),  x = "log_ATX_all_ug_org_mat", y = "log_predicted_gene_abundance",
                      rsquared = TRUE))
})
names(predfunctions_sfe_nt) <- unique(data_select$TM$functional_grouping)
lapply(names(predfunctions_sfe_nt), function(x) print(predfunctions_sfe_nt[[x]]$plot + ggtitle(x)))
# 1 and 3 are lightly significant (cobalamin & N-fix)

## (b) SFE-M TM
predfunctions_sfe_tm <- lapply(unique(data_select$TM$functional_grouping), function(f) {
  return(linear_model(data = data_select$TM %>% filter(site == "SFE-M") %>% 
                        filter(functional_grouping == f),  x = "log_ATX_all_ug_org_mat", y = "log_predicted_gene_abundance",
                      rsquared = TRUE))
})
names(predfunctions_sfe_tm) <- unique(data_select$TM$functional_grouping)
lapply(names(predfunctions_sfe_tm), function(x) print(predfunctions_sfe_tm[[x]]$plot + ggtitle(x)))
# 2 is lightly significant * (nitrification)

## (c) SFE-M TAC
predfunctions_sfe_tac <- lapply(unique(data_select$TM$functional_grouping), function(f) {
  return(linear_model(data = data_select$TAC %>% filter(site == "SFE-M") %>% 
                        filter(functional_grouping == f),  x = "log_ATX_all_ug_org_mat", y = "log_predicted_gene_abundance",
                      rsquared = TRUE))
})
names(predfunctions_sfe_tac) <- unique(data_select$TM$functional_grouping)
lapply(names(predfunctions_sfe_tac), function(x) print(predfunctions_sfe_tac[[x]]$plot + ggtitle(x)))
# nothing!

## (d) RUS NT
predfunctions_rus_nt <- lapply(unique(data_select$TM$functional_grouping), function(f) {
  return(linear_model(data = data_select$NT %>% filter(site == "RUS") %>% 
                        filter(functional_grouping == f),  x = "log_ATX_all_ug_org_mat", y = "log_predicted_gene_abundance",
                      rsquared = TRUE))
})
names(predfunctions_rus_nt) <- unique(data_select$TM$functional_grouping)
lapply(names(predfunctions_rus_nt), function(x) print(predfunctions_rus_nt[[x]]$plot + ggtitle(x)))
# nothing!

## (e) RUS TAC
predfunctions_rus_tac <- lapply(unique(data_select$TM$functional_grouping), function(f) {
  return(linear_model(data = data_select$TAC %>% filter(site == "RUS") %>% 
                        filter(functional_grouping == f),  x = "log_ATX_all_ug_org_mat", y = "log_predicted_gene_abundance",
                      rsquared = TRUE))
})
names(predfunctions_rus_tac) <- unique(data_select$TM$functional_grouping)
lapply(names(predfunctions_rus_tac), function(x) print(predfunctions_rus_tac[[x]]$plot + ggtitle(x)))
# nothing!

#### (4) NMDS plots ####

# will try two versions of data: (1) hellinger transform and (2) total predicted abundances
# for both all orthologs and our selected functional groups

## (a) entire functional profile (hellinger transformed)

# pivot wider & Hellinger-transform
data_all_wide_hellinger <- lapply(data_all, function(x){
  
  # pivot wider and then use decostand to relative
  y <- x %>% 
    # need to group same ko id's in same sample
    group_by(site, site_reach, ko_id, field_date, sample_type, atx_detected, atx_group) %>% 
    dplyr::summarize(total_abundance = sum(predicted_gene_abundance)) %>% 
    pivot_wider(names_from = "ko_id", values_from = "total_abundance", values_fill = 0)
  y[,7:ncol(y)] <- decostand(y[,7:ncol(y)], method = "hellinger")
  
  return(y)
})

# get NMDS for each dataframe (sqrt-transformed!)
NMDS_list_all_hel <- list()
set.seed(1)
for(i in c("NT", "TAC", "TM")) {
  for(s in c("SFE-M", "RUS")) {
    if(s == "RUS" & i == "TM") {
    } else {
      name = paste(s, i, sep = " ")
      NMDS_list_all_hel[[name]] <- getNMDSdata(data_all_wide_hellinger[[i]] %>% 
                                                   filter(site == s), 7, ASV = TRUE)
    }
  }
}
NMDS_plots_all_hel <- lapply(NMDS_list_all_hel, function(x)
                             makeNMDSplot(x, FALSE, FALSE, color = "atx_detected", shape = "atx_detected"))
lapply(NMDS_plots_all_hel, print)
# NT RUS is problematic again! just going to ignore it for now!

## (b) entire functional profile (total abundances)

# pivot wider
data_all_wide <- lapply(data_all, function(x){
  
  # pivot wider and then use decostand to relative
  y <- x %>% 
    # need to group same ko id's in same sample
    group_by(site, site_reach, ko_id, field_date, sample_type, atx_detected, atx_group) %>% 
    dplyr::summarize(total_abundance = sum(predicted_gene_abundance)) %>% 
    pivot_wider(names_from = "ko_id", values_from = "total_abundance", values_fill = 0)
  
  return(y)
})

# get NMDS for each dataframe
NMDS_list_all_wide <- list()
set.seed(1)
for(i in c("NT", "TAC", "TM")) {
  for(s in c("SFE-M", "RUS")) {
    if(s == "RUS" & i == "TM") {
    } else {
      name = paste(s, i, sep = " ")
      NMDS_list_all_wide[[name]] <- getNMDSdata(data_all_wide_hellinger[[i]] %>% 
                                                 filter(site == s), 7, ASV = TRUE)
    }
  }
}
NMDS_plots_all_wide <- lapply(NMDS_list_all_wide, function(x)
  makeNMDSplot(x, FALSE, FALSE, color = "atx_group", shape = "atx_group"))
lapply(NMDS_plots_all_wide, print)

## (c) select functions (hellinger transformed)

# pivot wider & Hellinger-transform
data_select_wide_hellinger <- lapply(data_select, function(x){
  
  # pivot wider and then use decostand to relative
  y <- x %>% 
    select(field_date, site, site_reach, functional_grouping, predicted_gene_abundance,
           atx_detected, atx_group) %>% 
    # need to group same ko id's in same sample
    pivot_wider(names_from = "functional_grouping", values_from = "predicted_gene_abundance", values_fill = 0)
  y[,6:ncol(y)] <- decostand(y[,6:ncol(y)], method = "hellinger")
  
  return(y)
})

# get NMDS for each dataframe (hellinger transformed!)
set.seed(1)
NMDS_list_select_hel <- list()
set.seed(1)
for(i in c("NT", "TAC", "TM")) {
  for(s in c("SFE-M", "RUS")) {
    if(s == "RUS" & i == "TM") {
    } else {
      name = paste(s, i, sep = " ")
      NMDS_list_select_hel[[name]] <- getNMDSdata(data_all_wide_hellinger[[i]] %>% 
                                                  filter(site == s), 7, ASV = TRUE)
    }
  }
}
NMDS_plots_select_hel <- lapply(NMDS_list_select_hel, function(x)
  makeNMDSplot(x, FALSE, FALSE, color = "atx_group", shape = "atx_group"))
lapply(NMDS_plots_select_hel, print)

## (d) select functions (total abundances!)

# pivot wider
data_select_wide <- lapply(data_select, function(x){
  
  # pivot wider and then use decostand to relative
  y <- x %>% 
    select(field_date, site, site_reach, functional_grouping, predicted_gene_abundance,
           atx_detected, atx_group) %>% 
    pivot_wider(names_from = "functional_grouping", values_from = "predicted_gene_abundance", values_fill = 0)
  
  return(y)
})

# get NMDS for each dataframe (raw abundances!)
NMDS_list_select_wide <- list()
set.seed(1)
for(i in c("NT", "TAC", "TM")) {
  for(s in c("SFE-M", "RUS")) {
    if(s == "RUS" & i == "TM") {
    } else {
      name = paste(s, i, sep = " ")
      NMDS_list_select_wide[[name]] <- getNMDSdata(data_all_wide_hellinger[[i]] %>% 
                                                    filter(site == s), 7, ASV = TRUE)
    }
  }
}
NMDS_plots_select_wide <- lapply(NMDS_list_select_wide, function(x)
  makeNMDSplot(x, FALSE, FALSE, color = "atx_group", shape = "atx_group"))
lapply(NMDS_plots_select_wide, print)

#### (5) PERMANOVA ####

# run PERMANOVAs (there is probably a more efficient set-up!
permanovas_all_hell_det <- list()
permanovas_all_total_det <- list()
permanovas_select_hell_det <- list()
permanovas_select_total_det <- list()
permanovas_all_hell_group <- list()
permanovas_all_total_group <- list()
permanovas_select_hell_group <- list()
permanovas_select_total_group <- list()
set.seed(1)
for(i in c("NT", "TAC", "TM")) {
  for(s in c("SFE-M", "RUS")) {
    if(s == "RUS" & i == "TM") {
    } else {
      name = paste(s, i, sep = " ")
      permanovas_all_hell_det[[name]] <- runPERMANOVA(data = data_all_wide_hellinger[[i]] %>% filter(site == s), 
                                              start_col = 7, group = (data_all_wide_hellinger[[i]] %>% 
                                                filter(site == s))$`atx_detected`, na.action = "na.omit")
      permanovas_all_total_det[[name]] <- runPERMANOVA(data = data_all_wide[[i]] %>% filter(site == s), 
                                                      start_col = 7, group = (data_all_wide[[i]] %>% 
                                                        filter(site == s))$`atx_detected`, na.action = "na.omit")
      permanovas_select_hell_det[[name]] <- runPERMANOVA(data = data_select_wide_hellinger[[i]] %>% filter(site == s), 
                                                      start_col = 7, group = (data_select_wide_hellinger[[i]] %>% 
                                                        filter(site == s))$`atx_detected`, na.action = "na.omit")
      permanovas_select_total_det[[name]] <- runPERMANOVA(data = data_select_wide[[i]] %>% filter(site == s), 
                                                      start_col = 7, group = (data_select_wide[[i]] %>% 
                                                        filter(site == s))$`atx_detected`, na.action = "na.omit")
      permanovas_all_hell_group[[name]] <- runPERMANOVA(data = data_all_wide_hellinger[[i]] %>% filter(site == s), 
                                                      start_col = 7, group = (data_all_wide_hellinger[[i]] %>% 
                                                        filter(site == s))$`atx_group`, na.action = "na.omit")
      permanovas_all_total_group[[name]] <- runPERMANOVA(data = data_all_wide[[i]] %>% filter(site == s), 
                                                       start_col = 7, group = (data_all_wide[[i]] %>% 
                                                         filter(site == s))$`atx_group`, na.action = "na.omit")
      permanovas_select_hell_group[[name]] <- runPERMANOVA(data = data_select_wide_hellinger[[i]] %>% filter(site == s), 
                                                         start_col = 7, group = (data_select_wide_hellinger[[i]] %>% 
                                                           filter(site == s))$`atx_group`, na.action = "na.omit")
      permanovas_select_total_group[[name]] <- runPERMANOVA(data = data_select_wide[[i]] %>% filter(site == s), 
                                                          start_col = 7, group = (data_select_wide[[i]] %>% 
                                                            filter(site == s))$`atx_group`, na.action = "na.omit")
      
    }
  }
}

# empty table for permanova outputs
p_table <- data.frame(test = NA,
                      group = NA,
                      sample_type = NA,
                      all_or_select = NA,
                      data_transformation = NA,
                      p_value = NA,
                      F_stat = NA)

# print and add test results to table
for(i in 1:length(permanovas_all_hell_group)) {
  
  # save stats to table
  ## detected versus non-detected first
  # permanovs for all hellinger transformed
  p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                       group = "detected",
                                       sample_type = names(permanovas_all_hell_det[i]),
                                       all_or_select = "all",
                                       data_transformation = "hellinger",
                                       p_value = permanovas_all_hell_det[[i]]$`Pr(>F)`[1],
                                       F_stat = permanovas_all_hell_det[[i]]$`F`[1]))
  
  # permanovas for square-rooted total abundances
  p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                       group = "detected",
                                       sample_type = names(permanovas_all_total_det[i]),
                                       all_or_select = "all",
                                       data_transformation = "total_abundances",
                                       p_value = permanovas_all_total_det[[i]]$`Pr(>F)`[1],
                                       F_stat = permanovas_all_total_det[[i]]$`F`[1]))
  
  # permanovas for select hellinger transformed
  p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                       group = "detected",
                                       sample_type = names(permanovas_select_hell_det[i]),
                                       all_or_select = "select",
                                       data_transformation = "hellinger",
                                       p_value = permanovas_select_hell_det[[i]]$`Pr(>F)`[1],
                                       F_stat = permanovas_select_hell_det[[i]]$`F`[1]))
  
  # permanovas for square-rooted total abundances
  p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                       group = "detected",
                                       sample_type = names(permanovas_select_total_det[i]),
                                       all_or_select = "select",
                                       data_transformation = "total_abundaces",
                                       p_value = permanovas_select_total_det[[i]]$`Pr(>F)`[1],
                                       F_stat = permanovas_select_total_det[[i]]$`F`[1]))
  
  ## atx groupings 
  p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                       group = "group",
                                       sample_type = names(permanovas_all_hell_group[i]),
                                       all_or_select = "all",
                                       data_transformation = "hellinger",
                                       p_value = permanovas_all_hell_group[[i]]$`Pr(>F)`[1],
                                       F_stat = permanovas_all_hell_group[[i]]$`F`[1]))
  
  # permanovas for square-rooted total abundances
  p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                       group = "group",
                                       sample_type = names(permanovas_all_total_group[i]),
                                       all_or_select = "all",
                                       data_transformation = "total_abundances",
                                       p_value = permanovas_all_total_group[[i]]$`Pr(>F)`[1],
                                       F_stat = permanovas_all_total_group[[i]]$`F`[1]))
  
  # permanovas for select hellinger transformed
  p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                       group = "group",
                                       sample_type = names(permanovas_select_hell_group[i]),
                                       all_or_select = "select",
                                       data_transformation = "hellinger",
                                       p_value = permanovas_select_hell_group[[i]]$`Pr(>F)`[1],
                                       F_stat = permanovas_select_hell_group[[i]]$`F`[1]))
  
  # permanovas for square-rooted total abundances
  p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                       group = "group",
                                       sample_type = names(permanovas_select_total_group[i]),
                                       all_or_select = "select",
                                       data_transformation = "total_abundaces",
                                       p_value = permanovas_select_total_group[[i]]$`Pr(>F)`[1],
                                       F_stat = permanovas_select_total_group[[i]]$`F`[1]))
}
view(p_table)
# Okay in order: 
# differences only for SFE NT w/ Hellinger transformation! between detected and grouping
# no predicted functional differences for mats!
