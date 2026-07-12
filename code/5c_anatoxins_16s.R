#### Comparing 16s molecular data with regard to anatoxin concentrations
### Jordan Zabrecky
## last edited: 07.11.2026

# This script examines how communities as identified by microscopy
# change with increasing anatoxin concentrations with PERMANOVA
# and PCoA. Note that as only one Salmon sample contained detectable 
# anatoxin, we omit the Salmon River samples from this analyses

# Note (7/11): for assemblage comparisons (PERMANOVA), we verified that there
# were no significant differences among reaches except for NT and TAC
# Russian River samples, found on script:
# "S5d_differences_among_reaches.R"
# For these, we added a strata argument to permanovas

#### (1) Loading libraries & data ####

# libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot"), require, character.only = T)

# read in files (data transformed in previous script, "4b_amongrivers_16s.R")
data <- lapply(list.files(path = "./data/molecular/transformed/", pattern = "16s_nochimera"),
               function(x) read.csv(paste("./data/molecular/transformed/", x, sep = "")))
names(data) <- c("NT", "TAC", "TM")

# also read in long for making bar plots
nt <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_FINAL.csv") %>% 
  filter(sample_type == "NT")
tm <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_TM_nomicro.csv")
tac <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_TAC_noanacyl.csv")

# add into list
data_long <- list(nt, tm, tac)
names(data_long) <- c("NT", "TM", "TAC")

# also read in diversity data
diversity <- lapply(list.files(path = "./data/molecular/shannon_diversity/", pattern = ".csv"),
                   function(x) read.csv(paste("./data/molecular/shannon_diversity/", x, sep = "")))
names(diversity) <- c("NT", "TAC", "TM")

# read in toxin data
atx <- read.csv("./data/field_and_lab/atx_w_categorical_groupings.csv")

# split
atx <- split(atx, atx$sample_type)

# join in atx w/ assemblage data
data <- lapply(names(atx), function(x) {
  left_join(atx[[x]], data[[x]], by = c("field_date", "site", "site_reach", "sample_type"))
})
names(data) <- c("NT", "TAC", "TM")
diversity <- lapply(names(atx), function(x) {
  left_join(diversity[[x]], atx[[x]], by = c("field_date", "site", "site_reach"))
})
names(diversity) <- c("NT", "TAC", "TM")
data_long <- lapply(names(atx), function(x) {
  left_join(data_long[[x]] %>% mutate(field_date = mdy(field_date)),
            atx[[x]] %>% mutate(field_date = ymd(field_date)),
            by = c("field_date", "site", "site_reach", "sample_type")) %>% 
    dplyr::rename(relative_abundance = picrust2_relative_abundance) %>% 
    mutate(phylum_class = paste(phylum, " - ", class))
})
names(data_long) <- c("NT", "TAC", "TM")

# load additional functions for analyses
source("./code/supplemental_code/S4b_community_analyses_func.R")
source("./code/supplemental_code/S4c_grouping_func.R")

# set start_column
start_col <- 10

# need to remove ATX data where we did not have 16s information
data$TAC <- data$TAC[-which(is.na(data$TAC$triplicate)),] # one missing Russian and one in SFE
data$NT <- data$NT[-which(is.na(data$NT$triplicate)),] # one NT Russian

#### (2) Relative Abundance Bar Plots ####

# add in broader data categories for plotting
# lastly, add in broader categories to data_long so we don't have several with barplots
# could do this quickly with forcats package, but since we are lumping, may make sense
# to have more customized categories
phylums <- lapply(data_long, function(x) x %>%
                    # calculate average abundance for each phylum
                    group_by(phylum) %>% 
                    dplyr::summarize(mean = mean(relative_abundance)) %>% 
                    # put into "Other" category if NA or % is less than #%
                    mutate(phylum_cat = case_when(is.na(phylum) ~ "Other",
                                                  mean < 0.03  ~ "Other",
                                                  TRUE ~ phylum)) %>% 
                    select(phylum, phylum_cat))
classes <- lapply(data_long, function(x) x %>%
                    # calculate average abundance for each (phylum - ) class
                    group_by(phylum_class) %>% 
                    dplyr::summarize(mean = mean(relative_abundance)) %>% 
                    # put into "Other" category if NA or % is less than #%
                    mutate(classes_cat = case_when(is.na(phylum_class) ~ "Other",
                                                  mean < 0.05  ~ "Other",
                                                  TRUE ~ phylum_class)) %>% 
                    select(phylum_class, classes_cat))
cyano_order <- lapply(data_long, function(x) x %>%
                        # filter for cyanobacteria phylum
                        filter(phylum == "Cyanobacteria") %>% 
                        # calculate average abundance for each order
                        group_by(order) %>% 
                        dplyr::summarize(mean = mean(relative_abundance)) %>% 
                        # put into "Other" category if NA or % is less than #%
                        mutate(cyano_order_cat = case_when(is.na(order) ~ "Other",
                                                      mean < 0.035  ~ "Other",
                                                      TRUE ~ order)) %>% 
                        select(order, cyano_order_cat))
cyano_genus <- lapply(data_long, function(x) x %>%
                        # filter for cyanobacteria phylum
                        filter(phylum == "Cyanobacteria") %>% 
                        # calculate average abundance for each genus
                        group_by(genus) %>% 
                        dplyr::summarize(mean = mean(relative_abundance)) %>% 
                        # put into "Other" category if NA or % is less than #%
                        mutate(cyano_genus_cat = case_when(is.na(genus) ~ "Other",
                                                      mean < 0.11  ~ "Other",
                                                      TRUE ~ genus)) %>% 
                        select(genus, cyano_genus_cat))
# add these categories into our long dataframes
data_long_broader <- data_long
for(i in 1:length(data_long)) {
  data_long_broader[[i]] <- left_join(data_long_broader[[i]], phylums[[i]], by = "phylum")
  data_long_broader[[i]] <- left_join(data_long_broader[[i]], classes[[i]], by = "phylum_class")
  data_long_broader[[i]] <- left_join(data_long_broader[[i]], cyano_order[[i]], by = "order")
  data_long_broader[[i]] <- left_join(data_long_broader[[i]], cyano_genus[[i]], by = "genus")
}

# put bar plots into lists
# (note if issues- something is masked; restart R)
barplot_phylum_plots <- lapply(data_long_broader, function(x) barplot(x, x = "atx_group", y = "relative_abundance",
                                                                      fill = "phylum_cat", facet_wrap = "site"))
barplot_class_plots <- lapply(data_long_broader, function(x) barplot(x, x = "atx_group", y = "relative_abundance",
                                                                     fill = "classes_cat") + facet_wrap(~site))
barplot_cyanoorder_plots <- lapply(data_long_broader, function(x) barplot(x %>% filter(phylum == "Cyanobacteria"),
                                                                          x = "atx_group", y = "relative_abundance",
                                                                          fill = "cyano_order_cat") + facet_wrap(~site))
barplot_cyanogenus_plots <- lapply(data_long_broader, function(x) barplot(x %>% filter(phylum == "Cyanobacteria"),
                                                                          x = "atx_group", y = "relative_abundance",
                                                                          fill = "cyano_genus_cat") + facet_wrap(~site))
 
# titles for plots
titles <- c("Non-Target Samples", 
            "Anabaena/Cylindrospermum Samples (excluding AC & GA)",
            "Microcoleus Samples (excluding M)")

# view plots (all bacteria categories)
for(i in 1:length(barplot_phylum_plots)) {
  print(barplot_phylum_plots[[i]] + labs(title = titles[i]))
  print(barplot_class_plots[[i]] + labs(title = titles[i]))
}
# note: NA is not present for NT!
# NT notes: seeming to be more difference for those with non-detections versus those with
# more cyanobacteria in Eel with non-detection which is interesting
# seems to be a mis-labeling
# not seeing as many differences with TAC and TM in SFE

#### (3) Alpha Diversity ####

# load script for linear models
source("./code/supplemental_code/S5c_linear_analyses.R")

# diversity probably does not need to be log-transformed!
hist(diversity$NT$shannon_diversity)

## (a) SFE-M NT
diversity_sfe_nt <- linear_model(data = diversity$NT %>% filter(site == "SFE-M"),
                                 x = "shannon_diversity", y = "log_ATX_all_ug_org_mat",
                                 rsquared = TRUE)
diversity_sfe_nt$plot
summary(diversity_sfe_nt$model) # relationship!

## (b) SFE-M TM
diversity_sfe_tm <- linear_model(data = diversity$TM %>% filter(site == "SFE-M"),
                                 x = "shannon_diversity", y = "log_ATX_all_ug_org_mat",
                                 rsquared = TRUE)
diversity_sfe_tm$plot 
summary(diversity_sfe_tm$model) # no relationship

## (c) SFE-M TAC
diversity_sfe_tac <- linear_model(data = diversity$TAC %>% filter(site == "SFE-M"),
                                 x = "shannon_diversity", y = "log_ATX_all_ug_org_mat",
                                 rsquared = TRUE)
diversity_sfe_tac$plot
summary(diversity_sfe_tac$model) # no relationship

## (d) RUS NT
diversity_rus_nt <- linear_model(data = diversity$NT %>% filter(site == "RUS"),
                                  x = "shannon_diversity", y = "log_ATX_all_ug_org_mat",
                                  rsquared = TRUE)
diversity_rus_nt$plot
summary(diversity_rus_nt$model) # no relationship

## (e) RUS TAC
diversity_rus_tac <- linear_model(data = diversity$TAC %>% filter(site == "RUS"),
                                 x = "shannon_diversity", y = "log_ATX_all_ug_org_mat",
                                 rsquared = TRUE)
diversity_rus_tac$plot # 0.13
summary(diversity_rus_tac$model) # no relationship

# higher anatoxins when diversity is overall higher in the non-target samples!

#### (4) PCoA ####

## (a) South Fork Eel River
set.seed(1)
PCoA_list_eel <- lapply(c("NT", "TM", "TAC"), function(x) getPCoAdata(data[[x]] %>% 
                                                                filter(site == "SFE-M"), start_col))
names(PCoA_list_eel) <- c("NT", "TM", "TAC")

# making plots
PCoA_plots_eel <- lapply(PCoA_list_eel, function(x) makePCoAplot(x, color = "atx_group", shape = "atx_group"))
lapply(PCoA_plots_eel, print)
# seems like difference between detected and non-detected, but not for variability
# TAC outlier annoying me but it is what it is

## (b) Russian River
set.seed(1)
PCoA_list_rus <- lapply(c("NT", "TAC"), function(x) getPCoAdata(data[[x]] %>% 
                                                                  filter(site == "RUS"), start_col))

# making plots
PCoA_plots_rus <- lapply(PCoA_list_rus, function(x) makePCoAplot(x, color = "atx_group", shape = "atx_group"))
lapply(PCoA_plots_rus, print)
# no difference here for TAC

#### (5) PERMANOVA ####

# separate out data by river

# create table to save results
p_table <-  data.frame(test = NA,
                       sample_type = NA,
                       river = NA,
                       atx = NA,
                       p_value = NA,
                       F_stat = NA)

## (a) Detected versus no Detected anatoxins?

set.seed(1)
for(s in c("SFE-M", "RUS")) {
  for(i in c("NT", "TM", "TAC")) {
    if(s == "RUS" & i == "TM") {
    } else {
      if(s != "RUS") {
        permanova = runPERMANOVA(data = data[[i]] %>% filter(site == s), start_col = start_col, 
                                 group = (data[[i]] %>% filter(site == s))$atx_detected)
      } else {
        permanova = runPERMANOVA(data = data[[i]] %>% filter(site == s), start_col = start_col, 
                                 strata =  (data[[i]] %>% filter(site == s))$site_reach,
                                 group = (data[[i]] %>% filter(site == s))$atx_detected)
      }
      
      p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                           sample_type = i,
                                           river = s,
                                           atx = "detected",
                                           p_value = permanova$`Pr(>F)`[1],
                                           F_stat = permanova$`F`[1]))
      
      permdisp = betadisper(vegdist((data[[i]] %>% filter(site == s))[,start_col:ncol(data[[i]] %>% filter(site == s))], 
                                    method = "bray"), (data[[i]] %>% filter(site == s))$atx_detected)
      test = adonis2(dist(permdisp$distances) ~ (data[[i]] %>% filter(site == s))$atx_detected)
      
      p_table <- rbind(p_table, data.frame(test = "PERMDISP",
                                           sample_type = i,
                                           river = s,
                                           atx = "detected",
                                           p_value = test$`Pr(>F)`[1],
                                           F_stat = test$`F`[1]))
    }
  }
}
# Difference with South Fork Eel but not for Russian River

# ATX groups (including non-detects!)
set.seed(1)
for(s in c("SFE-M", "RUS")) {
  for(i in c("NT", "TM", "TAC")) {
    if(s == "RUS" & i == "TM") {
    } else {
      if(s != "RUS") {
        permanova = runPERMANOVA(data = data[[i]] %>% filter(site == s), start_col = start_col, 
                                 group = (data[[i]] %>% filter(site == s))$atx_group)
      } else {
        permanova = runPERMANOVA(data = data[[i]] %>% filter(site == s), start_col = start_col, 
                                 strata =  (data[[i]] %>% filter(site == s))$site_reach,
                                 group = (data[[i]] %>% filter(site == s))$atx_group)
      }
    
      p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                           sample_type = i,
                                           river = s,
                                           atx = "atx_group",
                                           p_value = permanova$`Pr(>F)`[1],
                                           F_stat = permanova$`F`[1]))
      
      permdisp = betadisper(vegdist((data[[i]] %>% filter(site == s))[,start_col:ncol(data[[i]] %>% filter(site == s))], 
                                    method = "bray"), (data[[i]] %>% filter(site == s))$atx_group)
      test = adonis2(dist(permdisp$distances) ~ (data[[i]] %>% filter(site == s))$atx_group)
      
      p_table <- rbind(p_table, data.frame(test = "PERMDISP",
                                           sample_type = i,
                                           river = s,
                                           atx = "atx_group",
                                           p_value = test$`Pr(>F)`[1],
                                           F_stat = test$`F`[1]))
    }
  }
}

# ATX groups (including non-detects!)
set.seed(1)
for(s in c("SFE-M", "RUS")) {
  for(i in c("NT", "TM", "TAC")) {
    if(s == "RUS" & i == "TM") {
    } else {
      
      if(s != "RUS") {
        permanova = runPERMANOVA(data = data[[i]] %>% filter(site == s & atx_detected == "y"), start_col = start_col, 
                                 group = (data[[i]] %>% filter(site == s & atx_detected == "y"))$atx_group)
      } else {
        permanova = runPERMANOVA(data = data[[i]] %>% filter(site == s & atx_detected == "y"), start_col = start_col, 
                                 strata =  (data[[i]] %>% filter(site == s & atx_detected == "y"))$site_reach,
                                 group = (data[[i]] %>% filter(site == s & atx_detected == "y"))$atx_group)
      }
     
      p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                           sample_type = i,
                                           river = s,
                                           atx = "atx_group_no_nondetects",
                                           p_value = permanova$`Pr(>F)`[1],
                                           F_stat = permanova$`F`[1]))
      
      permdisp = betadisper(vegdist((data[[i]] %>% filter(site == s & atx_detected == "y"))[,start_col:ncol(data[[i]] %>% filter(site == s))], 
                                    method = "bray"), (data[[i]] %>% filter(site == s & atx_detected == "y"))$atx_group)
      test = adonis2(dist(permdisp$distances) ~ (data[[i]] %>% filter(site == s & atx_detected == "y"))$atx_group)
      
      p_table <- rbind(p_table, data.frame(test = "PERMDISP",
                                           sample_type = i,
                                           river = s,
                                           atx = "atx_group_no_nondetects",
                                           p_value = test$`Pr(>F)`[1],
                                           F_stat = test$`F`[1]))
    }
  }
}

# save PERMANOVA results
write.csv(p_table[-1,], "./data/PERMANOVA_results/Q2_molecular.csv", row.names = FALSE)
# again, too small of sample size for RUS (high versus low), but will save results anyways
