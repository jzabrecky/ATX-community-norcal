#### Comparing molecular predicted functional profiles with varying ATX concentrations
### Jordan Zabrecky
## last edited: 05.11.2026

# This code compares normalized select orthologs/functions predicted via PICRUSt2-SC,
# from NT, TM, and TAC samples with varying anatoxin concentrations for Q2.
# Data is analyzed using Kruskal-Wallis Tests, linear models, visualizations 
# for selected genes and PERMANOVA & NMDS for the full predicted functional profile

## SEE HOW TARGET GENES DIFFER BETWEEN DETECTS AND NON-DETECTS :)

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

# as decided in the supplemental script, "S4e_testing_data_sqrtformations_predgenes.R",
# we will square-root the predicted gene abundances to minimize impact of high gene counts
# so load in from transformed folder
data_all <- lapply(list.files(path = "./data/molecular/transformed/", pattern = "KO_all"),
                     function(x) read.csv(paste("./data/molecular/transformed/", x, sep = "")))
names(data_all) <- c("NT", "TM", "TAC")

# read in toxin data
atx <- read.csv("./data/field_and_lab/atx_w_categorical_groupings.csv")

# split
atx <- split(atx, atx$sample_type)

# join in w/ functional data
data_all <- lapply(names(data_all), function(x) left_join(data_all[[x]] %>% 
                                                            mutate(field_date = ymd(field_date)),
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
source("./code/supplemental_code/S5c_linear_analyses.R")

# set start col for predicted function matrix
start_col = 5
# note: atx data added to end, so will have to use end_col arguments

#### (2) Kruskal Wallis Tests ####

# make box plots!!!
boxplots <- lapply(data_select, function(x) {
  plot = ggplot(x, aes(x = site, y = predicted_gene_abundance, fill = site)) +
    geom_boxplot() +
    facet_wrap(~functional_grouping, scales = "free")
  print(plot)
  return(plot)
})

# make empty results table
kruskal_test_results = data.frame(function_groups = NA,
                                  sample_type = NA,
                                  site = NA,
                                  test = NA,
                                  kruskal_test = NA)

# empty list for boxplots
boxplots_detected <- list()
boxplots_group <- list()

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
                                                                       data = (data_select[[i]] %>% filter(functional_grouping == f & site == s)))$p.value))
       
        # atx groupings
        kruskal_test_results =  rbind(kruskal_test_results,
                                      data.frame(function_groups = f,
                                           sample_type = i,
                                           site = s,
                                           test = "atx_groups",
                                           kruskal_test = kruskal.test(atx_group~predicted_gene_abundance, 
                                                                       data = (data_select[[i]] %>% filter(functional_grouping == f & site == s)))$p.value))
      }
      
      # make boxplots
      boxplots_detected[[i]][[s]] <- ggplot(data_select[[i]] %>% filter(site == s) %>% na.omit(), 
                                             aes(x = atx_detected, y = predicted_gene_abundance, fill = atx_detected)) +
        geom_boxplot() +
        geom_jitter() +
        facet_wrap(~functional_grouping, scales = "free") +
        ggtitle(paste(s, i))
      
      boxplots_group[[i]][[s]] <- ggplot(data_select[[i]] %>% filter(site == s) %>% na.omit(), 
                                            aes(x = atx_group, y = predicted_gene_abundance, fill = atx_group)) +
        geom_boxplot() +
        geom_jitter() +
        facet_wrap(~functional_grouping, scales = "free") +
        ggtitle(paste(s, i))
      
      # display boxplots
      print(boxplots_detected[[i]][[s]])
      print(boxplots_group[[i]][[s]])
    }
  }
}

view(kruskal_test_results) # nothing significant

# real quick to make sure my code worked
test <- kruskal.test(atx_detected~predicted_gene_abundance, 
             data = data_select$NT %>% filter(functional_grouping == "cobalamin_B12" & site == "SFE-M"))
test # p-value is 0.4544 which matches table

# save results
write.csv(kruskal_test_results,
          "./data/kruskal_wallis_results/Q2_selectfunctions.csv", row.names = FALSE)

# curious which is more interesting; this or LM

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

# as decided in the supplemental script, "S4e_testing_data_sqrtformations_predgenes.R",
# we will square-root the predicted gene abundances to minimize impact of high gene counts

# get NMDS for each dataframe (sqrt-transformed!)
NMDS_list <- list()
set.seed(1)
for(i in c("NT", "TAC", "TM")) {
  for(s in c("SFE-M", "RUS")) {
    if(s == "RUS" & i == "TM") {
    } else {
      name = paste(s, i, sep = " ")
      NMDS_list[[name]] <- getNMDSdata(data_all[[i]] %>% 
                                                   filter(site == s), start_col,
                                       # need to set end_col as we have atx data at end
                                       end_col = (ncol(data_all[[i]] %>% 
                                                        filter(site == s)) - 4), ASV = TRUE)
    }
  }
}
NMDS_plots <- lapply(NMDS_list, function(x)
                             makeNMDSplot(x, FALSE, FALSE, color = "atx_detected", shape = "atx_detected"))
lapply(NMDS_plots, print)
# NT RUS is problematic again! just going to ignore it for now!

#### (5) PERMANOVA ####

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
      # remove rows with no atx information & set end_col
      data = data_all[[i]] %>% filter(site == s) %>% na.omit()
      end_col = (ncol(data_all[[i]] %>% filter(site == s)) - 4)
      
      permanova = runPERMANOVA(data = data, start_col, 
                               end_col = end_col,
                               group = data$atx_detected)
      
      p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                           sample_type = i,
                                           river = s,
                                           atx = "detected",
                                           p_value = permanova$`Pr(>F)`[1],
                                           F_stat = permanova$`F`[1]))
      
      permdisp = betadisper(vegdist(data[,start_col:end_col], 
                                    method = "bray"), data$atx_detected)
      test = adonis2(dist(permdisp$distances) ~ (data$atx_detected))
      
      p_table <- rbind(p_table, data.frame(test = "PERMDISP",
                                           sample_type = i,
                                           river = s,
                                           atx = "detected",
                                           p_value = test$`Pr(>F)`[1],
                                           F_stat = test$`F`[1]))
    }
  }
}
# Difference with South Fork Eel NT but nothing else!

# ATX groups (including non-detects!)
set.seed(1)
for(s in c("SFE-M", "RUS")) {
  for(i in c("NT", "TM", "TAC")) {
    if(s == "RUS" & i == "TM") {
    } else {
      # remove rows with no atx information & set end_col
      data = data_all[[i]] %>% filter(site == s) %>% na.omit()
      end_col = (ncol(data_all[[i]] %>% filter(site == s)) - 4)
      
      permanova = runPERMANOVA(data = data, start_col, 
                               end_col = end_col,
                               group = data$atx_group)
      
      p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                           sample_type = i,
                                           river = s,
                                           atx = "atx_group",
                                           p_value = permanova$`Pr(>F)`[1],
                                           F_stat = permanova$`F`[1]))
      
      permdisp = betadisper(vegdist(data[,start_col:end_col], 
                                    method = "bray"), data$atx_group)
      test = adonis2(dist(permdisp$distances) ~ (data$atx_group))
      
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
      # remove rows with no atx information & set end_col
      data = data_all[[i]] %>% filter(site == s & atx_detected == "y") %>% na.omit()
      end_col = (ncol(data_all[[i]] %>% filter(site == s)) - 4)
      
      permanova = runPERMANOVA(data = data, start_col, 
                               end_col = end_col,
                               group = data$atx_group)
      
      p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                           sample_type = i,
                                           river = s,
                                           atx = "atx_group_nonondetects",
                                           p_value = permanova$`Pr(>F)`[1],
                                           F_stat = permanova$`F`[1]))
      
      permdisp = betadisper(vegdist(data[,start_col:end_col], 
                                    method = "bray"), data$atx_group)
      test = adonis2(dist(permdisp$distances) ~ (data$atx_group))
      
      p_table <- rbind(p_table, data.frame(test = "PERMDISP",
                                           sample_type = i,
                                           river = s,
                                           atx = "atx_group_nonondetects",
                                           p_value = test$`Pr(>F)`[1],
                                           F_stat = test$`F`[1]))
    }
  }
}

# save PERMANOVA results
write.csv(p_table[-1,], "./data/PERMANOVA_results/Q2_predfunc.csv", row.names = FALSE)
