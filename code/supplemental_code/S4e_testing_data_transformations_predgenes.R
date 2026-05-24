#### Playing with data transformations for predicted functional groups
### Jordan Zabrecky
## last edited: 05.24.2026

# Do different data transformations changes our predicted functional group results?
# This script explores: (1) using untransformed and (2) square-root transformed data
# and also (3) relativizing that data (i.e., relative abundance of genes) and
# (4) hellinger-transform (square-root of relative abudance)

# In doing so, we also compare differences in results of targeted samples 
# with the targeted taxa remove, and, in the case of Anabaena samples (TAC), 
# how removing green algae--typically the substrate--alters results

# We look predominantly at PERMANOVAs and then visually assess potential 
# differences with NMDS plots

#### (1) Loading libraries & data ####

# libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot"), require, character.only = T)

# load in PICRUSt2 predicted functions
files_picrust <- list.files(path = "./data/molecular/", pattern = "predicted_KO_all")
data_long_picrust <- lapply(files_picrust, function(x) read.csv(paste("./data/molecular/", x, sep = "")) %>% 
                          mutate(field_date = mdy(field_date),
                                 month = month(field_date)) %>% # add month tag
                          relocate(month, .before = "field_date"))
names(data_long_picrust) <- files_picrust

# need to separate out tm, tac, and nt from first csv loaded
full <- data_long_picrust[[1]]
data_long_picrust[[1]] <- full %>% filter(sample_type == "NT")
data_long_picrust[[4]] <- full %>% filter(sample_type == "TM")
data_long_picrust[[5]] <- full %>% filter(sample_type == "TAC")
names(data_long_picrust)[c(1, 4, 5)] <- c("PICRUSt2_predicted_KO_all_nt",
                                      "PICRUSt2_predicted_KO_all_tm_w_M",
                                      "PICRUSt2_predicted_KO_all_tac_w_AC")

# need to remove an outlier for Russian TAC that is affecting NMDS plotting
data_long_picrust$PICRUSt2_predicted_KO_all_tac_w_AC <- 
  data_long_picrust$PICRUSt2_predicted_KO_all_tac_w_AC %>%
  filter(!(site_reach == "RUS-2" & field_date == as.Date("2022-07-20")))
data_long_picrust$PICRUSt2_predicted_KO_all_tac_noanacyl.csv <- 
  data_long_picrust$PICRUSt2_predicted_KO_all_tac_noanacyl.csv %>%
  filter(!(site_reach == "RUS-2" & field_date == as.Date("2022-07-20")))

# load in anatoxin data to do tests with Q3 questions
atx <- read.csv("./data/field_and_lab/atx_w_categorical_groupings.csv") %>% 
  mutate(field_date = ymd(field_date))

# lastly, load analysis functions
source("./code/supplemental_code/S4b_community_analyses_func.R")

#### (2) Data Transformations ####

# join atx and predicted function data
data_long_picrust = lapply(data_long_picrust, function(x) left_join(x, atx, by = c("field_date", "site", "site_reach",
                                                                                   "sample_type")))

# pivot wider
data <- lapply(data_long_picrust, function(x) {
  y <- x %>% 
    # put atx information at beginning
    #relocate(atx_detected, .before = "ko_id") %>% 
    #relocate(atx_group, .before = "ko_id") %>% 
    # remove atx columns we don't care about 
    dplyr::select(!c("ATX_all_ug_org_mat", "log_ATX_all_ug_org_mat")) %>% 
    # group separate reads of ko_id together %>% 
    group_by(site, site_reach, field_date, sample_type, atx_detected, atx_group, ko_id) %>% 
    dplyr::summarize(predicted_gene_abundance = sum(predicted_gene_abundance)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = ko_id, values_from = predicted_gene_abundance, values_fill = 0)
  return(y)
})

# set start_column (column where gene values begin)
start_col = 7

# square-root transform values
data_sqrttrans <- lapply(data, function(x) {
  x[,start_col:ncol(x)] <- sqrt(x[,start_col:ncol(x)])
  return(x)
})

# relativize values
data_relativized <- lapply(data, function(x) {
  x[,start_col:ncol(x)] <- decostand(x[,start_col:ncol(x)], method = "total")
  return(x)
})
# yes, this does right thing
rowSums(data_relativized$PICRUSt2_predicted_KO_all_tm_nomicro.csv[,7:ncol(data_relativized$PICRUSt2_predicted_KO_all_tm_nomicro.csv)])

# use rowsums to see if this relativized data
data_hellinger <- lapply(data, function(x) {
  x[,start_col:ncol(x)] <- decostand(x[,start_col:ncol(x)], method = "hellinger")
  return(x)
})

#### (3) Q1 PERMANOVAs ####

# create summary table
Q1_permanovas <- data.frame(data = NA,
                            transformation = NA,
                            significant = NA)

# run PERMANOVAs for Q1
for(i in 1:length(data)) {
  # make a temporary dataframe for all PERMANOVAs at index i
  set.seed(1)
  temp = data.frame(data = c(names(data)[i], names(data_sqrttrans)[i]),
                    transformation = c("none", "squareroot-transformation",
                                       "relative", "hellinger"),
                    significant = c(runPERMANOVA(data[[i]], start_col = start_col, group = data[[i]]$site)$`Pr(>F)`[1], 
                                    runPERMANOVA(data_sqrttrans[[i]], start_col = start_col, group = data_sqrttrans[[i]]$site)$`Pr(>F)`[1],
                                    runPERMANOVA(data_relativized[[i]], start_col = start_col, group = data_relativized[[i]]$site)$`Pr(>F)`[1],
                                    runPERMANOVA(data_hellinger[[i]], start_col = start_col, group = data_hellinger[[i]]$site)$`Pr(>F)`[1]))
  
  # add to existing dataframe
  Q1_permanovas <- rbind(Q1_permanovas, temp)
}
view(Q1_permanovas)
# significantly different for all NT
# gets slightly signifcant for AC without transformation, but the more we relativize 
# the less different it gets
# for TM relativizing makes it more significantly different which is... strange?
# with Microcoleus included, as expected, there is no difference
# maybe this is a dispersion "false flag" significance...

# test for dispersion in relativized TM data
PERMDISP = betadisper(vegdist(data_relativized$PICRUSt2_predicted_KO_all_tm_nomicro.csv[,start_col:ncol(data_relativized$PICRUSt2_predicted_KO_all_tm_nomicro.csv)], 
                              method = "bray"), data_relativized$PICRUSt2_predicted_KO_all_tm_nomicro.csv$site)
test = adonis2(dist(PERMDISP$distances) ~ data_relativized$PICRUSt2_predicted_KO_all_tm_nomicro.csv$site)
print(test) # significant dispersion which can influence PERMANOVA
# we shall see how it plots

# also want to look at TAC
PERMDISP = betadisper(vegdist(data$PICRUSt2_predicted_KO_all_tac_noanacyl.csv[,start_col:ncol(data$PICRUSt2_predicted_KO_all_tac_noanacyl.csv)], 
                              method = "bray"), data$PICRUSt2_predicted_KO_all_tac_noanacyl.csv$site)
test = adonis2(dist(PERMDISP$distances) ~ data$PICRUSt2_predicted_KO_all_tac_noanacyl.csv$site)
print(test) # yes significant

#### (4) Q1 NMDS ####

# get data for NMDS plots
set.seed(1)
NMDS_list <- list()
NMDS_list$`untransformed` <- lapply(data, function(x) getNMDSdata(x, start_col, ASV = TRUE))
NMDS_list$`sqrt_transformed` <- lapply(data_sqrttrans, function(x) getNMDSdata(x, start_col, ASV = TRUE))
NMDS_list$`relativized` <- lapply(data_relativized, function(x) getNMDSdata(x, start_col, ASV = TRUE))
NMDS_list$`hellinger` <- lapply(data_hellinger, function(x) getNMDSdata(x, start_col, ASV = TRUE))

# viewing plots against each other
for(i in 1:length(data_sqrttrans)) {
  plots <- list()
  for(j in 1:length(NMDS_list)) {
    plots[[j]] = makeNMDSplot(NMDS_list[[j]][[i]], FALSE, FALSE, shape = "month", color = "site") +
      labs(title = paste(names(data)[i], names(NMDS_list)[j], sep = " "))
  }
  print(plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol = 2))
}
# TM no Micro with relative does seem like it could be similar here
# so likely, they aren't that different
# NT looks more similar than I would think, how about dispersion there?
PERMDISP = betadisper(vegdist(data_sqrttrans$PICRUSt2_predicted_KO_all_nt[,start_col:ncol(data_sqrttrans$PICRUSt2_predicted_KO_all_nt)], 
                              method = "bray"), data_sqrttrans$PICRUSt2_predicted_KO_all_nt$site)
test = adonis2(dist(PERMDISP$distances) ~ data_sqrttrans$PICRUSt2_predicted_KO_all_nt$site)
print(test) # not significant dispersion so true differences

#### (5) Q3 PERMANOVAs ####

# create summary table
Q3_permanovas <- data.frame(data = NA,
                            atx = NA,
                            site = NA,
                            transformation = NA,
                            significant = NA)

# run PERMANOVAs for Q3
for(i in 1:length(data)) {
  
  ## (a) do all sample types for South Fork Eel
  set.seed(1)
  # make a temporary dataframe for all PERMANOVAs at index i for atx detected
  temp = data.frame(data = c(names(data)[i], names(data_sqrttrans)[i]),
                    atx = c("detected", "detected"),
                    site = "SFE-M",
                    transformation = c("none", "squareroot-transformation", "relativized", "hellinger"),
                    significant = c(runPERMANOVA(data[[i]] %>% filter(site == "SFE-M"), start_col = start_col, 
                                                 group = (data[[i]] %>%filter(site == "SFE-M"))$atx_detected, na.action = "na.omit")$`Pr(>F)`[1], 
                                    runPERMANOVA(data_sqrttrans[[i]] %>% filter(site == "SFE-M"), start_col = start_col, 
                                                 group = (data_sqrttrans[[i]] %>% filter(site == "SFE-M"))$atx_detected, na.action = "na.omit")$`Pr(>F)`[1],
                                    runPERMANOVA(data_relativized[[i]] %>% filter(site == "SFE-M"), start_col = start_col, 
                                                 group = (data_relativized[[i]] %>% filter(site == "SFE-M"))$atx_detected, na.action = "na.omit")$`Pr(>F)`[1],
                                    runPERMANOVA(data_hellinger[[i]] %>% filter(site == "SFE-M"), start_col = start_col, 
                                                 group = (data_hellinger[[i]] %>% filter(site == "SFE-M"))$atx_detected, na.action = "na.omit")$`Pr(>F)`[1]))
  
  # make a temporary dataframe for all PERMANOVAs at index i for atx grouping
  temp2 = data.frame(data = c(names(data)[i], names(data_sqrttrans)[i]),
                   atx = c("group", "group"),
                   site = "SFE-M",
                   transformation = c("none", "squareroot-transformation", "relativized", "hellinger"),
                   significant = c(runPERMANOVA(data[[i]] %>% filter(site == "SFE-M"), start_col = start_col, 
                                                group = (data[[i]] %>% filter(site == "SFE-M"))$atx_group, na.action = "na.omit")$`Pr(>F)`[1], 
                                   runPERMANOVA(data_sqrttrans[[i]] %>% filter(site == "SFE-M"), start_col = start_col, 
                                                group = (data_sqrttrans[[i]] %>% filter(site == "SFE-M"))$atx_group, na.action = "na.omit")$`Pr(>F)`[1],
                                   runPERMANOVA(data_relativized[[i]] %>% filter(site == "SFE-M"), start_col = start_col, 
                                                group = (data_relativized[[i]] %>% filter(site == "SFE-M"))$atx_group, na.action = "na.omit")$`Pr(>F)`[1],
                                   runPERMANOVA(data_hellinger[[i]] %>% filter(site == "SFE-M"), start_col = start_col, 
                                                group = (data_hellinger[[i]] %>% filter(site == "SFE-M"))$atx_group, na.action = "na.omit")$`Pr(>F)`[1]))
  
  # add to existing dataframe
  Q3_permanovas <- rbind(Q3_permanovas, temp, temp2)
  
  ## (b) only do TAC & NT for RUS
  if(!(names(data)[i] == "PICRUSt2_predicted_KO_all_tm_nomicro.csv" | 
       names(data)[i] == "PICRUSt2_predicted_KO_all_tm_w_M")) {
    
    # make a temporary dataframe for all PERMANOVAs at index i for atx detected
    temp3 = data.frame(data = c(names(data)[i], names(data_sqrttrans)[i]),
                      atx = c("detected", "detected"),
                      site = "RUS",
                      transformation = c("none", "squareroot-transformation", "relativized", "hellinger"),
                      significant = c(runPERMANOVA(data[[i]] %>% filter(site == "RUS"), start_col = start_col, 
                                                   group = (data[[i]] %>%filter(site == "RUS"))$atx_detected, na.action = "na.omit")$`Pr(>F)`[1], 
                                      runPERMANOVA(data_sqrttrans[[i]] %>% filter(site == "RUS"), start_col = start_col, 
                                                   group = (data_sqrttrans[[i]] %>% filter(site == "RUS"))$atx_detected, na.action = "na.omit")$`Pr(>F)`[1],
                                      runPERMANOVA(data_relativized[[i]] %>% filter(site == "RUS"), start_col = start_col, 
                                                   group = (data_relativized[[i]] %>% filter(site == "RUS"))$atx_detected, na.action = "na.omit")$`Pr(>F)`[1],
                                      runPERMANOVA(data_hellinger[[i]] %>% filter(site == "RUS"), start_col = start_col, 
                                                   group = (data_hellinger[[i]] %>% filter(site == "RUS"))$atx_detected, na.action = "na.omit")$`Pr(>F)`[1]))
    
    # make a temporary dataframe for all PERMANOVAs at index i for atx grouping
    temp4 = data.frame(data = c(names(data)[i], names(data_sqrttrans)[i]),
                       atx = c("group", "group"),
                       site = "RUS",
                       transformation = c("none", "squareroot-transformation", "relativized", "hellinger"),
                       significant = c(runPERMANOVA(data[[i]] %>% filter(site == "RUS"), start_col = start_col, 
                                                    group = (data[[i]] %>% filter(site == "RUS"))$atx_group, na.action = "na.omit")$`Pr(>F)`[1], 
                                       runPERMANOVA(data_sqrttrans[[i]] %>% filter(site == "RUS"), start_col = start_col, 
                                                    group = (data_sqrttrans[[i]] %>% filter(site == "RUS"))$atx_group, na.action = "na.omit")$`Pr(>F)`[1],
                                       runPERMANOVA(data_relativized[[i]] %>% filter(site == "RUS"), start_col = start_col, 
                                                    group = (data_relativized[[i]] %>% filter(site == "RUS"))$atx_group, na.action = "na.omit")$`Pr(>F)`[1],
                                       runPERMANOVA(data_hellinger[[i]] %>% filter(site == "RUS"), start_col = start_col, 
                                                    group = (data_hellinger[[i]] %>% filter(site == "RUS"))$atx_group, na.action = "na.omit")$`Pr(>F)`[1]))
    
    # add to existing dataframe
    Q3_permanovas <- rbind(Q3_permanovas, temp3, temp4)
  }
  
}
view(Q3_permanovas)
# SFE Results:
# difference for NT w/ detection and mostly difference (all except raw) for NT w/ group
# no difference for TAC or TM

# RUS Results:
# no difference in either

#### (6) Q3 NMDS ####

## (a) SFE

# make separate NMDS lists for ATX
NMDS_sfe_list <- list()
NMDS_sfe_list$`untransformed` <- lapply(data, function(x) getNMDSdata(x %>% filter(site == "SFE-M"), start_col, ASV = TRUE))
NMDS_sfe_list$`sqrt_transformed` <- lapply(data_sqrttrans, function(x) getNMDSdata(x %>% filter(site == "SFE-M"), start_col, ASV = TRUE))
NMDS_sfe_list$`relativized` <- lapply(data_relativized, function(x) getNMDSdata(x %>% filter(site == "SFE-M"), start_col, ASV = TRUE))
NMDS_sfe_list$`hellinger` <- lapply(data_hellinger, function(x) getNMDSdata(x %>% filter(site == "SFE-M"), start_col, ASV = TRUE))

# viewing plots against each other
for(i in 1:length(data_sqrttrans)) {
  plots <- list()
  for(j in 1:length(NMDS_sfe_list)) {
    plots[[j]] = makeNMDSplot(NMDS_sfe_list[[j]][[i]], FALSE, FALSE, shape = "atx_group", color = "atx_detected") +
      labs(title = paste(names(data)[i], names(NMDS_sfe_list)[j], sep = " "))
  }
  print(plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol = 2))
}
# visuals seem to support PERMANOVA across data transformations

## (b) RUS
NMDS_rus_list <- list()
NMDS_rus_list$`untransformed` <- lapply(list(data$PICRUSt2_predicted_KO_all_nt,
                                             data$PICRUSt2_predicted_KO_all_tac_noanacyl.csv,
                                             data$PICRUSt2_predicted_KO_all_tac_w_AC), function(x) getNMDSdata(x %>% filter(site == "RUS"), start_col, ASV = TRUE))
NMDS_rus_list$`sqrt_transformed` <- lapply(list(data_sqrttrans$PICRUSt2_predicted_KO_all_nt,
                                                data_sqrttrans$PICRUSt2_predicted_KO_all_tac_noanacyl.csv,
                                                data_sqrttrans$PICRUSt2_predicted_KO_all_tac_w_AC), function(x) getNMDSdata(x %>% filter(site == "RUS"), start_col, ASV = TRUE))
NMDS_rus_list$`relativized` <- lapply(list(data_relativized$PICRUSt2_predicted_KO_all_nt,
                                           data_relativized$PICRUSt2_predicted_KO_all_tac_noanacyl.csv,
                                           data_relativized$PICRUSt2_predicted_KO_all_tac_w_AC), function(x) getNMDSdata(x %>% filter(site == "RUS"), start_col, ASV = TRUE))
NMDS_rus_list$`hellinger` <- lapply(list(data_hellinger$PICRUSt2_predicted_KO_all_nt,
                                         data_hellinger$PICRUSt2_predicted_KO_all_tac_noanacyl.csv,
                                         data_hellinger$PICRUSt2_predicted_KO_all_tac_w_AC), function(x) getNMDSdata(x %>% filter(site == "RUS"), start_col, ASV = TRUE))
# with the relative data, seem to be having an issue!

# viewing plots against each other
for(i in c(1:3)) {
  names = c("NT", "TAC_noanacyl", "TAC_w_AC")
  plots <- list()
  for(j in 1:length(NMDS_rus_list)) {
    plots[[j]] = makeNMDSplot(NMDS_rus_list[[j]][[i]], FALSE, FALSE, shape = "atx_detected", color = "atx_group") +
      labs(title = paste(names[i], names(NMDS_rus_list)[j], sep = " "))
  }
  print(plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol = 2))
}
# visuals seem to support PERMANOVA which is nothing is grouping!

# decided to use square-root transform to preserve predicted abudances but limit impact of
# orthologs with very high counts!
# (plus have issue with relativized data)
