#### Playing with data transformations for predicted functional groups
### Jordan Zabrecky
## last edited: 05.15.2026

# Do different data transformations changes our predicted functional group results?
# This script explores: (1) using untransformed and (2) square-root transformed data

# In doing so, we also compare differences in results of targeted samples 
# with the targeted taxa remove, and, in the case of Anabaena samples (TAC), 
# how removing green algae--typically the substrate--alters results

# We look predominantly at PERMANOVAs and Species Indicator Analysis Results
# and then visually assess potential differences in NMDS

#### (1) Loading libraries & data ####

# set seed for reproducibility
set.seed(2025)

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
  # does this do relative abundance??
  # LEFT OFF HERE!!!!
  x[,start_col:ncol(x)] <- decostand(x[,start_col:ncol(x)], method = "total")
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
  temp = data.frame(data = c(names(data)[i], names(data_sqrttrans)[i]),
                    transformation = c("none", "squareroot-transformation"),
                    significant = c(runPERMANOVA(data[[i]], start_col = start_col, group = data[[i]]$site)$`Pr(>F)`[1], 
                                    runPERMANOVA(data_sqrttrans[[i]], start_col = start_col, group = data_sqrttrans[[i]]$site)$`Pr(>F)`[1]))
  
  # add to existing dataframe
  Q1_permanovas <- rbind(Q1_permanovas, temp)
}
view(Q1_permanovas)
# square-root transformation seems to make differences more significant, but not in a way
# that would significantly (lol) alter results

#### (4) Q1 NMDS ####

NMDS_list <- list()
NMDS_list$`untransformed` <- lapply(data, function(x) getNMDSdata(x, start_col, ASV = TRUE))
NMDS_list$`sqrt_transformed` <- lapply(data_sqrttrans, function(x) getNMDSdata(x, start_col, ASV = TRUE))

# viewing plots against each other
for(i in 1:length(data_sqrttrans)) {
  plots <- list()
  for(j in 1:length(NMDS_list)) {
    plots[[j]] = makeNMDSplot(NMDS_list[[j]][[i]], FALSE, FALSE, shape = "month", color = "site") +
      labs(title = paste(names(data)[i], names(NMDS_list)[j], sep = " "))
  }
  print(plot_grid(plots[[1]], plots[[2]], ncol = 1))
}

#### (5) Q3 PERMANOVAs ####

# create summary table
Q3_permanovas <- data.frame(data = NA,
                            atx = NA,
                            transformation = NA,
                            significant = NA)

# run PERMANOVAs for Q3
for(i in 1:length(data)) {
  # make a temporary dataframe for all PERMANOVAs at index i for atx detected
  temp = data.frame(data = c(names(data)[i], names(data_sqrttrans)[i]),
                    atx = c("detected", "detected"),
                    transformation = c("none", "squareroot-transformation"),
                    significant = c(runPERMANOVA(data[[i]], start_col = start_col, group = data[[i]]$atx_detected, na.action = "na.omit")$`Pr(>F)`[1], 
                                    runPERMANOVA(data_sqrttrans[[i]], start_col = start_col, group = data_sqrttrans[[i]]$atx_detected, na.action = "na.omit")$`Pr(>F)`[1]))
  
  # make a temporary dataframe for all PERMANOVAs at index i for atx grouping
  temp2 = data.frame(data = c(names(data)[i], names(data_sqrttrans)[i]),
                   atx = c("group", "group"),
                   transformation = c("none", "squareroot-transformation"),
                   significant = c(runPERMANOVA(data[[i]], start_col = start_col, group = data[[i]]$atx_group, na.action = "na.omit")$`Pr(>F)`[1], 
                                   runPERMANOVA(data_sqrttrans[[i]], start_col = start_col, group = data_sqrttrans[[i]]$atx_group, na.action = "na.omit")$`Pr(>F)`[1]))
  
  # add to existing dataframe
  Q3_permanovas <- rbind(Q3_permanovas, temp, temp2)
}
view(Q3_permanovas)
# all not significant!

# think about square root when comparing select predicted genes- probably wouldn't want to do this?





# microscopy data
files <- list.files(path = "./data/morphological/", pattern = ".csv")
data_wide <- lapply(files, function(x) read.csv(paste("./data/morphological/", x, sep = "")) %>% 
                      mutate(field_date = ymd(field_date),
                             month = month(field_date)) %>% # add month tag
                      relocate(month, .before = "field_date") %>% 
                      filter(year(field_date) == 2022)) # filter for 2022 only data
names(data_wide) <- files

# assemblage results?
# This script explores: (1) using un-transformed relative abundances
# (2) squareroot transformed relative abundances (Hellinger-transformation)
# (3) and Hellinger-transformed data with rare taxa removed

# In doing so, we also compare differences in results of targeted samples 
# with the targeted taxa remove, and, in the case of Anabaena samples (TAC), 
# how removing green algae--typically the substrate--alters results

# We look predominantly at PERMANOVAs and Species Indicator Analysis Results
# and then visually assess potential differences in NMDS
