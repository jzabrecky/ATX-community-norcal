#### Comparing morphologically-identified assemblages among rivers
### Jordan Zabrecky
## last edited: 04.21.2026

# This code compares microscopy data from NT, TM, and TAC samples
# across rivers to answer Q1. First data is transformed (sqrt).
# Data is analyzed using NMDS, PERMANOVA, and ISA. We also averaged across all samples
# from a river and created bar plots to visually compare average samples at each river

#### (1) Loading libraries & data ####

# libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot", "indicspecies"), require, character.only = T)

# we will use the NT data, and TM excluding Microcoleus, and TAC excluding Anabaena and Green Algae
nt <- read.csv("./data/morphological/nt_algalonly.csv")
tm <- read.csv("./data/morphological/tm_algalonly_nomicro.csv")
tac <- read.csv("./data/morphological/tac_algalonly_noanacylgreenalgae.csv")

# add into list
unaltered_data <- list(nt, tm, tac)
names(unaltered_data) <- c("nt", "tm", "tac")

# filter for data from year 2022 (for three river comparison)
unaltered_data <- lapply(unaltered_data, function(x) x %>% filter(year(ymd(field_date)) == 2022))

# set column where abundance data starts in dataframe
start_col <- 5

# remove any taxa where there are no observations at any river for each sample type
# note, need to add 4 to match the subset taken for colSums
no_taxa <- lapply(unaltered_data, function(x) colnames(x)[c((start_col - 1) 
                                                                      + which(colSums(x[,start_col:ncol(x)]) == 0))])
# all taxa recorded in NT, whereas there are quite a few for TM and TAC
# remove these taxa (that are not present in any samples) from the dataframes
data <- lapply(c(1:3), function(x) unaltered_data[[x]] %>% 
                 dplyr::select(!c(no_taxa[[x]])))
names(data) <- names(unaltered_data)

# create a longer version of the unaltered data for bar plots of relative abundances
data_longer <- lapply(data, 
                      function(x) pivot_longer(x, cols = all_of(c(start_col:ncol(x))), values_to = "percent",
                                               names_to = "taxa"))

# see our exploration with data transformation in another script, "S4a_testing_data_transformations.R"
# decided on square-root transformation on the relative abundances (Hellinger transformation)
for(i in 1:length(data)) {
  data[[i]][,start_col:ncol(data[[i]])] <- sqrt(data[[i]][,start_col:ncol(data[[i]])])
}

# save this data to use in future scripts (RUN ONCE)
#write.csv(data$nt, "./data/morphological/transformed/nt_algalonly_sqrttransformed.csv",
#          row.names = FALSE)
#write.csv(data$tm, "./data/morphological/transformed/tm_algalonly_nomicro_sqrttransformed.csv",
#          row.names = FALSE)
#write.csv(data$tac, "./data/morphological/transformed/tac_algalonly_noanacylgreenalgae_sqrttransformed.csv",
#          row.names = FALSE)

# add in broader group classification to longer dataframe
# load in functions from supplemental script
source("./code/supplemental_code/S4c_grouping_func.R")

# use functions
data_longer$tm <- target_broader(data_longer$tm)
data_longer$tac <- target_broader(data_longer$tac)
data_longer$nt <- nontarget_broader(data_longer$nt)

#### (2) Functions for Analyses ####

# load from supplemental script
source("./code/supplemental_code/S4b_community_analyses_func.R")

# summarize function
# @param data_long is relative abundance data in long format
# @param grouping is either "taxa" or "broader" groups
summarize_site <- function(data_long, grouping) {
  data = data_long %>% 
    dplyr::group_by(site, .data[[grouping]]) %>%
    dplyr::summarize(avg_percent = mean(percent)) %>% 
    arrange(-avg_percent) %>% 
    ungroup()
  
  return(data)
}

#### (3) Relative Abundance Bar Plots ####

# put bar plots into lists
barplot_taxa_plots <- lapply(data_longer, function(x) barplot(x, x = "site", y = "percent", fill = "taxa"))
barplot_broader_plots <- lapply(data_longer, function(x) barplot(x, x = "site", y  = "percent", fill = "broader"))

# titles for plots
titles <- c("Non-Target Samples", 
            "Microcoleus Samples (excluding M)",
            "Anabaena/Cylindrospermum Samples (excluding AC & GA)")

# view plots
for(i in 1:length(barplot_taxa_plots)) {
  print(barplot_taxa_plots[[i]] + labs(title = titles[i]))
  print(barplot_broader_plots[[i]] + labs(title = titles[i]))
}

#### (4) Q: What is dominant taxa across samples? ####

# get summaries for each boader group
summaries_broader <- lapply(data_longer, function(x) summarize_site(x, "broader"))
lapply(summaries_broader, function(x) head(x))
# RESULTS:
# NT: diatoms & Microcoleus for Salmon; Spirogyra & Cladorphora for South Fork Eel;
# and Spirogyra and diatoms other than Epithemia for Russian
# TM: diatoms other than Epithemia and green algae for Salmon; Epithemia,
# green algae, diatoms, abd other ATX producers (Geitlerinema/Leptolyngbya) for South Fork Eel
# TAC: diatoms other than Epithemia dominate for all three, then Epithemia

# get summaries for each taxa
summaries_taxa <- lapply(data_longer, function(x) summarize_site(x, "taxa"))
lapply(summaries_taxa, function(x) head(x))
# similar top groupings as above

#### (5) NMDS Plots ####

# get NMDS for each dataframe (sqrt-transformed!)
set.seed(1)
NMDS_list <- lapply(data, function(x) getNMDSdata(x, start_col))

# making plots
NMDS_plots <- lapply(NMDS_list, function(x) makeNMDSplot(x, TRUE, TRUE, 
                                                         color = "site", shape = "month"))

lapply(NMDS_plots, print)
# RESULT: TM and NT groups are visually distinct among rivers, but not for TAC

#### (6) Q: Are communities from each river significantly different? (PERMANOVA) ####

# empty table for permanova outputs
p_table <- data.frame(test = NA,
                      sample_type = NA,
                      p_value = NA,
                      F_stat = NA)

# run PERMANOVAs
set.seed(1)
permanovas <- lapply(data, function(x) runPERMANOVA(data = x, 
                                                    start_col = start_col, 
                                                    group = x$`site`))

# print and add test results to table
for(i in 1:length(permanovas)) {
  
  # print test results to console
  print(names(permanovas[i]))
  print(permanovas[[i]])
  
  # save stats to table
  p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                       sample_type = names(permanovas[i]),
                                       p_value = permanovas[[i]]$`Pr(>F)`[1],
                                       F_stat = permanovas[[i]]$`F`[1]))
}
# RESULTS: significant difference for TM and NT across rivers, but not TAC

# check dispersion to see if that influences results
for(i in 1:length(data)) {
  
  # run ANOVA
  set.seed(1) 
  anova = anova(betadisper(vegdist(data[[i]][,start_col:ncol(data[[i]])], method = "bray"), 
                           data[[i]]$site))
  
  # print results
  print(names(data)[i])
  print(anova)
  
  # add results table
  p_table <- rbind(p_table, data.frame(test = "PERMDISP",
                                       sample_type = names(data)[i],
                                       p_value = anova$`Pr(>F)`[1],
                                       F_stat = anova$`F value`[1]))
}
# TAC not significantly different, but TM and NT are
# however, based on https://www.youtube.com/watch?v=oLf0EpMJ4yA
# and his paper https://www.nature.com/articles/ismej20085
# this may not affect results of adonis2, especially if NMDS shows that groups are very far apart
# which we do see in our NMDS plots (With the exception maybe of the NT plot, but the centroids
# for those groups are different)

# save tests
write.csv(p_table[-1,], "./data/PERMANOVA_results/Q1_microscopy.csv", row.names = FALSE)

#### (7) Q: What explains these differences? Species Indicator Analyses ####

# note: if issues, another package may be masking the function unique() ???
# if so, restart R and run this script only

# run each separately:

# (i) NT
set.seed(1)
nt_test <- multipatt(data$nt[,start_col:ncol(data$nt)], data$nt$site, func = "r.g", control = how(nperm = 999))
summary(nt_test)
write.csv(nt_test$sign, "./data/ISA_results/Q1_nt_microscopy.csv")
# SAL: homoethrix, leptolyngbya, coccoids, unknown green algae, ulothrix
# SFE: cladophora, stauridium, nostoc, coelastrum, unknown, tetraedron, cosmarium, rivularia,
# ankistrodesmus, lacunastrum, aphanothece
# RUS: mougeotia
# RUS + SAL: diatoms, stigeoclonium
# RUS + SFE: spirogyra, epithemia, anabaena, scenedesmus, odeogonium, rhopalodia
# SAL + SFE: microcoleus

# (ii) TM
set.seed(2)
tm_test <- multipatt(data$tm[,start_col:ncol(data$tm)], data$tm$site, func = "r.g", control = how(nperm = 999))
summary(tm_test)
write.csv(tm_test$sign, "./data/ISA_results/Q1_tm_microscopy.csv")
# SAL: diatoms
# SFE: anabaena, epithemia, nostoc

# (iii) TAC
set.seed(3)
tac_test <- multipatt(data$tac[,start_col:ncol(data$tac)], data$tac$site, func = "r.g", control = how(nperm = 999))
summary(tac_test)
write.csv(tac_test$sign, "./data/ISA_results/Q1_tac_microscopy.csv")
# nothing!

#### (8) Misc. Q's ####

## (1) How many more taxa groups were identified in South Fork Eel samples than Salmon River samples?
lapply(data, function(x) specnumber(x[,start_col:ncol(x)], groups = x$`site`))
# NT: RUS 28, SAL 30, SFE-M 36
# TM: SAL 6, SFE-M 11
# TAC: RUS 10, SAL 6, SFE-M 13

## (2) For target taxa, what is the range of the % the target taxa is present?
# (plus green algae abundance for A/C samples)
tm_w_target <- read.csv("./data/morphological/tm_algalonly.csv")
min(tm_w_target$microcoleus) # 30.4%
max(tm_w_target$microcoleus) # 94.8%

tac_w_target <- read.csv("./data/morphological/tac_algalonly.csv")
min(tac_w_target$anabaena_and_cylindrospermum) # 13.2%
max(tac_w_target$anabaena_and_cylindrospermum) # 75.7%
min(tac_w_target$green_algae) # 2.7%
max(tac_w_target$green_algae) # 70.7%

## (3) If we group broader for NT samples, what is most abundant for each river? 
# five groups: diatom, spirogyra, cladophora, diazotrophic cyanos, non-diazo filamentous cyanos,
# other green algae, coccoidal cyanobacteria
# join diatoms, spirogyra, cladophora, diazotrophic cyanobacteria, filamentous cyanobacteria (non-diaztotrophic),
# other green algae... so five griyo
even_broader_NT <- data_longer$nt %>% 
  mutate(even_broader = case_when(broader == "Cladophora" ~ "Cladophora",
                                  broader == "Spirogyra" ~ "Spirogyra",
                                  broader == "Epithemia or Rhopalodia" ~ "Diatoms",
                                  broader == "Diatoms Other than Epithemia or Rhopalodia" ~ "Diatoms",
                                  broader == "Nostoc" ~ "Diazotrophic Cyanobacteria",
                                  broader == "Anabaena or Cylindrospermum" ~ "Diazotrophic Cyanobacteria",
                                  broader == "Other N-fixing Cyanobacteria" ~ "Diazotrophic Cyanobacteria",
                                  broader == "Unicellular Cyanobacteria" ~ "Coccoidal Cyanobacteria",
                                  broader == "Microcoleus" ~ "Other Filamentous Cyanobacteria",
                                  broader == "Other Green Algae" ~ "Other Green Algae",
                                 TRUE ~ broader)) %>% 
  # merge groups for total in each broader group (i.e., reduce rows)
  dplyr::group_by(site, site_reach, field_date, even_broader) %>% 
  dplyr::summarize(total = sum(percent)) %>% 
  dplyr::ungroup() %>% 
  # regroup to calculate average per site across all samples
  dplyr::group_by(site, even_broader) %>% 
  dplyr::summarize(mean = mean(total))
view(even_broader_NT)
