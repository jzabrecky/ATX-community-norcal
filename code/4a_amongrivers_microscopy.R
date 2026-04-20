#### Comparing morphologically-identified assemblages among rivers
### Jordan Zabrecky
## last edited: 04.16.2026

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
# green algae and other ATX producers for South Fork Eel
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

# run PERMANOVAs
set.seed(1)
permanovas <- lapply(data, function(x) runPERMANOVA(data = x, 
                                                    start_col = start_col, 
                                                    group = x$`site`))
lapply(permanovas, print)
# RESULTS: significant difference for TM and NT across rivers, but not TAC

# check dispersion to see if that influences results
for(i in 1:length(data)) {
  print(names(data)[i])
  set.seed(1)
  print(anova(betadisper(vegdist(data[[i]][,start_col:ncol(data[[i]])], method = "bray"), 
                         data[[i]]$site)))
}
# TAC not significantly different, but TM and NT are
# however, based on https://www.youtube.com/watch?v=oLf0EpMJ4yA
# and his paper https://www.nature.com/articles/ismej20085
# this may not affect results of adonis2, especially if NMDS shows that groups are very far apart
# which we do see in our NMDS plots (With the exception maybe of the NT plot, but the centroids
# for those groups are different)

#### (7) Q: What explains these differences? Loadings & Species Indicator Analyses ####

## (a) NMDS loadings
# note these are to be interpretted as supplemental, not explanatory
# get p-values from envfit in vegan ran earlier
envfit_pvalues <- lapply(NMDS_list, function(x) {
  newdata = as.data.frame(x$vs$vectors$pvals)
  newdata = rownames_to_column(newdata, var = "taxa")
  colnames(newdata) = c("taxa", "pval")
  newdata = newdata %>% arrange(pval) })

# view results
lapply(envfit_pvalues, function(x) head(x, 15))
# RESULT:
# NT: Cladophora, Epithemia, Homoethrix, Microcoleus, Mougeotia, Diatoms, Nostoc, 
# Oedogonium, Rhopalodia, Spirogyra, and Stigeoclonium are most influential
# TM: Diatoms, Epithemia, Nostoc, Green Algae, Coccoids, Geilerinema are most influential
# TAC: Epithemia, Microcoleus, Nostoc, Oscillatoria, Phormidium, Diatoms, Geitlerinema
# are most influential

## (b) Species Indicator Analyses

# note: if issues, another package may be masking the function unique() ???
# if so, restart R and run this script only

# run each separately:

# (i) NT
set.seed(1)
nt_test <- multipatt(data$nt[,start_col:ncol(data$nt)], data$nt$site, func = "r.g", control = how(nperm = 999))
summary(nt_test)
#write.csv(nt_test$sign, "./data/ISA_results/Q1_nt_microscopy.csv")
# SAL: homoethrix, leptolyngbya, coccoids, unknown green algae
# SFE: cladophora, stauridium, nostoc, coelastrum, unknown, tetraedron, cosmarium, rivularia,
# ankistrodesmus, lacunastrum
# RUS: mougeotia, phormidium
# RUS + SAL: diatoms, stigeoclonium
# RUS + SFE: spirogyra, epithemia, anabaena, scenedesmus, odeogonium, rhopalodia
# SAL + SFE: microcoleus

# (ii) TM
set.seed(1)
tm_test <- multipatt(data$tm[,start_col:ncol(data$tm)], data$tm$site, func = "r.g", control = how(nperm = 999))
summary(tm_test)
#write.csv(tm_test$sign, "./data/ISA_results/Q1_tm_microscopy.csv")
# SAL: diatoms
# SFE: anabaena, nostoc

# (iii) TAC
set.seed(1)
tac_test <- multipatt(data$tac[,start_col:ncol(data$tac)], data$tac$site, func = "r.g", control = how(nperm = 999))
summary(tac_test)
#write.csv(tac_test$sign, "./data/ISA_results/Q1_tac_microscopy.csv")
# RUS: phormidium, oscillatoria
# SFE: nodularia, microcoleus

#### (8) Misc. Q's ####

## How many more taxa groups were identified in South Fork Eel samples than Salmon River samples?
lapply(data, function(x) specnumber(x[,start_col:ncol(x)], groups = x$`site`))
# NT: RUS 30, SAL 32, SFE-M 37
# TM: SAL 7, SFE-M 12
# TAC: RUS 11, SAL 7, SFE-M 14

## Is Leptolyngbya/Geitlerinema present in all TM samples?
data$tm$leptolyngbya_geitlerinema
count(data$tm$leptolyngbya_geitlerinema > 0)
# present in 17 out of 23

## What about the presence of Geitlerinema in TAC samples?
data$tac$leptolyngbya_geitlerinema
count(data$tac$leptolyngbya_geitlerinema > 0)
# all of them which is crazy >1%!

## How about Microcoleus in TAC samples? 
## (particularly interested in Russian River where we did not obsere M. macroscopically)
data$tac$microcoleus
count(data$tac$microcoleus > 0 & data$tac$site == "RUS")
# true for 22/28
# 10 of those are russian river
count(data$tac$site == "RUS") # of 15 samples

## Finally, let's group a little farther for NT samples
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
                                  broader == "Microcoleus" ~ "Non-diaztrophic Cyanobacteria",
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

### MAY MOVE BELOW TO Q3

## Let's look at only other anatoxin associated taxa in all samples
# using list from Christensen & Khan et al. (2019): Anabaena, Aphanizomenon,
# Aphanothece, Arthospira, Cylindrospermopsis, Cylindrospermum, Gomphosphaeria,
# Limnothrix, Lyngbya, Microcystis, Nostoc, Oscillatoria, Phormidium/Microcoleus,
# Planktothrix, Planktolyngbia, Synechocystis, Psuedoanabaena,
# Raphidopsis, Tychonema
# noting that it doesn't have Geilerinema, so we should also include list
# of ATX producers from Wood et al. (2020) which adds: Fisherella, 
# Geitlerinema, Leptolyngbya, Microseira (formerly Lyngbya), planktothrix
lapply(data, function(x) colnames(x[,5:ncol(x)]))
atx_taxa_only <- lapply(data_longer, function(x) {
  
  # make dataframe with only taxa in list above 
  # (only writing what taxa we recorded from that list)
  df = x %>% 
    filter(taxa %in% c("aphanothece", "anabaena_and_cylindrospermum",
                       "other_coccoids", "geitlerinema", "leptolyngbya", 
                       "lyngbya", "other_coccoids", "nostoc",
                       "oscillatoria", "phormidium_unknown", "microcoleus")) %>% 
    mutate(sample_name = paste("(", month(field_date), "-", day(field_date), ") ", site_reach, sep = ""))
  
  # make bar plot (show each sample individually)
  plot <- ggplot(data = df, aes(x = sample_name, y = percent / 100, fill = taxa)) +
    geom_bar(stat = "identity") +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(x = NULL, y = "Relative Abundance") +
    facet_wrap(~site, scales = "free_x")
  print(plot) # view plot
    
  # return a list including dataframe, then plot
  return(list(df, plot))
})

# NOTE: may shove this to a later script