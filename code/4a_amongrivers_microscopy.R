#### Comparing microscopy data among rivers
### Jordan Zabrecky
## last edited: 01.05.2026

# This code compares microscopy data from NT, TM, and TAC samples
# across rivers to answer Q1. First data is transformed (sqrt).
# Data is analyzed using NMDS, PERMANOVA, and ISA. We also averaged across all samples
# from a river and created bar plots to visually compare average samples at each river

#### (1) Loading libraries & data ####

# set seed for reproducibility
set.seed(2025)

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

# see our exploration with data transformation in another script, "S4a_testing_data_transformations.R"
# decided on square-root transformation on the relative abundances (Hellinger transformation)
data <- unaltered_data
for(i in 1:length(data)) {
  data[[i]][,start_col:ncol(data[[i]])] <- sqrt(data[[i]][,start_col:ncol(data[[i]])])
}

# save this data to use in future scripts (RUN ONCE)
#write.csv(data$nt, "./data/morphological/transformed/nt_algalonly_sqrttransformed.csv",
#          row.names = FALSE)
#write.csv(data$tm, "./data/morphological/transformed/tm_algalonly_nomicro_sqrttransformed.csv",
#          row.names = FALSE)
#write.csv(data$tac, "./data/morphological/transformed/nt_algalonly_noanacylgreenalgae_sqrttransformed.csv",
#          row.names = FALSE)

# create a longer version of the unaltered data for bar plots of relative abundances
data_longer <- lapply(unaltered_data, 
                      function(x) pivot_longer(x, cols = all_of(c(start_col:ncol(x))), values_to = "percent",
                                                           names_to = "taxa"))

# add in broader group classification
# grouping for TM & TAC
for(i in 2:length(data_longer)) {
  data_longer[[i]] <- data_longer[[i]] %>% 
    mutate(broader = case_when(taxa == "lyngbya" | taxa == "nodularia" |  taxa == "calothrix" |
                                 taxa == "scytonema" | taxa == "gloeotrichia" ~ "Other N-fixing Cyanobacteria",
                               taxa == "nostoc" ~ "Nostoc",
                               taxa == "chroococcus" | taxa == "other_coccoids"
                               ~ "Unicellullar Cyanobacteria",
                               taxa == "anabaena_and_cylindrospermum" ~ "Anabaena or Cylindrospermum",
                               taxa == "e_diatoms" ~ "Epithemia",
                               taxa == "geitlerinema" ~ "Other Anatoxin-Associated Cyanobacteria",
                               taxa == "green_algae" ~ "Green Algae",
                               taxa == "oscillatoria" | taxa == "phormidium_unknown" |
                                 taxa == "leptolyngbya" | taxa == "homoeothrix"
                               ~ "Other Filamentous Cyanobacteria",
                               taxa == "microcoleus" ~ "Microcoleus",
                               taxa == "non_e_diatoms" ~ "Diatoms Other than Epithemia",
                               taxa == "unknown" ~ "Unknown"
    ))
}

# grouping for NT
data_longer$nt <- data_longer$nt %>% 
  mutate(broader = case_when(taxa == "lyngbya" | taxa == "nodularia" |  taxa == "calothrix" |
                               taxa == "scytonema" | taxa == "gloeotrichia" | taxa == "rivularia" |
                               taxa == "tolypothrix"
                             ~ "Other N-fixing Cyanobacteria",
                             taxa == "nostoc" ~ "Nostoc",
                             taxa == "chroococcus" | taxa == "other_coccoids" | taxa == "aphanothece"
                             ~ "Unicellullar Cyanobacteria",
                             taxa == "anabaena_and_cylindrospermum" ~ "Anabaena or Cylindrospermum",
                             taxa == "epithemia" ~ "Epithemia",
                             taxa == "geitlerinema" ~ "Other Anatoxin-Associated Cyanobacteria",
                             taxa == "oscillatoria" | taxa == "phormidium_unknown" |
                               taxa == "leptolyngbya" | taxa == "homoeothrix"
                             ~ "Other Filamentous Cyanobacteria",
                             taxa == "microcoleus" ~ "Microcoleus",
                             taxa == "non_e_r_diatoms" ~ "Diatoms Other than Epithemia or Rhopalodia",
                             taxa == "unknown" | taxa == "chantransia" | taxa == "euglenoid" |
                               taxa == "unknown_green_algae"
                             ~ "Other",
                             taxa == "ankistrodesmus" | taxa == "gloeocystis" | taxa == "lacunastrum" | 
                               taxa == "oocystis" | taxa == "pediastrum" | taxa == "scenedesmus_no_spines" |
                               taxa == "stauridium" | taxa == "tetraedron" | taxa == "coelastrum" |
                               taxa == "cosmarium" | taxa == "desmodesmus_spines" | taxa == "closterium"
                             ~ "Unicellular Green Algae",
                             taxa == "cladophora" ~ "Cladophora",
                             taxa == "mougeotia" | taxa == "ulothrix" | taxa == "zygnema" |
                               taxa == "stigeoclonium"
                             ~ "Other Filamentous Green Algae",
                             taxa == "oedogonium" ~ "Oedogonium",
                             taxa == "rhopalodia" ~ "Rhopalodia",
                             taxa == "spirogyra" ~ "Spirogyra"
  ))

#### (2) Functions for Analyses ####

# load from supplemental script
source("./code/supplemental_code/S4a_community_analyses_func.R")
source("./code/supplemental_code/S4c_barplot_func.R")

# summarize function
# @param data_long is relative abundance data in long format
# @param grouping is either "taxa" or "broader" groups
summarize_site <- function(data_long, grouping) {
  data = data_long %>% 
    dplyr::group_by(site, .data[[grouping]]) %>%
    dplyr::summarize(avg_percent = mean(percent)) %>% 
    arrange(-avg_percent)
  
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
NMDS_list <- lapply(data, function(x) getNMDSdata(x, start_col))

# making plots
NMDS_plots <- lapply(NMDS_list, function(x) makeNMDSplot(x, TRUE, TRUE, 
                                                         color = "site", shape = "month"))

lapply(NMDS_plots, print)
# RESULT: TM and NT groups are visually distinct among rivers, but not for TAC

#### (6) Q: Are communities from each river significantly different? (PERMANOVA) ####

# run PERMANOVAs
permanovas <- lapply(data, function(x) runPERMANOVA(x, start_col, x$`site`))
lapply(permanovas, print)
# RESULTS: significant difference for TM and NT across rivers, but not TAC

# check dispersion to see if that influences results
for(i in 1:length(data)) {
  print(names(data)[i])
  print(anova(betadisper(vegdist(data[[i]][,start_col:ncol(data[[i]])], method = "bray"), 
                         data[[i]]$site)))
}
# TAC not significantly different, but TM and NT are
# however, based on https://www.youtube.com/watch?v=oLf0EpMJ4yA
# and his paper https://www.nature.com/articles/ismej20085
# this may not affect results of adonis2, especially if NMDS shows that groups are very far apart
# which we do see in our NMDS plots (With the exception maybe of the NT plot, but the centroids
# for those groups are different)

# let's lastly compare the distance between centroids
# (mostly to compare for Q2 issue that came up)
centroid_distance <- lapply(NMDS_list, function(x) { 
  # calculate centroids
  centroids = x[[1]] %>% 
    dplyr::group_by(site) %>% 
    dplyr::summarize(axis1 = mean(NMDS1),
                     axis2 = mean(NMDS2)) %>% 
    ungroup()
  # calculate distances between centroids
  distances = dist(centroids[,2:3], method = "euclidean")
  return(mean(distances))}
  
)
lapply(names(centroid_distance), function(x) print(paste(x, ": ", centroid_distance[[x]], sep = "")))
# nt: 0.658672167241778
# tm: 0.684185197772008
# tac: 0.107665714801558

#### (7) Q: What explains these differences? Loading & Species Indicator Analyses ####

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
summary(multipatt(data$nt[,start_col:ncol(data$nt)], data$nt$site, func = "r.g", control = how(nperm = 999)))
# SAL: homoethrix, leptolyngbya, coccoids, unknown green algae
# SFE: cladophora, stauridium, nostoc, coelastrum, unknown, tetraedron, cosmarium, rivularia,
# ankistrodesmus, lacunastrum
# RUS: mougeotia, phormidium
# RUS + SAL: diatoms, stigeoclonium
# RUS + SFE: spirogyra, epithemia, anabaena, scenedesmus, odeogonium, rhopalodia
# SAL + SFE: microcoleus

# (ii) TM
summary(multipatt(data$tm[,start_col:ncol(data$tm)], data$tm$site, func = "r.g", control = how(nperm = 999)))
# SAL: diatoms
# SFE: anabaena, nostoc

# (iii) TAC
# omit single salmon sample
tac_sub <- data$tac %>% filter(site != "SAL")
summary(multipatt(tac_sub[,start_col:ncol(tac_sub)], tac_sub$site, func = "r.g", control = how(nperm = 999)))
# RUS: phormidium, oscillatoria
# SFE: nodularia, microcoleus

#### (8) Conclusions ####

## In conclusion, 
## come back & revisit this to rewrite this