#### Comparing microscopy data among rivers
### Jordan Zabrecky
## last edited: 12.17.2025

## This code compares microscopy data from NT, TM, and TAC samples
## across rivers to answer Q1. First data is transformed (sqrt). For reference
## we also complete analyses on untransformed data and data with rare taxa removed.
## Data is analyzed using NMDS, PERMANOVA, and ISA. We also averaged across all samples
## from a river and created bar plots to visually compare average samples at each river

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

# see our exploration with data transformation in another script, "S4a_testing_data_transformations.R"
# decided on square-root transformation on the relative abundances (Hellinger transformation)
data <- unaltered_data
for(i in 1:length(data)) {
  data[[i]][,5:ncol(data[[i]])] <- sqrt(data[[i]][,5:ncol(data[[i]])])
}

# save this data to use in future scripts (RUN ONCE)
#write.csv(data$nt, "./data/morphological/transformed/nt_algalonly_sqrttransformed.csv",
#          row.names = FALSE)
#write.csv(data$tm, "./data/morphological/transformed/nt_algalonly_nomicro_sqrttransformed.csv",
#          row.names = FALSE)
#write.csv(data$tac, "./data/morphological/transformed/nt_algalonly_noanacylgreenalgae_gsqrttransformed.csv",
#          row.names = FALSE)

# create a longer version of the unaltered data for bar plots of relative abundances
data_longer <- lapply(unaltered_data, 
                      function(x) pivot_longer(x, cols = c(5:ncol(x)), values_to = "percent",
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

#### (4) NMDS Plots ####

# get NMDS for each dataframe (sqrt-transformed!)
NMDS_list <- lapply(data, function(x) getNMDSdata(x, 5))

# making plots
NMDS_plots <- lapply(NMDS_list, function(x) makeNMDSplot(x, TRUE, TRUE))

# compare with non-transformed data or sqrt-transformed w/ rare taxa removed (get data and then plots)
NMDS_list_nontransformed <- lapply(data_wide, function(x) getNMDSdata(x))
NMDS_list_raretaxaremoved <- lapply(data_filtered_wide, function(x) getNMDSdata(x))
NMDS_plots_nontransformed <- lapply(NMDS_list_nontransformed, function(x) makeNMDSplot(x, TRUE, TRUE))
NMDS_plots_raretaxaremoved <- lapply(NMDS_list_raretaxaremoved, function(x) makeNMDSplot(x, TRUE, TRUE))

# viewing plots against each other
for(i in 1:length(NMDS_plots)) {
  print(plot_grid(NMDS_plots[[i]] + labs(title = titles[i]), 
            NMDS_plots_nontransformed[[i]] + labs(title = "non-transformed"),
            NMDS_plots_raretaxaremoved[[i]] + labs(title = "rare taxa removed"), ncol = 1))
}

#### (6) Q: What is dominant taxa across samples? ####

# get summaries for each taxa
summaries <- lapply(data, function(x) summarize_site(x))

#### (7) Q: Are communities from each river significantly different? (PERMANOVA) ####

# run PERMANOVAs (on square-root transformed, unaltered, and sqrt-transformed w/ rare taxa removed)
permanovas <- lapply(data_wide_sqrt, function(x) runPERMANOVA(x))
permanovas_nontransformed <- lapply(data_wide, function(x) runPERMANOVA(x))
permanovas_raretaxaremoved <- lapply(data_filtered_wide, function(x) runPERMANOVA(x))

# create summary table
permanova_summaries <- data.frame(data = NA,
                                  altering = NA,
                                  significant = NA)
for(i in 1:length(permanovas)) {
  # make dataframe for all PERMANOVAs at index i
  temp <- data.frame(data = c(names(permanovas)[i], names(permanovas_nontransformed)[i], names(permanovas_raretaxaremoved)[i]),
             altering = c("sqrt-transformed", "unaltered", "sqrt-transformed & rare taxa removed"),
             significant = c(permanovas[[i]]$`Pr(>F)`[1], permanovas_nontransformed[[i]]$`Pr(>F)`[1], permanovas_raretaxaremoved[[i]]$`Pr(>F)`[1]))
  
  # add to existing dataframe
  permanova_summaries <- rbind(permanova_summaries, temp)
}

view(permanova_summaries)
# TM & NT significant across all while TAC not significant across all!

# strata test
adonis2(vegdist(data_wide_sqrt$nt_algalonly.csv[,6:ncol(data_wide_sqrt$nt_algalonly.csv)], method = "bray") ~ site, 
        data = data_wide_sqrt$nt_algalonly.csv,
        strata = data_wide_sqrt$nt_algalonly.csv$field_date)

## checking beta dispersion
for(i in 1:length(data_wide_sqrt)) {
  print(anova(betadisper(vegdist(data_wide_sqrt[[i]][,6:ncol(data_wide_sqrt[[i]])], method = "bray"), 
                         data_wide_sqrt[[i]]$site)))
}
# there is a significant difference in beta-dispersion for all
# however, based on https://www.youtube.com/watch?v=oLf0EpMJ4yA
# and his paper https://www.nature.com/articles/ismej20085
# this may not affect results of adonis2, especially if NMDS shows that groups are very far apart
# which we definitely see for Microcoleus and for Non-target which are our significant samples via adonis2

#### (8) Q: What explains these differences? ####

# get p-values from envfit in vegan ran earlier
envfit_pvalues <- lapply(NMDS_list, function(x) as.data.frame(x$vs$vectors$pvals))

# view results from dataframes we care about
view(envfit_pvalues$tm_algalonly_nomicro.csv)
view(envfit_pvalues$tac_algalonly_noanacylgreenalgae.csv)
view(envfit_pvalues$nt_algalonly.csv)

#### (9) Conclusions ####

## In conclusion, 
## (1) Microcoleus samples were significantly different among rivers with those from the South Fork Eel 
## characterized by more Epithemia and Nitrogen-fixers (predominantly Anabaena & Nostoc)
## (2) Anabaena/Cylindrospermum were NOT significantly different among rivers, all rivers had samples with
## Epithemia and other anatoxin-associated cyanobacteria, particularly Geitlerinema
## (3) Non-target communities were significantly different among rivers with South Fork Eel having more
## Nostoc, Salmon having less green algae (could be due to sampling interruption)

# save sqrt-transformed values for future analyses:
for(i in 1:length(data_wide)) {
  write.csv(data_wide[[i]], paste("./data/morphological/transformed/", 
                                  str_replace(names(data_wide)[i], ".csv", ""),
                                  "_sqrttransformed.csv", sep = ""),
            row.names = FALSE)
}
