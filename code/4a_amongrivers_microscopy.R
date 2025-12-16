#### Comparing microscopy data among rivers
### Jordan Zabrecky
## last edited: 12.15.2025

## This code compares microscopy data from NT, TM, and TAC samples
## across rivers to answer Q1. First data is transformed (sqrt). For reference
## we also complete analyses on untransformed data and data with rare taxa removed.
## Data is analyzed using NMDS and ANOSIM. We also averaged across all samples
## from a river and created bar plots to visually compare average samples at each river

#### (1) Loading libraries & data ####

# set seed for reproducibility
set.seed(2025)

# libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot"), require, character.only = T)

# read in files (note: two target taxa csvs- one with target taxa included, other with it excluded)
# doing as a list rather than one csv as each has a different # of columns
files <- list.files(path = "./data/morphological/", pattern = ".csv")
data_wide <- lapply(files, function(x) read.csv(paste("./data/morphological/", x, sep = "")))
names(data_wide) <- files

#### (2) Data Transformation

## (a) pivot longer and sqrt-transform % values

# focus on 2022 data for the three river comparison
data_wide <- lapply(data_wide, function(x) x %>% mutate(year = year(ymd(field_date))) %>% filter(year == 2022) %>% 
                      relocate(year, .before = "field_date"))

# beforehand, curious if any taxa has 0 for all samples
# (adding 5 to match indexing that starts at 6 for colSums)
notaxa <- lapply(data_wide, function(x) colnames(x)[5 + c(which(colSums(x[,6:ncol(x)]) == 0))])

# remove these taxa (that are not present in any samples) from the dataframes
for(i in 1:length(data_wide)) {
  data_wide[[i]] <- data_wide[[i]] %>% 
    dplyr::select(!c(notaxa[[i]]))
}

# pivot longer
data <- lapply(data_wide, function(x) x %>% pivot_longer(cols = c(6:ncol(x)), names_to = "taxa",
                                                    values_to = "percent"))

# histogram of raw percent (highly right-skew)
lapply(data, function(x) hist(x$percent))

# square root transform for multivariate analyses since data is highly right-skewed
# also keeps values positive for Bray-Curtis
data <- lapply(data, function(x) x %>% 
                          mutate(sqrt_percent = sqrt(percent)))
data_wide_sqrt <- lapply(data, function(x)
  pivot_wider(x %>% select(!percent), names_from = taxa, values_from = sqrt_percent))

# updated histogram - still have zero-inflation, but otherwise closer to normally-distributed
lapply(data, function(x) hist(x$sqrt_percent, breaks=seq(0,10,l=20)))

## (b) make broader groups

# grouping for TM & TAC
for(i in 2:length(data)) {
  data[[i]] <- data[[i]] %>% 
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
data$nt_algalonly.csv <- data$nt_algalonly.csv %>% 
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

## (c) consider removing "rare" taxa

# consider removing taxa that compose less than 1% across all samples
rare_taxa <- lapply(data, function(x) x %>% 
    dplyr::group_by(taxa) %>% 
    dplyr::summarize(max = max(percent),
                   min = min(percent),
                   mean = mean(percent)) %>% 
    filter(max < 1))

# separate list of dataframes with rare taxa removed
data_filtered <- list()
for(i in 1:length(data)) {
  data_filtered[[i]] <- data[[i]] %>% 
    filter(! taxa %in% rare_taxa[[i]]$taxa)
}
names(data_filtered) <- names(data)

# histogram of data with "rare" taxa removed (better but still zero-inflation!)
lapply(data_filtered, function(x) hist(x$sqrt_percent,  breaks=seq(0, 10,l=20)))

# make a wider version
data_filtered_wide <- lapply(data, function(x) colnames(x))
data_filtered_wide <- lapply(data, function(x)
  pivot_wider(x %>% select(!c(percent, broader)), names_from = taxa, values_from = sqrt_percent))
  
#### (3) Functions for Analysis ####

# set universal plot theme
theme_set(theme_bw() + theme(panel.grid = element_blank(),
                             panel.border = element_rect(fill = NA, color = "black"),
                             legend.position = "right"))

# bar plot by site fill taxa (argument data is data in long format))
barplot_taxa <- function(data) {
  plot = ggplot(data = data, aes(x = site, y = percent, fill = taxa)) +
    geom_bar(position = "fill", stat = "identity")
  
  return(plot)
}

# bar plot by site fill broader groupings (argument data is data in long format)
barplot_broader <- function(data) {
  plot = ggplot(data = data, aes(x = site, y = percent, fill = broader)) +
    geom_bar(position = "fill", stat = "identity")
  
  return(plot)
}

# creates NMDS data point coordinates and loadings (argument data is data in wide format)
getNMDSdata <- function(data) {
  # use vegan to calculate NMDS distances
  nmds = metaMDS(as.matrix(data[,6:ncol(data)]),
                 distance = "bray",
                 trymax = 500,
                 autotransform = TRUE)
  # bind x & y positions to site information
  nmds_final = cbind(as.data.frame(scores(nmds, "sites")), 
        data %>% select(site_reach, site, field_date)) %>% 
    mutate(field_date = ymd(field_date),
           year = year(field_date),
           month = as.character(month(field_date)))
  # get loadings for taxa
  vs = envfit(nmds, as.matrix(data[,6:ncol(data)]), perm = 999)
  coord = as.data.frame(scores(vs, "vectors"))
  stress = nmds$stress
  
  # return a named list with both dataframes
  list <- list(nmds_final, vs, coord, stress)
  names(list) = c("nmds", "vs", "coord", "stress")
  return(list)
}

# make NMDS plots (without loadings; data is nmds data, loading is a TRUE/FALSE argument,
# and significant is a TRUE/FALSE argument)
makeNMDSplot <- function(data, loading, significant) {
  
  # separating out data to be able to easily call each
  nmds_data = data$nmds
  stress = data$stress
  loadings = data$coord
  pvalues = as.data.frame(data$vs$vectors$pvals)
  colnames(pvalues) = "pvalue"
  
  # make plot
  plot = ggplot(nmds_data, aes(x = NMDS1, y = NMDS2)) +
    geom_point(aes(color = site, shape = month), size = 4) +
    stat_ellipse(aes(color = site), type = "t", linetype = 2, size = 0.5) +
    scale_color_manual(values = c("SAL" = "#62a7f8",
                                  "SFE-M" = "#416f16", 
                                  "RUS" = "#bdb000")) +
    labs(subtitle = paste("Stress:", round(stress, 3)),
         x = "NMDS Axis 1",
         y = "NMDS Axis 2")
  
  # add in loadings
  if(loading) {
    
    if(significant) {
      loadings = cbind(loadings, pvalues) %>% 
        filter(pvalue < 0.05)
      
    }
    
    plot = plot + geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
                   data = loadings, size =1, alpha = 0.5, colour = "grey30") +
                  geom_text(data = loadings, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
                    fontface = "bold", label = rownames(loadings))
  }
  
  return(plot)
}

# run PERMANOVA test (argument data is for wide format)
runPERMANOVA <- function(data) {
  # create distance matrix based on Bray-Curtis distances
  dist_matrix = vegdist(data[,6:ncol(data)], method = "bray")
  
  # return PERMANOVA test results
  return(adonis2(dist_matrix ~ site, data = data))
}

# summarize taxa abundance across samples within a site (input: data long format)
summarize_site <- function(data) {
  summary = data %>% 
    dplyr::group_by(site, taxa) %>% 
    dplyr::summarize(mean = mean(percent))
}

#### (4) Bar Plots ####

# put bar plots into lists
barplot_taxa_plots <- lapply(data, function(x) barplot_taxa(x))
barplot_broader_plots <- lapply(data, function(x) barplot_broader(x))

# titles for plots
titles <- c("Non-Target Samples", "Anabaena/Cylindrospermum Samples (including)", 
            "Anabaena/Cylindrospermum Samples (excluding)", 
            "Anabaena/Cylindrospermum Samples (also excluding GA)",
            "Microcoleus Samples (including)",
            "Microcoleus Samples (excluding)")

# view plots
for(i in 1:length(barplot_taxa_plots)) {
  print(barplot_taxa_plots[[i]] + labs(title = titles[i]))
  print(barplot_broader_plots[[i]] + labs(title = titles[i]))
}

#### (5) NMDS ####

# get NMDS for each dataframe (sqrt-transformed!)
NMDS_list <- lapply(data_wide_sqrt, function(x) getNMDSdata(x))

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
# TM & NT significant across all while TAC not significant across all!http://127.0.0.1:39173/graphics/plot_zoom_png?width=1184&height=861

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
