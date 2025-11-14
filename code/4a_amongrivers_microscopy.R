#### Comparing microscopy data among rivers
### Jordan Zabrecky
## last edited: 11.12.2025

#### (1) Loading libraries & data ####

# set seed
set.seed(2025)

# libraries
lapply(c("tidyverse", "plyr", "vegan"), require, character.only = T)

# read in files (note: two target taxa csvs- one with target taxa included, other with it excluded)
# doing as a list rather than one csv as each has a different # of columns
files <- list.files(path = "./data/morphological/", pattern = ".csv")
data_wide <- lapply(files, function(x) read.csv(paste("./data/morphological/", x, sep = "")))
names(data_wide) <- files

#### (2) Data Transformation

## (a) pivot longer and log-transform % values

# focus on 2022 data for the three river comparison
data_wide <- lapply(data_wide, function(x) x %>% mutate(year = year(ymd(field_date))) %>% filter(year == 2022) %>% 
                      relocate(year, .before = "field_date"))

# beforehand, curious if any taxa has 0 for all samples
# (adding 5 to match indexing that starts at 6 for colSums)
notaxa <- lapply(data_wide, function(x) colnames(x)[5 + c(which(colSums(x[,6:ncol(x)]) == 0))])

# remove these taxa from the dataframes
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

# log-transform percent values adding a small amount for zero for both long & wide dataframes
data <- lapply(data, function(x) x %>% 
                          mutate(sqrt_percent = sqrt(percent)))
data_wide_sqrt <- lapply(data, function(x)
  pivot_wider(x %>% select(!percent), names_from = taxa, values_from = sqrt_percent))

# updated histogram - still have zero-inflation, but otherwise closer to normally-distributed
lapply(data, function(x) hist(x$sqrt_percent, breaks=seq(0,10,l=20)))

## (b) make broader groups

# grouping for TM & TAC
for(i in 2:5) {
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
data_filtered_wide <- lapply(data_filtered, function(x)
  pivot_wider(x %>% select(!percent), names_from = taxa, values_from = sqrt_percent))
  
#### (3) Functions for Analysis ####

# set universal plot theme
theme_set(theme_bw() + theme(panel.grid = element_blank(),
                             panel.border = element_rect(fill = NA, color = "black"),
                             legend.position = "right"))

# bar plot by site fill taxa
barplot_taxa <- function(data) {
  plot = ggplot(data = data, aes(x = site, y = percent, fill = taxa)) +
    geom_bar(position = "fill", stat = "identity")
  
  return(plot)
}

# bar plot by site fill broader groupings
barplot_broader <- function(data) {
  plot = ggplot(data = data, aes(x = site, y = percent, fill = broader)) +
    geom_bar(position = "fill", stat = "identity")
  
  return(plot)
}

getNMDSdata <- function(data) {
  # use vegan to calculate NMDS distances
  nmds = metaMDS(as.matrix(data_wide_sqrt$tm_algalonly.csv[,6:ncol(data_wide_sqrt$tm_algalonly.csv)]),
                 distance = "bray",
                 trymax = 500,
                 autotransform = TRUE)
  # bind x & y positions to site information
  nmds_final = cbind(as.data.frame(scores(tm_nmds, "sites")), 
        data_wide_sqrt$tm_algalonly.csv %>% select(site_reach, site, field_date)) %>% 
    mutate(field_date = ymd(field_date),
           year = year(field_date),
           month = as.character(month(field_date)))
  # get loadings for taxa
  
  return(list(nmds, nmds_final))
}

getloadings <- function()


#### (4) Bar Plots ####

# make empty lists to hold bar plots
barplot_taxa_list <- list()
barplot_broader_list <- list()

for(i in 1:length(data)) {
  barplot_taxa_list[[i]] <- barplot_taxa(data[[i]])
  barplot_broader_list[[i]] <- barplot_broader(data[[i]])
  #print(barplot_taxa_list[[i]])
  print(barplot_broader_list[[i]])
}

#### (4) NMDS

#### (3) TM Samples ####

## (a) average across site bar plot

## with Microcoleus included
# NMDS data
tm_nmds <- metaMDS(as.matrix(data_wide_sqrt$tm_algalonly.csv[,6:ncol(data_wide_sqrt$tm_algalonly.csv)]),
                     distance = "bray",
                     trymax = 500,
                     autotransform = TRUE)
tm_nmds_data <- cbind(as.data.frame(scores(tm_nmds, "sites")), 
                      data_wide_sqrt$tm_algalonly.csv %>% select(site_reach, site, field_date)) %>% 
  mutate(field_date = ymd(field_date),
         year = year(field_date),
         month = as.character(month(field_date)))

# get loadings
tm_vf <- envfit(tm_nmds, as.matrix(data_wide_sqrt$tm_algalonly.csv[,6:ncol(data_wide_sqrt$tm_algalonly.csv)]), perm = 999)
tm_coord <- as.data.frame(scores(tm_vf, "vectors")) * ordiArrowMul(tm_vf)

# plot NMDS
tm_nmds_plot <- ggplot(tm_nmds_data, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = site, shape = month), size = 4) +
  stat_ellipse(aes(color = site), type = "t", linetype = 2, size = 0.5) +
  scale_color_manual(values = c("SAL" = "#62a7f8",
                                "SFE-M" = "#416f16")) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = tm_coord, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = tm_coord, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(tm_coord)) +
  labs(title = "NMDS of Morphologically-Identified Microcoleus Mat Community",
       subtitle = paste("Stress:", round(tm_nmds$stress, 3)),
       x = "NMDS Axis 1",
       y = "NMDS Axis 2")
tm_nmds_plot

##  without Microcoleus included
# NMDS data
tm_nmds_nomicro <- metaMDS(as.matrix(data_wide_sqrt$tm_algalonly_nomicro.csv[,6:ncol(data_wide_sqrt$tm_algalonly_nomicro.csv)]),
                   distance = "bray",
                   trymax = 500,
                   autotransform = TRUE)
tm_nmds_data_nomicro <- cbind(as.data.frame(scores(tm_nmds, "sites")), 
                      data_wide_sqrt$tm_algalonly_nomicro.csv %>% select(site_reach, site, field_date)) %>% 
  mutate(field_date = ymd(field_date),
         year = year(field_date),
         month = as.character(month(field_date)))

# get loadings
tm_vf_nomicro <- envfit(tm_nmds_nomicro, as.matrix(data_wide_sqrt$tm_algalonly_nomicro.csv[,6:ncol(data_wide_sqrt$tm_algalonly_nomicro.csv)]), perm = 999)
tm_coord_nomicro <- as.data.frame(scores(tm_vf_nomicro, "vectors")) * ordiArrowMul(tm_vf_nomicro)

# plot NMDS
tm_nmds_plot_nomicro <- ggplot(tm_nmds_data_nomicro, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = site, shape = month), size = 4) +
  stat_ellipse(aes(color = site), type = "t", linetype = 2, size = 0.5) +
  scale_color_manual(values = c("SAL" = "#62a7f8",
                                "SFE-M" = "#416f16",
                                "SFE-SH" = "#a8ff82")) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = tm_coord_nomicro, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = tm_coord_nomicro, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(tm_coord_nomicro)) +
  labs(title = "NMDS of Morphologically-Identified Microcoleus Mat Community",
       subtitle = paste("Stress:", round(tm_nmds$stress, 3)),
       x = "NMDS Axis 1",
       y = "NMDS Axis 2")
tm_nmds_plot_nomicro

## (c) misc. questions

# Excluding Microcoleus, what are the dominant taxa across all samples?
tm_summary <- data$tm_algalonly_nomicro.csv %>% 
  dplyr::group_by(site, taxa) %>% 
  dplyr::summarize(mean = mean(percent))
# For Salmon TM: non-e diatoms, green algae, other coccoids
# For SfkEel TM: Epithemia, green algae, non-e diatoms, and has many other cyanos that SAL TM does not

# Are the communities from each river significantly different via PERMANOVA?
tm_dist <- vegdist(data_wide_sqrt$tm_algalonly.csv[6:ncol(data_wide_sqrt$tm_algalonly.csv)], method = "bray")
tm_permanova <- adonis2(tm_dist ~ site, data = data_wide_sqrt$tm_algalonly.csv)  # 0.001 ***
tm_dist_nomicro <- vegdist(data_wide_sqrt$tm_algalonly_nomicro.csv[6:ncol(data_wide_sqrt$tm_algalonly_nomicro.csv)], method = "bray")
tm_permanova_nomicro <- adonis2(tm_dist_nomicro ~ site, data = data_wide_sqrt$tm_algalonly_nomicro.csv)  # 0.001 ***
# Yes!

# How does the above change with rare taxa removed?
tm_dist_raretaxa_nomicro <- vegdist(data_filtered$tm_algalonly_nomicro.csv[6:ncol(data_filtered$tm_algalonly_nomicro.csv)], method = "bray")
tm_permanova_raretaxa_nomicro <- adonis2(tm_dist_raretaxa_nomicro ~ site, data = data_filtered$tm_algalonly_nomicro.csv)  # 0.001 ***

#### (4) TAC Samples ####

## (a) average across site bar plot

# with A/C included
TAC_bar_taxa <- ggplot(data = data$tac_algalonly.csv, aes(x = site, y = percent, fill = taxa)) +
  geom_bar(position = "fill", stat = "identity")
TAC_bar_taxa
TAC_bar_broad <- ggplot(data = data$tac_algalonly.csv, aes(x = site, y = percent, fill = broader)) +
  geom_bar(position = "fill", stat = "identity")
TAC_bar_broad

# without Microcoleus included
TM_bar_taxa_excl <- ggplot(data = data$tm_algalonly_nomicro.csv, aes(x = site, y = percent, fill = taxa)) +
  geom_bar(position = "fill", stat = "identity")
TM_bar_taxa_excl
TM_bar_broad_excl <- ggplot(data = data$tm_algalonly_nomicro.csv, aes(x = site, y = percent, fill = broader)) +
  geom_bar(position = "fill", stat = "identity")
TM_bar_broad_excl

## (b) NMDS

## with Microcoleus included
# NMDS data
tm_nmds <- metaMDS(as.matrix(data_wide_sqrt$tm_algalonly.csv[,6:ncol(data_wide_sqrt$tm_algalonly.csv)]),
                   distance = "bray",
                   trymax = 500,
                   autotransform = TRUE)
tm_nmds_data <- cbind(as.data.frame(scores(tm_nmds, "sites")), 
                      data_wide_sqrt$tm_algalonly.csv %>% select(site_reach, site, field_date)) %>% 
  mutate(field_date = ymd(field_date),
         year = year(field_date),
         month = as.character(month(field_date)))

# get loadings
tm_vf <- envfit(tm_nmds, as.matrix(data_wide_sqrt$tm_algalonly.csv[,6:ncol(data_wide_sqrt$tm_algalonly.csv)]), perm = 999)
tm_coord <- as.data.frame(scores(tm_vf, "vectors")) * ordiArrowMul(tm_vf)

# plot NMDS
tm_nmds_plot <- ggplot(tm_nmds_data, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = site, shape = month), size = 4) +
  stat_ellipse(aes(color = site), type = "t", linetype = 2, size = 0.5) +
  scale_color_manual(values = c("SAL" = "#62a7f8",
                                "SFE-M" = "#416f16")) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = tm_coord, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = tm_coord, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(tm_coord)) +
  labs(title = "NMDS of Morphologically-Identified Microcoleus Mat Community",
       subtitle = paste("Stress:", round(tm_nmds$stress, 3)),
       x = "NMDS Axis 1",
       y = "NMDS Axis 2")
tm_nmds_plot

##  without Microcoleus included
# NMDS data
tm_nmds_nomicro <- metaMDS(as.matrix(data_wide_sqrt$tm_algalonly_nomicro.csv[,6:ncol(data_wide_sqrt$tm_algalonly_nomicro.csv)]),
                           distance = "bray",
                           trymax = 500,
                           autotransform = TRUE)
tm_nmds_data_nomicro <- cbind(as.data.frame(scores(tm_nmds, "sites")), 
                              data_wide_sqrt$tm_algalonly_nomicro.csv %>% select(site_reach, site, field_date)) %>% 
  mutate(field_date = ymd(field_date),
         year = year(field_date),
         month = as.character(month(field_date)))

# get loadings
tm_vf_nomicro <- envfit(tm_nmds_nomicro, as.matrix(data_wide_sqrt$tm_algalonly_nomicro.csv[,6:ncol(data_wide_sqrt$tm_algalonly_nomicro.csv)]), perm = 999)
tm_coord_nomicro <- as.data.frame(scores(tm_vf_nomicro, "vectors")) * ordiArrowMul(tm_vf_nomicro)

# plot NMDS
tm_nmds_plot_nomicro <- ggplot(tm_nmds_data_nomicro, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = site, shape = month), size = 4) +
  stat_ellipse(aes(color = site), type = "t", linetype = 2, size = 0.5) +
  scale_color_manual(values = c("SAL" = "#62a7f8",
                                "SFE-M" = "#416f16",
                                "SFE-SH" = "#a8ff82")) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = tm_coord_nomicro, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = tm_coord_nomicro, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(tm_coord_nomicro)) +
  labs(title = "NMDS of Morphologically-Identified Microcoleus Mat Community",
       subtitle = paste("Stress:", round(tm_nmds$stress, 3)),
       x = "NMDS Axis 1",
       y = "NMDS Axis 2")
tm_nmds_plot_nomicro

## (c) misc. questions

# Excluding Microcoleus, what are the dominant taxa across all samples?
tm_summary <- data$tm_algalonly_nomicro.csv %>% 
  dplyr::group_by(site, taxa) %>% 
  dplyr::summarize(mean = mean(percent))
# For Salmon TM: non-e diatoms, green algae, other coccoids
# For SfkEel TM: Epithemia, green algae, non-e diatoms, and has many other cyanos that SAL TM does not

# Are the communities from each river significantly different via PERMANOVA?
tm_dist <- vegdist(data_wide_sqrt$tm_algalonly.csv[6:ncol(data_wide_sqrt$tm_algalonly.csv)], method = "bray")
tm_permanova <- adonis2(tm_dist ~ site, data = data_wide_sqrt$tm_algalonly.csv)  # 0.001 ***
tm_dist_nomicro <- vegdist(data_wide_sqrt$tm_algalonly_nomicro.csv[6:ncol(data_wide_sqrt$tm_algalonly_nomicro.csv)], method = "bray")
tm_permanova_nomicro <- adonis2(tm_dist ~ site, data = data_wide_sqrt$tm_algalonly_nomicro.csv)  # 0.001 ***
# Yes!


#### (5) NT Samples ####
