#### Comparing microscopy data among rivers
### Jordan Zabrecky
## last edited: 11.12.2025

#### (1) Loading libraries & data ####

# libraries
lapply(c("tidyverse", "plyr", "vegan"), require, character.only = T)

# read in files (note: two target taxa csvs- one with target taxa included, other with it excluded)
# doing as a list rather than one csv as each has a different # of columns
files <- list.files(path = "./data/morphological/", pattern = ".csv")
data_wide <- lapply(files, function(x) read.csv(paste("./data/morphological/", x, sep = "")))
names(data_wide) <- files

#### (2) Data Transformation

## (a) pivot longer and log-transform % values

# pivot longer
data <- lapply(data_wide, function(x) x %>% pivot_longer(cols = c(5:ncol(x)), names_to = "taxa",
                                                    values_to = "percent"))

# focus on 2022 data for the three river comparison
data <- lapply(data, function(x) x %>% mutate(year = year(ymd(field_date))) %>% filter(year == 2022) %>% 
                          relocate(year, .before = "field_date"))
data_wide <- lapply(data_wide, function(x) x %>% mutate(year = year(ymd(field_date))) %>% filter(year == 2022) %>% 
                 relocate(year, .before = "field_date"))

# histogram of raw percent (highly right-skew)
lapply(data, function(x) hist(x$percent))

# log-transforming for multivariate analyses since data is highly right-skewed
# minimum nonzero value across all dataframes is...
lapply(data, function(x) min(x$percent[x$percent != 0])) # 0.02079

# beforehand, curious if any taxa has 0 for all samples
notaxa <- lapply(data_wide, function(x) which(colSums(x[,6:ncol(x)]) == 0))
test <- notaxa$tac_algalonly_noanacyl.csv

## REMOVE TAXA WHERE THERE IS NONE
# no aphanothecetac_algalonly_noanacyl.csv# no aphanothece in microcoleus samples
data$tm_algalonly.csv <- data$tm_algalonly.csv %>% filter(taxa != "aphanothece")
data$tm_algalonly_nomicro.csv <- data$tm_algalonly_nomicro.csv %>% filter(taxa != "aphanothece")
data_wide$tm_algalonly.csv <- data_wide$tm_algalonly.csv %>% filter(taxa != "aphanothece")
data_wide$tm_algalonly_nomicro.csv <- data_wide$tm_algalonly_nomicro.csv %>% filter(taxa != "aphanothece")
 
# log-transform percent values adding a small amount for zero
data <- lapply(data, function(x) x %>% 
                          mutate(percent = case_when(percent == 0 ~ 0.01,
                                                     TRUE ~ percent),
                          log_percent = log(percent)))

# updated histogram - still have zero-inflation, but otherwise more normally-distributed
lapply(data, function(x) hist(x$log_percent,  breaks=seq(-10,10,l=20)))

## (b) make broader groups

# grouping for TM & TAC
for(i in 2:5) {
  data[[i]] <- data[[i]] %>% 
    mutate(broader = case_when(taxa == "calothrix" | taxa == "lyngbya" | taxa == "nodularia" | 
                       taxa == "rivularia" | taxa == "scytonema"| taxa == "tolypothrix" | 
                         taxa == "gloeotrichia"
                     ~ "Other N-fixing Cyanobacteria",
                     taxa == "nostoc" ~ "Nostoc",
                     taxa == "chroococcus" | taxa == "aphanothece" | taxa == "other_coccoids"
                     ~ "Unicellullar Cyanobacteria",
                     taxa == "anabaena_and_cylindrospermum" ~ "Anabaena or Cylindrospermum",
                     taxa == "e_diatoms" ~ "Epithemia",
                     taxa == "geitlerinema" ~ "Other Anatoxin-Associated Cyanobacteria",
                     taxa == "green_algae" ~ "Green Algae",
                     taxa == "homoeothrix" | taxa == "oscillatoria" | taxa == "phormidium_unknown" |
                       taxa == "leptolyngbya"
                     ~ "Other Filamentous Cyanobacteria",
                     taxa == "microcoleus" ~ "Microcoleus",
                     taxa == "non_e_diatoms" ~ "Diatoms Other than Epithemia",
                     taxa == "unknown" ~ "Unknown"
                     ))
}

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

# histogram of data with "rare" taxa removed (better but still some zero-inflation!)
lapply(data_filtered, function(x) hist(x$log_percent,  breaks=seq(-10,10,l=20)))

#### (3) TM Samples ####

# set universal plot theme
theme_set(theme_bw())

## (a) average across site bar plot

# with Microcoleus included
TM_bar_taxa <- ggplot(data = data$tm_algalonly.csv, aes(x = site, y = percent, fill = taxa)) +
  geom_bar(position = "fill", stat = "identity")
TM_bar_taxa
TM_bar_broad <- ggplot(data = data$tm_algalonly.csv, aes(x = site, y = percent, fill = broader)) +
  geom_bar(position = "fill", stat = "identity")
TM_bar_broad

# without Microcoleus included
TM_bar_taxa_excl <- ggplot(data = data$tm_algalonly_nomicro.csv, aes(x = site, y = percent, fill = taxa)) +
  geom_bar(position = "fill", stat = "identity")
TM_bar_taxa_excl
TM_bar_broad_excl <- ggplot(data = data$tm_algalonly_nomicro.csv, aes(x = site, y = percent, fill = broader)) +
  geom_bar(position = "fill", stat = "identity")
TM_bar_broad_excl
 
## (b) NMDS

# with Microcoleus included
tm_nmds <- metaMDS(as.matrix(data_wide$tm_algalonly.csv[,6:ncol(data_wide$tm_algalonly.csv)]),
                     distance = "bray",
                     trymax = 500,
                     autotransform = TRUE)
tm_nmds_data <- cbind(as.data.frame(scores(tm_nmds, "sites")), 
                   data_wide$tm_algalonly.csv %>% select(site_reach, site, field_date)) %>% 
  mutate(field_date = ymd(field_date),
         year = year(field_date),
         month = as.character(month(field_date)))

# NEED TO FINISH THIS
tm_nmds_plot <- ggplot(tm_nmds_data, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = site, shape = month), size = 4) +
  stat_ellipse(aes(color = site), type = "t", linetype = 2, size = 0.5) +
  scale_color_manual(values = c("SAL" = "#62a7f8",
                                "SFE-M" = "#416f16",
                                "SFE-SH" = "#a8ff82")) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "right") +
  labs(title = "NMDS of Morphologically-Identified Microcoleus Mat Community",
       subtitle = paste("Stress:", round(tm_nmds$stress, 3)),
       x = "NMDS Axis 1",
       y = "NMDS Axis 2")
tm_nmds_plot

## (c) misc. questions

# Excluding Microcoleus, what are the dominant taxa across all samples?
tm_summary <- data$tm_algalonly_nomicro.csv %>% 
  dplyr::group_by(site, taxa) %>% 
  dplyr::summarize(mean = mean(percent))

#### (4) TAC Samples ####

#### (5) NT Samples ####
