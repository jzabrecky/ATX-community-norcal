#### Comparing molecular 16s data among rivers
### Jordan Zabrecky
## last edited: 04.04.2026

# This code compares normalized 16s relative data from NT, TM, and TAC samples
# across rivers to answer Q1. First data is transformed (sqrt).
# Data is analyzed using NMDS and PERMANOVA. We also averaged across all samples
# from a river and created bar plots to visually compare average samples at each river

# TBD on keeping ISA results

#### (1) Loading libraries & data ####

# load libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot", "indicspecies"), require, character.only = T)

# we will use the NT data, and TM excluding Microcoleus, and TAC excluding Anabaena and Green Algae
nt <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_FINAL.csv") %>% 
  filter(sample_type == "NT")
tm <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_TM_nomicro.csv")
tac <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_TAC_noanacyl.csv")

# add into list
data_long <- list(nt, tm, tac)
names(data_long) <- c("nt", "tm", "tac")

# change format of field_data
data_long <- lapply(data_long, function(x) x <- x %>% 
                      mutate(field_date = mdy(field_date)))

# pivot wider using ASVs for community matrix
data <- lapply(data_long, function(x) {
  # subset columns we care about and pivot wider
  y = x %>% 
    select(site_reach, site, field_date, sample_type, triplicate, feature_ID, picrust2_relative_abundance) %>% 
    dplyr::rename(relative_abundance = picrust2_relative_abundance) %>% 
    pivot_wider(names_from = feature_ID, values_from = relative_abundance)
  # replace NA (indicating that ASV was not present in sample) with 0
  y[,6:ncol(y)][is.na(y[,6:ncol(y)])] = 0
  return(y)})

# set start col for community matrix
start_col <- 6  

# see our exploration with data transformation in another script, "S4a_testing_data_transformations.R"
# decided on square-root transformation on the relative abundances (Hellinger transformation)
for(i in 1:length(data)) {
  data[[i]][,start_col:ncol(data[[i]])] <- sqrt(data[[i]][,start_col:ncol(data[[i]])])
}

# save this data to use in future scripts (RUN ONCE)
#write.csv(data$nt, "./data/molecular/transformed/16s_nochimera_rarefied_95_NT_sqrttransformed.csv",
#          row.names = FALSE)
#write.csv(data$tm, "./data/molecular/transformed/16s_nochimera_rarefied_95_TM_nomicro_sqrttransformed.csv",
#          row.names = FALSE)
#write.csv(data$tac, "./data/molecular/transformed/16s_nochimera_rarefied_95_TAC_noanacyl_sqrttransformed.csv",
#          row.names = FALSE)

# lastly, add in broader categories to data_long so we don't have several with barplots
# could do this quickly with forcats package, but since we are lumping, may make sense
# to have more customized categories
data_long <- lapply(data_long, function(x) x %>% dplyr::rename(relative_abundance = picrust2_relative_abundance))
phylums <- lapply(data_long, function(x) x %>%
                    # calculate average abundance for each phylum
                    group_by(phylum) %>% 
                    dplyr::summarize(mean = mean(relative_abundance)) %>% 
                    # put into "Other" category if NA or % is less than #%
                    mutate(phylum_cat = case_when(is.na(phylum) ~ "Other",
                                                  mean < 0.03  ~ "Other",
                                                  TRUE ~ phylum)) %>% 
                    select(phylum, phylum_cat))
classes <- lapply(data_long, function(x) x %>%
                    # calculate average abundance for each class
                    group_by(class) %>% 
                    dplyr::summarize(mean = mean(relative_abundance)) %>% 
                    # put into "Other" category if NA or % is less than #%
                    mutate(classes_cat = case_when(is.na(class) ~ "Other",
                                                  mean < 0.05  ~ "Other",
                                                  TRUE ~ class)) %>% 
                    select(class, classes_cat))
cyano_order <- lapply(data_long, function(x) x %>%
                        # filter for cyanobacteria phylum
                        filter(phylum == "Cyanobacteria") %>% 
                        # calculate average abundance for each order
                        group_by(order) %>% 
                        dplyr::summarize(mean = mean(relative_abundance)) %>% 
                        # put into "Other" category if NA or % is less than #%
                        mutate(cyano_order_cat = case_when(is.na(order) ~ "Other",
                                                      mean < 0.035  ~ "Other",
                                                      TRUE ~ order)) %>% 
                        select(order, cyano_order_cat))
cyano_genus <- lapply(data_long, function(x) x %>%
                        # filter for cyanobacteria phylum
                        filter(phylum == "Cyanobacteria") %>% 
                        # calculate average abundance for each genus
                        group_by(genus) %>% 
                        dplyr::summarize(mean = mean(relative_abundance)) %>% 
                        # put into "Other" category if NA or % is less than #%
                        mutate(cyano_genus_cat = case_when(is.na(genus) ~ "Other",
                                                      mean < 0.11  ~ "Other",
                                                      TRUE ~ genus)) %>% 
                        select(genus, cyano_genus_cat))
# add these categories into our long dataframes
data_long_broader <- data_long
for(i in 1:length(data_long)) {
  data_long_broader[[i]] <- left_join(data_long_broader[[i]], phylums[[i]], by = "phylum")
  data_long_broader[[i]] <- left_join(data_long_broader[[i]], classes[[i]], by = "class")
  data_long_broader[[i]] <- left_join(data_long_broader[[i]], cyano_order[[i]], by = "order")
  data_long_broader[[i]] <- left_join(data_long_broader[[i]], cyano_genus[[i]], by = "genus")
}


#### (2) Functions for Analyses ####

# load from supplemental script
source("./code/supplemental_code/S4b_community_analyses_func.R")

# function to calculate diversity of ASVs
# @param data is dataframe of wide abundances
# @param start_col is column where abundance data starts
calc_diversity <- function(data, start_col) {
  vector = diversity(data[,start_col:ncol(data)])
}

#### (3) Relative Abundance Bar Plots ####

# put bar plots into lists
# (note if issues- something is masked; restart R)
barplot_phylum_plots <- lapply(data_long_broader, function(x) barplot(x, x = "site", y = "relative_abundance",
                                                              fill = "phylum_cat"))
barplot_class_plots <- lapply(data_long_broader, function(x) barplot(x, x = "site", y = "relative_abundance",
                                                              fill = "classes_cat"))
barplot_cyanoorder_plots <- lapply(data_long_broader, function(x) barplot(x %>% filter(phylum == "Cyanobacteria"),
                                                             x = "site", y = "relative_abundance",
                                                              fill = "cyano_order_cat"))
barplot_cyanogenus_plots <- lapply(data_long_broader, function(x) barplot(x %>% filter(phylum == "Cyanobacteria"),
                                                                  x = "site", y = "relative_abundance",
                                                                  fill = "cyano_genus_cat"))

# titles for plots
titles <- c("Non-Target Samples", 
            "Microcoleus Samples (excluding M)",
            "Anabaena/Cylindrospermum Samples (excluding AC & GA)")

# view plots (all bacteria categories)
for(i in 1:length(barplot_phylum_plots)) {
  print(barplot_phylum_plots[[i]] + labs(title = titles[i]))
  print(barplot_class_plots[[i]] + labs(title = titles[i]))
}
# could combine TAC categories more! will do so in figure script

# view plots (cyanobacteria categories)
for(i in 1:length(barplot_phylum_plots)) {
  print(barplot_cyanoorder_plots[[i]] + labs(title = titles[i]))
  print(barplot_cyanogenus_plots[[i]] + labs(title = titles[i]))
}

#### (4) Alpha Diversity Metrics ####

# calculate diversity for each dataframe
diversity <- lapply(data, function(x) {
  x = x %>% 
    mutate(shannon_diversity = calc_diversity(x, start_col)) %>% 
    select(field_date, site_reach, site, shannon_diversity)})

# plot diversity as boxplots
for(i in 1:length(diversity)) {
  boxplot = ggplot(data = diversity[[i]], aes(x = site, y = shannon_diversity, fill = site)) +
    geom_boxplot()
  print(boxplot)
}

# Does diversity differ across rivers?
set.seed(1)
lapply(diversity, function(x) kruskal.test(shannon_diversity~site, data = x))
# significantly different for TM (p = 0.02), close but not for NT (0.05), and not for TAC (p = 0.45)

# What is the mean of each?
means_medians <- lapply(diversity, function(x) x %>% 
                          dplyr::group_by(site) %>% 
                          dplyr::summarize(mean = mean(shannon_diversity),
                                           median = median(shannon_diversity)))
view(means_medians$nt)

# save diversity calculations (RUN ONCE)
#lapply(names(diversity), function(x) write.csv(diversity[[x]], 
#                                               paste("./data/molecular/shannon_diversity/", x, "_diversity.csv", sep = ""),
#                                               row.names = FALSE))

#### (5) NMDS Plots ####

# get NMDS for each dataframe (sqrt-transformed!)
set.seed(1) # set seed for reproducibility
NMDS_list <- lapply(data, function(x) getNMDSdata(x, start_col, ASV = TRUE))

# making plots
NMDS_plots <- lapply(NMDS_list, function(x) makeNMDSplot(x, FALSE, FALSE, 
                                                         color = "site", shape = "month"))

lapply(NMDS_plots, print)
# RESULT: TM and NT groups are visually distinct among rivers, but not for TAC

#### (6) Q: Are communities from each river significantly different? (PERMANOVA) ####

# run PERMANOVAs
set.seed(1)
permanovas <- lapply(data, function(x) runPERMANOVA(x, start_col, group = x$`site`))
lapply(permanovas, print)
# RESULTS: significant difference for all but not convinced about the TAC
# as visually they plotted on top of each other but had different dispersion

# check dispersion to see if that influences results
for(i in 1:length(data)) {
  set.seed(1)
  print(names(data)[i])
  print(anova(betadisper(vegdist(data[[i]][,start_col:ncol(data[[i]])], method = "bray"), 
                         data[[i]]$site)))
}
# dispersion not significantly different for NT, a little for TM (*) and very for
# TAC (***)

# Since, I am not convinced about TAC being significantly different visually, 
# let's compare centroid distance differences
set.seed(1)
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
# nt: 1.46816727428784
# tm: 1.58190801917234
# tac: 0.519095592038606
# TAC centroid distances are much closer on average, they are less than the distances
# for algal assemblages for TM and NT and but much more than TAC (which was 0.108)

## what if we remove the single salmon sample?
test <- data$tac %>% filter(site != "SAL")
set.seed(1)
print(anova(betadisper(vegdist(test[start_col:ncol(test)], method = "bray"), 
                       test$site)))
# still significantly different, but less so (**)
set.seed(1)
runPERMANOVA(test, start_col, group = test$`site`) # still very significant here (***)

#### (7) Q: What explains these differences? Species Indicator Analyses ####

# will likely not report on this because there is so much and its hard to parse
# what is important but here it is

## (a) phylums

phylums <- lapply(data_long, function(x) {
    # subset columns we care about and pivot wider
      y = x %>% 
        select(site_reach, site, field_date, sample_type, phylum, relative_abundance) %>% 
        dplyr::group_by(site_reach, site, field_date, sample_type, phylum) %>% 
        # remove multiples of phylums due to differing ASVs
        dplyr::summarize(relative_abundance = sum(relative_abundance)) %>% 
        pivot_wider(names_from = phylum, values_from = relative_abundance)
      # if there is a column NA, remove it
      if("NA" %in% colnames(y)) {
        y = y %>% 
          select(!c("NA"))
      }
      # replace NA (indicating that ASV was not present in sample) with 0
      y[,5:ncol(y)][is.na(y[,5:ncol(y)])] = 0
      # lastly, square-root transform the data
      y[,5:ncol(y)] <- sqrt(y[,5:ncol(y)])
      return(y)})

# (i) NT
summary(multipatt(phylums$nt[,5:ncol(phylums$nt)], phylums$nt$site, func = "r.g", control = how(nperm = 999)))
# lots of unique identified with RUS including ***: Actinobacteriota, Elusimicrobiota,
# Desulfobacteria, Planctomycetota, Verrucomicrobiota
# SAL: Deincoccota (***), Armatimonadota (**)
# RUS + SAL: Proteobacteria *
# SAL + SFE-M: Cyanobacteria **

# (ii) TM
summary(multipatt(phylums$tm[,5:ncol(phylums$tm)], phylums$tm$site, func = "r.g", control = how(nperm = 999)))
# only identified for SFE-M including **: Desulfobacteria, Cyanobacteria, Verrucomicrobiota,
# Chloroflexiota

# (iii) TAC
summary(multipatt(phylums$tac[,5:ncol(phylums$tac)], phylums$tac$site, func = "r.g", control = how(nperm = 999)))
# only identified for SAL: Fibrobacterota *, Spirochaetota * 

## (b) classes

classes <- lapply(data_long, function(x) {
  # subset columns we care about and pivot wider
  y = x %>% 
    select(site_reach, site, field_date, sample_type, class, relative_abundance) %>% 
    dplyr::group_by(site_reach, site, field_date, sample_type, class) %>% 
    # remove multiples of phylums due to differing ASVs
    dplyr::summarize(relative_abundance = sum(relative_abundance)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = class, values_from = relative_abundance)
  # if there is a column NA, remove it
  if("NA" %in% colnames(y)) {
    y = y %>% 
      select(!c("NA"))
  }
  # replace NA (indicating that ASV was not present in sample) with 0
  y[,5:ncol(y)][is.na(y[,5:ncol(y)])] = 0
  # lastly, square-root transform the data
  y[,5:ncol(y)] <- sqrt(y[,5:ncol(y)])
  return(y)})

# (i) NT
summary(multipatt(classes$nt[,5:ncol(classes$nt)], classes$nt$site, func = "r.g", control = how(nperm = 999)))
# lots for RUS (not listing)
# SAL: Fimbriimonadia ***, Deinococci ***
# SFE-M: Rhodothermia *
# RUS + SAL: Gammaproteobacteria **, Acidimicrobiia *
# RUS + SFE-M: including Desulfuromonadia **, Desulfovibrionia **
# SAL + SFE-M: Cyanobacteria ***, Gracilibacteria *

# (ii) TM
summary(multipatt(classes$tm[,5:ncol(classes$tm)], classes$tm$site, func = "r.g", control = how(nperm = 999)))
# SAL: Bacteroidia *
# SFE-M: many but no with ***

# (iii) TAC
summary(multipatt(classes$tac[,5:ncol(classes$tac)], classes$tac$site, func = "r.g", control = how(nperm = 999)))
# only for SAL including Vampirivibrionia ***
 
## (c) cyanobacteria orders (also will likely not discuss in paper)

cyano_orders <- lapply(data_long, function(x) {
  # subset columns we care about and pivot wider
  y = x %>% 
    filter(phylum == "Cyanobacteria") %>% 
    select(site_reach, site, field_date, sample_type, order, relative_abundance) %>% 
    dplyr::group_by(site_reach, site, field_date, sample_type, order) %>% 
    # remove multiples of phylums due to differing ASVs
    dplyr::summarize(relative_abundance = sum(relative_abundance)) %>% 
    pivot_wider(names_from = order, values_from = relative_abundance)
  # if there is a column NA, remove it
  if("NA" %in% colnames(y)) {
    y = y %>% 
      select(!c("NA"))
  }
  # replace NA (indicating that ASV was not present in sample) with 0
  y[,5:ncol(y)][is.na(y[,5:ncol(y)])] = 0
  # lastly, square-root transform the data
  y[,5:ncol(y)] <- sqrt(y[,5:ncol(y)])
  return(y)})

# (i) NT
summary(multipatt(cyano_orders$nt[,5:ncol(cyano_orders$nt)], cyano_orders$nt$site, func = "r.g", control = how(nperm = 999)))
# RUS: Obscuribacterales **
# SAL: all *** Gloeobacterales, Gomontiellales, Chroococcidiopsidales, Leptolynbyales
# SFE-M: Nodosineales ***, Choococcidiopsidaceae ***, Psuedanabaenales *
# RUS + SAL: Pleurocapsales **
# RUS + SFE-M: Chroococcales **, Synechococcales **
# SAL + SFE-M: Nostocales ***

# (ii) TM
summary(multipatt(cyano_orders$tm[,5:ncol(cyano_orders$tm)], cyano_orders$tm$site, func = "r.g", control = how(nperm = 999)))
# only identified for SFE-M: Pseudanabaenales **, Chroococcales ***, Synechococcales ***,
# Sericytochromatia **, Endosymbiotic Diazoplast **, Leptolynbyaceae *

# (iii) TAC
summary(multipatt(cyano_orders$tac[,5:ncol(cyano_orders$tac)], cyano_orders$tac$site, func = "r.g", control = how(nperm = 999)))
# only identified for SAL: Gastranaerophilales **, Vampirovibrionales **

## (d) Cyanobacteria Genera (Likely will not report due to high resolution of genus!)

cyano_genera <- lapply(data_long, function(x) {
  # subset columns we care about and pivot wider
  y = x %>% 
    filter(phylum == "Cyanobacteria") %>% 
    select(site_reach, site, field_date, sample_type, genus, relative_abundance) %>% 
    dplyr::group_by(site_reach, site, field_date, sample_type, genus) %>% 
    # remove multiples of genera due to differing ASVs
    dplyr::summarize(relative_abundance = sum(relative_abundance)) %>% 
    pivot_wider(names_from = genus, values_from = relative_abundance)
  # if there is a column NA, remove it
  if("NA" %in% colnames(y)) {
    y = y %>% 
      select(!c("NA"))
  }
  # replace NA (indicating that ASV was not present in sample) with 0
  y[,5:ncol(y)][is.na(y[,5:ncol(y)])] = 0
  # lastly, square-root transform the data
  y[,5:ncol(y)] <- sqrt(y[,5:ncol(y)])
  return(y)})

# (i) NT
summary(multipatt(cyano_genera$nt[,5:ncol(cyano_genera$nt)], cyano_genera$nt$site, func = "r.g", control = how(nperm = 999)))
# RUS: Chroococcopsis **, Vampirovibrio *
# SAL: all *** Calothrix, Gloeobacter, Cyanothece, Aliterella, Oscillatoria, Phormidium, 
# Chamaesiphon, and ** Pleurocapsa
# SFE-M: Synechoccus ***, Snowella ***, Nodosilinea ***, Arthronema ***, Tolypothrix **,
# Synechocystis **, Nodularia **, Cuspidothrix **, Limnothrix **, Sphaerospermopsis **,
# Merismopedia *, Gloeotrichia *, Aphanizomenon *
# RUS + SFE-M: Microcystis **, Cyanobium **, Geminocystis *
# SAL + SFE-M: Schizothrix **

# (ii) TM
summary(multipatt(cyano_genera$tm[,5:ncol(cyano_genera$tm)], cyano_genera$tm$site, func = "r.g", control = how(nperm = 999)))
# only identified for SFE-M: Psuedanabaena **, Cyanobium ***, Geminocystis *, Sericytochromatia **,
# Microcystis **, Limnothrix *, Endosymbiotic Diazoplast **, Leptolyngbya *, Nostoc *

# (iii) TAC
summary(multipatt(cyano_genera$tac[,5:ncol(cyano_genera$tac)], cyano_genera$tac$site, func = "r.g", control = how(nperm = 999)))
# only identified for SAL: Aphanizomenon *, Gastranerophilales **, Kamptonema *, Arthronema *,
# Vampirovibrio **

#### (8) Misc. Questions ####

## How present are other anatoxin associated taxa?
# using ATX taxa as identified in Christensen & Khan (2019) and Wood et al. (2020)
# using list from Christensen & Khan et al. (2019): Anabaena, Aphanizomenon,
# Aphanothece, Arthospira, Cylindrospermopsis, Cylindrospermum, Gomphosphaeria,
# Limnothrix, Lyngbya, Microcystis, Nostoc, Oscillatoria, Phormidium/Microcoleus,
# Planktothrix, Planktolyngbia, Synechocystis, Psuedoanabaena,
# Raphidopsis, Tychonema
# noting that it doesn't have Geilerinema, so we should also include list
# of ATX producers from Wood et al. (2020) which adds: Fisherella, 
# Geitlerinema, Leptolyngbya, Microseira (formerly Lyngbya)
atx_taxa_only <- lapply(data_long, function(x) {
  
  # make dataframe with only taxa in list above 
  # (only writing what taxa we recorded from that list)
  df = x %>% 
    filter(genus %in% c("Anabaena", "Aphanizomenon", "Cylindrospermum", "Cylindrospermopsis",
                        "Limnothrix", "Leptolyngbya", "Microcystis", "Nostoc", "Oscillatoria",
                        "Phormidium", "Microcoleus", "Planktothrix", "Synechocystis", 
                        "Geitlerinema", "Pseudanabaena"
    )) %>% 
    # searching the TAC and TM dataframes do not yield: aphanothece, arthospira, gomphosphaeria,
    # lyngbya, planktolyngbya, microseira
    mutate(sample_name = paste("(", month(field_date), "-", day(field_date), ") ", site_reach, sep = "")) %>% 
    select(sample_name, field_date, sample_type, site, site_reach, relative_abundance, order, genus) %>% 
    mutate(order_genus = paste(order, " - ", genus, sep = ""))
  
  # make bar plot (show each sample individually)
  plot <- ggplot(data = df, aes(x = sample_name, y = relative_abundance / 100, fill = genus)) +
    geom_bar(stat = "identity") +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(x = NULL, y = "Relative Abundance") +
    facet_wrap(~site, scales = "free_x")
  print(plot) # view plot
  
  # return a list including dataframe, then plot
  return(list(df, plot))
})

# some discrepancies here between and microscopy data
# could be limits of 16s with genus level or the poor resolution of database
# for more, see Dvorak et al. (2025)