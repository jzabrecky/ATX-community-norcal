#### Comparing microscopy data across time
### Jordan Zabrecky
## last edited: 12.15.2025

# This code compares microscopy data from NT, TM, and TAC samples
# across time to answer Q2

#### (1) Loading libraries & data ####

# set seed for reproducibility
set.seed(2025)

# libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot"), require, character.only = T)

# read in files (data transformed in previous script, "4a_amongrivers_microscopy.R")
data <- lapply(list.files(path = "./data/morphological/transformed/", pattern = ".csv"),
               function(x) read.csv(paste("./data/morphological/transformed/", x, sep = "")))
names(data) <- c("nt", "tac", "tm")

# also load un-transformed relative abundances and make it longer for bar plots through time
nt <- read.csv("./data/morphological/nt_algalonly.csv")
tm <- read.csv("./data/morphological/tm_algalonly_nomicro.csv")
tac <- read.csv("./data/morphological/tac_algalonly_noanacylgreenalgae.csv")

# add into list
unaltered_data <- list(nt, tm, tac)
names(unaltered_data) <- c("nt", "tm", "tac")

# finally, pivot longer for those bar plots and select only 2022 data
data_longer <- lapply(unaltered_data, 
                      function(x) x %>% pivot_longer(cols = c(5:ncol(x)), values_to = "percent",
                                               names_to = "taxa") %>% 
                        filter(year(ymd(field_date)) == 2022))

#### (2) Add  Columns for Sampling Event & Broader Taxa ####

# we have field dates for sampling distinguishing samples, however, let's change it
# to be the the number of sampling event (where x is data)
add_event_no <- function(data) {
  data %>% 
    mutate(field_date = ymd(field_date),
           month = month(field_date)) %>% 
    mutate(event_no = case_when((field_date >= ymd("2022-06-24") & field_date <= ymd("2022-06-29")) ~ 1,
                                (field_date >= ymd("2022-07-06") & field_date <= ymd("2022-07-14")) ~ 2,
                                (field_date >= ymd("2022-07-20") & field_date <= ymd("2022-07-28")) ~ 3,
                                (field_date >= ymd("2022-08-02") & field_date <= ymd("2022-08-10")) ~ 4,
                                (field_date >= ymd("2022-08-17") & field_date <= ymd("2022-08-23")) ~ 5,
                                (field_date >= ymd("2022-09-01") & field_date <= ymd("2022-09-06")) ~ 6,
                                (field_date >= ymd("2022-09-15") & field_date <= ymd("2022-09-22")) ~ 7)) %>% 
    relocate(event_no, .before = "field_date") %>% 
    relocate(month, .before = "field_date")
}

# apply to all dataframes
data <- lapply(data, add_event_no)
data_longer <- lapply(data_longer, add_event_no)

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

##### (3) Function for Analyses ####

# load from supplemental script
source("./code/supplemental_code/S4a_community_analyses_func.R")
source("./code/supplemental_code/S4c_barplot_func.R")

#### (4) Barplots through Time ####

barplot_taxa_plots <- lapply(data_longer, function(x) barplot(x, x = "event_no", y = "percent", 
                                                              fill = "taxa", facet_wrap = "site"))
barplot_broader_plots <- lapply(data_longer, function(x) barplot(x, x = "event_no", y = "percent", 
                                                              fill = "broader", facet_wrap = "site"))

# titles for plots
titles <- c("Non-Target Samples", 
            "Microcoleus Samples (excluding M)",
            "Anabaena/Cylindrospermum Samples (excluding AC & GA)")

# view plots
for(i in 1:length(barplot_taxa_plots)) {
  print(barplot_taxa_plots[[i]] + labs(title = titles[i]))
  print(barplot_broader_plots[[i]] + labs(title = titles[i]))
}

#### (5) Q: What are the dominant five taxa at each event for each river? ####

temporal_overview <- lapply(data_longer, function(x) {
  summary = x %>% 
    group_by(site, event_no, broader, taxa, sample_type) %>% 
    dplyr::summarize(mean_abun = mean(percent)) %>% 
    # pivot_wider(names_from = event_no, values_from = mean_abun) %>% 
    ungroup() %>% 
    dplyr::group_by(site, event_no) %>% 
    arrange(desc(mean_abun)) %>% 
    slice_head(n = 5) %>% 
    ungroup()
  
  #summary = summary %>% # get table
  #  mutate(ranking = rep(seq(1:5), nrow(summary) / 5)) %>% 
  #  select(!c(mean_abun, broader)) %>% 
  #  pivot_wider(names_from = event_no, values_from = taxa)
  #return(summary)
  
  # split to have separate colors for each plot
  summaries = split(summary, summary$site)
  for(i in 1:length(summaries)) {
    plot = ggplot(summaries[[i]], aes(x = event_no, y = mean_abun, color = taxa)) +
              geom_line() +
              geom_point() +
              ggtitle(paste(summary$sample_type[1], ": ", names(summaries)[i], sep = ""))
    print(plot)
  }
})
# RUS NT: spirogyra largely dominant after event 1, with decreasing amounts of non er diatoms, 
# epithemia in top 5 after event 2, cladophora, oedogonium, stigeoclonium, also present
# ana/cyl in top five last two events
# SAL NT: dominance of diatoms that decreases, increase in microcoleus 2 to 3, 
# leptolyngyba and homoeothrix also present in early samples with some cladophora
# other coccoids only in top 5 first event, and geitlernema for event 2, green algaes for event 5
# SFE-M NT: decrease in nostoc from beginning, and microcoleus, increase in spirogyra,
# anabaena only present in top 5 2 & 3 and microcoleus 1 & 2, presence of epithemia at 3 that decreases

# SAL TM: mostly diatoms, consistent across the two dates, also present are geitlerinema, 
# other coccoids, and green algae
# SFE-M TM: increase in green algae (likely due to sampling), increase in diatoms at 3 & 4 then 
# a decrease, coccoids present in later samples (3 to 7) and anabaena in top 5 for 2-4,
# geitlerinema consistently present (sans 2) diatoms fluctuate throughout

# RUS TAC: mostly diatoms, decrease in epithemia after 5, presence of phormidium 6+7
# oscillatoria in top 5 for 5, microcoleus only for top 5 at 3, geitlerinema was a consistent presence
# as well as other coccoid cyanobacteria
# SFE-M TAC: microcoleus in top 5 from 2-5, other coccoids and geitlerinema present steady throughout,
# e diatoms peak at 3 and inconsistent decrease, oscillatoria present at 7
# SAL TAC: this sample also has geitlerinema

#### (6) PERMANOVA ####

# set column where abundance data starts
start_col <- 7

# separate data out by site
data_river <- lapply(data, function(x) split(x, x$`site`))

## (a) with the strata argument for EVENT NO.
permanovas_event <- lapply(data, function(x) runPERMANOVA(x, start_col, x$`event_no`, strata = x$'site'))
lapply(permanovas_event, print)
# NT***, TAC*** are significant among sampling dates, but not TM
# TM must have been overshadowed by the two dates of Salmon samples (see below)

## (b) separated out by river

permanovas_event_NT_sep <- lapply(data_river$nt, function(x) runPERMANOVA(x, start_col, x$`event_no`))
lapply(permanovas_event_NT_sep, print)
# NT significantly different for SAL**, RUS**, and SFE-M***

permanovas_event_TM_sep <- lapply(data_river$tm, function(x) runPERMANOVA(x, start_col, x$`event_no`))
lapply(permanovas_event_TM_sep, print)
# TM significant for SFE-M** but not SAL but only 2 days
# this is likely influencing the results of strata above, thus we should use this result

permanovas_event_TAC_sep <- lapply(data_river$tac, function(x) runPERMANOVA(x, start_col, x$`event_no`))
lapply(permanovas_event_TAC_sep, print)
# TAC significant for SFE-M* and especially RUS*** (only one data point for SAL, so NA)

## (c) check dispersion in groups for event no.

for(i in 1:length(data_river$nt)) {
  print(names(data_river$nt)[i])
  print(anova(betadisper(vegdist(data_river$nt[[i]][,start_col:ncol(data_river$nt[[i]])], method = "bray"), 
                         data_river$nt[[i]]$event_no)))
}
# dispersion not significant for any NT (as expected!)

for(i in 1:length(data_river$tm)) {
  print(names(data_river$tm)[i])
  print(anova(betadisper(vegdist(data_river$tm[[i]][,start_col:ncol(data_river$tm[[i]])], method = "bray"), 
                         data_river$tm[[i]]$event_no)))
}
# also not for TM

# skip the salmon sample which is index 2
for(i in c(1,3)) {
  print(names(data_river$tac)[i])
  print(anova(betadisper(vegdist(data_river$tac[[i]][,start_col:ncol(data_river$tac[[i]])], method = "bray"), 
                         data_river$tac[[i]]$event_no)))
}
# also not for TAC!

#### (7) Species Indicator Analysis ####

## (a) NT
for(i in 1:length(data_river$nt)) {
  print(names(data_river$nt)[i])
  print(summary(multipatt(data_river$nt[[i]][,start_col:ncol(data_river$nt[[i]])],
                          data_river$nt[[i]]$event_no, func = "r.g", control = how(nperm = 999))))
}
# RUS: 2+3+4+5+6+7 geitlerinema
# SAL: 1 other coccoids, 2 geitlerinema, 7 green algae (spirogyra, mougeotia, epithemia, oedogonium)
# SFE-M: 1+2 nostoc, 2+3+4 anabaena & cyl, 3+4+5+6+7 epithemia, 2+5+6+7 spirogyra
# 3&7 geitlerinema

## (b) TM
for(i in 1:length(data_river$tm)) {
  print(names(data_river$tm)[i])
  print(summary(multipatt(data_river$tm[[i]][,start_col:ncol(data_river$tm[[i]])],
                          data_river$tm[[i]]$event_no, func = "r.g", control = how(nperm = 999))))
}
# no unique for SAL, but only 2 events
# SFE-M: 1+2 nostoc, 2+3+4 anabaena

## (c) TAC
# skipping SAL
for(i in c(1,3)) {
  print(names(data_river$tac)[i])
  print(summary(multipatt(data_river$tac[[i]][,start_col:ncol(data_river$tac[[i]])],
                          data_river$tac[[i]]$event_no, func = "r.g", control = how(nperm = 999))))
}
# RUS: 6 microcoleus, 6+7 phormidium unknown, 5+6+7 oscillatoria, 3+4+5+6 epithemia
# SFE-M:  2+5 nodularia, e_diatoms 2+3+5

#### STILL CURIOUS- which is changing the most??? ####
# would this just be distances of centroids from start to finish?
# read the Avolio paper (2019)
# would likely need all sample types we are comparing on the same NMDS
# e.g., does NT change more than TM, etc.

#### (8) Rate of community change ####

#### (8) Conclusions ####

# will write when revisiting

# GEITLERINEMA IS ALWAYS PRESENT.
# curious how consistently present it is in the NT samples