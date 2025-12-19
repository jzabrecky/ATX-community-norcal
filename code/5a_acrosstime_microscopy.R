#### Comparing microscopy data across time
### Jordan Zabrecky
## last edited: 12.15.2025

## This code compares microscopy data from NT, TM, and TAC samples
## across rivers to answer Q2

## to-do:
## ANOSIM w/ time stamp
## indicator species analsis? vs. SIMPER?
## look into how other papers are determining which is the species most responsible

## CONSIDER WHAT TO DO WITH 2023 DATA

# MAYBE MOVE NMDS FUNCTIONS TO SUPPLEMENTAL CODE FILES

# ALSO ADD IN MONTH TO THIS

#### (1) Loading libraries & data ####

# set seed for reproducibility
set.seed(2025)

# libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot"), require, character.only = T)

# read in files (data transformed in previous script, "4a_amongrivers_microscopy.R")
data <- lapply(list.files(path = "./data/morphological/transformed/", pattern = ".csv"),
               function(x) read.csv(paste("./data/morphological/transformed/", x, sep = "")))
names(data) <- c("nt", "tac", "tm")

# load un-transformed relative abundances and make it longer for bar plots through time
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

#### (3) PERMANOVA ####

# set column where abundance data starts
start_col <- 7

# separate data out by site
data_river <- lapply(data, function(x) split(x, x$`site`))

## (a) with the strata argument for EVENT NO.
permanovas_event <- lapply(data, function(x) runPERMANOVA(x, start_col, x$`event_no`, strata = x$'site'))
lapply(permanovas_event, print)
# NT, TAC are significant among sampling dates, but not TAC

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

## (c) check dispersion in groups
# LEFT OFF HERE

#### (4) CLUSTER ANALYSES!?!??!??!?!???!????!?!?

# maybe heat map stuff WOULD be cool but like idk right now :)

#### (3) test indicator species analysis

# Probably want to add this to script 5a WILL GO BACK AND DO THIS NEXT
# maybe try out both???
# can you use indicator species analyses on relative abundances?

# more curious how this will change over time

# does this need to be hellinger transformed and why?
library(indicspecies)

# work with NT data to start
NT <- data_wide$nt_algalonly_sqrttransformed.csv

# func = "r.g" is for relative abundances 
groups <- NT$site
test <- multipatt(NT[,7:ncol(NT)], groups, func = "r.g", control = how(nperm = 999))
summary(test)
