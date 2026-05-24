# to-do: compare species indicator analyses with SIMPER

# algal data (microscopy)

# libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot", "indicspecies"), require, character.only = T)

# read in files (data transformed in previous script, "4a_amongrivers_microscopy.R")
data <- lapply(list.files(path = "./data/morphological/transformed/", pattern = ".csv"),
               function(x) read.csv(paste("./data/morphological/transformed/", x, sep = "")))
names(data) <- c("nt", "tac", "tm")

# test across rivers Q1

start_col = 5

simper_river <- lapply(data, function(x) simper(x[,start_col:ncol(x)], x$site))
indic_river <- lapply(data, function(x) multipatt(x[,start_col:ncol(x)], x$site,
                                                  func = "r.g", control = how(nperm = 999)))

summary(simper$tm)
# SIMPER: epithemia and non-epithemia diatoms (doesn't identify groups)
summary(indic_river$tm)
# indicspecies: diatoms for SAL, e diatoms and anabaena and nostoc*

summary(simper$tac)
# SIMPER: microcoleus and scytonema for RUS vs. SFE-M, nothing for others
summary(indic_river$tac)
# also none

summary(simper$nt)
# spirogyra, cladophora, microcoleus, non diatoms, nostoc, homeothrix, leptolygbya, anabaena
# geitlerinema
summary(indic_river$nt)


# okay definitely prefer Indicator species but let's see what SIMPER says for 
# anatoxin groupings

env_sep <- read.csv("./data/field_and_lab/env_tox_standardized_byriver.csv")
data_sep <- lapply(data, function(x) left_join(x, env_sep, by = c("site", "site_reach", "field_date")))
sample_types <- c("tac", "nt", "tm")

# TM 22 was a "not sure if this is microcoleus sample" and there was not enough material
# to analyze for anatoxins
data_tog$tm <- data_tog$tm[-22,]
data_sep$tm <- data_sep$tm[-22,]

summary(simper(data_sep$tm[,start_col:ncol(data$tm)], data_sep$tm$TM_atx_category))
# medium vs. high anabaena & oscillatoria
# leptolyngbya & oscillatoria in low vs. high
summary(multipatt(data_sep$tm[,start_col:ncol(data$tm)], data_sep$tm$TM_atx_category,
        func = "r.g", control = how(nperm = 999)))
# for high oscillatoria and leptolyngbya
# likely different results when looking within a single river!