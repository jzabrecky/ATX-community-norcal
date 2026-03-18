#### Comparing microscopy data with regard to anatoxin concentrations
### Jordan Zabrecky
## last edited: 03.12.2026

# This script examines how communities as identified by microscopy
# change with increasing anatoxin concentrations with PERMANOVA, NMDS,
# and ISA

#### (1) Loading libraries & data ####

# libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot", 
         "indicspecies"), require, character.only = T)

# read in relative abundance files (data transformed in previous script, "4a_amongrivers_microscopy.R")
data <- lapply(list.files(path = "./data/morphological/transformed/", pattern = ".csv"),
               function(x) read.csv(paste("./data/morphological/transformed/", x, sep = "")))
sample_types <- c("nt", "tac", "tm")
names(data) <- sample_types

# read in environmental covariates & toxin data
atx <- read.csv("./data/field_and_lab/atx_w_categorical_groupings.csv")
# as a reminder, NAs mean no sample was taken!!! so this is different 
# from 0's or non-detects when a sample was present, but no anatoxins were detected

# join data and anatoxin data
data_tog <- lapply(data, function(x) left_join(x, atx, by = c("field_date", "site_reach", "site")))

#### (2) How do samples change with varying anatoxin concentrations among all rivers? ####

# load community analyses functions from other script
source("./code/supplemental_code/S4b_community_analyses_func.R")

# set start column for community data
start_col <- 5

## (a) PERMANOVAs

## (i) log-anatoxins
set.seed(1)
lapply(sample_types, function(x) {
  # set column name for anatoxin based on sample type
  atx_col = paste("log_", str_to_upper(x), "_ATX_all_ug_orgmat_g", sep = "")
  
  # PERMANOVA test
  print(paste(x, "PERMANOVA"))
  print(runPERMANOVA(data_tog[[x]], start_col, end_col = ncol(data[[x]]), 
               group = as.vector(data_tog[[x]][atx_col])[[1]], na.action = "na.omit"))
  
  # BETADISPER test
  print(paste(x, "BETADISPER"))
  print(anova(betadisper(vegdist(data_tog[[x]][,start_col:ncol(data[[x]])], method = "bray"), 
                   as.vector(data_tog[[x]][atx_col])[[1]])))
  
  return()
})
# NT, *** permanova, * betadisper
# TAC: * permanova, 0.484 betadisper
# TM: ** permanova, ** betadisper

# (ii) anatoxin groupings
set.seed(1)
lapply(sample_types, function(x) {
  # set column name for anatoxin based on sample type
  atx_col = paste(str_to_upper(x), "_atx_category", sep = "")
  
  # PERMANOVA test
  print(paste(x, "PERMANOVA"))
  print(runPERMANOVA(data_tog[[x]], start_col, end_col = ncol(data[[x]]), 
                     group = as.vector(data_tog[[x]][atx_col])[[1]], na.action = "na.omit"))
  
  # BETADISPER test
  print(paste(x, "BETADISPER"))
  print(anova(betadisper(vegdist(data_tog[[x]][,start_col:ncol(data[[x]])], method = "bray"), 
                         as.vector(data_tog[[x]][atx_col])[[1]])))
  
  return()
})
# NT: *, *
# TAC: *, not
# TM: *, ***

# (iii) anatoxin binary
set.seed(1)
lapply(sample_types, function(x) {
  # set column name for anatoxin based on sample type
  atx_col = paste(str_to_upper(x), "_atx_detected", sep = "")
  
  # PERMANOVA test
  print(paste(x, "PERMANOVA"))
  print(runPERMANOVA(data_tog[[x]], start_col, end_col = ncol(data[[x]]), 
                     group = as.vector(data_tog[[x]][atx_col])[[1]], na.action = "na.omit"))
  
  # BETADISPER test
  print(paste(x, "BETADISPER"))
  print(anova(betadisper(vegdist(data_tog[[x]][,start_col:ncol(data[[x]])], method = "bray"), 
                         as.vector(data_tog[[x]][atx_col])[[1]])))
  
  return()
})
# NT: **, *
# TAC: not, not
# TM: *, ***

## (b) ISA

# (i) atx categories
set.seed(1)
lapply(sample_types, function(x) {
  # set column name for anatoxin based on sample type
  atx_col = paste(str_to_upper(x), "_atx_category", sep = "")
  
  # remove columns where atx_col is na
  if(x == "nt" | x == "tm") { 
    # tm also included here to remove sample from 8-23-2022 as we did not have enough to analyze
    temp = data_tog[[x]][-which(is.na(data_tog[[x]][atx_col])),]
  } else {
    temp = data_tog[[x]]
  }
  
  # run ISA
  print(paste(x, "ISA"))
  summary(multipatt(temp[,start_col:ncol(data[[x]])], as.vector(temp[atx_col])[[1]]), 
                               func = "r.g", control = how(nperm = 999))
})
# NT: assortment of unicellular green algae and rivularia, rhopalodia & anabaena for high & low
# TAC: microcoleus associated with high and low samples
# TM: geitlerinema associated with high and low samples

# (ii) atx detected
set.seed(1)
lapply(sample_types, function(x) {
  # set column name for anatoxin based on sample type
  atx_col = paste(str_to_upper(x), "_atx_detected", sep = "")
  
  # remove columns where atx_col is na
  if(x == "nt" | x == "tm") { 
    # tm also included here to remove sample from 8-23-2022 as we did not have enough to analyze
    temp = data_tog[[x]][-which(is.na(data_tog[[x]][atx_col])),]
  } else {
    temp = data_tog[[x]]
  }
  
  # run ISA
  print(paste(x, "ISA"))
  summary(multipatt(temp[,start_col:ncol(data[[x]])], as.vector(temp[atx_col])[[1]]), 
          func = "r.g", control = how(nperm = 999))
})
# weirder results here:
# NT: detected- cosmarium
# TAC: detected- microcoleus
# TM: none


## (c) NMDS

# get NMDS data
NMDS_data <- lapply(sample_types, function(x) getNMDSdata(data_tog[[x]], start_col, 
                    end_col = ncol(data[[x]])))
names(NMDS_data) <- sample_types

# (i) with atx category
lapply(sample_types, function(x) {
  
  # set column name for anatoxin based on sample type
  atx_col = paste(str_to_upper(x), "_atx_category", sep = "")
  
  makeNMDSplot(NMDS_data[[x]], FALSE, FALSE, color = atx_col, shape = "site")
})
# visually distinct for TM (note only one high sample for TM with current designation)
# NT obviously a mess with all rivers in one graph

# (ii) with atx detected
lapply(sample_types, function(x) {
  
  # set column name for anatoxin based on sample type
  atx_col = paste(str_to_upper(x), "_atx_detected", sep = "")
  
  makeNMDSplot(NMDS_data[[x]], FALSE, FALSE, color = atx_col, shape = "site")
})
# looks more visually distinct for TM

#### (3) How do samples change with varying anatoxin concentrations within a river? ####

# split data based on river!
data_river <- lapply(data_tog, function(x) split(x, x$site))

## (a) PERMANOVA

## (i) log-anatoxins
set.seed(1)
lapply(sample_types, function(x) {
  
  lapply(names(data_river[[x]]), function(y) {
    if(y == "SAL" & x == "tac") {
      print("no test")
    } else {
      # set column name for anatoxin based on sample type
      atx_col = paste("log_", str_to_upper(x), "_ATX_all_ug_orgmat_g", sep = "")
      
      # PERMANOVA test
      print(paste(y, x, "PERMANOVA"))
      print(runPERMANOVA(data_river[[x]][[y]], start_col, end_col = ncol(data[[x]]), 
                         group = as.vector(data_river[[x]][[y]][atx_col])[[1]], na.action = "na.omit"))
      
      # BETADISPER test
      print(paste(y, x, "BETADISPER"))
      print(anova(betadisper(vegdist(data_river[[x]][[y]][,start_col:ncol(data[[x]])], method = "bray"), 
                             as.vector(data_river[[x]][[y]][atx_col])[[1]])))
    
      return()
    }
  })
})
# RUS:
# nt not, *
# tac *, not
# SAL:
# nt not, not
# tm not, not
# tac no test
# SFE:
# nt *, not
# tm **, **
# tac: ***, ***

## (ii) atx category
set.seed(1)
lapply(sample_types, function(x) {
  
  lapply(names(data_river[[x]]), function(y) {
    if(y == "SAL" & x == "tac") {
      print("no test")
    } else {
      
      # set column name for anatoxin based on sample type
      atx_col = paste(str_to_upper(x), "_atx_category", sep = "")
      
      # PERMANOVA test
      print(paste(y, x, "PERMANOVA"))
      print(runPERMANOVA(data_river[[x]][[y]], start_col, end_col = ncol(data[[x]]), 
                         group = as.vector(data_river[[x]][[y]][atx_col])[[1]], na.action = "na.omit"))
      
      # BETADISPER test
      print(paste(y, x, "BETADISPER"))
      print(anova(betadisper(vegdist(data_river[[x]][[y]][,start_col:ncol(data[[x]])], method = "bray"), 
                             as.vector(data_river[[x]][[y]][atx_col])[[1]])))
      
      return()
    }
  })
})
# RUS:
# nt not, not
# tac not, not
# SAL:
# nt no test (for whatever reason was not working) -- SHOULD BE FIXED NOW
# will first double-check grouping with Joanna
# tm not, not
# tac no test
# SFE:
# nt not, not
# tm **, *
# tac: **, not

## (iii) atx detected
set.seed(1)
lapply(sample_types, function(x) {
  
  lapply(names(data_river[[x]]), function(y) {
    if(y == "SAL" & x == "tac") {
      print("no test")
    } else {
      
      # set column name for anatoxin based on sample type
      atx_col = paste(str_to_upper(x), "_atx_detected", sep = "")
      
      # PERMANOVA test
      print(paste(y, x, "PERMANOVA"))
      print(runPERMANOVA(data_river[[x]][[y]], start_col, end_col = ncol(data[[x]]), 
                         group = as.vector(data_river[[x]][[y]][atx_col])[[1]], na.action = "na.omit"))
      
      # BETADISPER test
      print(paste(y, x, "BETADISPER"))
      print(anova(betadisper(vegdist(data_river[[x]][[y]][,start_col:ncol(data[[x]])], method = "bray"), 
                             as.vector(data_river[[x]][[y]][atx_col])[[1]])))
      
      return()
    }
  })
})
#RUS:
# nt not, not
# tac not, not
# SAL:
# nt not not
# tm not, not
# tac *, not
# SFE:
# nt **, not
# tm **, *
# tac: ***, not

## (b) ISA

## (i) anatoxin groupings
set.seed(1)
lapply(sample_types, function(x) {
  lapply(names(data_river[[x]]), function(y) {
    if(y == "SAL" & x == "tac") {
      print("no test")
    } else {
      # set column name for anatoxin based on sample type
      atx_col = paste(str_to_upper(x), "_atx_category", sep = "")
      
      # remove columns where atx_col is na
      if(x == "nt" | (x == "tm" & y == "SFE-M")) { 
        # tm also included here to remove sample from 8-23-2022 as we did not have enough to analyze
        temp = data_river[[x]][[y]][-which(is.na(data_river[[x]][[y]][atx_col])),]
      } else {
        temp = data_river[[x]][[y]]
      }
      
      # run ISA
      print(paste(x, y, "ISA"))
      summary(multipatt(temp[,start_col:ncol(data[[x]])], as.vector(temp[atx_col])[[1]]), 
              func = "r.g", control = how(nperm = 999))
    }
  })
})
# RUS:
# nt (low) phormidium & oscillatoria
# tac (low) phormidium
# SAL:
# nt nothing
# tm not, not
# tac no test
# SFE:
# nt (high & low) rophalodia, anabaena_and_cylindrospermum, leptolyngbya
# tm (high & none) nostoc, (high) leptolyngbya
# tac: (high & low) microcoleus, (none) nostoc

## (ii) anatoxin binary
set.seed(1)
lapply(sample_types, function(x) {
  lapply(names(data_river[[x]]), function(y) {
    if(y == "SAL" & x == "tac") {
      print("no test")
    } else {
      # set column name for anatoxin based on sample type
      atx_col = paste(str_to_upper(x), "_atx_detected", sep = "")
      
      # remove columns where atx_col is na
      if(x == "nt" | (x == "tm" & y == "SFE-M")) { 
        # tm also included here to remove sample from 8-23-2022 as we did not have enough to analyze
        temp = data_river[[x]][[y]][-which(is.na(data_river[[x]][[y]][atx_col])),]
      } else {
        temp = data_river[[x]][[y]]
      }
      
      # run ISA
      print(paste(x, y, "ISA"))
      summary(multipatt(temp[,start_col:ncol(data[[x]])], as.vector(temp[atx_col])[[1]]), 
              func = "r.g", control = how(nperm = 999))
    }
  })
})
# RUS:
# nt (detected) phormidium & oscillatoria
# tac (detected) phormidium
# SAL:
# nt nothing
# tm not, not
# tac no test
# SFE:
# nt (detected) rophalodia
# tm (none) nostoc
# tac: (detected) microcoleus, (none) nostoc


## (c) NMDS

# get NMDS data
NMDS_data_river <- lapply(sample_types, function(x) {
    temp_list = lapply(names(data_river[[x]]), function(y) {
                      if((y == "SFE-M") | (y == "RUS" & x != "tm")) {
                        return(getNMDSdata(data_river[[x]][[y]], start_col, 
                                    end_col = ncol(data[[x]])))
                      } else {
                        return(NA)
                      }
  })
  names(temp_list) = names(data_river[[x]])
  return(temp_list)
})
names(NMDS_data_river) = sample_types

# (i) with atx category
lapply(sample_types, function(x) {
  lapply(names(data_river[[x]]), function(y) {
    
    if(class(NMDS_data_river[[x]][[y]]) == "list") {
  
      # set column name for anatoxin based on sample type
      atx_col = paste(str_to_upper(x), "_atx_category", sep = "")
      
      makeNMDSplot(NMDS_data_river[[x]][[y]], FALSE, FALSE, color = atx_col, shape = "site")
    }
  })
})
# visually distinct for TM (note only one high sample for TM with current designation)
# NT obviously a mess with all rivers in one graph

#### (4) How continuously present are other anatoxin producers? ####

# copy code from 4a