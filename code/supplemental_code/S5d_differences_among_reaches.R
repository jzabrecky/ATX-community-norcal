#### Testing for difference among reaches
### Jordan Zabrecky
## last edited: 07.11.2026

# Here, we test to see if there are differences among assemblages 
# to make sure we are properly address Q2. In the case that there ARE,
# we will use the 'strata' argument for PERMANOVA

#### (1) Loading libraries & data ####

# libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot"), require, character.only = T)

# source community functions
source("./code/supplemental_code/S4b_community_analyses_func.R")

# name sample types
sample_types <- c("NT", "TAC", "TM")

# read in transformed microscopy data
microscopy_data <- lapply(list.files(path = "./data/morphological/transformed/", pattern = ".csv"),
                    function(x) read.csv(paste("./data/morphological/transformed/", x, sep = "")) %>% 
                      mutate(month = month(field_date)) %>% 
                      relocate(month, .after = field_date))
names(microscopy_data) <- sample_types

# read in transformed 16s data
molec_data <- lapply(list.files(path = "./data/molecular/transformed/", pattern = "16s_nochimera"),
               function(x) read.csv(paste("./data/molecular/transformed/", x, sep = ""))  %>% 
                 mutate(month = month(field_date)) %>% 
                 relocate(month, .after = field_date))
names(molec_data) <- sample_types

# read in predicted functional content
func_data <- lapply(list.files(path = "./data/molecular/transformed/", pattern = "PICRUSt2"),
                    function(x) read.csv(paste("./data/molecular/transformed/", x, sep = ""))  %>% 
                      mutate(month = month(field_date)) %>% 
                      relocate(month, .after = field_date))
names(func_data) <- sample_types

# split for each river
microscopy_river <- lapply(microscopy_data, function(x) split(x, x$site))
molec_river <- lapply(molec_data, function(x) split(x, x$site))
func_river <- lapply(func_data, function(x) split(x, x$site))

# assign start columns
micro_start <- 6
molec_start <- 7
func_start <- 6

# make p table to hold significant tests
p_table <-  data.frame(test = NA,
                       sample_type = NA,
                       river = NA,
                       data_type = NA,
                       p_value = NA,
                       F_stat = NA)

#### (2) Do they differ? Microscopy ####

# run through all samples
set.seed(1)
for(s in c("SFE-M", "RUS")) {
  for(i in c("NT", "TM", "TAC")) {
    if(s == "RUS" & i == "TM") {
    } else {
      permanova = runPERMANOVA(data = microscopy_river[[i]][[s]], start_col = micro_start, 
                               group = microscopy_river[[i]][[s]]$site_reach)
      
      p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                           sample_type = i,
                                           river = s,
                                           data_type = "microscopy",
                                           p_value = permanova$`Pr(>F)`[1],
                                           F_stat = permanova$`F`[1]))
      
      permdisp = betadisper(vegdist(microscopy_river[[i]][[s]][,micro_start:ncol(microscopy_river[[i]][[s]])], 
                                    method = "bray"), microscopy_river[[i]][[s]]$site_reach)
      test = adonis2(dist(permdisp$distances) ~ microscopy_river[[i]][[s]]$site_reach)
      
      p_table <- rbind(p_table, data.frame(test = "PERMDISP",
                                           sample_type = i,
                                           river = s,
                                           data_type = "microscopy",
                                           p_value = test$`Pr(>F)`[1],
                                           F_stat = test$`F`[1]))
    }
  }
}

#view(p_table %>% filter(data_type == "microscopy"))
# no significant PERMANOVA difference!

## making PCoAplots

# SFE Microscopy
set.seed(1)
PCoA_list_sfe <- lapply(sample_types, function(x) getPCoAdata(microscopy_river[[x]][["SFE-M"]], micro_start))
PCoA_plots_sfe <- lapply(PCoA_list_sfe, function(x) makePCoAplot(x, color = "site_reach", shape = "site_reach",
                                                                 stat_ellipse = FALSE))
names(PCoA_plots_sfe) <- sample_types
lapply(sample_types, function(x) print(PCoA_plots_sfe[[x]] + ggtitle(paste("SFE-M", x))))

# RUS Microscopy
set.seed(1)
PCoA_list_rus <- lapply(c("NT", "TAC"), function(x) getPCoAdata(microscopy_river[[x]][["RUS"]], micro_start))
PCoA_plots_rus <- lapply(PCoA_list_rus, function(x) makePCoAplot(x, color = "site_reach", shape = "site_reach",
                                                                 stat_ellipse = FALSE))
names(PCoA_plots_rus) <- c("NT", "TAC")
lapply(sample_types, function(x) print(PCoA_plots_rus[[x]] + ggtitle(paste("RUS", x))))

# SAL Microscopy
set.seed(1)
PCoA_list_sal <- lapply(c("NT", "TM"), function(x) getPCoAdata(microscopy_river[[x]][["SAL"]], micro_start))
PCoA_plots_sal <- lapply(PCoA_list_sal, function(x) makePCoAplot(x, color = "site_reach", shape = "site_reach",
                                                                 stat_ellipse = FALSE))
names(PCoA_plots_sal) <- c("NT", "TM")
lapply(sample_types, function(x) print(PCoA_plots_sal[[x]] + ggtitle(paste("SAL", x))))

#### (3) Do they differ? Molecular ####

# run through all samples
set.seed(1)
for(s in c("SFE-M", "RUS")) {
  for(i in c("NT", "TM", "TAC")) {
    if(s == "RUS" & i == "TM") {
    } else {
      permanova = runPERMANOVA(data = molec_river[[i]][[s]], start_col = molec_start, 
                               group = molec_river[[i]][[s]]$site_reach)
      
      p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                           sample_type = i,
                                           river = s,
                                           data_type = "molecular",
                                           p_value = permanova$`Pr(>F)`[1],
                                           F_stat = permanova$`F`[1]))
      
      permdisp = betadisper(vegdist(molec_river[[i]][[s]][,molec_start:ncol(molec_river[[i]][[s]])], 
                                    method = "bray"), molec_river[[i]][[s]]$site_reach)
      test = adonis2(dist(permdisp$distances) ~ molec_river[[i]][[s]]$site_reach)
      
      p_table <- rbind(p_table, data.frame(test = "PERMDISP",
                                           sample_type = i,
                                           river = s,
                                           data_type = "molecular",
                                           p_value = test$`Pr(>F)`[1],
                                           F_stat = test$`F`[1]))
    }
  }
}

#view(p_table %>% filter(data_type == "molecular"))
# Russian Molecular samples differ <- use strata argument for these?

## Make PCoA plots

# SFE-Molecular
set.seed(1)
PCoA_list_sfe_molec <- lapply(sample_types, function(x) getPCoAdata(molec_river[[x]][["SFE-M"]], molec_start))
PCoA_plots_sfe_molec <- lapply(PCoA_list_sfe_molec, function(x) makePCoAplot(x, color = "site_reach", shape = "site_reach",
                                                                             stat_ellipse = FALSE))
names(PCoA_plots_sfe_molec) <- sample_types
lapply(sample_types, function(x) print(PCoA_plots_sfe_molec[[x]] + ggtitle(paste("SFE-M", x))))

# RUS- Molecular
set.seed(1)
PCoA_list_rus_molec <- lapply(c("NT", "TAC"), function(x) getPCoAdata(molec_river[[x]][["RUS"]], molec_start))
PCoA_plots_rus_molec <- lapply(PCoA_list_rus_molec, function(x) makePCoAplot(x, color = "site_reach", shape = "site_reach",
                                                                             stat_ellipse = FALSE))
names(PCoA_plots_rus_molec) <- c("NT", "TAC")
lapply(sample_types, function(x) print(PCoA_plots_rus_molec[[x]] + ggtitle(paste("RUS", x))))

# SAL- Molecular
set.seed(1)
PCoA_list_sal_molec <- lapply(c("NT", "TM"), function(x) getPCoAdata(molec_river[[x]][["SAL"]], molec_start))
PCoA_plots_sal_molec <- lapply(PCoA_list_sal_molec, function(x) makePCoAplot(x, color = "site_reach", shape = "site_reach",
                                                                             stat_ellipse = FALSE))
names(PCoA_plots_sal_molec) <- c("NT", "TM")
lapply(sample_types, function(x) print(PCoA_plots_sal_molec[[x]] + ggtitle(paste("SAL", x))))

#### (4) Do they differ? Functional ####

# run through all samples
set.seed(1)
for(s in c("SFE-M", "RUS")) {
  for(i in c("NT", "TM", "TAC")) {
    if(s == "RUS" & i == "TM") {
    } else {
      permanova = runPERMANOVA(data = func_river[[i]][[s]], start_col = func_start, 
                               group = func_river[[i]][[s]]$site_reach)
      
      p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                           sample_type = i,
                                           river = s,
                                           data_type = "functional",
                                           p_value = permanova$`Pr(>F)`[1],
                                           F_stat = permanova$`F`[1]))
      
      permdisp = betadisper(vegdist(func_river[[i]][[s]][,func_start:ncol(func_river[[i]][[s]])], 
                                    method = "bray"), func_river[[i]][[s]]$site_reach)
      test = adonis2(dist(permdisp$distances) ~ func_river[[i]][[s]]$site_reach)
      
      p_table <- rbind(p_table, data.frame(test = "PERMDISP",
                                           sample_type = i,
                                           river = s,
                                           data_type = "functional",
                                           p_value = test$`Pr(>F)`[1],
                                           F_stat = test$`F`[1]))
    }
  }
}

# view(p_table %>% filter(data_type == "functional"))
# Russian TAC interestingly?


## Make PCoA plots

# SFE-Functional
set.seed(1)
PCoA_list_sfe_func <- lapply(sample_types, function(x) getPCoAdata(func_river[[x]][["SFE-M"]], func_start))
PCoA_plots_sfe_func <- lapply(PCoA_list_sfe_func, function(x) makePCoAplot(x, color = "site_reach", shape = "site_reach",
                                                                           stat_ellipse = FALSE))
names(PCoA_plots_sfe_func) <- sample_types
lapply(sample_types, function(x) print(PCoA_plots_sfe_func[[x]] + ggtitle(paste("SFE-M", x))))

# RUS- Functional
set.seed(1)
PCoA_list_rus_func <- lapply(c("NT", "TAC"), function(x) getPCoAdata(func_river[[x]][["RUS"]], func_start))
PCoA_plots_rus_func <- lapply(PCoA_list_rus_func, function(x) makePCoAplot(x, color = "site_reach", shape = "site_reach",
                                                                           stat_ellipse = FALSE))
names(PCoA_plots_rus_func) <- c("NT", "TAC")
lapply(sample_types, function(x) print(PCoA_plots_rus_func[[x]] + ggtitle(paste("RUS", x))))

# SAL- Functional
set.seed(1)
PCoA_list_sal_func <- lapply(c("NT", "TM"), function(x) getPCoAdata(func_river[[x]][["SAL"]], func_start))
PCoA_plots_sal_func <- lapply(PCoA_list_sal_func, function(x) makePCoAplot(x, color = "site_reach", shape = "site_reach",
                                                                           stat_ellipse = FALSE))
names(PCoA_plots_sal_func) <- c("NT", "TM")
lapply(sample_types, function(x) print(PCoA_plots_sal_func[[x]] + ggtitle(paste("SAL", x))))
