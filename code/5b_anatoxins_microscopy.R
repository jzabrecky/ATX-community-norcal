#### Comparing microscopy data with regard to anatoxin concentrations
### Jordan Zabrecky
## last edited: 06.01.2026

# This script examines how communities as identified by microscopy
# change with increasing anatoxin concentrations with PERMANOVA, PCoA,
# and ISA. Note that as only one Salmon sample contained detectable 
# anatoxin, we omit the Salmon River samples from this analyses

# Note (7/11): for assemblage comparisons (PERMANOVA), we verified that there
# were no significant differences among reaches, found on script:
# "S5d_differences_among_reaches.R"

#### (1) Loading libraries & data ####

# libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot", 
         "indicspecies"), require, character.only = T)

# read in relative abundance files (data transformed in previous script, "4a_amongrivers_microscopy.R")
abun_data <- lapply(list.files(path = "./data/morphological/transformed/", pattern = ".csv"),
               function(x) read.csv(paste("./data/morphological/transformed/", x, sep = "")))
sample_types <- c("NT", "TAC", "TM")
names(abun_data) <- sample_types

# read in toxin data
atx <- read.csv("./data/field_and_lab/atx_w_categorical_groupings.csv")

# split
atx <- split(atx, atx$sample_type)

# join in data
data <- lapply(names(atx), function(x) {
  left_join(atx[[x]], abun_data[[x]], by = c("field_date", "site", "site_reach", "sample_type"))
})
names(data) <- sample_types

# set start_column
start_col <- 9

# lastly, make long 
data_long <- lapply(data, function(x) {
  return(pivot_longer(x, cols = c(9:ncol(x)), values_to = "percent",
               names_to = "taxa"))
})
names(data_long) <- sample_types

# source code for analyses
source("./code/supplemental_code/S4b_community_analyses_func.R")
source("./code/supplemental_code/S4c_grouping_func.R")

# make broader categories for plotting
data_long$TM <- target_broader(data_long$TM)
data_long$TAC <- target_broader(data_long$TAC)
data_long$NT <- nontarget_broader(data_long$NT)

# lastly, want to analyze rivers separately
# as we know the assemblages often differ among rivers!
data_river <- lapply(data, function(x) split(x, x$`site`))

#### (2) Relative Abundance Bar Plots ####

# put bar plots into lists
barplot_taxa_plots <- lapply(data_long, function(x) barplot(x, x = "atx_group", y = "percent", fill = "taxa"))
barplot_broader_plots <- lapply(data_long, function(x) barplot(x, x = "atx_group", y  = "percent", fill = "broader"))

# titles for plots
titles <- c("Non-Target Samples", 
            "Microcoleus Samples (excluding M)",
            "Anabaena/Cylindrospermum Samples (excluding AC & GA)")

# view plots
for(i in 1:length(barplot_taxa_plots)) {
  print(barplot_taxa_plots[[i]] + labs(title = titles[i]) + facet_wrap(~site))
  print(barplot_broader_plots[[i]] + labs(title = titles[i]) + facet_wrap(~site))
}
# feel like the one main noticeable thing is Nostoc in none samples for target
# would not think this is causation but rather a function of time/succession!

#### (4) PCoA Plots ####

# get PCoA for each dataframe (sqrt-transformed!)

## (a) South Fork Eel River
set.seed(1)
PCoA_list_eel <- lapply(sample_types, function(x) getPCoAdata(data_river[[x]]$`SFE-M`, start_col))

# making plots
PCoA_plots_eel <- lapply(PCoA_list_eel, function(x) makePCoAplot(x, color = "atx_group", shape = "atx_group"))
lapply(PCoA_plots_eel, print)
# seems like difference between detected and non-detected, but less for variability

## (b) Russian River
set.seed(1)
PCoA_list_rus <- lapply(c("TAC", "NT"), function(x) getPCoAdata(data_river[[x]]$`RUS`, start_col))

# making plots
PCoA_plots_rus <- lapply(PCoA_list_rus, function(x) makePCoAplot(x, color = "atx_group", shape = "atx_group"))
lapply(PCoA_plots_rus, print)
# less of a difference observable here

#### (5) Do communities differ with varying ATX concentrations? (PERMANOVA) ####

# create table to save results
p_table <-  data.frame(test = NA,
                       sample_type = NA,
                       river = NA,
                       atx = NA,
                       p_value = NA,
                       F_stat = NA)

## (a) Detected versus no Detected anatoxins?

set.seed(1)
for(s in c("SFE-M", "RUS")) {
  for(i in sample_types) {
    if(s == "RUS" & i == "TM") {
    } else {
      permanova = runPERMANOVA(data = data_river[[i]][[s]], start_col = start_col, 
                                 group = data_river[[i]][[s]]$`atx_detected`)
      
      p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                           sample_type = i,
                                           river = s,
                                           atx = "detected",
                                           p_value = permanova$`Pr(>F)`[1],
                                           F_stat = permanova$`F`[1]))
      
      permdisp = betadisper(vegdist(data_river[[i]][[s]][,start_col:ncol(data_river[[i]][[s]])], 
                                    method = "bray"), data_river[[i]][[s]]$atx_detected)
      test = adonis2(dist(permdisp$distances) ~ data_river[[i]][[s]]$atx_detected)
      
      p_table <- rbind(p_table, data.frame(test = "PERMDISP",
                                           sample_type = i,
                                           river = s,
                                           atx = "detected",
                                           p_value = test$`Pr(>F)`[1],
                                           F_stat = test$`F`[1]))
    }
  }
}
# Difference with South Fork Eel but not for Russian River

set.seed(1)
for(s in c("SFE-M", "RUS")) {
  for(i in sample_types) {
    if(s == "RUS" & i == "TM") {
    } else {
      permanova = runPERMANOVA(data = data_river[[i]][[s]], start_col = start_col, 
                               group = data_river[[i]][[s]]$`atx_group`)
      
      p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                           sample_type = i,
                                           river = s,
                                           atx = "atx_group",
                                           p_value = permanova$`Pr(>F)`[1],
                                           F_stat = permanova$`F`[1]))
      
      permdisp = betadisper(vegdist(data_river[[i]][[s]][,start_col:ncol(data_river[[i]][[s]])], 
                                    method = "bray"), data_river[[i]][[s]]$atx_group)
      test = adonis2(dist(permdisp$distances) ~ data_river[[i]][[s]]$atx_group)
      
      p_table <- rbind(p_table, data.frame(test = "PERMDISP",
                                           sample_type = i,
                                           river = s,
                                           atx = "atx_group",
                                           p_value = test$`Pr(>F)`[1],
                                           F_stat = test$`F`[1]))
    }
  }
}
# Differences in groups for SFE-M and TAC RUS (but barely)
# the second point is interesting

# How about if we remove non-detects from the group?
set.seed(1)
for(s in c("SFE-M", "RUS")) {
  for(i in sample_types) {
    if(s == "RUS" & i == "TM") {
    } else {
      permanova = runPERMANOVA(data = data_river[[i]][[s]] %>% filter(atx_detected == "y"), start_col = start_col, 
                               group = (data_river[[i]][[s]] %>% filter(atx_detected == "y"))$`atx_group`)
      
      p_table <- rbind(p_table, data.frame(test = "PERMANOVA",
                                           sample_type = i,
                                           river = s,
                                           atx = "atx_group_no_nondetect",
                                           p_value = permanova$`Pr(>F)`[1],
                                           F_stat = permanova$`F`[1]))
      
      permdisp = betadisper(vegdist((data_river[[i]][[s]] %>% filter(atx_detected == "y"))[,start_col:ncol(data_river[[i]][[s]])], 
                                    method = "bray"), (data_river[[i]][[s]] %>% filter(atx_detected == "y"))$atx_group)
      test = adonis2(dist(permdisp$distances) ~ (data_river[[i]][[s]] %>% filter(atx_detected == "y"))$atx_group)
      
      p_table <- rbind(p_table, data.frame(test = "PERMDISP",
                                           sample_type = i,
                                           river = s,
                                           atx = "atx_group_no_nondetect",
                                           p_value = test$`Pr(>F)`[1],
                                           F_stat = test$`F`[1]))
    }
  }
}
# no significant differences with non_detects removed
# also seems like there are not enough samples for the Russian River so we will skip that
# for the analyses

# save PERMANOVA results
write.csv(p_table[-1,], "./data/PERMANOVA_results/Q2_microscopy.csv", row.names = FALSE)

#### (6) What taxa may be indicative of ATX groups? ####

# just going to do each separately

## (a) SFE NT

# detected vs non-detected
set.seed(1)
eel_nt_test_det <- multipatt(data_river$NT$`SFE-M`[,start_col:ncol(data_river$NT$`SFE-M`)], 
                     data_river$NT$`SFE-M`$atx_detected, func = "r.g", control = how(nperm = 999))
summary(eel_nt_test_det)
# detected: rophalodia, epithemia, anabaena
write.csv(eel_nt_test_det$sign, "./data/ISA_results/Q2_nt_microscopy_SFE_detects.csv")

# atx groups (when detected)
set.seed(1)
eel_nt_test_group <- multipatt((data_river$NT$`SFE-M` %>% filter(atx_detected == "y"))[,start_col:ncol(data_river$NT$`SFE-M`)], 
                               (data_river$NT$`SFE-M` %>% filter(atx_detected == "y"))$atx_group, func = "r.g", control = how(nperm = 999))
summary(eel_nt_test_group)
write.csv(eel_nt_test_group$sign, "./data/ISA_results/Q2_nt_microscopy_SFE_atx_groups.csv")
# high & medium: epithemia

## (b) SFE TM

# detected versus non-detected
set.seed(1)
eel_tm_test_det <- multipatt(data_river$TM$`SFE-M`[,start_col:ncol(data_river$TM$`SFE-M`)], 
                             data_river$TM$`SFE-M`$atx_detected, func = "r.g", control = how(nperm = 999))
summary(eel_tm_test_det)
# detected: other coccoids & geitlerinema; not detected: nostoc
write.csv(eel_tm_test_det$sign, "./data/ISA_results/Q2_tm_microscopy_SFE_detects.csv")

# atx groups (when detected)
set.seed(1)
eel_tm_test_group <- multipatt((data_river$TM$`SFE-M` %>% filter(atx_detected == "y"))[,start_col:ncol(data_river$TM$`SFE-M`)], 
                                (data_river$TM$`SFE-M` %>% filter(atx_detected == "y"))$atx_group, func = "r.g", control = how(nperm = 999))
summary(eel_tm_test_group)
write.csv(eel_tm_test_group$sign, "./data/ISA_results/Q2_tm_microscopy_SFE_atx_groups.csv")
# nothing!

## (c) SFE TAC

# detected versus non-detected
set.seed(1)
eel_tac_test_det <- multipatt(data_river$TAC$`SFE-M`[,start_col:ncol(data_river$TAC$`SFE-M`)], 
                             data_river$TAC$`SFE-M`$atx_detected, func = "r.g", control = how(nperm = 999))
summary(eel_tac_test_det)
# nostoc in low detection
write.csv(eel_tac_test_det$sign, "./data/ISA_results/Q2_tac_microscopy_SFE_detects.csv")

# atx groups (when detected)
set.seed(1)
eel_tac_test_group <- multipatt((data_river$TAC$`SFE-M` %>% filter(atx_detected == "y"))[,start_col:ncol(data_river$TAC$`SFE-M`)], 
                                (data_river$TAC$`SFE-M` %>% filter(atx_detected == "y"))$atx_group, func = "r.g", control = how(nperm = 999))
summary(eel_tac_test_group)
write.csv(eel_tac_test_group$sign, "./data/ISA_results/Q2_tac_microscopy_SFE_atx_groups.csv")
# low & medium: nodularia

## (d) RUS NT

# detected versus non-detected
set.seed(1)
rus_nt_test_det <- multipatt(data_river$NT$`RUS`[,start_col:ncol(data_river$NT$`RUS`)], 
                              data_river$NT$`RUS`$atx_detected, func = "r.g", control = how(nperm = 999))
summary(rus_nt_test_det)
# detected: phormidium_unknown, oscillatoria; not-detected: scenedesmus
write.csv(rus_nt_test_det$sign, "./data/ISA_results/Q2_nt_microscopy_RUS_detects.csv")

# atx groups (when detected)
set.seed(1)
rus_nt_test_group <- multipatt((data_river$NT$`RUS` %>% filter(atx_detected == "y"))[,start_col:ncol(data_river$NT$`RUS`)], 
                                (data_river$NT$`RUS` %>% filter(atx_detected == "y"))$atx_group, func = "r.g", control = how(nperm = 999))
summary(rus_nt_test_group)
write.csv(rus_nt_test_group$sign, "./data/ISA_results/Q2_nt_microscopy_RUS_atx_group.csv")
# low: non_e_r_diatoms

## (e) RUS TAC

# detects versus non-detects
set.seed(1)
rus_tac_test_det <- multipatt(data_river$TAC$`RUS`[,start_col:ncol(data_river$TAC$`RUS`)], 
                              data_river$TAC$`RUS`$atx_detected, func = "r.g", control = how(nperm = 999))
summary(rus_tac_test_det)
# nothing!
write.csv(rus_tac_test_det$sign, "./data/ISA_results/Q2_tac_microscopy_RUS_detects.csv")

# atx groups (when detected)
set.seed(1)
rus_tac_test_group <- multipatt((data_river$TAC$`RUS` %>% filter(atx_detected == "y"))[,start_col:ncol(data_river$TAC$`RUS`)], 
                                (data_river$TAC$`RUS` %>% filter(atx_detected == "y"))$atx_group, func = "r.g", control = how(nperm = 999))
summary(rus_tac_test_group)
write.csv(rus_tac_test_group$sign, "./data/ISA_results/Q2_tac_microscopy_RUS_atx_groups.csv")
# nothing!

# dropped the low versus high for Russian samples due to lower sample sizes but saving results just in case