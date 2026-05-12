#### Comparing microscopy data with regard to anatoxin concentrations
### Jordan Zabrecky
## last edited: 05.11.2026

# This script examines how communities as identified by microscopy
# change with increasing anatoxin concentrations with PERMANOVA, NMDS,
# and ISA. Note that as only one Salmon sample contained detectable 
# anatoxin, we omit the Salmon River samples from this analyses

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

#### (4) NMDS Plots ####

# get NMDS for each dataframe (sqrt-transformed!)

## (a) South Fork Eel River
set.seed(1)
NMDS_list_eel <- lapply(sample_types, function(x) getNMDSdata(data_river[[x]]$`SFE-M`, start_col))

# making plots
NMDS_plots_eel <- lapply(NMDS_list_eel, function(x) makeNMDSplot(x, TRUE, TRUE, 
                                                         color = "atx_group", shape = "atx_group"))

lapply(NMDS_plots_eel, print)
# seems like difference between detected and non-detected, but not for variability

## (b) Russian River
set.seed(1)
NMDS_list_rus <- lapply(c("TAC", "NT"), function(x) getNMDSdata(data_river[[x]]$`RUS`, start_col))

# making plots
NMDS_plots_rus <- lapply(NMDS_list_rus, function(x) makeNMDSplot(x, TRUE, TRUE, 
                                                                 color = "atx_group", shape = "atx_group"))

lapply(NMDS_plots_rus, print)
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

# save PERMANOVA results
write.csv(p_table[-1,], "./data/PERMANOVA_results/Q3_microscopy.csv", row.names = FALSE)

#### (6) What taxa may be indicative of ATX groups? ####

# just going to do each separately

## (a) SFE NT
set.seed(1)
eel_nt_test_det <- multipatt(data_river$NT$`SFE-M`[,start_col:ncol(data_river$NT$`SFE-M`)], 
                     data_river$NT$`SFE-M`$atx_detected, func = "r.g", control = how(nperm = 999))
summary(eel_nt_test_det)
eel_nt_test_group <- multipatt(data_river$NT$`SFE-M`[,start_col:ncol(data_river$NT$`SFE-M`)], 
                             data_river$NT$`SFE-M`$atx_group, func = "r.g", control = how(nperm = 999))
summary(eel_nt_test_group)
write.csv(eel_nt_test_group$sign, "./data/ISA_results/Q3_nt_microscopy_SFE.csv")
# detected: rophalodia, epithemia, anabaena
# high & medium: epithemia, anabaena
# high & medium & low: rophalodia

## (b) SFE TM
set.seed(1)
eel_tm_test_det <- multipatt(data_river$TM$`SFE-M`[,start_col:ncol(data_river$TM$`SFE-M`)], 
                             data_river$TM$`SFE-M`$atx_detected, func = "r.g", control = how(nperm = 999))
summary(eel_tm_test_det)
eel_tm_test_group <- multipatt(data_river$TM$`SFE-M`[,start_col:ncol(data_river$TM$`SFE-M`)], 
                               data_river$TM$`SFE-M`$atx_group, func = "r.g", control = how(nperm = 999))
summary(eel_tm_test_group)
write.csv(eel_tm_test_group$sign, "./data/ISA_results/Q3_tm_microscopy_SFE.csv")
# detected: other coccoids & geitlerinema; not detected: nostoc
# none: nostoc, medium: other coccoids

## (c) SFE TAC
set.seed(1)
eel_tac_test_det <- multipatt(data_river$TAC$`SFE-M`[,start_col:ncol(data_river$TAC$`SFE-M`)], 
                             data_river$TAC$`SFE-M`$atx_detected, func = "r.g", control = how(nperm = 999))
summary(eel_tac_test_det)
eel_tac_test_group <- multipatt(data_river$TAC$`SFE-M`[,start_col:ncol(data_river$TAC$`SFE-M`)], 
                               data_river$TAC$`SFE-M`$atx_group, func = "r.g", control = how(nperm = 999))
summary(eel_tac_test_group)
write.csv(eel_tac_test_group$sign, "./data/ISA_results/Q3_tac_microscopy_SFE.csv")
# not detected: nostoc
# low & medium: nodularia

## (d) RUS NT
set.seed(1)
rus_nt_test_det <- multipatt(data_river$NT$`RUS`[,start_col:ncol(data_river$NT$`RUS`)], 
                              data_river$NT$`RUS`$atx_detected, func = "r.g", control = how(nperm = 999))
summary(rus_nt_test_det)
rus_nt_test_group <- multipatt(data_river$NT$`RUS`[,start_col:ncol(data_river$NT$`RUS`)], 
                                data_river$NT$`RUS`$atx_group, func = "r.g", control = how(nperm = 999))
summary(rus_nt_test_group)
write.csv(rus_nt_test_group$sign, "./data/ISA_results/Q3_nt_microscopy_RUS.csv")
# detected: phormidium_unknown, oscillatoria; not-detected: scenedesmus
# high: miscellaneous oscillatoriales

## (e) RUS TAC
set.seed(1)
rus_tac_test_det <- multipatt(data_river$TAC$`RUS`[,start_col:ncol(data_river$TAC$`RUS`)], 
                              data_river$TAC$`RUS`$atx_detected, func = "r.g", control = how(nperm = 999))
summary(rus_tac_test_det)
rus_tac_test_group <- multipatt(data_river$TAC$`RUS`[,start_col:ncol(data_river$TAC$`RUS`)], 
                                data_river$TAC$`RUS`$atx_group, func = "r.g", control = how(nperm = 999))
summary(rus_tac_test_group)
write.csv(rus_tac_test_group$sign, "./data/ISA_results/Q3_tac_microscopy_RUS.csv")
# miscellaneous oscillatoriales for high; nothing for detected versus not detected
