#### Barplots for individual samples
### Jordan Zabrecky
## last edited: 02.10.2026

# This script creates bar plots for each individual sample rather than aggregating
# by time, river, etc. One supplemental figure will be created for each
# sample type

# TO-DO: move other categories to the bottom
# working on microbial relative abundance plots

#### (1) Loading libraries & data ####

# libraries
lapply(c("tidyverse", "plyr", "ggtext"), require, character.only = T)

# read in untransformed relative abundances for microscopy
nt_mic <- read.csv("./data/morphological/nt_algalonly.csv")
tm_mic <- read.csv("./data/morphological/tm_algalonly_nomicro.csv")
tac_mic <- read.csv("./data/morphological/tac_algalonly_noanacylgreenalgae.csv")

# add into list
microscopy_wide <- list(nt_mic, tm_mic, tac_mic)
names(microscopy_wide) <- c("nt", "tm", "tac")

# create a longer version of the unaltered data for bar plots of relative abundances
# and filter for 2022 data only
microscopy <- lapply(microscopy_wide, 
                      function(x) x = x %>% pivot_longer(cols = all_of(c(5:ncol(x))), values_to = "percent",
                                               names_to = "taxa") %>% 
                       filter(year(field_date) == 2022))

# read in untransformed microbial/16s data
nt_16s <- read.csv("./data/molecular/16s_nochimera_rarefied_95_FINAL.csv") %>% 
  filter(sample_type == "NT")
tm_16s <- read.csv("./data/molecular/16s_nochimera_rarefied_95_TM_nomicro.csv")
tac_16s <- read.csv("./data/molecular/16s_nochimera_rarefied_95_TAC_noanacyl.csv")

# add into list
molecular <- list(nt_16s, tm_16s, tac_16s)
names(molecular) <- c("nt", "tm", "tac")

# get barplot functions
source("./code/supplemental_code/S4b_community_analyses_func.R")
source("./code/supplemental_code/S4c_grouping_func.R")

# add in site reach w/ date
microscopy <- lapply(microscopy, function(x) 
  x = x %>% mutate(date_site_reach = paste(field_date, ", ", site_reach, sep = "")) %>% 
                     relocate(date_site_reach, .before = "site_reach"))
molecular <- lapply(molecular, function(x) 
  x = x %>% mutate(date_site_reach = paste(field_date, ", ", site_reach, sep = "")) %>% 
                     relocate(date_site_reach, .before = "site_reach"))

#### (2) Morphology Bar Plots ####

# add broader group categories
microscopy$tm <- target_broader(microscopy$tm)
microscopy$tac <- target_broader(microscopy$tac)
microscopy$nt <- nontarget_broader(microscopy$nt)

# add in universal plot themes
theme_set(theme_bw() + theme(strip.background = element_blank(),
            axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
            plot.title = element_text(hjust = 0.5), legend.text = element_markdown(),
            text = element_text(size = 11), strip.text = element_text(size = 12)))

## (a) TM (plotting each separately due to different categories)
tm_morpho <- barplot(microscopy$tm, x = "date_site_reach", y = "percent", fill = "broader", facet_wrap = "site") +
  labs(x = NULL, y = "Relative Abundance", title = "Morphological Data") +
  scale_fill_discrete(name = NULL, 
                      palette = "Set3",
                      labels = c("*Anabaena* or *Cylindrospermum*", "Diatoms (other than *Epithemia*)",
                                 "*Epithemia*", "*Geitlerinema*", "Green Algae", "*Nostoc*", "Other Filamentous Cyanobacteria",
                                 "Other N-fixing Cyanobacteria", "Unicellular Cyanobacteria")) +
  facet_wrap(~site, scales = "free_x",
             labeller = as_labeller(c(`SAL` = "Salmon River", `SFE-M` = "South Fork Eel River")))
tm_morpho

## (b) TAC
tac_morpho <- barplot(microscopy$tac, x = "date_site_reach", y = "percent", fill = "broader", facet_wrap = "site") +
  labs(x = NULL, y = "Relative Abundance", title = "Morphological Data") +
  scale_fill_discrete(name = NULL,
                      palette = "Set3",
                      labels = c("Diatoms (other than *Epithemia*)", "*Epithemia*",
                                 "*Geitlerinema*", "*Microcoleus*", "*Nostoc*", "Other Filamentous Cyanobacteria",
                                 "Other N-fixing Cyanobacteria", "Unicellular Cyanobacteria", "Unknown")) +
  facet_wrap(~site, scales = "free_x",
             labeller = as_labeller(c(`SAL` = "Salmon River", `SFE-M` = "South Fork Eel River",
                                      `RUS` = "Russian River")))
tac_morpho

## (c) NT
nt_morpho <- barplot(microscopy$nt, x = "date_site_reach", y = "percent", fill = "broader", facet_wrap = "site") +
  labs(x = NULL, y = "Relative Abundance", title = "Morphological Data") +
  scale_fill_discrete(name = NULL,
                      palette = "Set3",
                      labels = c("*Anabaena* or *Cylindrospermum*", "*Cladophora*",
                                 "Diatoms (other than *Epithemia* or *Rophalodia*)",
                                 "*Epithemia* or *Rhopalodia*", "*Microcoleus*",
                                 "Miscellaneous Other", "*Nostoc*", "Other Filamentous Cyanobacteria",
                                 "Other Green Algae", "Other N-fixing Cyanobacteria", "*Spirogyra*",
                                 "Unicellular Cyanobacteria")) +
  facet_wrap(~site, scales = "free_x",
             labeller = as_labeller(c(`SAL` = "Salmon River", `SFE-M` = "South Fork Eel River",
                                      `RUS` = "Russian River")))
nt_morpho

# if deciding how to further group things, view summary of groups
#summary <- microscopy$nt %>% 
#  dplyr::group_by(site, broader) %>% 
#  dplyr::summarize(mean = mean(percent))
#view(summary)

#### (4) 16s/Microbial Barplots -- Phylum & Class ####

# pair phylum and class names
molecular <- lapply(molecular, function(x) {
  y <- x %>% 
    mutate(phylum_class = paste(phylum, " - ", class))
  
  return(y)
})
test <- microbial_grouping(molecular$nt, "phylum_class", 0.05)

barplot(test, x = "date_site_reach", y = "relative_abundance", fill = "broader")

# want to maybe always add "other" to the bottom
# see what other papers do