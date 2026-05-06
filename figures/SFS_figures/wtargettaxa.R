


# libraries
lapply(c("tidyverse", "plyr", "ggtext", "cowplot"), require, character.only = T)

# read in untransformed relative abundances for microscopy
tm_mic <- read.csv("./data/morphological/tm_algalonly.csv")
tac_mic <- read.csv("./data/morphological/tac_algalonly.csv")

# add into list
microscopy_wide <- list(tm_mic, tac_mic)
names(microscopy_wide) <- c("tm", "tac")

# create a longer version of the unaltered data for bar plots of relative abundances
# and filter for 2022 data only
microscopy <- lapply(microscopy_wide, 
                     function(x) x = x %>% pivot_longer(cols = all_of(c(5:ncol(x))), values_to = "percent",
                                                        names_to = "taxa") %>% 
                       filter(year(field_date) == 2022))

# get barplot functions
source("./code/supplemental_code/S4b_community_analyses_func.R")
source("./code/supplemental_code/S4c_grouping_func.R")

# add in site reach w/ date
microscopy <- lapply(microscopy, function(x) 
  x = x %>% mutate(date_site_reach = paste(field_date, ", ", site_reach, sep = "")) %>% 
    relocate(date_site_reach, .before = "site_reach"))

microscopy$tm <- target_broader(microscopy$tm)
microscopy$tac <- target_broader(microscopy$tac)


theme_set(theme_bw() + theme(strip.background = element_blank(),
                             axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
                             plot.title = element_text(hjust = 0.5), legend.text = element_markdown(),
                             text = element_text(size = 10), strip.text = element_text(size = 10)))

# custom palette
palette <- c("#FBF6B0", "#C5BD53", "#777122", "#C2DFFF", "#5E9DE0", "#205288", 
             "#C0ED96", "#7AB048", "#3D631A", "#CBC5F6", "#8A80CF", "#61389E")
# otherwise use built-in palette, "Set3"

# color for other or unknown
end_color <- "lightgray"

## (a) TM (plotting each separately due to different categories)
tm_morpho <- barplot(microscopy$tm %>% filter(broader != "Unknown"), x = "date_site_reach", y = "percent", fill = "broader", facet_wrap = "site") +
  labs(x = NULL, y = "Relative Abundance", title = "Morphological Data") +
  scale_fill_discrete(name = NULL, 
                      palette = palette,
                      labels = c("*Anabaena* or *Cylindrospermum*", "Diatoms (other than *Epithemia*)",
                                 "*Epithemia*", "*Geitlerinema*/*Leptolygnbya*", "Green Algae", "*Microcoleus*", 
                                 "*Nostoc*", "Other Filamentous Cyanobacteria",
                                 "Other N-fixing Cyanobacteria", "Unicellular Cyanobacteria")) +
  facet_wrap(~site, scales = "free_x",
             labeller = as_labeller(c(`SAL` = "Salmon River", `SFE-M` = "South Fork Eel River")))
tm_morpho

## (b) TAC
tac_morpho <- barplot(microscopy$tac %>% filter(site == "SFE-M" | site == "RUS"), x = "date_site_reach", y = "percent", fill = "broader", facet_wrap = "site") +
  labs(x = NULL, y = "Relative Abundance", title = "Morphological Data") +
  scale_fill_discrete(name = NULL,
                      palette = c(palette[1:(length(unique(microscopy$tac$broader)) - 1)], end_color),
                      labels = c("*Anabaena* or *Cylindrospermum*", "Diatoms (other than *Epithemia*)", "*Epithemia*",
                                 "*Geitlerinema*/*Lyngbya*", "Green Algae", "*Microcoleus*", "*Nostoc*", "Other Filamentous Cyanobacteria",
                                 "Other N-fixing Cyanobacteria", "Unicellular Cyanobacteria", "Unknown")) +
  facet_wrap(~site, scales = "free_x",
             labeller = as_labeller(c(`SAL` = "Salmon River", `SFE-M` = "South Fork Eel River",
                                      `RUS` = "Russian River")))
tac_morpho
