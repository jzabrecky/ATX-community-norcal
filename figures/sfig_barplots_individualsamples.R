#### Barplots for individual samples
### Jordan Zabrecky
## last edited: 02.12.2026

# This script creates bar plots for each individual sample rather than aggregating
# by time, river, etc. One supplemental figure will be created for each
# sample type

#### (1) Loading libraries & data ####

# libraries
lapply(c("tidyverse", "plyr", "ggtext", "cowplot"), require, character.only = T)

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
            text = element_text(size = 10), strip.text = element_text(size = 10)))

# custom palette
palette <- c("#FBF6B0", "#C5BD53", "#777122", "#C2DFFF", "#5E9DE0", "#205288", 
             "#C0ED96", "#7AB048", "#3D631A", "#CBC5F6", "#8A80CF", "#61389E")
# otherwise use built-in palette, "Set3"

# color for other or unknown
end_color <- "lightgray"

# no unknown in TM samples- remove that category
microscopy$tm <- microscopy$tm %>% 
  filter(broader != "Unknown")

# put "misc. other" at the end of NT
microscopy$nt <- microscopy$nt %>% 
  # keep alphabetical otherwise
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(microscopy$nt$broader))[-which(str_sort(unique(microscopy$nt$broader)) == "Misc. Other")],
                                 "Misc. Other")))

## (a) TM (plotting each separately due to different categories)
tm_morpho <- barplot(microscopy$tm, x = "date_site_reach", y = "percent", fill = "broader", facet_wrap = "site") +
  labs(x = NULL, y = "Relative Abundance", title = "Morphological Data") +
  scale_fill_discrete(name = NULL, 
                      palette = palette,
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
                      palette = c(palette[1:(length(unique(microscopy$tac$broader)) - 1)], end_color),
                      labels = c("Diatoms (other than *Epithemia*)", "*Epithemia*",
                                 "*Geitlerinema*", "*Microcoleus*", "*Nostoc*", "Other Filamentous Cyanobacteria",
                                 "Other N-fixing Cyanobacteria", "Unicellular Cyanobacteria", "Unknown")) +
  facet_wrap(~site, scales = "free_x",
             labeller = as_labeller(c(`SAL` = "Salmon River", `SFE-M` = "South Fork Eel River",
                                      `RUS` = "Russian River")))
tac_morpho

## (c) NT
nt_morpho <- barplot(microscopy$nt, x = "date_site_reach", y = "percent", fill = "broader_factor", facet_wrap = "site") +
  labs(x = NULL, y = "Relative Abundance", title = "Morphological Data") +
  scale_fill_discrete(name = NULL,
                      palette = c(palette[-length(palette)], end_color),
                      labels = c("*Anabaena* or *Cylindrospermum*", "*Cladophora*",
                                 "Diatoms (other than *Epithemia* or *Rophalodia*)",
                                 "*Epithemia* or *Rhopalodia*", "*Microcoleus*", "*Nostoc*", "Other Filamentous Cyanobacteria",
                                 "Other Green Algae", "Other N-fixing Cyanobacteria", "*Spirogyra*",
                                 "Unicellular Cyanobacteria", "Miscellaneous Other")) +
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

## (a) TM

# add broader grouping
phylum_class_tm <- microbial_grouping(molecular$tm, "phylum_class", 0.04) # pick threshold based on ~12 groupings

# factor so other is at the end
phylum_class_tm <- phylum_class_tm %>% 
  mutate(broader_factor = factor(broader,
         levels = c(str_sort(unique(phylum_class_tm$broader))[-which(str_sort(unique(phylum_class_tm$broader)) == "Other")],
                    "Other")))

# plot
tm_molec<- barplot(phylum_class_tm, x = "date_site_reach", y = "relative_abundance", 
                   fill = "broader_factor", facet_wrap = "site") +
  labs(x = NULL, y = "Relative Abundance", title = "Molecular Data") +
  scale_fill_discrete(name = "Phylum - Class", 
                      # customizing a bit further so Proteobacteria all have same hue (purple)
                      palette = c(palette[-c(8)], end_color)) +
  facet_wrap(~site, scales = "free_x",
             labeller = as_labeller(c(`SAL` = "Salmon River", `SFE-M` = "South Fork Eel River")))
tm_molec

## (b) TAC

# add broader grouping
phylum_class_tac <- microbial_grouping(molecular$tac, "phylum_class", 0.08) # pick threshold based on ~12 groupings

# factor so other is at the end
phylum_class_tac <- phylum_class_tac %>% 
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(phylum_class_tac$broader))[-which(str_sort(unique(phylum_class_tac$broader)) == "Other")],
                                            "Other")))

# plot
tac_molec <- barplot(phylum_class_tac, x = "date_site_reach", y = "relative_abundance", 
                   fill = "broader_factor", facet_wrap = "site") +
  labs(x = NULL, y = "Relative Abundance", title = "Molecular Data") +
  scale_fill_discrete(name = "Phylum - Class",
                      # weirdness because I don't want light purple next to light gray
                      palette = c(palette[c(1:(length(unique(phylum_class_tac$broader_factor))-2),
                                            (length(unique(phylum_class_tac$broader_factor))))], end_color)) +
  facet_wrap(~site, scales = "free_x",
             labeller = as_labeller(c(`SAL` = "Salmon River", `SFE-M` = "South Fork Eel River",
                                      `RUS` = "Russian River")))
tac_molec

## (c) NT

# add broader grouping
phylum_class_nt <- microbial_grouping(molecular$nt, "phylum_class", 0.045) # pick threshold based on ~12 groupings

# factor so other is at the end
phylum_class_nt <- phylum_class_nt %>% 
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(phylum_class_nt$broader))[-which(str_sort(unique(phylum_class_nt$broader)) == "Other")],
                                            "Other")))

# plot
nt_molec <- barplot(phylum_class_nt, x = "date_site_reach", y = "relative_abundance", 
                     fill = "broader_factor", facet_wrap = "site") +
  labs(x = NULL, y = "Relative Abundance", title = "Molecular Data") +
  scale_fill_discrete(name = "Phylum - Class",
                      # weirdness because I don't want light purple next to light gray
                      palette = c(palette, end_color)) +
  facet_wrap(~site, scales = "free_x",
             labeller = as_labeller(c(`SAL` = "Salmon River", `SFE-M` = "South Fork Eel River",
                                      `RUS` = "Russian River")))
nt_molec

##### (3) 16S Cyanobacteria Order & Genus ####

# pair phylum and class names
molecular <- lapply(molecular, function(x) {
  y <- x %>% 
    mutate(order_genus = paste(order, " - *", genus, "*", sep = ""))
  
  return(y)
})

## (a) TM

# add broader grouping
cyano_order_genus_tm <- microbial_grouping(molecular$tm, "order_genus", 0.12, TRUE) # pick threshold based on ~12 groupings

# factor so other is at the end
cyano_order_genus_tm <- cyano_order_genus_tm %>% 
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(cyano_order_genus_tm$broader))[-which(str_sort(unique(cyano_order_genus_tm$broader)) == "Other")],
                                            "Other")))

# plot
tm_molec_cyano <- barplot(cyano_order_genus_tm, x = "date_site_reach", y = "relative_abundance", 
                   fill = "broader_factor", facet_wrap = "site") +
  labs(x = NULL, y = "Relative Abundance", title = "Molecular Data") +
  scale_fill_discrete(name = "Phylum - Class", 
                      palette = c(palette, end_color)) +
  facet_wrap(~site, scales = "free_x",
             labeller = as_labeller(c(`SAL` = "Salmon River", `SFE-M` = "South Fork Eel River")))
tm_molec_cyano
# deciding that there are too much other so do each separately with own scale!

## (a.1) TM- SFE only
cyano_order_genus_tm_sfe <- microbial_grouping(molecular$tm %>% filter(site == "SFE-M"), "order_genus", 0.12, TRUE) # pick threshold based on ~12 groupings
cyano_order_genus_tm_sfe <- cyano_order_genus_tm_sfe %>% 
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(cyano_order_genus_tm_sfe$broader))[-which(str_sort(unique(cyano_order_genus_tm_sfe$broader)) == "Other")],
                                            "Other")))
tm_molec_cyano_sfe <- barplot(cyano_order_genus_tm_sfe, x = "date_site_reach", y = "relative_abundance", 
                              fill = "broader_factor") +
  labs(x = NULL, y = "Relative Abundance", title = "South Fork Eel River TM") +
  scale_fill_discrete(name = "Order - *Genus*", 
                      palette = c(palette, end_color))
tm_molec_cyano_sfe

## (a.2) TM- SAL only
cyano_order_genus_tm_sal <- microbial_grouping(molecular$tm %>% filter(site == "SAL"), "order_genus", 0.005, TRUE) # pick threshold based on ~12 groupings
cyano_order_genus_tm_sal <- cyano_order_genus_tm_sal %>% 
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(cyano_order_genus_tm_sal$broader))[-which(str_sort(unique(cyano_order_genus_tm_sal$broader)) == "Other")],
                                            "Other")))
tm_molec_cyano_sal <- barplot(cyano_order_genus_tm_sal, x = "date_site_reach", y = "relative_abundance", 
                              fill = "broader_factor") +
  labs(x = NULL, y = "Relative Abundance", title = "Salmon River TM") +
  scale_fill_discrete(name = "Order - *Genus*", 
                      palette = c(palette, end_color))
tm_molec_cyano_sal

## (b) TAC

# add broader grouping
cyano_order_genus_tac <- microbial_grouping(molecular$tac, "order_genus", 0.40, TRUE) # pick threshold based on ~12 groupings

# factor so other is at the end
cyano_order_genus_tac <- cyano_order_genus_tac %>% 
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(cyano_order_genus_tac$broader))[-which(str_sort(unique(cyano_order_genus_tac$broader)) == "Other")],
                                            "Other")))

# plot
tac_molec_cyano <- barplot(cyano_order_genus_tac, x = "date_site_reach", y = "relative_abundance", 
                          fill = "broader_factor", facet_wrap = "site") +
  labs(x = NULL, y = "Relative Abundance", title = "Molecular Data") +
  scale_fill_discrete(name = "Phylum - Class", 
                      palette = c(palette, end_color)) +
  facet_wrap(~site, scales = "free_x",
             labeller = as_labeller(c(`SAL` = "Salmon River", `SFE-M` = "South Fork Eel River",
                                      `RUS` = "Russian River")))
tac_molec_cyano
# deciding that there are too much other so do each separately with own scale!

## (b.1) TAC- SFE only
cyano_order_genus_tac_sfe <- microbial_grouping(molecular$tac %>% filter(site == "SFE-M"), "order_genus", 0.15, TRUE) # pick threshold based on ~12 groupings
cyano_order_genus_tac_sfe <- cyano_order_genus_tac_sfe %>% 
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(cyano_order_genus_tac_sfe$broader))[-which(str_sort(unique(cyano_order_genus_tac_sfe$broader)) == "Other")],
                                            "Other")))
tac_molec_cyano_sfe <- barplot(cyano_order_genus_tac_sfe, x = "date_site_reach", y = "relative_abundance", 
                              fill = "broader_factor") +
  labs(x = NULL, y = "Relative Abundance", title = "South Fork Eel River TAC") +
  scale_fill_discrete(name = "Order - *Genus*", 
                      palette = c(palette, end_color))
tac_molec_cyano_sfe

## (b.2) TAC- RUS only
cyano_order_genus_tac_rus <- microbial_grouping(molecular$tac %>% filter(site == "RUS"), "order_genus", 0.213, TRUE) # pick threshold based on ~12 groupings
cyano_order_genus_tac_rus <- cyano_order_genus_tac_rus %>% 
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(cyano_order_genus_tac_rus$broader))[-which(str_sort(unique(cyano_order_genus_tac_rus$broader)) == "Other")],
                                            "Other")))
tac_molec_cyano_rus <- barplot(cyano_order_genus_tac_rus, x = "date_site_reach", y = "relative_abundance", 
                               fill = "broader_factor") +
  labs(x = NULL, y = "Relative Abundance", title = "Russian River TAC") +
  scale_fill_discrete(name = "Order - *Genus*", 
                      palette = c(palette, end_color))
tac_molec_cyano_rus

## (b.3) TAC- SAL only
cyano_order_genus_tac_sal <- microbial_grouping(molecular$tac %>% filter(site == "SAL"), "order_genus", 0.09, TRUE) # pick threshold based on ~12 groupings
cyano_order_genus_tac_sal <- cyano_order_genus_tac_sal %>% 
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(cyano_order_genus_tac_sal$broader))[-which(str_sort(unique(cyano_order_genus_tac_sal$broader)) == "Other")],
                                            "Other")))
tac_molec_cyano_sal <- barplot(cyano_order_genus_tac_sal, x = "date_site_reach", y = "relative_abundance", 
                               fill = "broader_factor") +
  labs(x = NULL, y = "Relative Abundance", title = "Russian River TAC") +
  scale_fill_discrete(name = "Order - *Genus*", 
                      palette = c(palette[1:(length(unique(cyano_order_genus_tac_sal$broader))-1)], end_color))
tac_molec_cyano_sal

## (c) NT

# add broader grouping
cyano_order_genus_nt <- microbial_grouping(molecular$nt, "order_genus", 0.15, TRUE) # pick threshold based on ~12 groupings

# factor so other is at the end
cyano_order_genus_nt <- cyano_order_genus_nt %>% 
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(cyano_order_genus_nt$broader))[-which(str_sort(unique(cyano_order_genus_nt$broader)) == "Other")],
                                            "Other")))

# plot
nt_molec_cyano <- barplot(cyano_order_genus_nt, x = "date_site_reach", y = "relative_abundance", 
                           fill = "broader_factor", facet_wrap = "site") +
  labs(x = NULL, y = "Relative Abundance", title = "Molecular Data") +
  scale_fill_discrete(name = "Phylum - Class", 
                      palette = c(palette, end_color)) +
  facet_wrap(~site, scales = "free_x",
             labeller = as_labeller(c(`SAL` = "Salmon River", `SFE-M` = "South Fork Eel River",
                                      `RUS` = "Russian River")))
nt_molec_cyano
# deciding that there are too much other so do each separately with own scale!

## (b.1) NT- SFE only
cyano_order_genus_nt_sfe <- microbial_grouping(molecular$nt %>% filter(site == "SFE-M"), "order_genus", 0.15, TRUE) # pick threshold based on ~12 groupings
cyano_order_genus_nt_sfe <- cyano_order_genus_nt_sfe %>% 
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(cyano_order_genus_nt_sfe$broader))[-which(str_sort(unique(cyano_order_genus_nt_sfe$broader)) == "Other")],
                                            "Other")))
nt_molec_cyano_sfe <- barplot(cyano_order_genus_tac_sfe, x = "date_site_reach", y = "relative_abundance", 
                               fill = "broader_factor") +
  labs(x = NULL, y = "Relative Abundance", title = "South Fork Eel River NT") +
  scale_fill_discrete(name = "Order - *Genus*", 
                      palette = c(palette, end_color))
nt_molec_cyano_sfe

## (b.2) NT- RUS only
cyano_order_genus_nt_rus <- microbial_grouping(molecular$nt %>% filter(site == "RUS"), "order_genus", 0.0822, TRUE) # pick threshold based on ~12 groupings
cyano_order_genus_nt_rus <- cyano_order_genus_nt_rus %>% 
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(cyano_order_genus_nt_rus$broader))[-which(str_sort(unique(cyano_order_genus_nt_rus$broader)) == "Other")],
                                            "Other")))
nt_molec_cyano_rus <- barplot(cyano_order_genus_nt_rus, x = "date_site_reach", y = "relative_abundance", 
                               fill = "broader_factor") +
  labs(x = NULL, y = "Relative Abundance", title = "Russian River NT") +
  scale_fill_discrete(name = "Order - *Genus*", 
                      palette = c(palette, end_color))
nt_molec_cyano_rus

## (b.3) NT- SAL only
cyano_order_genus_nt_sal <- microbial_grouping(molecular$nt %>% filter(site == "SAL"), "order_genus", 0.13, TRUE) # pick threshold based on ~12 groupings
cyano_order_genus_nt_sal <- cyano_order_genus_nt_sal %>% 
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(cyano_order_genus_nt_sal$broader))[-which(str_sort(unique(cyano_order_genus_nt_sal$broader)) == "Other")],
                                            "Other")))
nt_molec_cyano_sal <- barplot(cyano_order_genus_nt_sal, x = "date_site_reach", y = "relative_abundance", 
                               fill = "broader_factor") +
  labs(x = NULL, y = "Relative Abundance", title = "Salmon River NT") +
  scale_fill_discrete(name = "Order - *Genus*", 
                      palette = c(palette[1:(length(unique(cyano_order_genus_nt_sal$broader))-1)], end_color))
nt_molec_cyano_sal
#### (4) Functional Groups TBD ####

#### (5) Putting All Plots Together ####

# function to make legend smaller
small_legend <- function(plot) {
  new_plot = plot + 
    theme(legend.key.size = unit(0.5, 'cm'),
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'),
          legend.title = element_text(size=10),
          legend.text = element_markdown())
  
  return(new_plot)
}

## (a) TM
tm_plots <- plot_grid(tm_morpho, tm_molec, ncol = 1, align = "hv")
tm_plots

# save
ggsave("./figures/tiff_files/sfig_tm_relativeabundances.tiff", 
       width = 18, height = 25, unit = "cm", dpi = 600)

# do cyanobacteria order and genus separately
tm_cyano <- plot_grid(tm_molec_cyano_sal, tm_molec_cyano_sfe, ncol = 1, align = "hv")
tm_cyano

# save
ggsave("./figures/tiff_files/sfig_tm_cyano_relativeabundances.tiff", 
       width = 18, height = 25, unit = "cm", dpi = 600)

## (b) TAC
tac_plots <- plot_grid(tac_morpho, tac_molec, ncol = 1, align = "hv")
tac_plots

# save
ggsave("./figures/tiff_files/sfig_tac_relativeabundances.tiff", 
       width = 18, height = 25, unit = "cm", dpi = 600)

# do cyanobacteria order and genus separately
tac_cyano <- plot_grid(small_legend(tac_molec_cyano_rus), 
                       small_legend(tac_molec_cyano_sal), 
                       small_legend(tac_molec_cyano_sfe), ncol = 1, align = "hv")
tac_cyano
# not perfect but will manually move the russian legend down

# save
ggsave("./figures/tiff_files/sfig_tac_cyano_relativeabundances.tiff", 
       width = 18, height = 25, unit = "cm", dpi = 600)

## (c) NT
nt_plots <- plot_grid(nt_morpho, nt_molec, ncol = 1, align = "hv")
nt_plots

# save
ggsave("./figures/tiff_files/sfig_nt_relativeabundances.tiff", 
       width = 18, height = 25, unit = "cm", dpi = 600)

# do cyanobacteria order and genus separately
nt_cyano <- plot_grid(small_legend(nt_molec_cyano_rus), 
                       small_legend(nt_molec_cyano_sal), 
                       small_legend(nt_molec_cyano_sfe), ncol = 1, align = "hv")
nt_cyano
# not perfect but will manually move the russian legend down

# save
ggsave("./figures/tiff_files/sfig_nt_cyano_relativeabundances.tiff", 
       width = 18, height = 25, unit = "cm", dpi = 600)
