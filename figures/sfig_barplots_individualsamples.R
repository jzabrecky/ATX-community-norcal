#### Barplots for individual samples
### Jordan Zabrecky
## last edited: 07.06.2026

# This script creates bar plots for each individual sample rather than aggregating
# by time, river, etc. One supplemental figure will be created for each
# sample type

#### (1) Loading data & libraries ####

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
nt_16s <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_FINAL.csv") %>% 
  filter(sample_type == "NT")
tm_16s <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_TM_nomicro.csv")
tac_16s <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_TAC_noanacyl.csv")

# add into list
molecular <- list(nt_16s, tm_16s, tac_16s)
names(molecular) <- c("nt", "tm", "tac")

# lastly, load in anatoxin data to add '*' to samples with anatoxin detected
atx <- read.csv("./data/field_and_lab/atx_w_categorical_groupings.csv")

# get barplot functions
source("./code/supplemental_code/S4b_community_analyses_func.R")
source("./code/supplemental_code/S4c_grouping_func.R")

# add in atx data, site reach w/ date, and site labels for facet wrapping
microscopy <- lapply(microscopy, function(x)
  x = left_join(x, atx %>% select(field_date, sample_type, site_reach, atx_detected),
                by = c("field_date", "sample_type", "site_reach")) %>% 
    # originally had separate labels, but hard to get what I want, so will do it manually
    mutate(date_site_reach = case_when(atx_detected == "y" ~ paste(field_date, ", ", site_reach, sep = ""),
                                       TRUE ~ paste(field_date, ", ", site_reach, sep = ""))) %>% 
    relocate(date_site_reach, .before = "site_reach") %>% 
    mutate(site_labels = case_when(site == "SFE-M" ~ "South Fork<br>Eel River",
                                   site == "SAL" ~ "Salmon<br>River",
                                   site == "RUS" ~ "Russian<br>River")))
molecular <- lapply(molecular, function(x)
  x = left_join(x %>% mutate(field_date = mdy(field_date)), 
                atx %>% select(field_date, sample_type, site_reach, atx_detected) %>% 
                  mutate(field_date = ymd(field_date)),
                by = c("field_date", "sample_type", "site_reach")) %>% 
    # originally had separate labels, but hard to get what I want, so will do it manually
    mutate(date_site_reach = case_when(atx_detected == "y" ~ paste(field_date, ", ", site_reach, sep = ""),
                                       TRUE ~ paste(field_date, ", ", site_reach, sep = ""))) %>% 
    relocate(date_site_reach, .before = "site_reach") %>% 
    mutate(site_labels = case_when(site == "SFE-M" ~ "South Fork<br>Eel River",
                                   site == "SAL" ~ "Salmon<br>River",
                                   site == "RUS" ~ "Russian<br>River")))

#### (2) Morphology Bar Plots ####

# add in universal plot themes
theme_set(theme_bw() + theme(strip.background = element_blank(),
                             axis.text.x = element_text(size = 7, angle = 60, vjust = 1, hjust=1),
                             plot.title = element_text(hjust = 0.5), legend.text = element_markdown(),
                             text = element_text(size = 7), legend.title = element_markdown(size = 8),
                             strip.text = element_markdown(size = 8)))

# custom palette
palette <- c("#FBF6B0", "#C5BD53", "#777122", "#C2DFFF", "#5E9DE0", "#205288", 
             "#C0ED96", "#7AB048", "#3D631A", "#CBC5F6", "#8A80CF", "#61389E")

# color for other or unknown
end_color <- "lightgray"

# add broader group categories for microscopy
microscopy$tm <- target_broader(microscopy$tm)
microscopy$tac <- target_broader(microscopy$tac)
microscopy$nt <- nontarget_broader(microscopy$nt)

# factor those groups for desired order
microscopy <- lapply(microscopy, function(x) {
  x <- x %>% 
    mutate(figure_groups = case_when(broader == "Nostoc" ~ "Other N-fixing Cyanobacteria",
                                     taxa == "leptolyngbya_and_geitlerinema" ~ "*Leptolyngbya* / *Geitlerinema*",
                                     broader == "Unknown" ~ "Other",
                                     broader == "Misc. Other" ~ "Other",
                                     taxa == "rhopalodia" ~ "*Rhopalodia*",
                                     taxa == "epithemia" ~ "*Epithemia*",
                                     taxa == "e_diatoms" ~ "*Epithemia*",
                                     taxa == "non_e_r_diatoms" ~ "Diatoms (other than<br>*Epithemia* or *Rhopalodia*)",
                                     taxa == "non_e_diatoms" ~ "Diatoms Other than *Epithemia*",
                                     taxa == "cladophora" ~ "*Cladophora*",
                                     taxa == "spirogyra" ~ "*Spirogyra*",
                                     taxa == "anabaena_and_cylindrospermum" ~ "*Anabaena* / *Cylindrospermum*",
                                     taxa == "microcoleus" ~ "*Microcoleus*",
                                     TRUE ~ broader)) %>% 
    mutate(figure_groups_factored = factor(figure_groups, 
                                           levels = c("Diatoms (other than<br>*Epithemia* or *Rhopalodia*)", 
                                                      "Diatoms Other than *Epithemia*", "*Epithemia*",
                                                      "*Rhopalodia*", "Other N-fixing Cyanobacteria",
                                                      "Other Filamentous Cyanobacteria",
                                                      "Unicellular Cyanobacteria", "*Cladophora*", 
                                                      "*Spirogyra*", "Other Green Algae", "Green Algae",
                                                      "*Anabaena* / *Cylindrospermum*", "*Microcoleus*",
                                                      "*Leptolyngbya* / *Geitlerinema*", "Unknown", "Other")))
  return(x)
})

## (a) NT
nt_micro <- barplot(data = microscopy$nt,  x = "date_site_reach", y  = "percent", fill = "figure_groups_factored") +
  scale_fill_discrete("Taxa Group", palette = c(palette, end_color)) +
  labs(x = NULL, y = "Relative Abundance") +
  facet_wrap(~site_labels, scales = "free") +
  theme(axis.text.x = element_markdown(size = 7, color = "#333333"))
nt_micro

## (b) TM (removing taxa group Other as it has no abundance)
tm_micro <- barplot(data = microscopy$tm %>% filter(figure_groups != "Other"),  x = "date_site_reach", y  = "percent", fill = "figure_groups_factored") +
  scale_fill_discrete("Taxa Group", palette = c(palette[c(1:2, 4:6, 8, 10, 12)], end_color)) +
  labs(x = NULL, y = "Relative Abundance") +
  facet_wrap(~site_labels, scales = "free") +
  theme(axis.text.x = element_markdown(size = 7, color = "#333333"))
tm_micro

## (c) TAC
tac_micro <- barplot(data = microscopy$tac,  x = "date_site_reach", y  = "percent", fill = "figure_groups_factored") +
  scale_fill_discrete("Taxa Group", palette = c(palette[c(1:2, 4:6, 11:12)], end_color)) +
  labs(x = NULL, y = "Relative Abundance") +
  facet_wrap(~site_labels, scales = "free") +
  theme(axis.text.x = element_markdown(size = 7, color = "#333333"))
tac_micro

#### (3) Molecular Bar Plots ####

# create phylum-class column, change relative abundance col name, and add site labels
molecular <- lapply(molecular, function(x) x %>% mutate(phylum_class = paste(phylum, " - ", class)) %>% 
                      dplyr::rename(relative_abundance = picrust2_relative_abundance) %>% 
                      mutate(site_labels))

## (a) NT

# add in broader categories (separately for each taxa)
molecular_broader_nt <- microbial_grouping(molecular$nt , "phylum_class", 0.047) %>% 
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(broader))[-which(str_sort(unique(broader)) == "Other")],
                                            "Other")))

# make figure
nt_molec <- barplot(data = molecular_broader_nt,  x = "date_site_reach", y  = "relative_abundance", fill = "broader_factor") +
  scale_fill_discrete("Phylum - Class", palette = c(palette, end_color)) +
  labs(x = NULL, y = "Relative Abundance") +
  facet_wrap(~site_labels, scales = "free") +
  theme(axis.text.x = element_markdown(size = 7, color = "#333333"))
nt_molec

## (b) TM

# add in broader categories (separately for each taxa)
molecular_broader_tm <- microbial_grouping(molecular$tm , "phylum_class", 0.15) %>% 
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(broader))[-which(str_sort(unique(broader)) == "Other")],
                                            "Other")))

# make figure
tm_molec <- barplot(data = molecular_broader_tm,  x = "date_site_reach", y  = "relative_abundance", fill = "broader_factor") +
  scale_fill_discrete("Phylum - Class", palette = c(palette, end_color)) +
  labs(x = NULL, y = "Relative Abundance") +
  facet_wrap(~site_labels, scales = "free") +
  theme(axis.text.x = element_markdown(size = 7, color = "#333333"))
tm_molec

## (c) TAC

# add in broader categories (separately for each taxa)
molecular_broader_tac <- microbial_grouping(molecular$tac , "phylum_class", 0.075) %>% 
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(broader))[-which(str_sort(unique(broader)) == "Other")],
                                            "Other")))

# make figure
tac_molec <- barplot(data = molecular_broader_tac,  x = "date_site_reach", y  = "relative_abundance", fill = "broader_factor") +
  scale_fill_discrete("Phylum - Class", palette = c(palette, end_color)) +
  labs(x = NULL, y = "Relative Abundance") +
  facet_wrap(~site_labels, scales = "free") +
  theme(axis.text.x = element_markdown(size = 7, color = "#333333"))
tac_molec

#### (4) Putting Together Large Final Figure ####

## (a) nt
nt_fig <- plot_grid(nt_micro + theme(axis.text.x = element_markdown(size = 5),
                                     legend.key.size = unit(0.5, 'cm'),),
                    nt_molec + theme(axis.text.x = element_markdown(size = 5),
                                     legend.key.size = unit(0.5, 'cm'),), ncol = 1)
nt_fig

ggsave("./figures/tiff_files/sfig_nt_each_sample.tiff", dpi = 600,
       width=18, height=22, unit="cm")

## (b) tm
tm_fig <- plot_grid(tm_micro, tm_molec, ncol = 1)
tm_fig

ggsave("./figures/tiff_files/sfig_tm_each_sample.tiff", dpi = 600,
       width=18, height=22, unit="cm")

## (c) tac
tac_fig <- plot_grid(tac_micro, tac_molec, ncol = 1)
tac_fig

ggsave("./figures/tiff_files/sfig_tac_each_sample.tiff", dpi = 600,
       width=18, height=22, unit="cm")

# will put in asterisk manually above each bar plot as I am struggling here :)
