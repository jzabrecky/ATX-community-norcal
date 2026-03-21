#### Main figure to show differences in morphologically-identified assemblages among rivers
### Jordan Zabrecky
## last edited: 03.20.2026

# This script creates a main figure to show differences in morphologically-identified
# assemblages for all sample types including an (NMDS) and relative abundance barplots

# One figure is made for NT samples showing (a) relative abundance plots and (b) NMDS
# Another figure is made for TM and TAC samples showing (a) relative abundance plots WITH
# target taxa included, (b) relative abundance WITHOUT target taxa included (and green
# algae in the case of Anabaena), and (c) NMDS plots

#### (1) Loading libraries & data ####

# load from analysis script
source("./code/4a_amongrivers_microscopy.R")

# need to also load data with target taxa included
tm_w_m <- read.csv("./data/morphological/tm_algalonly.csv") %>% 
  pivot_longer(5:ncol(.), names_to = "taxa", values_to = "percent") %>% 
  filter(percent != 0 & year(field_date) == 2022)
tac_w_ac_g <- read.csv("./data/morphological/tac_algalonly.csv") %>% 
  pivot_longer(5:ncol(.), names_to = "taxa", values_to = "percent") %>% 
  filter(percent != 0 & year(field_date) == 2022)

# load additional libraries
lapply(c("cowplot", "ggtext"), require, character.only = T)

# set universal plot theme
theme_set(theme_bw() + theme(strip.background = element_blank(),
                             plot.title = element_text(hjust = 0.5), legend.text = element_markdown(),
                             text = element_text(size = 8), strip.text = element_text(size = 8),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),))

# custom palette
palette <- c("#FBF6B0", "#C5BD53", "#777122", "#C2DFFF", "#5E9DE0", "#205288", 
             "#C0ED96", "#7AB048", "#3D631A", "#CBC5F6", "#8A80CF", "#61389E")

# color for other or unknown
end_color <- "lightgray"

#### (2) Creating Individual Plots ####

# creating string vector to iterate through sample types
sample_types = c("nt", "tac", "tm")

## (a) NMDS

fig_a <- list()
for(i in sample_types) {
  fig_a[[i]] = makeNMDSplot(NMDS_list[[i]], FALSE, FALSE,
                            color = "site", shape = "site") +
    theme(legend.position = "none")
}
lapply(fig_a, print)

# add site labels (to not be abbreviations)

## (b) algal taxa bar plot

# create dataframes that are even broader to match NT and T as much as possible:
# overall aim:
# yellow 1: diatoms exluding Epithemia (and Ropalodia for NT)
# yellow 2: Epithemia
# yellow 3: Rhopalodia (NT only)
# blue 1: Diazotrophic Filamentous Cyanobacteria
# blue 2: Other Filamentous Cyanobacteria
# blue 3: Unicellular Cyanobacteria
# green 1: Cladophora (NT only)
# green 2: Spirogyra (NT only)
# green 3: other green algae (NT) green algae (T)
# purple 1: anabaena
# purple 2: microcoleus
# grey: other
figure_data <- lapply(data_longer, function(x) {
  y = x %>% 
    # move Nostoc to other N-fixing cyanobacteria
    mutate(figure_groups = case_when(broader == "Nostoc" ~ "Other N-fixing Cyanobacteria",
                                     taxa == "geitlerinema" ~ "Geitlerinema",
                                     broader == "Unknown" ~ "Other",
                                     broader == "Misc. Other" ~ "Other",
                                     taxa == "rhopalodia" ~ "Rhopalodia",
                                     taxa == "epithemia" ~ "Epithemia",
                                     TRUE ~ broader),
           figure_groups_factored = factor(figure_groups, 
                                           levels = c("Diatoms Other than Epithemia or Rhopalodia", 
                                                      "Diatoms Other than Epithemia", "Epithemia",
                                                      "Rhopalodia", "Other N-fixing Cyanobacteria",
                                                      "Other Filamentous Cyanobacteria",
                                                      "Unicellular Cyanobacteria", "Cladophora", 
                                                      "Spirogyra", "Other Green Algae", "Green Algae",
                                                      "Anabaena or Cylindrospermum", "Microcoleus",
                                                      "Geitlerinema", "Unknown", "Other"))) 
  # return new dataframe
  return(y)
})

# empty list
fig_b <- list()

# set attributes based on sample type
colors <- list(c(1:12), c(1:2, 4:6, 11:12), c(1:2, 4:6, 8, 10, 12)) # list for number of fill colors to end before "other" category
labels <- list(c("Diatoms (other than *Epithemia* or *Rhopalodia*)", "*Epithemia*",
               "*Rhopalodia*", "Other N-Fixing Cyanobacteria", "Other Filamentous Cyanobacteria",
               "Unicellular Cyanobacteria", "*Cladophora*", "*Spirogyra*", "Other Green Algae",
               "*Anabaena* and *Cylindrospermum*", "*Microcoleus*", "*Geitlerinema*", "Other"),
               c("Diatoms (other than *Epithemia*)", "*Epithemia*", "Other N-Fixing Cyanobacteria",
                 "Other Filamentous Cyanobacteria", "Unicellular Cyanobacteria", "*Microcoleus*", 
                 "*Geitlerinema*", "Other"),
               c("Diatoms (other than *Epithemia*)", "*Epithemia*", "Other N-Fixing Cyanobacteria",
                 "Other Filamentous Cyanobacteria", "Unicellular Cyanobacteria", "Green Algae",
                 "*Anabaena* and*Cylindrospermum*", "*Geitlerinema*"))
site_labels <- list(c("Russian<br>River", "Salmon<br>River", "South Fork<br>Eel River"),
                    c("Russian<br>River", "Salmon<br>River", "South Fork<br>Eel River"),
                    c("Salmon<br>River", "South Fork<br>Eel River"))
names(colors) <- sample_types
names(labels) <- sample_types
names(site_labels) <- sample_types

# make plots
for(i in sample_types) {
 fig_b[[i]] <- barplot(data = figure_data[[i]],
                       x = "site", y  = "percent", fill = "figure_groups_factored") +
   scale_fill_discrete("Taxa Group", palette = c(palette[colors[[i]]], end_color), 
                       labels = labels[[i]]) +
   labs(x = NULL, y = "Relative Abundance") +
   scale_x_discrete(labels = site_labels[[i]]) +
   theme(axis.text.x = element_markdown(size = 7))
}
lapply(fig_b, print)

## (c) lastly need target taxa with that taxa incorporated

# add factoring to target w/ target taxa dataframes
target_w_target_taxa <- lapply(list(tm_w_m, tac_w_ac_g), function(x) {
  
  # run broader groups function from previous script to this data
  y = target_broader(x)
  
  # add groupings to match barplots
  y = y %>% 
    # move Nostoc to other N-fixing cyanobacteria
    mutate(figure_groups = case_when(broader == "Nostoc" ~ "Other N-fixing Cyanobacteria",
                                     taxa == "geitlerinema" ~ "Geitlerinema",
                                     broader == "Unknown" ~ "Other",
                                     broader == "Misc. Other" ~ "Other",
                                     taxa == "rhopalodia" ~ "Rhopalodia",
                                     taxa == "epithemia" ~ "Epithemia",
                                     TRUE ~ broader),
           figure_groups_factored = factor(figure_groups, 
                                           levels = c("Diatoms Other than Epithemia or Rhopalodia", 
                                                      "Diatoms Other than Epithemia", "Epithemia",
                                                      "Rhopalodia", "Other N-fixing Cyanobacteria",
                                                      "Other Filamentous Cyanobacteria",
                                                      "Unicellular Cyanobacteria", "Cladophora", 
                                                      "Spirogyra", "Other Green Algae", "Green Algae",
                                                      "Anabaena or Cylindrospermum", "Microcoleus",
                                                      "Geitlerinema", "Unknown", "Other"))) 
  # return new dataframe
  return(y)
})
names(target_w_target_taxa) <- c("tm", "tac")

# make empty list
fig_c <- list()

# tm
fig_c[["tm"]] <- barplot(data = target_w_target_taxa$tm,  x = "site", y  = "percent", fill = "figure_groups_factored") +
  scale_fill_discrete("Taxa Group", palette = c(palette[c(1:2, 4:6, 8, 10:12)], end_color), 
                      labels = c("Diatoms (other than *Epithemia*)", "*Epithemia*", "Other N-Fixing Cyanobacteria",
                                 "Other Filamentous Cyanobacteria", "Unicellular Cyanobacteria", "Green Algae",
                                 "*Microcoleus*", "*Anabaena* and*Cylindrospermum*", "*Geitlerinema*")) +
  labs(x = NULL, y = "Relative Abundance") +
  scale_x_discrete(labels = c("Salmon<br>River", "South Fork<br>Eel River")) +
  theme(axis.text.x = element_markdown(size = 7))

# tac
fig_c[["tac"]] <- barplot(data = target_w_target_taxa$tac,  x = "site", y  = "percent", fill = "figure_groups_factored") +
  scale_fill_discrete("Taxa Group", palette = c(palette[c(1:2, 4:6, 8, 10:12)], end_color), 
                      labels = c("Diatoms (other than *Epithemia*)", "*Epithemia*", "Other N-Fixing Cyanobacteria",
                                 "Other Filamentous Cyanobacteria", "Unicellular Cyanobacteria", "Green Algae",
                                 "*Microcoleus*", "*Anabaena* and*Cylindrospermum*", "*Geitlerinema*", "Other")) +
  labs(x = NULL, y = "Relative Abundance") +
  scale_x_discrete(labels = c("Russian<br>River", "Salmon<br>River", "South Fork<br>Eel River")) +
  theme(axis.text.x = element_markdown(size = 7))

# show plots
lapply(fig_c, print)

#### (3) Putting Figures Together ####

## (a) NT figure

# put figure together
nt_figure <- plot_grid(fig_b$nt + theme(legend.position = "none") + labs(y = "Relative Abundance"),
                       fig_a$nt + labs(x = "NMDS1", y = "NMDS2"), align = "hv", ncol = 1)
nt_figure

# save
ggsave("./figures/tiff_files/Q1_nt_microscopy_figure.tiff", dpi = 600, 
       width=8.5, height=12, unit="cm")

## (b) T figure

# target taxa
t_figure <- plot_grid(fig_c$tm + theme(legend.position = "none") + labs(y = NULL),
                      fig_c$tac + theme(legend.position = "none") + labs(y = NULL),
                      fig_b$tm + theme(legend.position = "none") + labs(y = NULL),
                      fig_b$tac+ theme(legend.position = "none") + labs(x = NULL, y = NULL),
                      fig_a$tm + labs(y = NULL, x = NULL), 
                      fig_a$tac + labs(y = NULL, x = NULL),
                      align = "hv", ncol = 2)
t_figure

# save
ggsave("./figures/tiff_files/Q1_t_microscopy_figure.tiff", dpi = 600, 
       width=17.6, height=18, unit="cm")

## (c) legends (will add in manually in inkscape)

# create dummy figure for legend with all groups for targeted 
# (as TAC is missing AC and TM is missing M)
targetlegend_data <- data.frame(taxa = c("Diatoms(other than *Epithemia* )", "*Epithemia*", "Other N-Fixing Cyanobacteria",
                                         "Other Filamentous Cyanobacteria", "Unicellular  Cyanobacteria", "Green Algae",
                                         "*Anabaena* and*Cylindrospermum*", "*Microcoleus*", "*Geitlerinema*", "Other"),
                       abundance = rep(100 / 10, 10),
                       site = rep('river', 10))
targetlegend_plot <- barplot(data = targetlegend_data, 
                             x = "site", y  = "abundance", fill = "taxa") +
  scale_fill_discrete("Taxa Group", palette = c(palette[c(1:2, 4:6, 8, 10:12)], end_color),
                      labels = targetlegend_data$taxa) +
  labs(x = NULL, y = "Relative Abundance")

# visually confirm that these match what is in the ta & tm plots
targetlegend_plot
fig_b$tac
fig_b$tm

# putting legend stuff together
figure_legend <- plot_grid(fig_a$nt + scale_shape_discrete(labels = c("Russian  River", "Salmon  River", "South  Fork  Eel  River")) +
                             scale_color_discrete(labels = c("Russian  River", "Salmon  River", "South  Fork  Eel  River"),
                                                  palette = c("#bdb000", "#62a7f8", "#416f16")) +
                             theme(legend.position = "bottom", legend.title = element_blank()),
                           fig_b$nt + theme(legend.position = "bottom", legend.title = element_blank()) +
                             theme(legend.key.size = unit(0.3, 'cm'),
                                   legend.key.height = unit(0.3, 'cm'),
                                   legend.key.width = unit(0.3, 'cm')),
                           targetlegend_plot + theme(legend.position = "bottom") +
                             theme(legend.key.size = unit(0.3, 'cm'),
                                   legend.key.height = unit(0.3, 'cm'),
                                   legend.key.width = unit(0.3, 'cm')),
                           ncol = 1)
figure_legend

# save
ggsave("./figures/tiff_files/Q1_microscopy_legend.tiff", dpi = 600,
       width=18, height=13, unit="cm")
