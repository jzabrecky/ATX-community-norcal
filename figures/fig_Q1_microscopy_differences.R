#### Main figure to show differences in morphologically-identified assemblages among rivers
### Jordan Zabrecky
## last edited: 03.06.2026

# This script creates a main figure to show differences in morphologically-identified
# assemblages for all sample types including an (NMDS) and relative abundance barplo

# POSSIBLE TO-DO: grouping taxa better for NT and TM/TAC, 
# putting PERMANOVA & PERMDISP results

#### (1) Loading libraries & data ####

# load from analysis script
source("./code/4a_amongrivers_microscopy.R")

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
labels <- list(c("Diatoms (other than<br>*Epithemia* or *Rhopalodia*)", "*Epithemia*",
               "*Rhopalodia*", "Other N-Fixing Cyanobacteria", "Other Filamentous Cyanobacteria",
               "Unicellular Cyanobacteria", "*Cladophora*", "*Spirogyra*", "Other Green Algae",
               "*Anabaena* and<br>*Cylindrospermum*", "*Microcoleus*", "*Geitlerinema*", "Other"),
               c("Diatoms (other than *Epithemia*)", "*Epithemia*", "Other N-Fixing Cyanobacteria",
                 "Other Filamentous Cyanobacteria", "Unicellular Cyanobacteria", "*Microcoleus*", 
                 "*Geitlerinema*", "Other"),
               c("Diatoms (other than *Epithemia*)", "*Epithemia*", "Other N-Fixing Cyanobacteria",
                 "Other Filamentous Cyanobacteria", "Unicellular Cyanobacteria", "Green Algae",
                 "*Anabaena* and<br>*Cylindrospermum*", "*Geitlerinema*"))
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
   theme(axis.text.x = element_markdown(size = 7, angle = 60, vjust = 1, hjust=1))
}
lapply(fig_b, print)

#### (3) Putting Figures Together ####

## (a) figure

# will create top row and bottom row first separately, then put together for final figure
figure <- plot_grid(fig_a$nt + labs(x = NULL, y = NULL), 
                    fig_a$tm + labs(x = NULL, y = NULL),
                    fig_a$tac + labs(x = NULL, y = NULL),
                    fig_b$nt + theme(legend.position = "none") + labs(y = NULL), 
                    fig_b$tm + theme(legend.position = "none") + labs(y = NULL), 
                    fig_b$tac + theme(legend.position = "none") + labs(y = NULL), 
                    align = "hv")
figure

# save
ggsave("./figures/tiff_files/Q1_microscopy_figure.tiff", dpi = 600, 
       width=17.75, height=13, unit="cm")

## (b) legend (will add in manually in inkscape)

# create dummy figure for legend with all groups for targeted 
# (as TAC is missing AC and TM is missing M)
targetlegend_data <- data.frame(taxa = c("Diatoms<br>(other than *Epithemia* )", "*Epithemia*", "Other N-Fixing Cyanobacteria",
                                         "Other Filamentous Cyanobacteria", "Unicellular  Cyanobacteria", "Green Algae",
                                         "*Anabaena* and<br>*Cylindrospermum*", "*Microcoleus*", "*Geitlerinema*", "Other"),
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
                             theme(legend.key.size = unit(0.5, 'cm'),
                                   legend.key.height = unit(0.5, 'cm'),
                                   legend.key.width = unit(0.5, 'cm')),
                           targetlegend_plot + theme(legend.position = "bottom") +
                             theme(legend.key.size = unit(0.5, 'cm'),
                                   legend.key.height = unit(0.5, 'cm'),
                                   legend.key.width = unit(0.5, 'cm')),
                           ncol = 1)
figure_legend

# save
ggsave("./figures/tiff_files/Q1_microscopy_legend.tiff", dpi = 600,
       width=18, height=13, unit="cm")
