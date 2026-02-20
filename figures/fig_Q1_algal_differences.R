#### Main figure to show differences in morphologically-identified assemblages among rivers
### Jordan Zabrecky
## last edited: 01.15.2026

# This script creates a main figure to show differences in morphologically-identified
# assemblages for all sample types including an (NMDS) and relative abundance barplot

#### (1) Loading libraries & data ####

# load from analysis script
source("./code/4a_amongrivers_microscopy.R")

# load additional libraries
lapply(c("cowplot", "ggtext"), require, character.only = T)

# set universal plot theme
theme_set(theme_bw() + theme(strip.background = element_blank(),
                             axis.text.x = element_markdown(size = 10, angle = 60, vjust = 1, hjust=1),
                             plot.title = element_text(hjust = 0.5), legend.text = element_markdown(),
                             text = element_text(size = 10), strip.text = element_text(size = 10)))

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

# empty list
fig_b <- list()

# set attributes based on sample type
colors <- list(c(1:11), c(2:5, 7:10), c(1:4, 6:10)) # list for number of fill colors to end before "other" category
# note: for NT, need to refactor as desired label order is not alphabetical
labels <- list(c("*Anabaena* or *Cylindrospermum*", "*Cladophora*",
                 "Diatoms (other than<br>*Epithemia* or *Rophalodia*)",
                 "*Epithemia* or *Rhopalodia*", "*Microcoleus*", "*Nostoc*", "Other Filamentous Cyanobacteria",
                 "Other Green Algae", "Other N-fixing Cyanobacteria", "*Spirogyra*",
                 "Unicellular Cyanobacteria", "dummy test"),
               c("Diatoms (other than *Epithemia*)", "*Epithemia*",
                 "*Geitlerinema*", "*Microcoleus*", "*Nostoc*", "Other Filamentous Cyanobacteria",
                 "Other N-fixing Cyanobacteria", "Unicellular Cyanobacteria", "Unknown"),
               c("*Anabaena* or *Cylindrospermum*", "Diatoms (other than *Epithemia*)",
                 "*Epithemia*", "*Geitlerinema*", "Green Algae", "*Nostoc*", "Other Filamentous Cyanobacteria",
                 "Other N-fixing Cyanobacteria", "Unicellular Cyanobacteria"))
site_labels <- list(c("Russian River", "Salmon River", "South Fork<br>Eel River"),
                    c("Russian River", "Salmon River", "South Fork<br>Eel River"),
                    c("Salmon River", "South Fork<br>Eel River"))
names(colors) <- sample_types
names(labels) <- sample_types
names(site_labels) <- sample_types

# make TM and TAC plots (need to refactor for NT)
for(i in c("tm", "tac")) {
 fig_b[[i]] <- barplot_broader_plots[[i]] +
   scale_fill_discrete("Taxa Group", palette = c(palette[colors[[i]]], end_color),
                       labels = labels[[i]]) +
   labs(x = NULL, y = "Relative Abundance") +
   scale_x_discrete(labels = site_labels[[i]])
}
# change x and Y labels, customize x axis labels

# NT plot- use factor to get desired taxa group order
fig_b$nt <- data_longer$nt %>% 
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(data_longer$nt$broader))[-which(str_sort(unique(data_longer$nt$broader)) == "Misc. Other")],
                                            "Misc. Other"))) %>% 
  barplot(x = "site", y  = "percent", fill = "broader_factor") +
  scale_fill_discrete("Taxa Group", palette = c(palette[initial_num["nt"]:col_num["nt"]], end_color),
                      labels = labels["nt"]) +
  labs(x = NULL, y = "Relative Abundance") +
  scale_x_discrete(labels = site_labels[[i]])

test <- data_longer$nt %>% 
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(data_longer$nt$broader))[-which(str_sort(unique(data_longer$nt$broader)) == "Misc. Other")],
                                            "Misc. Other")))
barplot(test, x = "site", y  = "percent", fill = "broader_factor") +
  scale_fill_discrete("Taxa Group", palette = c(palette[1:col_num["nt"]], end_color),
                      labels = labels["nt"]) +
  labs(x = NULL, y = "Relative Abundance") +
  scale_x_discrete(labels = site_labels[[i]])

lapply(fig_b, print)
# put in line with supplemental figure
# why text markdown not working

#### (3) Putting Figures Together ####

# will create top row and bottom row first separately, then put together for final figure
figure <- cowplot::plot_grid(fig_a$nt, fig_a$tm, fig_a$tac,
                             fig_b$nt, fig_b$tm, fig_b$tac)
figure
# consider putting PERMANOVA and/or PERMDISP results on top of NMDS!
# also add in titles for title space!

# maybe consider legends on the bottom?
# also permanova and permdisp results?

# save when complete :)