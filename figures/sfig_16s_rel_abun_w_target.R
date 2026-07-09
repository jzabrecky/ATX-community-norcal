#### Supplementary figure showing relative abundance of 16s reads with target taxa included
### Jordan Zabrecky
## last edited: 07.03.2026

# This script creates a supplemental figure showing the relative abundance of 
# the target molecularly-identified assemblages including the respective
# target taxa

#### (1) Loading libraries & data ####

# source barplot function & broader microbial grouping function
source("./code/supplemental_code/S4b_community_analyses_func.R")
source("./code/supplemental_code/S4c_grouping_func.R")

# need to also load data with target taxa included
molecular_data <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_FINAL.csv") %>% 
  dplyr::rename(relative_abundance = picrust2_relative_abundance)

# assign smaller abundance phylums to other 
molecular_data_grouped <- microbial_grouping(molecular_data, "phylum", 0.041) %>% 
  # have some reads that are just under bacteria
  mutate(broader = case_when(is.na(broader) ~ "Other",
                             broader == "Abditibacteriota" ~ "Other", 
                             TRUE ~ broader)) %>% 
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(broader))[-which(str_sort(unique(broader)) == "Other")],
                                            "Other")),
         site_label = case_when(site == "SFE-M" ~ "South Fork<br>Eel River",
                                site == "SAL" ~ "Salmon<br>River",
                                site == "RUS" ~ "Russian<br>River"))

# want to match it to figure 5
unique(molecular_data_grouped$broader)

# the extra category we have is Abditibacteriota.... reassign to other (doing this up above
# so it's included in the factor)

# separate out for each target taxa
tm <- molecular_data_grouped %>% filter(sample_type == "TM")
tac <- molecular_data_grouped %>% filter(sample_type == "TAC")

# load additional libraries
lapply(c("cowplot", "ggtext"), require, character.only = T)

#### (2) Making Figure ####

# set universal plot theme
theme_set(theme_bw() + theme(strip.background = element_blank(),
                             plot.title = element_text(hjust = 0.5), legend.text = element_markdown(),
                             text = element_text(size = 8), strip.text = element_text(size = 8),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# custom palette
palette <- c("#FBF6B0", "#C5BD53", "#777122", "#C2DFFF", "#5E9DE0", "#205288", 
             "#C0ED96", "#7AB048", "#3D631A", "#CBC5F6", "#8A80CF", "#61389E")

# color for other or unknown
end_color <- "lightgray"

# TM barplot
tm_barplot <- barplot(data = tm,  x = "site_label", y  = "relative_abundance", fill = "broader_factor") +
  scale_fill_discrete("Taxa Group", palette = c(palette[-length(palette)], end_color)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_markdown(size = 7, color = "#333333"))
tm_barplot

# TAC barplot
tac_barplot <- barplot(data = tac,  x = "site_label", y  = "relative_abundance", fill = "broader_factor") +
  scale_fill_discrete("Taxa Group", palette = c(palette[-length(palette)], end_color)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_markdown(size = 7, color = "#333333"))
tac_barplot

# put together into one plot
figure <- plot_grid(tm_barplot + theme(legend.position = "none"), 
                    tac_barplot + theme(legend.position = "none"))
figure

# save!
ggsave("./figures/tiff_files/sfig_16s_rel_abun_w_target.tiff", dpi = 500,
       width=14.25, height=11, unit="cm")

# for legend on the right
tac_barplot

# save!
ggsave("./figures/tiff_files/sfig_16s_rel_abun_w_target_legend.tiff", dpi = 500,
       width=15.5, height=11, unit="cm")
