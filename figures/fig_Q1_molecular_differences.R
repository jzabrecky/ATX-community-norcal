#### Main figure to show differences in morphologically-identified assemblages among rivers
### Jordan Zabrecky
## last edited: 05.07.2026

# This script creates a main figure to show differences in assemblages associated 
# with nontarget, microcoleus, and anabaena samples

# THINGS TO CONSIDER: 
# - removing Russian outlier sample
# - separating classes for certain types then the rest as phylums?
# - change NMDS axes!

#### (1) Loading libraries & data ####

# load from analysis script
source("./code/4b_amongrivers_16s.R")

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

## (b) classes bar plot

# putting together data so they share the same ~12 Phylum-Classes for the figure!
all <- rbind(data_long$nt, data_long$tm, data_long$tac)

# use function 
source("./code/supplemental_code/S4c_grouping_func.R")
all_phylumclass <- microbial_grouping(all, "phylum_class", 0.065) %>% 
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(broader))[-which(str_sort(unique(broader)) == "Other")],
                                            "Other")))

# split dataframe back up to plot separately
all_split <- split(all_phylumclass, all_phylumclass$sample_type)

## make figures
fig_b <- list()
for(i in sample_types) {
  fig_b[[i]] = barplot(all_split[[str_to_upper(i)]], x = "site", y = "relative_abundance",
                       fill = "broader_factor") +
    scale_fill_discrete(palette = c(palette[-length(palette)], end_color)) +
    labs(x = NULL, y = NULL) #+
    #theme(legend.position = "none")
}
lapply(fig_b, print)

## okay, real quick let's see a phylum version

# putting together data so they share the same ~12 Phylum-Classes for the figure!
all_phylum <- microbial_grouping(all, "phylum", 0.05) %>% 
  # add in case if NA
  mutate(broader = case_when(is.na(phylum) ~ "Other",
                                    TRUE ~ broader)) %>% 
  # then factor
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(broader))[-which(str_sort(unique(broader)) == "Other")],
                                            "Other")))

# split dataframe back up to plot separately
all_p_split <- split(all_phylum, all_phylum$sample_type)

## make figures
fig_b_ver2 <- list()
for(i in sample_types) {
  fig_b_ver2[[i]] = barplot(all_p_split[[str_to_upper(i)]], x = "site", y = "relative_abundance",
                       fill = "broader_factor") +
    scale_fill_discrete(palette = c(palette[-length(palette)], end_color)) +
    labs(x = NULL, y = NULL) #+
  #theme(legend.position = "none")
}
lapply(fig_b_ver2, print)

# just to visually compare between the two
for(i in sample_types) {
  print(plot_grid(fig_b[[i]], fig_b_ver2[[i]]))
}

# Eh maybe just put class one in supplemental!

## (c) diversity box plot

fig_c <- list()
for(i in sample_types) {
  fig_c[[i]] = ggplot(data = diversity[[i]], aes(x = site, y = shannon_diversity, fill = site)) +
    geom_boxplot() +
    scale_fill_manual(values = c("SAL" = "#62a7f8",
                                 "SFE-M" = "#416f16",
                                 "RUS" = "#bdb000")) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none")
}
lapply(fig_c, print)

#### (3) Putting Figure Together ####

# put figure together
final <-  plot_grid(fig_a$nt, fig_a$tm, fig_a$tac, fig_b$nt, fig_b$tm, fig_b$tac,
                    fig_c$nt, fig_c$tm, fig_c$tac, align = "hv")
final

# save
ggsave("./figures/tiff_files/Q1_molecular.tiff", dpi = 600,
       width=18, height=20, unit="cm")

# saving sfs stuff
plot_grid(fig_a$tac, fig_a$tm)
ggsave("./figures/SFS_figures/target_molecular.tiff", dpi = 600,
       width=17.6, height=6, unit="cm")
plot_grid(fig_c$tac, fig_c$tm)
ggsave("./figures/SFS_figures/target_molecular_diversity.tiff", dpi = 600,
       width=17.6, height=6, unit="cm")
