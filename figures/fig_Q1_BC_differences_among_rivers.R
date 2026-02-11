#### Main figure to show differences in assemblages associated with target taxa among rivers
### Jordan Zabrecky
## last edited: 01.15.2026

# This script creates a main figure to show differences in assemblages associated 
# with Microcoleus and Anabaena among rivers sampled

#### (1) Loading libraries & data ####

# load from analysis script
source("./code/4a_amongrivers_microscopy.R")

# rename data to not get overrun
data_morphological <- data
data_longer_morphological <- data_longer
NMDS_list_morphological <- NMDS_list
NMDS_plots_morphological <- NMDS_plots

# load from second analysis script
source("./code/4b_amongrivers_16s.R")

# TBD- load from PICRUSt analysis

# load additional libraries
lapply(c("viridisLite", "cowplot"), require, character.only = T)

# set universal plot theme (just text size for now)
theme_set(theme_bw() + theme(text = element_text(size = 10)))

#### (2) Creating Individual Plots ####

# list of desired figures:
# a- NMDS algal
# b- algal taxa groups (bar plot)
# c- NMDS bacterial
# d- alpha diversity (scatterplot)
# e- venn diagram
# f- functional groups (bar plot)

# creating string vector to iterate through sample types
sample_types = c("nt", "tac", "tm")

## (a) NMDS algal

fig_a <- list()
for(i in sample_types) {
  fig_a[[i]] = makeNMDSplot(NMDS_list_morphological[[i]], FALSE, FALSE,
                            color = "site", shape = "site")
}
lapply(fig_a, print)

## (b) algal taxa bar plot

fig_b <- list()
col_num <- list(14, 9, 10) # list for number of fill colors
names(col_num) <- sample_types
for(i in sample_types) {
 fig_b[[i]] <- barplot_broader_plots[[i]] +
   scale_fill_discrete(palette = viridis(col_num[[i]], direction = -1))
}
lapply(fig_b, print)

# TO-DO: maybe want to reorder or rechoose groups / merge some together
# unknown needs to be removed for TM (no taxa contain unknown)

## (c) NMDS bacterial

fig_c <- list()
for(i in sample_types) {
  fig_c[[i]] <- makeNMDSplot(NMDS_list[[i]], FALSE, FALSE,
                             color = "site", shape = "site")
}
lapply(fig_c, print)

## (d) alpha diversity (shannon diversity index)

fig_d <- list()
d_colors <- list(c("#ebdf38", "#62a7f8", "#416f16"),
                 c("#ebdf38", "#62a7f8", "#416f16"),
                 c("#62a7f8", "#416f16"))
names(d_colors) <- sample_types
for(i in sample_types){
  fig_d[[i]] <- ggplot(data = diversity[[i]], aes(x = site, y = shannon_diversity, 
                                                  fill = site)) +
    geom_boxplot(alpha = 0.5) +
    geom_jitter(aes(color = site), alpha = 0.9, size = 2) +
    scale_fill_discrete(palette = d_colors[[i]]) +
    scale_color_discrete(palette = d_colors[[i]])
}
lapply(fig_d, print)

## (e) ASV venn diagram

# figures created in other script:
# "fig_Q1_ASV_venn_diagrams.R"
# will put everything together in Inkscape, but for now will use the previous plots as 
# a place holder

## (f) functional group bar plot

# TO-DO: to be created after using PiCRUST

#### (3) Putting Figures Together ####

# will create top row and bottom row first separately, then put together for final figure
final_figure <- list()
for(i in sample_types) {
  top_row = plot_grid(fig_a[[i]] + theme(legend.position = "none"),
                      fig_c[[i]] + theme(legend.position = "none"), 
                      # will add legends for the above separately, as they are small 
                      # and can fit into a different space!
                       fig_d[[i]] + theme(legend.position = "none"), 
                      nrow = 1, rel_widths = c(1.5, 1.5, 1))
  bottom_row = plot_grid(fig_b[[i]] + theme(legend.position = "right"), 
                        fig_b[[i]] + theme(legend.position = "right"), 
                        fig_d[[i]] + theme(legend.position = "none"), 
                        nrow = 1, rel_widths = c(1.5, 1.5, 1))
  final_figure[[i]] <- plot_grid(top_row, bottom_row, 
                                 ncol = 1, rel_heights = c(1.25, 1))
}
lapply(final_figure, print)

# consider putting PERMANOVA and/or PERMDISP results on top of NMDS!
# also add in titles for title space!

# save when complete :)