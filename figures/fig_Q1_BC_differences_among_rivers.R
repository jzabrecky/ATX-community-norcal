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
lapply(c("cowplot"), require, character.only = T)

# set universal plot theme (just text size for now)
theme_set(theme_bw() + theme(text = element_text(size = 10)))

# custom palette
palette <- c("#FBF6B0", "#C5BD53", "#777122", "#C2DFFF", "#5E9DE0", "#205288", 
             "#C0ED96", "#7AB048", "#3D631A", "#CBC5F6", "#8A80CF", "#61389E")

#### (2) Creating Individual Plots ####

# list of desired figures:
# a- NMDS algal
# b- NMDS bacterial
# c- NMDS functional
# d- algal taxa groups (bar plot)
# e- bacterial phylum class (bar plot)
# f- functional groups (bar plot)
# h- alpha diversity (scatterplot)
# i- venn diagram

# creating string vector to iterate through sample types
sample_types = c("nt", "tac", "tm")

## (a) NMDS algal

fig_a <- list()
for(i in sample_types) {
  fig_a[[i]] = makeNMDSplot(NMDS_list_morphological[[i]], FALSE, FALSE,
                            color = "site", shape = "site")
}
lapply(fig_a, print)


## (b) NMDS bacterial

fig_b <- list()
for(i in sample_types) {
  fig_b[[i]] <- makeNMDSplot(NMDS_list[[i]], FALSE, FALSE,
                             color = "site", shape = "site")
}
lapply(fig_b, print)

## (c) NMDS functional - TBD

## (d) algal taxa bar plot

fig_d <- list()
col_num <- list(14, 9, 10) # list for number of fill colors
names(col_num) <- sample_types
for(i in sample_types) {
 fig_d[[i]] <- barplot_broader_plots[[i]] +
   scale_fill_discrete(palette = palette)
}
lapply(fig_d, print)

# put in line with supplemental figure 

## (e) bacterial taxa bar plot

fig_e <- list()
for(i in sample_types) {
  fig_e[[i]] <- barplot_phylum_plots[[i]] + 
    scale_fill_discrete(palette = palette)
}
lapply(fig_e, print)

## (f) functional taxa TBD

## (e) alpha diversity (shannon diversity index)

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

## (f) ASV venn diagram

# figures created in other script:
# "fig_Q1_ASV_venn_diagrams.R"
# will put everything together in Inkscape, but for now will use the previous plots as 
# a place holder

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
                                 ncol = 1, rel_heights = c(1.25, 1),
                                 scale = 0.95)
}
lapply(final_figure, print)

# consider putting PERMANOVA and/or PERMDISP results on top of NMDS!
# also add in titles for title space!

# save when complete :)