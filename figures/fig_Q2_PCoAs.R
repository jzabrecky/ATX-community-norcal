#### Main figure for PCoA for sample types among atx concentrations
### Jordan Zabrecky
## last updated: 06.05.2026

#### (1) Loading libraries & data ####

# source microscopy script
source("./code/5b_anatoxins_microscopy.R")

# save
microscopy_PCoA_list_eel <- PCoA_list_eel
microscopy_PCoA_list_rus <- PCoA_list_rus

# remove unneeded in environment
rm(list = setdiff(ls(), c("microscopy_PCoA_list_eel", "microscopy_PCoA_list_rus")))

# source 16s/molecular script
source("./code/5c_anatoxins_16s.R")

# save 
molecular_PCoA_list_eel <- PCoA_list_eel
molecular_PCoA_list_rus <- PCoA_list_rus

# remove unneeded in environment
rm(list = setdiff(ls(), c("microscopy_PCoA_list_eel", "microscopy_PCoA_list_rus",
                          "molecular_PCoA_list_eel", "molecular_PCoA_list_rus")))

# source functional script
source("./code/5d_anatoxins_functions.R")

# save
functional_PCoA_lists <- PCoA_list

# remove unneeded in environment
rm(list = setdiff(ls(), c("microscopy_PCoA_list_eel", "microscopy_PCoA_list_rus",
                          "molecular_PCoA_list_eel", "molecular_PCoA_list_rus",
                          "functional_PCoA_lists")))

# lastly, reload community analyses functions
source("./code/supplemental_code/S4b_community_analyses_func.R")

# set universal plot theme
theme_set(theme_bw() + theme(strip.background = element_blank(),
                             plot.title = element_text(hjust = 0.5), legend.text = element_markdown(),
                             text = element_text(size = 7), strip.text = element_text(size = 7),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
                             axis.text.x = element_text(margin = margin(t = 0.005, unit = "cm")),
                             axis.text.y = element_text(margin = margin(t = 0.005, unit = "cm")),
                             axis.title.y = element_text(margin = margin(t = 0.005, unit = "cm"))))

#### (2) Making plots - opt 1 upper & lower quantile coloring ####

## (a) Microscopy
microscopy_PCoA_plot_eel <- lapply(microscopy_PCoA_list_eel, function(x) 
  makePCoAplot(x, color = "atx_group", shape = "atx_group", stat_ellipse = FALSE))
microscopy_PCoA_plot_rus <- lapply(microscopy_PCoA_list_rus, function(x) 
  makePCoAplot(x, color = "atx_group", shape = "atx_group", stat_ellipse = FALSE))
# looking at PCoA data, first is NT then TAC
names(microscopy_PCoA_plot_eel) <- c("NT", "TAC", "TM")
names(microscopy_PCoA_plot_rus) <- c("NT", "TAC")

## (b) Molecular
molecular_PCoA_plot_eel <- lapply(molecular_PCoA_list_eel, function(x) 
  makePCoAplot(x, color = "atx_group", shape = "atx_group", stat_ellipse = FALSE))
molecular_PCoA_plot_rus <- lapply(molecular_PCoA_list_rus, function(x) 
  makePCoAplot(x, color = "atx_group", shape = "atx_group", stat_ellipse = FALSE))
# looking at PCoA data, first is NT then TAC
names(molecular_PCoA_plot_rus) <- c("NT", "TAC")

## (c) Functional
functional_PCoA_plot <- lapply(functional_PCoA_lists, function(x) 
  makePCoAplot(x, color = "atx_group", shape = "atx_group", stat_ellipse = FALSE))

## (d) plot all together

# plot all with plot grid
final <- plot_grid(microscopy_PCoA_plot_eel$NT + theme(legend.position = "none") +
                     stat_ellipse(aes(color = atx_detected), type = "t", linetype = 2, linewidth = 0.5), 
                   microscopy_PCoA_plot_eel$TM + theme(legend.position = "none") +
                     stat_ellipse(aes(color = atx_detected), type = "t", linetype = 2, linewidth = 0.5),
                   microscopy_PCoA_plot_eel$TAC + theme(legend.position = "none") +
                     stat_ellipse(aes(color = atx_detected), type = "t", linetype = 2, linewidth = 0.5),
                   microscopy_PCoA_plot_rus$NT + theme(legend.position = "none"),
                   microscopy_PCoA_plot_rus$TAC + theme(legend.position = "none"),
                   molecular_PCoA_plot_eel$NT + theme(legend.position = "none") +
                     stat_ellipse(aes(color = atx_detected), type = "t", linetype = 2, linewidth = 0.5), 
                   molecular_PCoA_plot_eel$TM + theme(legend.position = "none") +
                     stat_ellipse(aes(color = atx_detected), type = "t", linetype = 2, linewidth = 0.5),
                   molecular_PCoA_plot_eel$TAC + theme(legend.position = "none") + 
                     stat_ellipse(aes(color = atx_detected), type = "t", linetype = 2, linewidth = 0.5),
                   molecular_PCoA_plot_rus$NT + theme(legend.position = "none"),
                   molecular_PCoA_plot_rus$TAC + theme(legend.position = "none"),
                   functional_PCoA_plot$`SFE-M NT` + theme(legend.position = "none"), 
                   functional_PCoA_plot$`SFE-M TM` + theme(legend.position = "none"),
                   functional_PCoA_plot$`SFE-M TAC` + theme(legend.position = "none"),
                   functional_PCoA_plot$`RUS NT` + theme(legend.position = "none"),
                   functional_PCoA_plot$`RUS TAC` + theme(legend.position = "none"),
                   align = "hv", nrow = 3)
final

# save figure
ggsave("./figures/tiff_files/Q2_atx_pcoas_opt1.tiff", dpi = 600,
       width=17.25, height=10.75, unit="cm")

# save legend
microscopy_PCoA_plot_eel$NT + theme(legend.position = "right")
ggsave("./figures/tiff_files/Q2_atx_pcoas_opt1_legend.tiff", dpi = 600,
       width=18, height=11, unit="cm")

#### (2) Making plots - opt 2 coloring based on ATX intervals ####

## (a) Microscopy
microscopy_PCoA_plot_eel_2 <- lapply(microscopy_PCoA_list_eel, function(x) 
  makePCoAplot(x, color = "atx_group_interval", shape = "atx_group", stat_ellipse = FALSE))
microscopy_PCoA_plot_rus_2 <- lapply(microscopy_PCoA_list_rus, function(x) 
  makePCoAplot(x, color = "atx_group_interval", shape = "atx_group", stat_ellipse = FALSE))
# looking at PCoA data, first is NT then TAC
names(microscopy_PCoA_plot_eel_2) <- c("NT", "TAC", "TM")
names(microscopy_PCoA_plot_rus_2) <- c("NT", "TAC")

## (b) Molecular
molecular_PCoA_plot_eel_2 <- lapply(molecular_PCoA_list_eel, function(x) 
  makePCoAplot(x, color = "atx_group_interval", shape = "atx_group", stat_ellipse = FALSE))
molecular_PCoA_plot_rus_2 <- lapply(molecular_PCoA_list_rus, function(x) 
  makePCoAplot(x, color = "atx_group_interval", shape = "atx_group", stat_ellipse = FALSE))
# looking at PCoA data, first is NT then TAC
names(molecular_PCoA_plot_rus_2) <- c("NT", "TAC")

## (c) Functional
functional_PCoA_plot_2 <- lapply(functional_PCoA_lists, function(x) 
  makePCoAplot(x, color = "atx_group_interval", shape = "atx_group", stat_ellipse = FALSE))

## (d) plot all together

# plot all with plot grid
final2 <- plot_grid(microscopy_PCoA_plot_eel_2$NT + theme(legend.position = "none") +
                     stat_ellipse(aes(color = atx_detected), type = "t", linetype = 2, linewidth = 0.5), 
                   microscopy_PCoA_plot_eel_2$TM + theme(legend.position = "none") +
                     stat_ellipse(aes(color = atx_detected), type = "t", linetype = 2, linewidth = 0.5),
                   microscopy_PCoA_plot_eel_2$TAC + theme(legend.position = "none") +
                     stat_ellipse(aes(color = atx_detected), type = "t", linetype = 2, linewidth = 0.5),
                   microscopy_PCoA_plot_rus_2$NT + theme(legend.position = "none"),
                   microscopy_PCoA_plot_rus_2$TAC + theme(legend.position = "none"),
                   molecular_PCoA_plot_eel_2$NT + theme(legend.position = "none") +
                     stat_ellipse(aes(color = atx_detected), type = "t", linetype = 2, linewidth = 0.5), 
                   molecular_PCoA_plot_eel_2$TM + theme(legend.position = "none") +
                     stat_ellipse(aes(color = atx_detected), type = "t", linetype = 2, linewidth = 0.5),
                   molecular_PCoA_plot_eel_2$TAC + theme(legend.position = "none") + 
                     stat_ellipse(aes(color = atx_detected), type = "t", linetype = 2, linewidth = 0.5),
                   molecular_PCoA_plot_rus_2$NT + theme(legend.position = "none"),
                   molecular_PCoA_plot_rus_2$TAC + theme(legend.position = "none"),
                   functional_PCoA_plot_2$`SFE-M NT` + theme(legend.position = "none"), 
                   functional_PCoA_plot_2$`SFE-M TM` + theme(legend.position = "none"),
                   functional_PCoA_plot_2$`SFE-M TAC` + theme(legend.position = "none"),
                   functional_PCoA_plot_2$`RUS NT` + theme(legend.position = "none"),
                   functional_PCoA_plot_2$`RUS TAC` + theme(legend.position = "none"),
                   align = "hv", nrow = 3)
final2

# save figure
ggsave("./figures/tiff_files/Q2_atx_pcoas_opt2.tiff", dpi = 600,
       width=17, height=10.5, unit="cm")

# save legend
microscopy_PCoA_plot_eel_2$NT + theme(legend.position = "right",  legend.margin = margin(0, 0, 0, 0),
                                      legend.spacing.x = unit(0, "mm"),
                                      legend.spacing.y = unit(0, "mm"))
ggsave("./figures/tiff_files/Q2_atx_pcoas_opt2_legend.tiff", dpi = 600,
       width=18, height=5, unit="cm")
