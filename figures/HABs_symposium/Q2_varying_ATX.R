####

# source supplemental figure on anatoxin concentration
source("./figures/sfig_atx_concentrations.R")

# customize plot formatting


set_theme(theme_bw() + theme(legend.position = "none",
                             axis.text.x = element_text(size = 22), axis.text.y = element_text(size = 22),
                             axis.title = element_text(size = 24),  panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(), 
                             panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
                             axis.ticks = element_line(linewidth = 1.5), axis.ticks.length=unit(.25, "cm")))

# ATX concentrations plot
plot[["layers"]]$geom_jitter$aes_params$size <- 6
plot + theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=2.5),
             axis.text.x = element_markdown(size = 24))

# save!
ggsave("./figures/HABs_symposium/atx.conc.png", dpi = 500,
       width=10, height=5, unit="in")

# get legend for sample type
legend + theme()

rm(list = ls())

# PCoA Plots

source("./figures/fig_Q2_PCoAs.R")

microscopy_PCoA_plot_eel_2$TAC[["layers"]]$geom_point$aes_params$size <- 6
microscopy_PCoA_plot_eel_2$TM[["layers"]]$geom_point$aes_params$size <- 6
molecular_PCoA_plot_eel_2$TAC[["layers"]]$geom_point$aes_params$size <- 6
molecular_PCoA_plot_eel_2$TM[["layers"]]$geom_point$aes_params$size <- 6
microscopy_PCoA_plot_eel_2$TAC[["layers"]]$geom_point$aes_params$stroke <- 2
microscopy_PCoA_plot_eel_2$TM[["layers"]]$geom_point$aes_params$stroke <- 2
molecular_PCoA_plot_eel_2$TAC[["layers"]]$geom_point$aes_params$stroke <- 2
molecular_PCoA_plot_eel_2$TM[["layers"]]$geom_point$aes_params$stroke <- 2


# set universal theme for text size like other script
set_theme(theme_bw() + theme(legend.position = "none",
                             axis.text.x = element_text(size = 22), axis.text.y = element_text(size = 22),
                             axis.title = element_text(size = 24),  panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(), 
                             panel.border = element_rect(colour = "black", fill=NA, linewidth=.5),
                             axis.ticks = element_line(linewidth = 1.5), axis.ticks.length=unit(.25, "cm")))

## microscopy plot
microscopy_pcoa <- plot_grid(microscopy_PCoA_plot_eel_2$TM +
                               stat_ellipse(aes(color = atx_detected), type = "t", 
                                            linetype = 2, linewidth = 1.5) +
                               scale_shape_manual(values = c("low" = 19,
                                                             "high" = 19,
                                                             "none" = 4)) +
                               theme(legend.position = "none"),
                             microscopy_PCoA_plot_eel_2$TAC +
                               stat_ellipse(aes(color = atx_detected), type = "t", 
                                            linetype = 2, linewidth = 1.5) +
                               scale_shape_manual(values = c("low" = 19,
                                                             "high" = 19,
                                                             "none" = 4)) +
                               theme(legend.position = "none"), 
                             nrow = 1)
microscopy_pcoa

# save plot!
ggsave("./figures/HABs_symposium/algal_pcoa_atx.png", dpi = 500,
       width=11.3, height=5, unit="in")

molecular_pcoa <- plot_grid(molecular_PCoA_plot_eel_2$TM + 
                              stat_ellipse(aes(color = atx_detected), type = "t", 
                                           linetype = 2, linewidth = 1.5) +
                              scale_shape_manual(values = c("low" = 19,
                                                            "high" = 19,
                                                            "none" = 4)) +
                              theme(legend.position = "none"),
                            molecular_PCoA_plot_eel_2$TAC + 
                              stat_ellipse(aes(color = atx_detected), type = "t", 
                                           linetype = 2, linewidth = 1.5) +
                              scale_shape_manual(values = c("low" = 19,
                                                            "high" = 19,
                                                            "none" = 4)) +
                              theme(legend.position = "none"),
                            nrow = 1)
molecular_pcoa

# save plot!
ggsave("./figures/HABs_symposium/bacterial_pcoa_atx.png", dpi = 500,
       width=11.3, height=5, unit="in")

# legend
molecular_PCoA_plot_eel_2$TM +
  scale_shape_manual(values = c("low" = 19,
                                "high" = 19,
                                "none" = 4)) +
  theme(legend.position = "bottom")

# svae plot!
ggsave("./figures/HABs_symposium/atx_legend.png", dpi = 500,
       width = 11.3, height = 5, unit = "in")
