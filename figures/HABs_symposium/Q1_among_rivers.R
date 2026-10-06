#####

# Format: PCoA Microcoleus (algal, bacterial)

# souce algal code
source("./figures/fig_Q1_microscopy_differences.R")

# cannot change ellipse- it just adds a new one so will have to 
# manually manipulate it in the object
# make preferred state
dummy_plot <- ggplot() + 
  stat_ellipse(aes(color = site), linewidth = 1.5, linetype = 2)
new_ellipse <- dummy_plot[["layers"]]$stat_ellipse

# replace stat ellipse with preferred aesthetics
fig_a$tm[["layers"]]$stat_ellipse <- new_ellipse
fig_a$tac[["layers"]]$stat_ellipse <- new_ellipse

# save figures that I care about
tm_microscopy_pcoa <- fig_a$tm
tac_microscopy_pcoa <- fig_a$tac

# maybe also bar plots TBD
tm_microscopy_bar <- fig_b$tm
tac_microscopy_bar <- fig_b$tac
# nt as well in case we want to use for legend
nt_microscopy_bar <- fig_b$nt

# remove everything else
rm(list=setdiff(ls(), c("tm_microscopy_pcoa", "tac_microscopy_pcoa",
                "tm_microscopy_bar", "tac_microscopy_bar", "nt_microscopy_bar", "new_ellipse")))

# source bacterial code
source("./figures/fig_Q1_molecular_differences.R")

# replace stat ellipse with preferred aesthetics
fig_a$tm[["layers"]]$stat_ellipse <- new_ellipse
fig_a$tac[["layers"]]$stat_ellipse <- new_ellipse

tm_molecular_pcoa <- fig_a$tm
tac_molecular_pcoa <- fig_a$tac

tm_molecular_diversity <- fig_c$tm
tac_molecular_diversity <- fig_c$tac

# remove everything else
rm(list=setdiff(ls(), c("tm_microscopy_pcoa", "tac_microscopy_pcoa",
                        "tm_microscopy_bar", "tac_microscopy_bar", "nt_microscopy_bar",
                        "tm_molecular_pcoa", "tac_molecular_pcoa",
                        "tm_molecular_diversity", "tac_molecular_diversity")))

#### putting together figure for poster

set_theme(theme_bw() + theme(legend.position = "none",
                axis.text.x = element_text(size = 22), axis.text.y = element_text(size = 22),
                axis.title = element_text(size = 24),  panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), 
                panel.border = element_rect(colour = "black", fill=NA, linewidth=1.5),
                axis.ticks = element_line(linewidth = 1.5), axis.ticks.length=unit(.25, "cm")))
tm_microscopy_pcoa

## Algal PCoA
algal_pcoa <- plot_grid(tm_microscopy_pcoa + coord_flip(clip = "off") +
                          geom_point(aes(color = site, shape = month), size = 6) +
                          stat_ellipse(aes(color = site), linewidth = 1.5, linetype = 2),
                        tac_microscopy_pcoa + coord_flip(clip = "off") +
                          geom_point(aes(color = site, shape = month), size = 6) +
                          stat_ellipse(aes(color = site), linewidth = 1.5, linetype = 2),
                        ncol = 2, align = "hv")
algal_pcoa

# save!
ggsave("./figures/HABs_symposium/algal_pcoa_rivers.png", dpi = 500,
       width=11.3, height=5, unit="in")

# get legend!
legend <- tac_microscopy_pcoa + coord_flip(clip = "off") +
  geom_point(aes(color = site, shape = month), size = 6) +
  stat_ellipse(aes(color = site)) +
  theme(legend.position = "bottom") +
  guides(
    color = guide_legend(override.aes = list(size = 5)) # Much bigger points in legend
  ) +
  theme(
    legend.key.size = unit(1.5, "cm"), 
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 20) 
  )
legend

# save legend!
ggsave("./figures/HABs_symposium/pcoa_rivers_legend_yellow.png",  dpi = 500,
       width=11.3, height=5, unit="in")


## Bacterial PCoA
bacteria_pcoa <- plot_grid(tm_molecular_pcoa + coord_flip(clip = "off") +
                             geom_point(aes(color = site, shape = month), size = 6) +
                             stat_ellipse(aes(color = site), linewidth = 1.5, linetype = 2),
                           tac_molecular_pcoa + coord_flip(clip = "off") +
                             geom_point(aes(color = site, shape = month), size = 6) +
                             stat_ellipse(aes(color = site), linewidth = 1.5, linetype = 2),
                           ncol = 2, align = "hv")
bacteria_pcoa

# save!
ggsave("./figures/HABs_symposium/bacteria_pcoa_rivers.png", dpi = 500,
       width=11.3, height=5, unit="in")

## Diversity metric

# change size of jitter
tac_molecular_diversity[["layers"]]$geom_jitter$aes_params$size <- 6
tm_molecular_diversity[["layers"]]$geom_jitter$aes_params$size <- 6
tac_molecular_diversity[["layers"]]$geom_jitter$aes_params$alpha <- 0.9
tm_molecular_diversity[["layers"]]$geom_jitter$aes_params$alpha <- 0.9
tac_molecular_diversity[["layers"]]$stat_summary$aes_params["linewidth"] <- 1
tm_molecular_diversity[["layers"]]$stat_summary$aes_params["linewidth"] <- 1
tac_molecular_diversity[["layers"]]$stat_summary$aes_params$color <- "gray"
tm_molecular_diversity[["layers"]]$stat_summary$aes_params$color <- "gray"

diversity <- plot_grid(tm_molecular_diversity + 
                         labs(y = "Shannon Diversity") +
                         theme(axis.text.x = element_markdown(size = 23)) + 
                         theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 3.25)),
                      tac_molecular_diversity + 
                        labs(y = "Shannon Diversity") +
                         theme(axis.text.x = element_markdown(size = 23)) + 
                        theme(panel.border = element_rect(colour = "black", fill=NA, linewidth= 3.25)),
                      ncol = 2, align = "hv", scale = 0.98)
diversity

# save!
ggsave("./figures/HABs_symposium/diversity_rivers.png", dpi = 500,
       width=11, height=5, unit="in")

# barplots
barplots <- plot_grid(tm_microscopy_bar + theme(axis.text.x = element_markdown(size = 22)) +
                        theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 3.25)),
                      tac_microscopy_bar + theme(axis.text.x = element_markdown(size = 22))  + 
                        theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 3.25)), ncol = 2, align = "hv",
                      scale = 1)
barplots


# save!
ggsave("./figures/HABs_symposium/barplots_rivers.png", dpi = 500,
       width=11, height=5, unit="in")

# bar plots legend
barplots_legend <- tm_microscopy_bar + theme(legend.position = "bottom")  +
  theme(
    legend.key.size = unit(1.25, "cm"), 
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15) 
  )
barplots_legend

# save!
ggsave("./figures/HABs_symposium/barplots_rivers_legend.png", dpi = 500,
       width=20, height=5, unit="in")

barplots_legend <- nt_microscopy_bar + theme(legend.position = "right")  +
  theme(
    legend.key.size = unit(1.25, "cm"), 
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15) 
  )
barplots_legend

# save!
ggsave("./figures/HABs_symposium/barplots_rivers_legend.png", dpi = 500,
       width=20, height=8, unit="in")

#### OLD ####

## Microcoleus !
microcoleus <- plot_grid(tm_microscopy_pcoa + coord_flip(clip = "off") +
                           geom_point(aes(color = site, shape = month, size = 15)) +
                           stat_ellipse(aes(color = site), linewidth = 1.5, linetype = 2),
                        tm_molecular_pcoa + coord_flip(clip = "off") +
                          geom_point(aes(color = site, shape = month, size = 15)), 
                                  ncol = 1, align = "hv")
microcoleus

# save!
ggsave("./figures/HABs_symposium/microcoleus_pcoa_rivers.png", dpi = 500,
       width=6.2, height=10, unit="in")

## Anabaena!
anabaena <- plot_grid(tac_microscopy_pcoa + coord_flip(clip = "off") +
                        geom_point(aes(color = site, shape = month, size = 15)),
                      NA,
                      tac_molecular_pcoa + coord_flip(clip = "off") +
                        geom_point(aes(color = site, shape = month, size = 15)),
                      tac_molecular_diversity + 
                        labs(y = "Shannon Diversity") +
                        theme(axis.text.x = element_markdown(size = 25)) + 
                        geom_boxplot(linewidth = 1, size = 6, aes(color = site)) + 
                        scale_color_manual(values = c("SAL" = "#2f6199",
                                                      "SFE-M" = "#224007",
                                                      "RUS" = "#696100")),
                      ncol = 2, align = "hv", scale = 0.95)
anabaena

# save!
ggsave("./figures/HABs_symposium/anabaena_pcoa_rivers.png", dpi = 500,
       width=12, height=10, unit="in")

# perfect