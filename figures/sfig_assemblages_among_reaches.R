#### Supplemental figure showing assemblage differences among rivers 
### Jordan Zabrecky
## 07.11.26

# This script creates a supplemental figure showing the PCoAs
# and PERMANOVA results from supplemental script,
# "S5d_differences_among_reaches.R"

#### (1) Loading libraries & data ####

# source supplemental script
source("./code/supplemental_code/S5d_differences_among_reaches.R")

# loading additional libraries
lapply(c("cowplot"), require, character.only = T)

# add label to p value table
p_table <- p_table %>% 
  mutate(label = case_when(p_value < 0.01 ~ paste("list(italic(p) == ", format(round(p_value, 3), nsmall = 2), ", italic(F) == ", 
                       format(round(F_stat, 2)), ")",
                       sep = ""),
                       TRUE ~ paste("list(italic(p) == ", format(round(p_value, 2), nsmall = 2), ", italic(F) == ", 
                                    format(round(F_stat, 2)), ")",
                                    sep = "")))

#### (2) Enhancing Individual Figures ####

# set universal plot theme
theme_set(theme_bw() + theme(strip.background = element_blank(),
                             plot.title = element_text(hjust = 0.5), legend.text = element_markdown(),
                             text = element_text(size = 7), strip.text = element_text(size = 7),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
                             axis.text.x = element_text(margin = margin(t = 0.005, unit = "cm")),
                             axis.text.y = element_text(margin = margin(t = 0.005, unit = "cm")),
                             axis.title.y = element_text(margin = margin(t = 0.005, unit = "cm"))))
# sfe palette
sfe_palette <- c("#a8e371", "#71a83d", "#3e6915")
rus_palette <- c("#e8db2e", "#ab9f00", "#635c00")

# enhance plots!

## (a) Microscopy PCoAs

# SFE NT
sfe_NT_micro <- PCoA_plots_sfe$NT +
  scale_color_discrete(palette = sfe_palette) +
  theme(legend.position = "none") + 
  annotate("text", x = -.09, y = .38,
           label = paste(p_table$label[which(p_table$sample_type == "NT" &
                                               p_table$river == "SFE-M" &
                                               p_table$test == "PERMANOVA" &
                                               p_table$data_type == "microscopy")]),
           parse = TRUE, color = "#3e6915", size = 7/.pt) +
  coord_cartesian(ylim = c(min(PCoA_list_sfe[[1]]$pcoa$pcoa2),
                           max(PCoA_list_sfe[[1]]$pcoa$pcoa2) + 0.1))
sfe_NT_micro

# SFE TM
sfe_TM_micro <- PCoA_plots_sfe$TM +
  scale_color_discrete(palette = sfe_palette) +
  theme(legend.position = "none") + 
  annotate("text", x = -.05, y = .35,
           label = paste(p_table$label[which(p_table$sample_type == "TM" &
                                               p_table$river == "SFE-M" &
                                               p_table$test == "PERMANOVA" &
                                               p_table$data_type == "microscopy")]),
           parse = TRUE, color = "#3e6915", size = 7/.pt) +
  coord_cartesian(ylim = c(min(PCoA_list_sfe[[3]]$pcoa$pcoa2),
                           max(PCoA_list_sfe[[3]]$pcoa$pcoa2) + 0.1))
sfe_TM_micro

# SFE TAC
sfe_TAC_micro <- PCoA_plots_sfe$TAC +
  scale_color_discrete(palette = sfe_palette) +
  theme(legend.position = "none") + 
  annotate("text", x = -.06, y = .15,
           label = paste(p_table$label[which(p_table$sample_type == "TAC" &
                                               p_table$river == "SFE-M" &
                                               p_table$test == "PERMANOVA" &
                                               p_table$data_type == "microscopy")]),
           parse = TRUE, color = "#3e6915", size = 7/.pt) +
  coord_cartesian(ylim = c(min(PCoA_list_sfe[[2]]$pcoa$pcoa2),
                           max(PCoA_list_sfe[[2]]$pcoa$pcoa2) + 0.07))
sfe_TAC_micro

# RUS NT
rus_NT_micro <- PCoA_plots_rus$NT +
  scale_color_discrete(palette = rus_palette) +
  theme(legend.position = "none") + 
  annotate("text", x = 0.02, y = .35,
           label = paste(p_table$label[which(p_table$sample_type == "NT" &
                                               p_table$river == "RUS" &
                                               p_table$test == "PERMANOVA" &
                                               p_table$data_type == "microscopy")]),
           parse = TRUE, color = "#635c00", size = 7/.pt) +
  coord_cartesian(ylim = c(min(PCoA_list_rus[[1]]$pcoa$pcoa2),
                           max(PCoA_list_rus[[1]]$pcoa$pcoa2) + 0.13))
rus_NT_micro

# RUS TAC
rus_TAC_micro <- PCoA_plots_rus$TAC +
  scale_color_discrete(palette = rus_palette) +
  theme(legend.position = "none") + 
  annotate("text", x = .01, y = .22,
           label = paste(p_table$label[which(p_table$sample_type == "TAC" &
                                               p_table$river == "RUS" &
                                               p_table$test == "PERMANOVA" &
                                               p_table$data_type == "microscopy")]),
           parse = TRUE, color = "#635c00", size = 7/.pt) +
  coord_cartesian(ylim = c(min(PCoA_list_rus[[2]]$pcoa$pcoa2),
                           max(PCoA_list_rus[[2]]$pcoa$pcoa2) + 0.06))
rus_TAC_micro

## (b) Molecular/16s PCoAs!

# SFE NT
sfe_NT_molec <- PCoA_plots_sfe_molec$NT +
  scale_color_discrete(palette = sfe_palette) +
  theme(legend.position = "none") + 
  annotate("text", x = .07, y = .4,
           label = paste(p_table$label[which(p_table$sample_type == "NT" &
                                               p_table$river == "SFE-M" &
                                               p_table$test == "PERMANOVA" &
                                               p_table$data_type == "molecular")]),
           parse = TRUE, color = "#3e6915", size = 7/.pt) +
  coord_cartesian(ylim = c(min(PCoA_list_sfe_molec[[1]]$pcoa$pcoa2),
                           max(PCoA_list_sfe_molec[[1]]$pcoa$pcoa2) + 0.12))
sfe_NT_molec

# SFE TM
sfe_TM_molec <- PCoA_plots_sfe_molec$TM +
  scale_color_discrete(palette = sfe_palette) +
  theme(legend.position = "none") + 
  annotate("text", x = -.10, y = .38,
           label = paste(p_table$label[which(p_table$sample_type == "TM" &
                                               p_table$river == "SFE-M" &
                                               p_table$test == "PERMANOVA" &
                                               p_table$data_type == "molecular")]),
           parse = TRUE, color = "#3e6915", size = 7/.pt) +
  coord_cartesian(ylim = c(min(PCoA_list_sfe_molec[[3]]$pcoa$pcoa2),
                           max(PCoA_list_sfe_molec[[3]]$pcoa$pcoa2) + 0.15))
sfe_TM_molec

# SFE TAC
sfe_TAC_molec <- PCoA_plots_sfe_molec$TAC +
  scale_color_discrete(palette = sfe_palette) +
  theme(legend.position = "none") + 
  annotate("text", x = -.07, y = .64,
           label = paste(p_table$label[which(p_table$sample_type == "TAC" &
                                               p_table$river == "SFE-M" &
                                               p_table$test == "PERMANOVA" &
                                               p_table$data_type == "molecular")]),
           parse = TRUE, color = "#3e6915", size = 7/.pt) +
  coord_cartesian(ylim = c(min(PCoA_list_sfe_molec[[2]]$pcoa$pcoa2),
                           max(PCoA_list_sfe_molec[[2]]$pcoa$pcoa2) + 0.15))
sfe_TAC_molec

# RUS NT (note: omitting outlier)
mod_rus_nt_molec <- PCoA_list_rus_molec[[1]]
mod_rus_nt_molec$pcoa <- mod_rus_nt_molec$pcoa %>% 
  filter(!(site_reach == "RUS-2" & field_date == "2022-08-17")) %>% 
  filter(!(site_reach == "RUS-1S" & field_date == "2022-09-01"))

rus_NT_molec <- makePCoAplot(mod_rus_nt_molec, shape = "site_reach", color = "site_reach") +
  scale_color_discrete(palette = rus_palette) +
  theme(legend.position = "none") + 
  annotate("text", x = 0.03, y = .36,
           label = paste(p_table$label[which(p_table$sample_type == "NT" &
                                               p_table$river == "RUS" &
                                               p_table$test == "PERMANOVA" &
                                               p_table$data_type == "molecular")]),
           parse = TRUE, color = "#635c00", size = 7/.pt) +
  stat_ellipse(aes(color = site_reach), type = "t", linetype = 2, linewidth = 0.5) +
  coord_cartesian(ylim = c(min(mod_rus_nt_molec$pcoa$pcoa2) - 0.1,
                           max(mod_rus_nt_molec$pcoa$pcoa2) + 0.25))
rus_NT_molec

# RUS TAC
rus_TAC_molec <- PCoA_plots_rus_molec$TAC +
  scale_color_discrete(palette = rus_palette) +
  theme(legend.position = "none") + 
  annotate("text", x = 0.07, y = .85,
           label = paste(p_table$label[which(p_table$sample_type == "TAC" &
                                               p_table$river == "RUS" &
                                               p_table$test == "PERMANOVA" &
                                               p_table$data_type == "molecular")]),
           parse = TRUE, color = "#635c00", size = 7/.pt) +
  stat_ellipse(aes(color = site_reach), type = "t", linetype = 2, linewidth = 0.5) +
  coord_cartesian(ylim = c(min(PCoA_list_rus_molec[[2]]$pcoa$pcoa2) - 0.4,
                           max(PCoA_list_rus_molec[[2]]$pcoa$pcoa2) + 0.53))
rus_TAC_molec

## (c) Functional PCoAs!

# SFE NT
sfe_NT_func <- PCoA_plots_sfe_func$NT +
  scale_color_discrete(palette = sfe_palette) +
  theme(legend.position = "none") + 
  annotate("text", x = .03, y = .13,
           label = paste(p_table$label[which(p_table$sample_type == "NT" &
                                               p_table$river == "SFE-M" &
                                               p_table$test == "PERMANOVA" &
                                               p_table$data_type == "functional")]),
           parse = TRUE, color = "#3e6915", size = 7/.pt) +
  coord_cartesian(ylim = c(min(PCoA_list_sfe_func[[1]]$pcoa$pcoa2),
                           max(PCoA_list_sfe_func[[1]]$pcoa$pcoa2) + 0.08))
sfe_NT_func

# SFE TM
sfe_TM_func <- PCoA_plots_sfe_func$TM +
  scale_color_discrete(palette = sfe_palette) +
  theme(legend.position = "none") + 
  annotate("text", x = -.06, y = .12,
           label = paste(p_table$label[which(p_table$sample_type == "TM" &
                                               p_table$river == "SFE-M" &
                                               p_table$test == "PERMANOVA" &
                                               p_table$data_type == "functional")]),
           parse = TRUE, color = "#3e6915", size = 7/.pt) +
  coord_cartesian(ylim = c(min(PCoA_list_sfe_func[[3]]$pcoa$pcoa2),
                           max(PCoA_list_sfe_func[[3]]$pcoa$pcoa2) + 0.05))
sfe_TM_func

# SFE TAC
sfe_TAC_func <- PCoA_plots_sfe_func$TAC +
  scale_color_discrete(palette = sfe_palette) +
  theme(legend.position = "none") + 
  annotate("text", x = -.02, y = .14,
           label = paste(p_table$label[which(p_table$sample_type == "TAC" &
                                               p_table$river == "SFE-M" &
                                               p_table$test == "PERMANOVA" &
                                               p_table$data_type == "functional")]),
           parse = TRUE, color = "#3e6915", size = 7/.pt) +
  coord_cartesian(ylim = c(min(PCoA_list_sfe_func[[2]]$pcoa$pcoa2),
                           max(PCoA_list_sfe_func[[2]]$pcoa$pcoa2) + 0.03))
sfe_TAC_func

# RUS NT
rus_NT_func <- PCoA_plots_rus_func$NT +
  scale_color_discrete(palette = rus_palette) +
  theme(legend.position = "none") + 
  annotate("text", x = -.06, y = .2,
           label = paste(p_table$label[which(p_table$sample_type == "NT" &
                                               p_table$river == "RUS" &
                                               p_table$test == "PERMANOVA" &
                                               p_table$data_type == "functional")]),
           parse = TRUE, color = "#635c00", size = 7/.pt) +
  coord_cartesian(ylim = c(min(PCoA_list_rus_func[[1]]$pcoa$pcoa2),
                           max(PCoA_list_rus_func[[1]]$pcoa$pcoa2) + 0.05))
rus_NT_func

# RUS TAC
rus_TAC_func <- PCoA_plots_rus_func$TAC +
  scale_color_discrete(palette = rus_palette) +
  theme(legend.position = "none") + 
  annotate("text", x = -.09, y = .51,
           label = paste(p_table$label[which(p_table$sample_type == "TAC" &
                                               p_table$river == "RUS" &
                                               p_table$test == "PERMANOVA" &
                                               p_table$data_type == "functional")]),
           parse = TRUE, color = "#635c00", size = 7/.pt) +
  stat_ellipse(aes(color = site_reach), type = "t", linetype = 2, linewidth = 0.5) +
  coord_cartesian(ylim = c(min(PCoA_list_rus_func[[2]]$pcoa$pcoa2) - .25,
                           max(PCoA_list_rus_func[[2]]$pcoa$pcoa2) + 0.32))
rus_TAC_func

##### (3) Making Main Figure ####

final <- plot_grid(sfe_NT_micro, sfe_TM_micro, sfe_TAC_micro, rus_NT_micro, rus_TAC_micro,
                   sfe_NT_molec, sfe_TM_molec, sfe_TAC_molec, rus_NT_molec, rus_TAC_molec,
                   sfe_NT_func, sfe_TM_func, sfe_TAC_func, rus_NT_func, rus_TAC_func,
                   nrow = 3, align = "hv")
final

ggsave("./figures/tiff_files/sfig_assemblages_among_reach.tiff", dpi = 500,
       width=17, height=10.5, unit="cm")

