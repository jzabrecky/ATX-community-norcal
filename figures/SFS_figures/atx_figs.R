## SFS ATX Figures

# source 
source("./code/6b_anatoxins_microscopy.R")

## (a) Morphological NMDS

# TM
NMDS_list_eel[[3]]
NMDS_plots_eel[[3]]
NMDS_plots_eel[[2]]
plot <- makeNMDSplot(NMDS_list_eel[[3]], TRUE, FALSE, color = "atx_group",
             shape = "atx_group") +
  scale_color_discrete(palette = rev(c("#949494", "#7ae06c", "#3f9633", "#13570a")))
plot
plot <- makeNMDSplot(NMDS_list_eel[[2]], TRUE, FALSE, color = "atx_group",
                     shape = "atx_group") +
  scale_color_discrete(palette = rev(c("#949494", "#7ae06c", "#3f9633", "#13570a")))
plot

stress = NMDS_list_eel[[3]]$stress
plot = ggplot(NMDS_list_eel[[3]]$nmds, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = atx_group, shape = atx_group), size = 2) +
  labs(subtitle = paste("Stress:", round(stress, 3)),
       x = "NMDS Axis 1",
       y = "NMDS Axis 2") +
  theme(plot.subtitle = element_text(hjust=1, vjust=0.5)) +
  scale_color_discrete(palette = rev(c("#949494", "#7ae06c", "#3f9633", "#13570a"))) +
  theme(legend.position = "bottom")
plot

# TAC
NMDS_plots_eel[[2]]


stress = NMDS_list_eel[[2]]$stress
plot2 = ggplot(NMDS_list_eel[[2]]$nmds, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = atx_group, shape = atx_group), size = 2) +
  labs(subtitle = paste("Stress:", round(stress, 3)),
       x = "NMDS Axis 1",
       y = "NMDS Axis 2") +
  theme(plot.subtitle = element_text(hjust=1, vjust=0.5)) +
  scale_color_discrete(palette = rev(c("#949494", "#7ae06c", "#3f9633", "#13570a"))) +
  theme(legend.position = "bottom")
plot2

both <- plot_grid(plot, plot2, nrow = 1)
both

ggsave("./figures/SFS_figures/atx_morpho.tiff", dpi = 600,
       width=17.6, height=8, unit="cm")

plot + theme(legend.position = "right")

ggsave("./figures/SFS_figures/atx_legend.tiff", dpi = 600,
       width=17.6, height=8, unit="cm")

#### EMPTY ENVIRONMENT :)

source("./code/6c_anatoxins_16s.R")

stress = NMDS_list_eel$TM$stress
plot3 = ggplot(NMDS_list_eel$TM$nmds, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = atx_group, shape = atx_group), size = 2) +
  labs(subtitle = paste("Stress:", round(stress, 3)),
       x = "NMDS Axis 1",
       y = "NMDS Axis 2") +
  theme(plot.subtitle = element_text(hjust=1, vjust=0.5)) +
  scale_color_discrete(palette = rev(c("#949494", "#7ae06c", "#3f9633", "#13570a"))) +
  theme(legend.position = "bottom")
plot3

stress = sfe_tac_NMDS_list$stress
plot4 = ggplot(sfe_tac_NMDS_list$nmds, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = atx_group, shape = atx_group), size = 2) +
  labs(subtitle = paste("Stress:", round(stress, 3)),
       x = "NMDS Axis 1",
       y = "NMDS Axis 2") +
  theme(plot.subtitle = element_text(hjust=1, vjust=0.5)) +
  scale_color_discrete(palette = rev(c("#949494", "#7ae06c", "#3f9633", "#13570a"))) +
  theme(legend.position = "bottom")
plot4


both <- plot_grid(plot3, plot4, nrow = 1)
both

ggsave("./figures/SFS_figures/atx_molecular.tiff", dpi = 600,
       width=17.6, height=8, unit="cm")

#### EMPTY ENVIRONMENT AGAIN :)

source("./code/6d_anatoxins_functions.R")

stress = NMDS_list_all_wide$`SFE-M TM`$stress
plot5 = ggplot(NMDS_list_all_wide$`SFE-M TM`$nmds, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = atx_group, shape = atx_group), size = 2) +
  labs(subtitle = paste("Stress:", round(stress, 3)),
       x = "NMDS Axis 1",
       y = "NMDS Axis 2") +
  theme(plot.subtitle = element_text(hjust=1, vjust=0.5)) +
  scale_color_discrete(palette = rev(c("#949494", "#7ae06c", "#3f9633", "#13570a"))) +
  theme(legend.position = "bottom")
plot5

stress = NMDS_list_all_wide$`SFE-M TAC`$stress
plot6 = ggplot(NMDS_list_all_wide$`SFE-M TAC`$nmds, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = atx_group, shape = atx_group), size = 2) +
  labs(subtitle = paste("Stress:", round(stress, 3)),
       x = "NMDS Axis 1",
       y = "NMDS Axis 2") +
  theme(plot.subtitle = element_text(hjust=1, vjust=0.5)) +
  scale_color_discrete(palette = rev(c("#949494", "#7ae06c", "#3f9633", "#13570a"))) +
  theme(legend.position = "bottom")
plot6

both <- plot_grid(plot5, plot6, nrow = 1)
both

ggsave("./figures/SFS_figures/atx_functional.tiff", dpi = 600,
       width=17.6, height=8, unit="cm")

eel_select <- rbind(data_select$TM %>% filter(site == "SFE-M"), data_select$TAC %>% filter(site == "SFE-M"))

plot = ggplot(data = eel_select, aes(x = log_predicted_gene_abundance, y = log_ATX_all_ug_org_mat)) +
  geom_smooth(aes(color = sample_type, fill = sample_type), method = "lm") +
  facet_wrap(~functional_grouping, scales = "free") + 
  scale_color_discrete(palette = c("#3f9633", "#6d4275")) +
  scale_fill_discrete(palette = c("#acf0a3", "#c79bcf")) +
  geom_point(aes(shape = sample_type, color = sample_type), size = 2, alpha = 0.7) +
  theme(strip.background = element_blank())
plot                    

ggsave("./figures/SFS_figures/atx_vs_gene.tiff", dpi = 600,
       width=17.6, height=8, unit="cm")


plot = ggplot(data = eel_select, aes(x = predicted_gene_abundance, y = log_ATX_all_ug_org_mat)) +
  geom_smooth(aes(color = sample_type, fill = sample_type), method = "lm") +
  facet_wrap(~functional_grouping, scales = "free") + 
  scale_color_discrete(palette = c("#3f9633", "#6d4275")) +
  scale_fill_discrete(palette = c("#acf0a3", "#c79bcf")) +
  geom_point(aes(shape = sample_type, color = sample_type), size = 2, alpha = 0.7) +
  theme(strip.background = element_blank())
plot 

eel_select_split <- split(eel_select, eel_select$sample_type)
tac_models <- lapply(unique(eel_select$functional_grouping), function(x) {
  model = lm(log_ATX_all_ug_org_mat ~ log_predicted_gene_abundance, 
     data = eel_select_split$TAC %>% filter(functional_grouping == x))
  print(x)
  print(round(summary(model)$coefficients[2,4], 2))
  return(model)
})
eel_select_split <- split(eel_select, eel_select$sample_type)
tm_models <- lapply(unique(eel_select$functional_grouping), function(x) {
  model = lm(log_ATX_all_ug_org_mat ~ log_predicted_gene_abundance, 
             data = eel_select_split$TM %>% filter(functional_grouping == x))
  print(x)
  print(round(summary(model)$coefficients[2,4], 2))
  return(model)
})

#### REMOVE ENVIRONMENT ONE MORE TIME :)

atx <- read.csv("./data/field_and_lab/atx_w_categorical_groupings.csv")

theme_set(theme_bw() + theme(strip.background = element_blank(),
                             plot.title = element_text(hjust = 0.5), legend.text = element_markdown(),
                             text = element_text(size = 8), strip.text = element_text(size = 8),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),))


plot = ggplot(data = atx %>% filter(sample_type != "NT"), aes(x = site, y = ATX_all_ug_org_mat, fill = site)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(alpha = 0.5, aes(fill = site), color = "black", shape = 21) + 
  scale_fill_manual(values = c("SAL" = "#62a7f8",
                               "SFE-M" = "#416f16",
                               "RUS" = "#bdb000")) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none")
plot

ggsave("./figures/SFS_figures/atx_rivers.tiff", dpi = 600,
       width=8, height=6, unit="cm")
