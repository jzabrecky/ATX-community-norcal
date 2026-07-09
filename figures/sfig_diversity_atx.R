#### Supplemental figures comparing ATX concentrations vs. Shannon Diversity
### Jordan Zabrecky
## last edited: 06.29.2026

# This script creates a supplemental figure showing the linear models
# between ATX concentrations and Shannon Diversity Indeces for 16s samples

#### (1) Loading libraries & data ####

# source from previous script
source("./code/5c_anatoxins_16s.R")

#### (2) 'lm' Test Results ####

# make tabel of lm results for plot text
model_info <- data.frame(plot_factor = factor(c("SFE-MNT", "RUSNT", "SFE-MTM", "RUSTM", "SFE-MTAC", "RUSTAC"), 
                                              levels = c("SFE-MNT", "RUSNT", "SFE-MTM", "RUSTM", "SFE-MTAC", "RUSTAC")),
                         sample_type = c("NT_text", "NT_text", "TM_text", "TM_text", "TAC_text", "TAC_text"),
                         r_squared = c(summary(diversity_sfe_nt$model)$r.squared,
                                       summary(diversity_rus_nt$model)$r.squared,
                                       summary(diversity_sfe_tm$model)$r.squared,
                                       summary(diversity_sfe_tm$model)$r.squared,
                                       summary(diversity_sfe_tac$model)$r.squared,
                                       summary(diversity_rus_tac$model)$r.squared),
                         p_value = c(summary(diversity_sfe_nt$model)$coefficients[2,4],
                                     summary(diversity_rus_nt$model)$coefficients[2,4],
                                     summary(diversity_sfe_tm$model)$coefficients[2,4],
                                     summary(diversity_sfe_tm$model)$coefficients[2,4],
                                     summary(diversity_sfe_tac$model)$coefficients[2,4],
                                     summary(diversity_rus_tac$model)$coefficients[2,4]),
                         x= c(min(diversity$NT$shannon_diversity[which(diversity$NT$site == "SFE-M")]) + 1/3.5,
                              min(diversity$NT$shannon_diversity[which(diversity$NT$site == "RUS")]) + 1/5,
                              min(diversity$TM$shannon_diversity[which(diversity$TM$site == "SFE-M")]) - 1/1.6,
                              min(diversity$TM$shannon_diversity[which(diversity$TM$site == "SFE-M")]) - 1/1.6,
                              min(diversity$TAC$shannon_diversity[which(diversity$TAC$site == "SFE-M")]) - 1/4,
                              min(diversity$TAC$shannon_diversity[which(diversity$TAC$site == "RUS")])) + 1,
                         y = c(4.5, 1.4, 5.0, 5.0, 4.6, 1.4)) %>%
  # add label, include asterisk based on p-value
  mutate(r_squared_mod = case_when(r_squared <0.01 ~ round(r_squared, 3),
                                   TRUE ~ round(r_squared, 2)),
         label = case_when(p_value > 0.05 ~ paste("*r<sup>2</sup>* = ", r_squared_mod, ", *p* = ", format(round(p_value, 2), nsmall = 2),
                                                  sep = ""),
                           # the single value is ** <0.05 and <0.01 but >0.001
                           TRUE ~ paste("*r<sup>2</sup>* = ", r_squared_mod, ", *p* = ", format(round(p_value, 3), nsmall = 2), "**",
                                        sep = "")))

# line for significant models only
sig_models <- model_info %>% 
  filter(p_value < 0.05) # SFE-M NT

#### (3) Making Plots ####

# set universal plot theme
theme_set(theme_bw() + theme(strip.background = element_blank(),
                             plot.title = element_text(hjust = 0.5), legend.text = element_markdown(),
                             text = element_text(size = 10),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# put data together all in one dataframe
data <- rbind(diversity$NT, diversity$TAC, diversity$TM)

# add in dummy sample to get a fake plot in TM space for RUS
data <- rbind((data[which(data$sample_type == "TM"),]) %>% 
                mutate(site = case_when(site == "SFE-M" ~ "RUS")), data)

# edit data frame
data <- data %>% 
  # remove Salmon River samples
  filter(site != "SAL") %>%
  # remove NAs with no matching ATX
  na.omit() %>% 
  # make factor for desired plot order (doing wrapping so combining to keep text on one level)
  mutate(site_sample_type = paste(site, sample_type, sep = ""),
         plot_factor = factor(site_sample_type, levels = c("SFE-MNT", "RUSNT", "SFE-MTM", "RUSTM", "SFE-MTAC", "RUSTAC")))

# dataset with only significant models
sig_model_data <- data %>% 
  filter(sample_type == "NT" & site == "SFE-M")

# make plot
plot <- ggplot(data = data, aes(y = log_ATX_all_ug_org_mat, x = shannon_diversity)) +
  geom_smooth(data = sig_model_data, aes(color = sample_type, fill = sample_type), method = "lm", alpha = 0.25) +
  geom_point(aes(color = sample_type, shape = sample_type), size = 2, alpha = 0.5) + 
  scale_color_manual(values = c("NT" = "#c26f11",
                                "TM" = "#6d4275",
                                "TAC" = "#247319",
                                "NT_text" = "#8c4d06",
                                "TM_text" = "#48234f",
                                "TAC_text" = "#247319")) +
  scale_fill_manual(values = c("NT" = "#f2b979",
                               "TM" = "#c79bcf",
                               "TAC" = "#acf0a3")) +
  scale_shape_manual(values = c("NT" = 15,
                                "TM" = 17,
                                "TAC" = 16)) +
  geom_richtext(
    data = model_info,
    mapping = aes(x = x, y = y, label = label, color = sample_type), size= 3,
    fill = NA, label.color = NA) + # add text.color based on sample type???!?!?!
  facet_wrap(~plot_factor, scales = "free", ncol = 2) +
  labs(x = "Shannon Diversity", y = "log(&mu;g anatoxins g <sup>-1</sup> OM)") +
  theme(legend.position = "none", axis.title.y = element_markdown(), strip.text = element_blank())
plot

# save plot!
ggsave("./figures/tiff_files/sfig_atx_vs_diversity.tiff", dpi = 500,
       width=17, height=14, unit="cm")
