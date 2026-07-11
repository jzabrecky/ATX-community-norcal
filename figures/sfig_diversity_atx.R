#### Supplemental figures comparing ATX concentrations vs. Shannon Diversity
### Jordan Zabrecky
## last edited: 07.11.2026

# This script creates a supplemental figure showing the linear models
# between ATX concentrations and Shannon Diversity Indeces for 16s samples
# Two versions are made: (1) without standardization for reach 
# (within each sample type) and (2)  with standardization within each reach 
# and sample type. We go with option 2 as it shows that nitrification in TM
# and cobalamin in NT are not signficant with (appropriate) standardization

#### (1) Loading libraries & data ####

# source from previous script
source("./code/5c_anatoxins_16s.R")

#### (2) Version without standardization ####

# make tabel of lm results for plot text
model_info <- data.frame(plot_factor = factor(c("SFE-MNT", "RUSNT", "SFE-MTM", "SFE-MTAC", "RUSTAC"), 
                                              levels = c("SFE-MNT", "RUSNT", "SFE-MTM", "SFE-MTAC", "RUSTAC")),
                         sample_type = c("NT_text", "NT_text", "TM_text", "TAC_text", "TAC_text"),
                         r_squared = c(summary(diversity_sfe_nt$model)$r.squared,
                                       summary(diversity_rus_nt$model)$r.squared,
                                       summary(diversity_sfe_tm$model)$r.squared,
                                       summary(diversity_sfe_tac$model)$r.squared,
                                       summary(diversity_rus_tac$model)$r.squared),
                         p_value = c(summary(diversity_sfe_nt$model)$coefficients[2,4],
                                     summary(diversity_rus_nt$model)$coefficients[2,4],
                                     summary(diversity_sfe_tm$model)$coefficients[2,4],
                                     summary(diversity_sfe_tac$model)$coefficients[2,4],
                                     summary(diversity_rus_tac$model)$coefficients[2,4]),
                         x= c(min(diversity$NT$shannon_diversity[which(diversity$NT$site == "SFE-M")]) + 1/3.5,
                              min(diversity$NT$shannon_diversity[which(diversity$NT$site == "RUS")]) + 1/5,
                              min(diversity$TM$shannon_diversity[which(diversity$TM$site == "SFE-M")]) - 1/1.6,
                              min(diversity$TAC$shannon_diversity[which(diversity$TAC$site == "SFE-M")]) - 1/4,
                              min(diversity$TAC$shannon_diversity[which(diversity$TAC$site == "RUS")])) + 1,
                         y = c(4.5, 1.4, 5.0, 4.6, 1.4)) %>%
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

# edit data frame
data <- data %>% 
  # remove Salmon River samples
  filter(site != "SAL") %>%
  # remove NAs with no matching ATX
  na.omit() %>% 
  # make factor for desired plot order (doing wrapping so combining to keep text on one level)
  mutate(site_sample_type = paste(site, sample_type, sep = ""),
         plot_factor = factor(site_sample_type, levels = c("SFE-MNT", "SFE-MTM", "SFE-MTAC", "RUSNT","RUSTAC")))

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
    fill = NA, label.color = NA) + 
  facet_wrap(~plot_factor, scales = "free", ncol = 1) +
  labs(x = "Shannon Diversity", y = "log(&mu;g anatoxins g <sup>-1</sup> OM)") +
  theme(legend.position = "none", axis.title.y = element_markdown(), strip.text = element_blank())
plot

#### (4) Standardized by Reach ####

# make a second dataframe for standardization
data_stnd <- data

# standardize within reach and sample type
data_stnd[,c("stnd_shannon_diversity", "stnd_log_ATX_all_ug_org_mat")] <- 
  apply(data_stnd[,c("shannon_diversity", "log_ATX_all_ug_org_mat")], MARGIN = 2,
        function(x) ave(x, data_stnd$site_reach, data_stnd$sample_type, FUN = scale))

# test on one reach & sample type to make sure that was done correctly
test <- data %>% filter(site_reach == "SFE-M-1S" & sample_type == "TM")
test[,c("stnd_shannon_diversity", "stnd_log_ATX_all_ug_org_mat")] <- 
  scale(test[,c("shannon_diversity", "log_ATX_all_ug_org_mat")])
all.equal(test, data_stnd %>% filter(site_reach == "SFE-M-1S" & sample_type == "TM")) # we good!

# apparently RUS-2 never has anatoxin detection for TAC! those can all be 0 which is center
# as there is no variability
data_stnd$stnd_log_ATX_all_ug_org_mat[which(data_stnd$sample_type == "TAC" &
                                              data_stnd$site_reach == "RUS-2")] <- 0

# remove old dataframes to make sure we call the right things
rm(data, model_info, diversity_sfe_nt, diversity_rus_tac, diversity_rus_nt,
   diversity_sfe_tm, diversity_sfe_tac, sig_models, sig_model_data)

# make new linear models using standardized data
stnd_diversity_sfe_nt <- linear_model(data_stnd %>% filter(site == "SFE-M" & sample_type == "NT"),
                                      x = "stnd_shannon_diversity", y = "stnd_log_ATX_all_ug_org_mat")
stnd_diversity_sfe_tm <- linear_model(data_stnd %>% filter(site == "SFE-M" & sample_type == "TM"),
                                      x = "stnd_shannon_diversity", y = "stnd_log_ATX_all_ug_org_mat")
stnd_diversity_sfe_tac <- linear_model(data_stnd %>% filter(site == "SFE-M" & sample_type == "TAC"),
                                      x = "stnd_shannon_diversity", y = "stnd_log_ATX_all_ug_org_mat")
stnd_diversity_rus_nt <- linear_model(data_stnd %>% filter(site == "RUS" & sample_type == "NT"),
                                      x = "stnd_shannon_diversity", y = "stnd_log_ATX_all_ug_org_mat")
stnd_diversity_rus_tac <- linear_model(data_stnd %>% filter(site == "RUS" & sample_type == "TAC"),
                                       x = "stnd_shannon_diversity", y = "stnd_log_ATX_all_ug_org_mat")

# make new model info table with linear models using standardized data
model_info_stnd <- data.frame(plot_factor = factor(c("SFE-MNT", "RUSNT", "SFE-MTM", "SFE-MTAC", "RUSTAC"), 
                                              levels = c("SFE-MNT", "RUSNT", "SFE-MTM", "SFE-MTAC", "RUSTAC")),
                         sample_type = c("NT_text", "NT_text", "TM_text", "TAC_text", "TAC_text"),
                         r_squared = c(summary(stnd_diversity_sfe_nt$model)$r.squared,
                                       summary(stnd_diversity_rus_nt$model)$r.squared,
                                       summary(stnd_diversity_sfe_tm$model)$r.squared,
                                       summary(stnd_diversity_sfe_tac$model)$r.squared,
                                       summary(stnd_diversity_rus_tac$model)$r.squared),
                         p_value = c(summary(stnd_diversity_sfe_nt$model)$coefficients[2,4],
                                     summary(stnd_diversity_rus_nt$model)$coefficients[2,4],
                                     summary(stnd_diversity_sfe_tm$model)$coefficients[2,4],
                                     summary(stnd_diversity_sfe_tac$model)$coefficients[2,4],
                                     summary(stnd_diversity_rus_tac$model)$coefficients[2,4]),
                         x= c(-1.5, -1.1, -1.75, -1.5, -1.4),
                         y = c(1.25, 1.5, 0.75, 1, 1)) %>%
  # add label, include asterisk based on p-value
  mutate(r_squared_mod = case_when(r_squared <0.01 ~ round(r_squared, 3),
                                   TRUE ~ round(r_squared, 2)),
         label = paste("*r<sup>2</sup>* = ", r_squared_mod, ", *p* = ", format(round(p_value, 2), nsmall = 2),
                                                  sep = ""))

# dataset with only significant models
stnd_sig_model_data <- data_stnd %>% 
  filter(sample_type == "NT" & site == "SFE-M")

# make plot
plot2 <- ggplot(data = data_stnd, aes(y = stnd_log_ATX_all_ug_org_mat, x = stnd_shannon_diversity)) +
  geom_smooth(data = stnd_sig_model_data, aes(color = sample_type, fill = sample_type), method = "lm", alpha = 0.25) +
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
    data = model_info_stnd,
    mapping = aes(x = x, y = y, label = label, color = sample_type), size= 3,
    fill = NA, label.color = NA) + 
  facet_wrap(~plot_factor, scales = "free", ncol = 1) +
  labs(x = "standardized Shannon Diversity", y = "standardized log(&mu;g anatoxins g <sup>-1</sup> OM)") +
  theme(legend.position = "none", axis.title.y = element_markdown(), strip.text = element_blank())
plot2

# save plot!
ggsave("./figures/tiff_files/sfig_atx_vs_diversity.tiff", dpi = 500,
       width=12, height=16, unit="cm")
