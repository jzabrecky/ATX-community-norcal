#### Main figure comparing ATX concentrations to predicted gene counts of selected orthologs
### Jordan Zabrecky
## last edited: 07.03.2026

# This script creates a main figure comparing ATX concentrations versus predicted
# gene counts of selected orthologs for South Fork Eel samples. Version with 
# Russian samples is relegated to the supplemental figures

#### (1) Loading libraries & data ####

# source from previous script
source("./code/5d_anatoxins_functions.R")

# load additional library
library(scales)

# put all data together in one dataframe
all <- rbind(data_select$TM , data_select$TAC,
             data_select$NT) %>% 
  mutate(target_nontarget = case_when(sample_type == "NT" ~ "nontarget",
                                      TRUE ~ "target"))

# get linear model r-squared and p-values
model_info = data.frame(sample_type = NA,
                        functional_grouping = NA,
                        site = NA,
                        r_squared = NA,
                        p_value = NA)
for(i in c("SFE-M", "RUS")) {
  for(j in unique(all$sample_type)) {
    for(k in unique(all$functional_grouping)) {
      if(!(i == "RUS" & j == "TM")) {
        model = linear_model((all %>% filter(sample_type == j & functional_grouping == k & site == i)),
                             x = "log_ATX_all_ug_org_mat", y = "log_predicted_gene_abundance")
        model_info = rbind(model_info,
                           data.frame(sample_type = j,
                                      functional_grouping = k,
                                      site = i,
                                      r_squared = summary(model$model)$r.squared,
                                      p_value = summary(model$model)$coefficients[2,4]))

      }
    }
  }
}

# add min & max values for plotting
model_info <- left_join(model_info[-1,], all %>% 
  dplyr::group_by(site, functional_grouping, sample_type) %>% 
  dplyr::summarize(x_min = min(log_predicted_gene_abundance),
                   x_max = max(log_predicted_gene_abundance),
                   y_max = max(log_ATX_all_ug_org_mat),
                   y_min = max(log_ATX_all_ug_org_mat)),
                   by = c("site", "sample_type", "functional_grouping"))

# make labels & add nontarget/target column
model_info <- model_info %>% 
  # different rounding for r-squared
  mutate(r_squared_mod = case_when(r_squared <0.01 ~ round(r_squared, 3),
                                   TRUE ~ round(r_squared, 2))) %>% 
  mutate(label = case_when(p_value > 0.05 ~ paste("*r<sup>2</sup>* = ", r_squared_mod, ", *p* = ", format(round(p_value, 2), nsmall = 2),
                                            sep = ""),
                           TRUE ~ paste("*r<sup>2</sup>* = ", r_squared_mod, ", *p* = ", format(round(p_value, 2), nsmall = 2), "*",
                                       sep = ""))) %>% 
  mutate(target_nontarget = case_when(sample_type == "NT" ~ "nontarget",
                                      TRUE ~ "target"))

# only plot lines when model is significant
sig_models <- model_info %>% 
  filter(p_value < 0.05) # only three models are significant

# dataset with only significant models
sig_model_data <- all %>% 
  filter((sample_type == "TM" & functional_grouping == "nitrification" & site == "SFE-M") |
           (sample_type == "NT" & functional_grouping == "nitrogen_fixation" & site == "SFE-M") |
           (sample_type == "NT" & functional_grouping == "cobalamin_B12" & site == "SFE-M"))

# set universal plot theme
theme_set(theme_bw() + theme(strip.background = element_blank(),
                             plot.title = element_text(hjust = 0.5), legend.text = element_markdown(),
                             text = element_text(size = 8), strip.text = element_text(size = 8),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#### (2) Making Main Figure (South Fork Eel) ####

# nontarget
subplot1_1 <- ggplot(data = all %>% filter(sample_type == "NT" & site == "SFE-M"), 
                     aes(x = log_predicted_gene_abundance, y = log_ATX_all_ug_org_mat)) +
  geom_smooth(data = sig_model_data %>% filter(sample_type == "NT", site == "SFE-M"),
              fill = "#f2b979", color = "#c26f11", method = "lm", alpha = 0.25) +
  facet_grid(target_nontarget~functional_grouping, scales = "free") + 
  geom_richtext(
    data = model_info %>% filter(sample_type == "NT" & site == "SFE-M"),
    mapping = aes(x = x_min + ((x_max - x_min) * .44), y = -5, label = label), size= 2.4,
    fill = NA, label.color = NA, text.color = "#8c4d06") +
  geom_point(color = "#c26f11", size = 2, alpha = 0.5, shape = 15) +
  scale_x_continuous(breaks = breaks_pretty(n = 3)) +
  scale_y_continuous(breaks = breaks_pretty(n = 4), limit = c(-5.4, 5.5)) + 
  theme(strip.background = element_blank(), legend.position = "bottom")
subplot1_1

subplot1_2 <- ggplot(data = all %>% filter(sample_type == "TAC" & site == "SFE-M"), 
                     aes(x = log_predicted_gene_abundance, y = log_ATX_all_ug_org_mat)) +
  geom_smooth(data = sig_model_data %>% filter(sample_type == "TAC", site == "SFE-M"),
              color = "#3f9633", fill = "#acf0a3", method = "lm", alpha = 0.25) +
  facet_grid(target_nontarget~functional_grouping, scales = "free") + 
  geom_point(shape = 16, color = "#3f9633", size = 2, alpha = 0.5) +
  geom_richtext(
    data = model_info %>% filter(sample_type == "TAC" & site == "SFE-M"),
    mapping = aes(x = x_min + ((x_max - x_min) * .42), y = -5.5, label = label), size= 2.4,
    fill = NA, label.color = NA, text.color = "#247319") +
  scale_x_continuous(breaks = breaks_pretty(n = 3)) +
  scale_y_continuous(breaks = breaks_pretty(n = 4), limit = c(-6, 6)) + 
  theme(strip.background = element_blank(), legend.position = "bottom")
subplot1_2

subplot1_3 <- ggplot(data = all %>% filter(sample_type == "TM" & site == "SFE-M"), 
                     aes(x = log_predicted_gene_abundance, y = log_ATX_all_ug_org_mat)) +
  geom_smooth(data = sig_model_data %>% filter(sample_type == "TM", site == "SFE-M"),
              color = "#6d4275", fill = "#c79bcf", method = "lm", alpha = 0.25) +
  facet_grid(target_nontarget~functional_grouping, scales = "free") + 
  geom_point(shape = 17, color = "#6d4275", size = 2, alpha = 0.5) +
  geom_richtext(
    data = model_info %>% filter(sample_type == "TM" & site == "SFE-M"),
    mapping = aes(x = x_min + ((x_max - x_min) * .42), y = -5.8, label = label), size= 2.4,
    fill = NA, label.color = NA, text.color = "#48234f") +
  scale_x_continuous(breaks = breaks_pretty(n = 3)) +
  theme(strip.background = element_blank(), legend.position = "bottom")
subplot1_3

plot_sfe <- plot_grid(subplot1_1, subplot1_3, subplot1_2, align = "hv", ncol = 1)
plot_sfe

# save figure
ggsave("./figures/tiff_files/Q2_atx_v_gene.tiff", dpi = 600,
       width=18, height=12, unit="cm")

#### (3) Making Supplemental Figure (Russian River) ####

# nontarget
subplot2_1 <- ggplot(data = all %>% filter(sample_type == "NT" & site == "RUS"), 
                     aes(x = log_predicted_gene_abundance, y = log_ATX_all_ug_org_mat)) +
  geom_smooth(data = sig_model_data %>% filter(sample_type == "NT" & site == "RUS"),
              fill = "#f2b979", color = "#c26f11", method = "lm", alpha = 0.25) +
  facet_grid(target_nontarget~functional_grouping, scales = "free") + 
  geom_richtext(
    data = model_info %>% filter(sample_type == "NT" & site == "RUS"),
    mapping = aes(x = x_min + ((x_max - x_min) * .42), y = -5, label = label), size= 2.4,
    fill = NA, label.color = NA, text.color = "#8c4d06") +
  geom_point(color = "#c26f11", size = 2, alpha = 0.5, shape = 15) +
  scale_x_continuous(breaks = breaks_pretty(n = 3)) +
  theme(strip.background = element_blank(), legend.position = "bottom")
subplot2_1

subplot2_2 <- ggplot(data = all %>% filter(sample_type == "TAC" & site == "RUS"), 
                     aes(x = log_predicted_gene_abundance, y = log_ATX_all_ug_org_mat)) +
  geom_smooth(data = sig_model_data %>% filter(sample_type == "TAC" & site == "RUS"),
              color = "#3f9633", fill = "#acf0a3", method = "lm", alpha = 0.25) +
  facet_grid(target_nontarget~functional_grouping, scales = "free") + 
  geom_point(shape = 16, color = "#3f9633", size = 2, alpha = 0.5) +
  geom_richtext(
    data = model_info %>% filter(sample_type == "TAC" & site == "RUS"),
    mapping = aes(x = x_min + ((x_max - x_min) * .44), y = -4.8, label = label), size= 2.4,
    fill = NA, label.color = NA, text.color = "#247319") +
  scale_x_continuous(breaks = breaks_pretty(n = 3)) +
  scale_y_continuous(breaks = breaks_pretty(n = 4), limit = c(-5, 2)) + 
  theme(strip.background = element_blank(), legend.position = "bottom")
subplot2_2

# put plot together
plot_rus <- plot_grid(subplot2_1, subplot2_2,align = "hv", ncol = 1)
plot_rus

# save figure
ggsave("./figures/tiff_files/sfig_rus_atx_v_gene.tiff", dpi = 600,
       width=18, height=8.5, unit="cm")
