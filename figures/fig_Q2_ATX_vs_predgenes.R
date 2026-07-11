#### Main figure comparing ATX concentrations to predicted gene counts of selected orthologs
### Jordan Zabrecky
## last edited: 07.11.2026

# This script creates a main figure comparing ATX concentrations versus predicted
# gene counts of selected orthologs for South Fork Eel samples. Version with 
# Russian samples is relegated to the supplemental figures. Two versions are made
# (1) without standardization for reach (within each sample type) and (2) 
# with standardization within each reach and sample type. We go with option 2
# as it shows that nitrification in TM and cobalamin in NT are not signficant
# with (appropriate) standardization

#### (1) Loading libraries & data ####

# source from previous script
source("./code/5d_anatoxins_functions.R")

# load additional library
library(scales)

# put all data together in one dataframe
all <- rbind(data_select$TM , data_select$TAC,
             data_select$NT)

# pivot wider for standardization
all_wider <- all %>% 
  select(site, site_reach, field_date, sample_type, log_ATX_all_ug_org_mat, functional_grouping, log_predicted_gene_abundance) %>% 
  pivot_wider(names_from = "functional_grouping", values_from = "log_predicted_gene_abundance", 
              # we had two samples that didn't have any nitrification genes, so sub in zero
              # (TM SFE-M 1S and SFE-M 3 8-10-22)
              # can just fill in with zero as the minimum log predicted gene abundance is 0.69 for the true value of 2
              values_fill = 0) 

# pivot back longer so zeroes are preserved
all <- all_wider %>% 
  pivot_longer(cols = c("cobalamin_B12", "nitrification",
                        "nitrogen_fixation", "phosphatase_transporters", 
                        "pyridoxal", "thiamine"), 
               values_to = "log_predicted_gene_abundance", names_to = "functional_grouping")

#### (2) No Reach Standardization ####

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
                             y = "log_ATX_all_ug_org_mat", x = "log_predicted_gene_abundance")
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
                                       sep = "")))

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

## SFE plots

# nontarget
subplot1_1 <- ggplot(data = all %>% filter(sample_type == "NT" & site == "SFE-M"), 
                     aes(x = log_predicted_gene_abundance, y = log_ATX_all_ug_org_mat)) +
  geom_smooth(data = sig_model_data %>% filter(sample_type == "NT", site == "SFE-M"),
              fill = "#f2b979", color = "#c26f11", method = "lm", alpha = 0.25) +
  facet_grid(~functional_grouping, scales = "free") + 
  geom_richtext(
    data = model_info %>% filter(sample_type == "NT" & site == "SFE-M"),
    mapping = aes(x = x_min + ((x_max - x_min) * .44), y = -5.3, label = label), size= 2.4,
    fill = NA, label.color = NA, text.color = "#8c4d06") +
  geom_point(aes(shape = site_reach), color = "#c26f11", size = 2, alpha = 0.5) +
  scale_x_continuous(breaks = breaks_pretty(n = 3)) +
  scale_y_continuous(breaks = breaks_pretty(n = 4), limit = c(-5.4, 5.5)) + 
  theme(strip.background = element_blank(), legend.position = "bottom")
subplot1_1

subplot1_2 <- ggplot(data = all %>% filter(sample_type == "TAC" & site == "SFE-M"), 
                     aes(x = log_predicted_gene_abundance, y = log_ATX_all_ug_org_mat)) +
  geom_smooth(data = sig_model_data %>% filter(sample_type == "TAC", site == "SFE-M"),
              color = "#3f9633", fill = "#acf0a3", method = "lm", alpha = 0.25) +
  facet_grid(~functional_grouping, scales = "free") + 
  geom_point(aes(shape = site_reach), color = "#3f9633", size = 2, alpha = 0.5) +
  geom_richtext(
    data = model_info %>% filter(sample_type == "TAC" & site == "SFE-M"),
    mapping = aes(x = x_min + ((x_max - x_min) * .44), y = -5.5, label = label), size= 2.4,
    fill = NA, label.color = NA, text.color = "#247319") +
  scale_x_continuous(breaks = breaks_pretty(n = 3)) +
  scale_y_continuous(breaks = breaks_pretty(n = 4), limit = c(-6, 6)) + 
  theme(strip.background = element_blank(), legend.position = "bottom")
subplot1_2

subplot1_3 <- ggplot(data = all %>% filter(sample_type == "TM" & site == "SFE-M"), 
                     aes(x = log_predicted_gene_abundance, y = log_ATX_all_ug_org_mat)) +
  geom_smooth(data = sig_model_data %>% filter(sample_type == "TM", site == "SFE-M"),
              color = "#6d4275", fill = "#c79bcf", method = "lm", alpha = 0.25) +
  facet_grid(~functional_grouping, scales = "free") + 
  geom_point(aes(shape = site_reach), color = "#6d4275", size = 2, alpha = 0.5) +
  geom_richtext(
    data = model_info %>% filter(sample_type == "TM" & site == "SFE-M"),
    mapping = aes(x = x_min + ((x_max - x_min) * .42), y = -5.8, label = label), size= 2.4,
    fill = NA, label.color = NA, text.color = "#48234f") +
  scale_x_continuous(breaks = breaks_pretty(n = 3)) +
  theme(strip.background = element_blank(), legend.position = "bottom")
subplot1_3

## RUS plots

# nontarget
subplot2_1 <- ggplot(data = all %>% filter(sample_type == "NT" & site == "RUS"), 
                     aes(x = log_predicted_gene_abundance, y = log_ATX_all_ug_org_mat)) +
  geom_smooth(data = sig_model_data %>% filter(sample_type == "NT" & site == "RUS"),
              fill = "#f2b979", color = "#c26f11", method = "lm", alpha = 0.25) +
  facet_grid(~functional_grouping, scales = "free") + 
  geom_richtext(
    data = model_info %>% filter(sample_type == "NT" & site == "RUS"),
    mapping = aes(x = x_min + ((x_max - x_min) * .44), y = -4.5, label = label), size= 2.4,
    fill = NA, label.color = NA, text.color = "#8c4d06") +
  geom_point(aes(shape = site_reach), color = "#c26f11", size = 2, alpha = 0.5) +
  scale_x_continuous(breaks = breaks_pretty(n = 3)) +
  theme(strip.background = element_blank(), legend.position = "bottom")
subplot2_1

subplot2_2 <- ggplot(data = all %>% filter(sample_type == "TAC" & site == "RUS"), 
                     aes(x = log_predicted_gene_abundance, y = log_ATX_all_ug_org_mat)) +
  geom_smooth(data = sig_model_data %>% filter(sample_type == "TAC" & site == "RUS"),
              color = "#3f9633", fill = "#acf0a3", method = "lm", alpha = 0.25) +
  facet_grid(~functional_grouping, scales = "free") + 
  geom_point(aes(shape = site_reach), color = "#3f9633", size = 2, alpha = 0.5) +
  geom_richtext(
    data = model_info %>% filter(sample_type == "TAC" & site == "RUS"),
    mapping = aes(x = x_min + ((x_max - x_min) * .44), y = -4.5, label = label), size= 2.4,
    fill = NA, label.color = NA, text.color = "#247319") +
  scale_x_continuous(breaks = breaks_pretty(n = 3)) +
  scale_y_continuous(breaks = breaks_pretty(n = 4)) + 
  theme(strip.background = element_blank(), legend.position = "bottom")
subplot2_2

#### (3) Z-Score test ####

# to make sure our results above aren't results of Simpson's paradox,
# standardize genes and anatoxin concentrations by reach

# make duplicate dataframe
all_standardized <- all_wider

# split by sample type to standardize within sample type

# standardize based on site reach
all_standardized[,c("log_ATX_all_ug_org_mat", "cobalamin_B12", "nitrification",
                    "nitrogen_fixation", "phosphatase_transporters", "pyridoxal", "thiamine")] <- 
  apply(all_standardized[,c("log_ATX_all_ug_org_mat", "cobalamin_B12", "nitrification",
                            "nitrogen_fixation", "phosphatase_transporters", "pyridoxal", "thiamine")], MARGIN = 2,
        function(x) ave(x, all_standardized$site_reach, all_standardized$sample_type, FUN = scale))
# note- have NaNs but they are only for Salmon River which is okay

# test on one reach & sample type to make sure that was done correctly
test <- all_wider %>% filter(site_reach == "SFE-M-1S" & sample_type == "TM")
test[,c("log_ATX_all_ug_org_mat", "cobalamin_B12", "nitrification",
        "nitrogen_fixation", "phosphatase_transporters", "pyridoxal", "thiamine")] <- 
  scale(test[,c("log_ATX_all_ug_org_mat", "cobalamin_B12", "nitrification",
                "nitrogen_fixation", "phosphatase_transporters", "pyridoxal", "thiamine")])
all.equal(test, all_standardized %>% filter(site_reach == "SFE-M-1S" & sample_type == "TM")) # we good!

# pivot standardized data back to longer
all_standardized_long <- all_standardized %>% 
  pivot_longer(cols = c("cobalamin_B12", "nitrification",
                        "nitrogen_fixation", "phosphatase_transporters", "pyridoxal", "thiamine"),
               names_to = "functional_grouping", values_to = "loggenes_standardized_by_reach") %>% 
  dplyr::rename(logatx_standardized_by_reach = log_ATX_all_ug_org_mat)

# remove other data to double-check that we aren't calling it
rm(model_info, all, sig_models, sig_model_data, subplot1_1, subplot1_2, subplot1_3,
   subplot2_1, subplot2_2)

# get linear model r-squared and p-values
model_info_stnd = data.frame(sample_type = NA,
                        functional_grouping = NA,
                        site = NA,
                        r_squared = NA,
                        p_value = NA)
for(i in c("SFE-M", "RUS")) {
  for(j in unique(all_standardized_long$sample_type)) {
    for(k in unique(all_standardized_long$functional_grouping)) {
      if(!(i == "RUS" & j == "TM")) {
        model = linear_model((all_standardized_long %>% filter(sample_type == j & functional_grouping == k & site == i)),
                             x = "loggenes_standardized_by_reach", y = "logatx_standardized_by_reach")
        model_info_stnd = rbind(model_info_stnd,
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
model_info_stnd <- left_join(model_info_stnd[-1,], all_standardized_long %>% 
                          dplyr::group_by(site, functional_grouping, sample_type) %>% 
                          dplyr::summarize(x_min = min(loggenes_standardized_by_reach),
                                           x_max = max(loggenes_standardized_by_reach),
                                           y_max = max(logatx_standardized_by_reach),
                                           y_min = max(logatx_standardized_by_reach)),
                        by = c("site", "sample_type", "functional_grouping"))

# make labels=
model_info_stnd <- model_info_stnd %>% 
  # different rounding for r-squared
  mutate(r_squared_mod = case_when(r_squared <0.01 ~ round(r_squared, 3),
                                   TRUE ~ round(r_squared, 2))) %>% 
  mutate(label = case_when(p_value > 0.05 ~ paste("*r<sup>2</sup>* = ", r_squared_mod, ", *p* = ", format(round(p_value, 2), nsmall = 2),
                                                  sep = ""),
                           TRUE ~ paste("*r<sup>2</sup>* = ", r_squared_mod, ", *p* = ", format(round(p_value, 2), nsmall = 2), "*",
                                        sep = "")))

# only plot lines when model is significant
sig_models_stnd <- model_info_stnd %>% 
  filter(p_value < 0.05) # only n-fix is significant

# dataset with only significant models
sig_model_data <- all_standardized_long %>% 
  filter((sample_type == "NT" & functional_grouping == "nitrogen_fixation" & site == "SFE-M"))

## SFE plots

# nontarget
subplot3_1 <- ggplot(data = all_standardized_long %>% filter(sample_type == "NT" & site == "SFE-M"), 
                     aes(x = loggenes_standardized_by_reach, y = logatx_standardized_by_reach)) +
  geom_smooth(data = sig_model_data %>% filter(sample_type == "NT", site == "SFE-M"),
              fill = "#f2b979", color = "#c26f11", method = "lm", alpha = 0.25) +
  facet_grid(~functional_grouping, scales = "free") + 
  geom_richtext(
    data = model_info_stnd %>% filter(sample_type == "NT" & site == "SFE-M"),
    mapping = aes(x = x_min + ((x_max - x_min) * .44), y = -2.6, label = label), size= 2.4,
    fill = NA, label.color = NA, text.color = "#8c4d06") +
  geom_point(aes(shape = site_reach), color = "#c26f11", size = 2, alpha = 0.5) +
  scale_x_continuous(breaks = breaks_pretty(n = 3)) +
  scale_y_continuous(breaks = breaks_pretty(n = 4)) + 
  theme(strip.background = element_blank(), legend.position = "bottom")
subplot3_1

subplot3_2 <- ggplot(data = all_standardized_long %>% filter(sample_type == "TAC" & site == "SFE-M"), 
                     aes(x = loggenes_standardized_by_reach, y = logatx_standardized_by_reach)) +
  geom_smooth(data = sig_model_data %>% filter(sample_type == "TAC", site == "SFE-M"),
              color = "#3f9633", fill = "#acf0a3", method = "lm", alpha = 0.25) +
  facet_grid(~functional_grouping, scales = "free") + 
  geom_point(aes(shape = site_reach), color = "#3f9633", size = 2, alpha = 0.5) +
  geom_richtext(
    data = model_info_stnd %>% filter(sample_type == "TAC" & site == "SFE-M"),
    mapping = aes(x = x_min + ((x_max - x_min) * .44), y = -1.5, label = label), size= 2.4,
    fill = NA, label.color = NA, text.color = "#247319") +
  scale_x_continuous(breaks = breaks_pretty(n = 3)) +
  scale_y_continuous(breaks = breaks_pretty(n = 4)) + 
  theme(strip.background = element_blank(), legend.position = "bottom")
subplot3_2

subplot3_3 <- ggplot(data = all_standardized_long %>% filter(sample_type == "TM" & site == "SFE-M"), 
                     aes(x = loggenes_standardized_by_reach, y = logatx_standardized_by_reach)) +
  geom_smooth(data = sig_model_data %>% filter(sample_type == "TM", site == "SFE-M"),
              color = "#6d4275", fill = "#c79bcf", method = "lm", alpha = 0.25) +
  facet_grid(~functional_grouping, scales = "free") + 
  geom_point(aes(shape = site_reach), color = "#6d4275", size = 2, alpha = 0.5) +
  geom_richtext(
    data = model_info_stnd %>% filter(sample_type == "TM" & site == "SFE-M"),
    mapping = aes(x = x_min + ((x_max - x_min) * .42), y = -1.9, label = label), size= 2.4,
    fill = NA, label.color = NA, text.color = "#48234f") +
  scale_x_continuous(breaks = breaks_pretty(n = 3)) +
  theme(strip.background = element_blank(), legend.position = "bottom")
subplot3_3

## Russian

subplot4_1 <- ggplot(data = all_standardized_long %>% filter(sample_type == "NT" & site == "RUS"), 
                     aes(x = loggenes_standardized_by_reach, y = logatx_standardized_by_reach)) +
  geom_smooth(data = sig_model_data %>% filter(sample_type == "NT" & site == "RUS"),
              fill = "#f2b979", color = "#c26f11", method = "lm", alpha = 0.25) +
  facet_grid(~functional_grouping, scales = "free") + 
  geom_richtext(
    data = model_info_stnd %>% filter(sample_type == "NT" & site == "RUS"),
    mapping = aes(x = x_min + ((x_max - x_min) * .44), y = -1.6, label = label), size= 2.4,
    fill = NA, label.color = NA, text.color = "#8c4d06") +
  geom_point(aes(shape = site_reach), color = "#c26f11", size = 2, alpha = 0.5) +
  scale_x_continuous(breaks = breaks_pretty(n = 3)) +
  theme(strip.background = element_blank(), legend.position = "bottom")
subplot4_1

subplot4_2 <- ggplot(data = all_standardized_long %>% filter(sample_type == "TAC" & site == "RUS"), 
                     aes(x = loggenes_standardized_by_reach, y = logatx_standardized_by_reach)) +
  geom_smooth(data = sig_model_data %>% filter(sample_type == "TAC" & site == "RUS"),
              color = "#3f9633", fill = "#acf0a3", method = "lm", alpha = 0.25) +
  facet_grid(~functional_grouping, scales = "free") + 
  geom_point(aes(shape = site_reach), color = "#3f9633", size = 2, alpha = 0.5) +
  geom_richtext(
    data = model_info_stnd %>% filter(sample_type == "TAC" & site == "RUS"),
    mapping = aes(x = x_min + ((x_max - x_min) * .44), y = -1.5, label = label), size= 2.4,
    fill = NA, label.color = NA, text.color = "#247319") +
  scale_x_continuous(breaks = breaks_pretty(n = 3)) +
  scale_y_continuous(breaks = breaks_pretty(n = 4)) + 
  theme(strip.background = element_blank(), legend.position = "bottom")
subplot4_2

#### (4) Putting Together Final Plots ####

# going with z-score'd version as two things are no longer significant in that version

# main figure for SFE!
plot_sfe <- plot_grid(subplot3_1 + theme(legend.position = "none"), 
                      subplot3_3 + theme(legend.position = "none"), 
                      subplot3_2 + theme(legend.position = "none"), align = "hv", ncol = 1)
plot_sfe

# save figure
ggsave("./figures/tiff_files/Q2_atx_v_gene.tiff", dpi = 500,
       width=17.5, height=12, unit="cm")

# sfig for RUS!
plot_rus <- plot_grid(subplot4_1 + theme(legend.position = "none"),
                      subplot4_2 + theme(legend.position = "none"),
                      align = "hv", ncol = 1)
plot_rus

# save figure !
ggsave("./figures/tiff_files/sfig_rus_atx_v_gene.tiff", dpi = 500,
       width=17.5, height=8, unit="cm")

# for legend in dark gray
legend_points <- ggplot(data = all %>% filter(sample_type == "TAC" & site == "SFE-M"), 
                        aes(x = log_predicted_gene_abundance, y = log_ATX_all_ug_org_mat)) +
  facet_grid(~functional_grouping, scales = "free") + 
  geom_point(aes(shape = site_reach), color = "black", size = 2, alpha = 0.5) +
  theme(strip.background = element_blank(), legend.position = "bottom")
legend_points

# save figure
ggsave("./figures/tiff_files/Q2_atx_v_gene_sitereach_legend.tiff", dpi = 500,
       width = 18, height = 4, unit = "cm")
