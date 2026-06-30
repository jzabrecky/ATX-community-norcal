#### Supplemental figure of boxplots of ATX concentrations
### Jordan Zabrecky
## last edited: 06.26.2026

# This script creates a supplemental figure showing the range of anatoxin 
# concentrations observed in each river

#### (1) Loading libraries & data ####

# loading libraries
lapply(c("tidyverse", "ggtext"), require, character.only = T)

# loading data
atx <- read.csv("./data/field_and_lab/atx_w_categorical_groupings.csv") %>% 
  # remove NT samples
  filter(sample_type != "NT")

#### (2) Summary of samples ####

# make summary of samples with detected anatoxin
atx_sum <- atx %>% 
  dplyr::group_by(site) %>% 
  dplyr::summarize(total = length(ATX_all_ug_org_mat),
                   detections = length(which(ATX_all_ug_org_mat > 0)))

# split based on site for easy indexing
atx_sum <- split(atx_sum, atx_sum$site)

#### (3) Making plot ####

# set universal plot theme
theme_set(theme_bw() + theme(strip.background = element_blank(),
                             plot.title = element_text(hjust = 0.5), legend.text = element_markdown(),
                             text = element_text(size = 10), strip.text = element_text(size = 10),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# set seed for consistent jitter :)
set.seed(6)
plot = ggplot(data = atx, aes(x = site, y = ATX_all_ug_org_mat, fill = site)) +
  stat_summary(fun = median, show.legend = FALSE, geom = "crossbar", color = "#424242") +
  geom_jitter(alpha = 0.5, aes(fill = site, shape = sample_type), color = "black", size = 2) + 
  scale_fill_manual(values = c("SAL" = "#81bbfc",
                               "SFE-M" = "#416f16",
                               "RUS" = "#ab9f00")) +
  scale_x_discrete(labels = c(paste("Russian<br>River<br><br>(N=", atx_sum$RUS$detections, " of ",
                                    atx_sum$RUS$total, ")", sep = ""), 
                              paste("Salmon<br>River<br><br>(N=", atx_sum$SAL$detections, " of ",
                                    atx_sum$SAL$total, ")", sep = ""),
                              paste("South Fork<br>Eel River<br><br>(N=", atx_sum$`SFE-M`$detections, " of ",
                                    atx_sum$`SFE-M`$total, ")", sep = ""))) +
  scale_shape_manual(values = c(22, 23)) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(1, 10, 100, 200, 300)) +
  labs(x = NULL, y = "&mu;g anatoxins g <sup>-1</sup> OM") +
  theme(axis.text.x = element_markdown(size = 9, color = "#333333"), axis.title.y = element_markdown()) + 
  theme(legend.position = "none")
plot

# save
ggsave("./figures/tiff_files/sfig_atx_rivers.tiff", dpi = 600,
       width = 11, height = 9, unit="cm")

# get legend
plot = ggplot(data = atx, aes(x = site, y = ATX_all_ug_org_mat, fill = site)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(alpha = 0.5, aes(fill = site, shape = sample_type), color = "black", size = 2) + 
  scale_fill_manual(values = c("SAL" = "#81bbfc",
                               "SFE-M" = "#416f16",
                               "RUS" = "#ab9f00")) +
  scale_x_discrete(labels = c("Russian River", "Salmon River", "South Fork<br>Eel River")) +
  scale_shape_manual(values = c(22, 23)) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  labs(x = NULL, y = expression(paste(mu, "g anatoxins per g OM"), sep = "")) +
  theme(legend.position = "right")
plot

# save
ggsave("./figures/tiff_files/sfig_atx_legend.tiff", dpi = 600,
       width = 11, height = 9, unit="cm")
