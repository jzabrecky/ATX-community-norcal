#### Main figure to show differences in morphologically-identified assemblages among rivers
### Jordan Zabrecky
## last edited: 06.02.2026

# This script creates a main figure to show differences in assemblages associated 
# with nontarget, microcoleus, and anabaena samples

#### (1) Loading libraries & data ####

# load from analysis script
source("./code/4b_amongrivers_16s.R")

# load additional libraries
lapply(c("cowplot", "ggtext"), require, character.only = T)

# set universal plot theme
theme_set(theme_bw() + theme(strip.background = element_blank(),
                             plot.title = element_text(hjust = 0.5), legend.text = element_markdown(),
                             text = element_text(size = 8), strip.text = element_text(size = 8),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),))

# custom palette
palette <- c("#FBF6B0", "#C5BD53", "#777122", "#C2DFFF", "#5E9DE0", "#205288", 
             "#C0ED96", "#7AB048", "#3D631A", "#CBC5F6", "#8A80CF", "#61389E")

# color for other or unknown
end_color <- "lightgray"

#### (2) Creating Individual Plots ####

# creating string vector to iterate through sample types
sample_types = c("nt", "tac", "tm")

## (a) PCoA

fig_a <- list()
for(i in sample_types) {
  fig_a[[i]] = makePCoAplot(PCoA_list[[i]], color = "site", shape = "month") +
    theme(legend.position = "none")
}
lapply(fig_a, print)

## (b) classes bar plot

# putting together data so they share the same ~12 Phylum-Classes for the figure!
all <- rbind(data_long$nt, data_long$tm, data_long$tac)

# use function 
source("./code/supplemental_code/S4c_grouping_func.R")
all_phylumclass <- microbial_grouping(all, "phylum_class", 0.065) %>% 
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(broader))[-which(str_sort(unique(broader)) == "Other")],
                                            "Other")),
         site_label = case_when(site == "SFE-M" ~ "South Fork<br>Eel River",
                                site == "SAL" ~ "Salmon<br>River",
                                site == "RUS" ~ "Russian<br>River"))

# split dataframe back up to plot separately
all_split <- split(all_phylumclass, all_phylumclass$sample_type)

## make figures
fig_b <- list()
for(i in sample_types) {
  fig_b[[i]] = barplot(all_split[[str_to_upper(i)]], x = "site_label", y = "relative_abundance",
                       fill = "broader_factor") +
    scale_fill_discrete(palette = c(palette[-length(palette)], end_color)) +
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_markdown(size = 7, color = "#333333"))
}
lapply(fig_b, print)

## okay, real quick let's see a phylum version

# putting together data so they share the same ~12 Phylum-Classes for the figure!
all_phylum <- microbial_grouping(all, "phylum", 0.05) %>% 
  # add in case if NA
  mutate(broader = case_when(is.na(phylum) ~ "Other",
                                    TRUE ~ broader),
         site_label = case_when(site == "SFE-M" ~ "South Fork<br>Eel River",
                               site == "SAL" ~ "Salmon<br>River",
                               site == "RUS" ~ "Russian<br>River")) %>% 
  # then factor
  mutate(broader_factor = factor(broader,
                                 levels = c(str_sort(unique(broader))[-which(str_sort(unique(broader)) == "Other")],
                                            "Other")))

# split dataframe back up to plot separately
all_p_split <- split(all_phylum, all_phylum$sample_type)

## make figures
fig_b_ver2 <- list()
for(i in sample_types) {
  fig_b_ver2[[i]] = barplot(all_p_split[[str_to_upper(i)]], x = "site_label", y = "relative_abundance",
                       fill = "broader_factor") +
    scale_fill_discrete(palette = c(palette[-length(palette)], end_color)) +
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_markdown(size = 7, color = "#333333"))
  #theme(legend.position = "none")
}
lapply(fig_b_ver2, print)

# just to visually compare between the two
for(i in sample_types) {
  print(plot_grid(fig_b[[i]], fig_b_ver2[[i]]))
}

## (c) diversity box plot

# add site label to diversity dataframe
diversity <- lapply(diversity, function(x) x %>% 
                      mutate(site_label = case_when(site == "SFE-M" ~ "South Fork<br>Eel River",
                                                    site == "SAL" ~ "Salmon<br>River",
                                                    site == "RUS" ~ "Russian<br>River")))

fig_c <- list()
for(i in sample_types) {
  fig_c[[i]] = ggplot(data = diversity[[i]], aes(x = site_label, y = shannon_diversity, fill = site)) +
    geom_boxplot(lwd = 0.3) +
    scale_fill_manual(values = c("SAL" = "#81bbfc",
                                 "SFE-M" = "#416f16",
                                 "RUS" = "#ab9f00")) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none", axis.text.x = element_markdown(size = 7, color = "#333333"))
}
lapply(fig_c, print)

#### (3) Putting Figure Together ####

# put figure together
final <-  plot_grid(fig_b_ver2$nt + theme(legend.position = "none"), 
                    fig_b_ver2$tm + theme(legend.position = "none"),
                    fig_b_ver2$tac + theme(legend.position = "none"),
                    fig_a$nt, fig_a$tm, fig_a$tac,
                    fig_c$nt, fig_c$tm, fig_c$tac,  align = "hv")
final

# save
ggsave("./figures/tiff_files/Q1_molecular.tiff", dpi = 600,
       width=18, height=16, unit="cm")
