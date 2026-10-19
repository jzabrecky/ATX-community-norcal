#### Relative abundance of other ATX producers in TM and TAC samples
### Jordan Zabrecky
## last edited: 10.16.2026

# This code produces figures to answer Q3 asking what is the full suite of
# potential anatoxin-producing taxa within target samples?

#### (1) Loading data ####

# source data from morphological script
source("./code/6a_otheratx_microscopy.R")

# rename 
microscopy_atx_taxa <- atx_taxa_only

# source data from molecular script
source("./code/6b_otheratx_molecular.R")

# rename for clarification
molecular_atx_taxa <- atx_taxa_only

# load additional libraries
lapply(c("ggtext", "cowplot"), require, character.only = T)

#### (2) Make figure ####

# set universal plot theme
theme_set(theme_bw() + theme(strip.background = element_blank(),
                             axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
                             plot.title = element_text(hjust = 0.5), legend.text = element_markdown(),
                             text = element_text(size = 10), strip.text = element_text(size = 10),))

# custom palette
palette <- c("#FBF6B0", "#C5BD53", "#777122", "#C2DFFF", "#5E9DE0", "#205288", 
             "#C0ED96", "#7AB048", "#3D631A", "#CBC5F6", "#8A80CF", "#61389E")

## (a) microscopy

# create taxa labels
microscopy_atx_taxa$tm$percent_samples <- microscopy_atx_taxa$tm$percent_samples %>% 
  mutate(taxa_label = case_when(taxa == "anabaena_and_cylindrospermum" ~ "*Anabaena* /<br>*Cylindrospermum*",
                                taxa == "leptolyngbya_and_geitlerinema" ~ "*Leptolyngbya* /<br>*Geitlerinema*",
                                taxa == "miscellaneous_oscillatoriales" ~ "Miscellaneous<br>Oscillatoriales",
                                taxa == "oscillatoria" ~ "*Oscillatoria*",
                                taxa == "phormidium_unknown" ~ "*Phormidium*"))

# plot microcoleus
present_morpho_tm <- ggplot(data = microscopy_atx_taxa$tm$percent_samples, aes(x = reorder(taxa_label, -percent_samples), y = percent_samples, fill = taxa)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none") +
  scale_fill_discrete(palette = c(palette[7], palette[2], palette[4], palette[6], palette[5])) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_markdown(angle = 60, vjust = 1, hjust=1, size = 9)) +
  scale_y_continuous(limits = c(0,100))
present_morpho_tm

# create taxa labels
microscopy_atx_taxa$tac$percent_samples <- microscopy_atx_taxa$tac$percent_samples %>% 
  mutate(taxa_label = case_when(taxa == "leptolyngbya_and_geitlerinema" ~ "*Leptolyngbya* /<br>*Geitlerinema*",
                                taxa == "miscellaneous_oscillatoriales" ~ "Miscellaneous<br>Oscillatoriales",
                                taxa == "oscillatoria" ~ "*Oscillatoria*",
                                taxa == "phormidium_unknown" ~ "*Phormidium*",
                                taxa == "microcoleus" ~ "*Microcoleus*"))

# plot anabaena
present_morpho_tac <- ggplot(data = microscopy_atx_taxa$tac$percent_samples, aes(x = reorder(taxa_label, -percent_samples), y = percent_samples, fill = taxa)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none") +
  scale_fill_discrete(palette = c(palette[2], palette[3], palette[4], palette[6], palette[5])) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_markdown(angle = 60, vjust = 1, hjust=1, size = 9)) +
  scale_y_continuous(limits = c(0,100))
present_morpho_tac

## (b) molecular

# mutate for italicized genus
molecular_atx_taxa$tac$percent_samples <- molecular_atx_taxa$tac$percent_samples %>% 
  mutate(genus_italicized = paste("*", genus, "*", sep = ""))
molecular_atx_taxa$tm$percent_samples <- molecular_atx_taxa$tm$percent_samples %>% 
  mutate(genus_italicized = paste("*", genus, "*", sep = ""))

# plot microcoleus
present_micro_tm <- ggplot(data = molecular_atx_taxa$tm$percent_samples, aes(x = reorder(genus_italicized, -percent_samples), y = percent_samples, fill = genus)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none") +
  scale_fill_discrete(palette = c(palette[7], palette[8], palette[9], palette[7], palette[2], palette[6], palette[5], 
                                  palette[12], palette[1], palette[10])) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_markdown(angle = 60, vjust = 1, hjust=1, size = 9)) +
  scale_y_continuous(limits = c(0,100))
present_micro_tm

# plot anabaena
present_micro_tac <- ggplot(data = molecular_atx_taxa$tac$percent_samples, aes(x = reorder(genus_italicized, -percent_samples), y = percent_samples, fill = genus)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none") +
  scale_fill_discrete(palette = c(palette[8], palette[9], palette[2], palette[11], palette[3], palette[6], palette[5], palette[12], 
                                  palette[1], palette[10])) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_markdown(angle = 60, vjust = 1, hjust=1, size = 9)) +
  scale_y_continuous(limits = c(0,100))
present_micro_tac

#### (3) Putting Figures Together ####

# put together plots
final <- plot_grid(present_morpho_tm, present_morpho_tac, present_micro_tm, present_micro_tac,
                   align = "hv", ncol = 2) +
  theme(plot.background = element_rect(fill = "white", color = "white"))
final

# save
ggsave("./figures/tiff_files/Q3_atx_taxa.tiff", dpi = 500,
       width=17, height=13, unit="cm")
