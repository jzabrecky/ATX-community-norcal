#### Relative abundance of other ATX producers in TM and TAC samples
### Jordan Zabrecky
## last edited: 05.06.2026

# This code produces figures to answer Q4 asking what is the full suite of
# anatoxin associated taxa within target samples?

#### (1) Loading data ####

# source data from morphological script
source("./code/7a_otheratx_microscopy.R")

# rename 
microscopy_atx_taxa <- atx_taxa_only

# source data from molecular script
source("./code/7b_otheratx_molecular.R")

# rename for clarification
molecular_atx_taxa <- atx_taxa_only

# load additional libraries
lapply(c("ggtext", "cowplot"), require, character.only = T)

#### OLD BELOW

#### (1) Morphological data ####

# source data from morphological script
source("./code/4a_amongrivers_microscopy.R")

# load additional libraries
lapply(c("ggtext", "cowplot"), require, character.only = T)

#### (2) Make figure ####

# set universal plot theme
theme_set(theme_bw() + theme(strip.background = element_blank(),
                             axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
                             plot.title = element_text(hjust = 0.5), legend.text = element_markdown(),
                             text = element_text(size = 10), strip.text = element_text(size = 10)))

# custom palette
palette <- c("#FBF6B0", "#C5BD53", "#777122", "#C2DFFF", "#5E9DE0", "#205288", 
             "#C0ED96", "#7AB048", "#3D631A", "#CBC5F6", "#8A80CF", "#61389E")

## (a) microscopy

# plot microcoleus
present_morpho_tm <- ggplot(data = microscopy_atx_taxa$tm$percent_samples, aes(x = reorder(taxa, -percent_samples), y = percent_samples, fill = taxa)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none") +
  scale_fill_discrete(palette = c("#C0ED96", "#C5BD53", "#C2DFFF", "#7AB048", "#205288", "#CBC5F6", "#5E9DE0")) +
  labs(x = NULL, y = NULL) +
  scale_x_discrete(labels = c("Other Coccoidal<br>Cyanobacteria", "*Leptolyngbya*<br>and *Geitlerinema*", 
                               "*Anabaena* and<br>*Cylindrospermum*", "*Nostoc*", 
                                "Miscellaneous<br>Oscillatoriales", "*Oscillatoria*", "*Phormidium*")) +
  theme(axis.text.x = element_markdown(angle = 60, vjust = 1, hjust=1, size = 9)) +
  scale_y_continuous(limits = c(0,100))
present_morpho_tm

# plot anabaena
present_morpho_tac <- ggplot(data = microscopy_atx_taxa$tac$percent_samples, aes(x = reorder(taxa, -percent_samples), y = percent_samples, fill = taxa)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none") +
  scale_fill_discrete(palette = c("#C5BD53", "#777122", "#C2DFFF", "#7AB048", "#205288", "#CBC5F6", "#5E9DE0")) +
  labs(x = NULL, y = NULL) +
  scale_x_discrete(labels = c("*Geitlerinema* and<br>*Leptolyngbya*", "Other Coccoidal<br>Cyanobacteria", "*Microcoleus*", 
                              "*Oscillatoria*", "*Nostoc*", "*Phormidium*", "Miscellaneous<br>Oscillatoriales")) +
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
  scale_fill_discrete(palette = c("#C0ED96", "#3D631A", "#C0ED96", "#C0ED96", "#C5BD53", "#C5BD53", "#8A80CF", 
                                  "#CBC5F6", "#7AB048", "#205288", "#5E9DE0", "#61389E", "#FBF6B0", "#CBC5F6")) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_markdown(angle = 60, vjust = 1, hjust=1, size = 9)) +
  scale_y_continuous(limits = c(0,100))
present_micro_tm

# plot anabaena
present_micro_tac <- ggplot(data = molecular_atx_taxa$tac$percent_samples, aes(x = reorder(genus, -percent_samples), y = percent_samples, fill = genus)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none") +
  scale_fill_discrete(palette = c("#3D631A", "#C5BD53", "#C5BD53", "#8A80CF", "#777122", "#CBC5F6", "#7AB048", "#205288", 
                                  "#5E9DE0", "#61389E", "#FBF6B0", "#CBC5F6")) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_markdown(angle = 60, vjust = 1, hjust=1, size = 9)) +
  scale_y_continuous(limits = c(0,100))
present_micro_tac

#### (3) Putting Figures Together ####

# put together plots
final <- plot_grid(present_morpho_tm, present_morpho_tac, present_micro_tm, present_micro_tac,
                   align = "hv", ncol = 2)
final

# save!
