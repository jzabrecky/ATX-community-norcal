#### Relative abundance of other ATX producers in TM and TAC samples
### Jordan Zabrecky
## last edited: 02.18.2026

## ~ to be written - note that I may shove other ATX producers to Q3 code scripts

#### (1) Morphological data ####

# source data from morphological script
source("./code/4a_amongrivers_microscopy.R")

# load additional libraries
lapply(c("ggtext", "cowplot"), require, character.only = T)

# set universal plot theme
theme_set(theme_bw() + theme(strip.background = element_blank(),
                             axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
                             plot.title = element_text(hjust = 0.5), legend.text = element_markdown(),
                             text = element_text(size = 10), strip.text = element_text(size = 10)))

# custom palette
palette <- c("#FBF6B0", "#C5BD53", "#777122", "#C2DFFF", "#5E9DE0", "#205288", 
             "#C0ED96", "#7AB048", "#3D631A", "#CBC5F6", "#8A80CF", "#61389E")

## (a) showing relative abundance of individual samples

# put TAC and TM together
morpho_bothtaxa <- bind_rows(atx_taxa_only$tm[[1]], atx_taxa_only$tac[[1]])

# add factor for custom rows
morpho_bothtaxa <- morpho_bothtaxa %>% 
  mutate(site_factor = factor(site, levels = c("SFE-M", "SAL", "RUS")),
         type_factor = factor(sample_type, levels = c("TM", "TAC")),
         taxa_factor = factor(taxa, levels = c(str_sort(unique(morpho_bothtaxa$taxa))[-which(str_sort(unique(morpho_bothtaxa$taxa)) == "other_coccoids")],
                                               "other_coccoids")))

# will merge these together in InkScape
morpho1 <- ggplot(data = morpho_bothtaxa, 
       aes(x = sample_name, y = percent / 100, fill = taxa_factor)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = NULL, y = "Relative Abundance") +
  facet_grid(type_factor~site_factor, scales = "free_x") + 
  scale_fill_discrete(name = NULL, 
                      palette = palette, 
                      labels = c("*Anabaena* or *Cylindrospermum*", "*Aphanothece*", "*Geitlerinema*",
                                 "*Leptolyngbya*", "*Lyngbya*", "*Microcoleus*", "*Nostoc*", 
                                 "*Oscillatoria*", "*Phormidium*", "Other Coccoidal Cyanobacteria"))
morpho1

morpho2 <- ggplot(data = morpho_bothtaxa, 
                  aes(x = sample_name, y = percent / 100, fill = taxa_factor)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = NULL, y = "Relative Abundance") +
  facet_grid(sample_type~site_factor, scales = "free_x") + 
  scale_fill_discrete(name = NULL, 
                      palette = palette, 
                      labels = c("*Anabaena* or *Cylindrospermum*", "*Aphanothece*", "*Geitlerinema*",
                                 "*Leptolyngbya*", "*Lyngbya*", "*Microcoleus*", "*Nostoc*", 
                                 "*Oscillatoria*", "*Phormidium*", "Other Coccoidal Cyanobacteria"))
morpho2

## (b) showing percent of samples with each taxa present

percent_present <- lapply(atx_taxa_only, function(x) {
  # access dataframe
  y = x[[1]]
  
  # get number of samples
  n_samples = length(unique(y$sample_name))
  
  # calculate number of observations for each taxa
  y = y %>% 
    # remove not present taxa
    filter(percent != 0) %>% 
    dplyr::group_by(taxa) %>% 
    dplyr::summarize(count = length(taxa)) %>% 
    ungroup() %>% 
    mutate(percent_samples = count / n_samples)
  
  return(y)
})

# TM barplot

# plot
present_morpho_tm <- ggplot(data = percent_present$tm, aes(x = reorder(taxa, -percent_samples), y = percent_samples, fill = taxa)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none") +
  scale_fill_discrete(palette = palette) +
  labs(x = NULL, y = "Percent of Samples Present") +
  scale_x_discrete(labels = c("Other Coccoidal<br>Cyanobacteria", "*Geitlerinema*", "*Anabaena* and<br>*Cylindrospermum*",
                              "*Nostoc*", "*Leptolyngbya*", "*Lyngbya*", "*Oscillatoria*", "*Phormidium*")) +
  theme(axis.text.x = element_markdown(angle = 60, vjust = 1, hjust=1, size = 9))
present_morpho_tm

# TAC barplot
present_morpho_tac <- ggplot(data = percent_present$tac, aes(x = reorder(taxa, -percent_samples), y = percent_samples, fill = taxa)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none") +
  scale_fill_discrete(palette = palette[c(2, 3, 4, 9, 5, 6, 7, 8)]) +
  labs(x = NULL, y = "Percent of Samples Present") +
  scale_x_discrete(labels = c("*Geitlerinema*", "Other Coccoidal<br>Cyanobacteria", "*Microcoleus*", 
                              "*Oscillatoria*", "*Nostoc*", "*Phormidium*", "*Lyngbya*", "*Leptolyngbya*")) +
  theme(axis.text.x = element_markdown(angle = 60, vjust = 1, hjust=1, size = 9))
present_morpho_tac

#### (2) Molecular Data ####

# source from molecular script
source("./code/4b_amongrivers_16s.R")

# need to reset theme-ing
# set universal plot theme
theme_set(theme_bw() + theme(strip.background = element_blank(),
                             axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
                             plot.title = element_text(hjust = 0.5), legend.text = element_markdown(),
                             text = element_text(size = 10), strip.text = element_text(size = 10)))

# want to merge taxa with least amounts together for a total of 12 taxa only
genus_amounts <- bind_rows(atx_taxa_only$tm[[1]], atx_taxa_only$tac[[1]]) %>% 
  dplyr::group_by(order_genus) %>% 
  dplyr::summarize(mean = mean(relative_abundance))
view(genus_amounts)
# merge microcystis and synechocystis as they are both coccoidal bacteria
# merge pseudanabaena and limnothrix as they are same order
# merge oscillatoria, phormidium, and plankthrix
molec_groupped <- lapply(atx_taxa_only, function(x) {
  # get dataframe
  y = x[[1]] %>%
    # need to group ASVs of same genus
    group_by(sample_name, site, genus) %>% 
    dplyr::summarize(relative_abundance = sum(relative_abundance)) %>% 
    # need to add microsys
    pivot_wider(names_from = "genus", values_from = "relative_abundance", 
                values_fill = 0) %>% 
    mutate(microcystis_and_synechocystis = Microcystis + Synechocystis,
           pseudanabaena_and_limnothrix = Pseudanabaena + Limnothrix,
           oscillatoria_phormidium_planktothrix = Oscillatoria + Phormidium + Planktothrix) %>% 
    select(!c(Microcystis, Synechocystis, Pseudanabaena, Limnothrix, 
              Oscillatoria, Phormidium, Planktothrix)) %>% 
    pivot_longer(!c("sample_name", "site"), names_to = "group", values_to = "relative_abundance")
}) 

# put TM and TAC on one dataframe
molec_bothtaxa <- rbind(molec_groupped$tm %>% mutate(sample_type = "TM"),
                       molec_groupped$tac %>% mutate(sample_type = "TAC"))

molecular1 <- ggplot(data = molec_bothtaxa, 
                     aes(x = sample_name, y = relative_abundance / 100, fill = group)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = NULL, y = "Relative Abundance") +
  facet_grid(sample_type~site, scales = "free") + 
  scale_fill_discrete(name = NULL)
molecular1

## REVISIT THIS. Hitchiker's Guide states that most assignments are wrong