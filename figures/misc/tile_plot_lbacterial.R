# test tile plots using bacterial classes

#### MOST ABUNDANT ####

# we will use the NT data, and TM excluding Microcoleus, and TAC excluding Anabaena and Green Algae
nt <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_FINAL.csv") %>% 
  filter(sample_type == "NT") %>% 
  mutate(field_date = mdy(field_date),
         phylum_class = paste(phylum, " - ", class))
tm <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_TM_nomicro.csv") %>%
  mutate(field_date = mdy(field_date),
         phylum_class = paste(phylum, " - ", class))
tac <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_TAC_noanacyl.csv") %>% 
  mutate(field_date = mdy(field_date),
         phylum_class = paste(phylum, " - ", class))

# figure out most abundant classes
abun <- nt %>% 
  dplyr::group_by(phylum_class) %>% 
  dplyr::summarize(mean = mean(picrust2_relative_abundance))
# Cyanobacteria, Alphaproteobacteria, Verrucomicrobiota, Gammaproteobacteria, Armatimonadia, Deinococci, Desulfobulbia
abun <- tm %>% 
  dplyr::group_by(phylum_class) %>% 
  dplyr::summarize(mean = mean(picrust2_relative_abundance))
# Gemmatimonadetes, Bacteroidia, Verrucomicrobiae
abun <- tac %>% 
  dplyr::group_by(phylum_class) %>% 
  dplyr::summarize(mean = mean(picrust2_relative_abundance))
# Actinobacteria

classes_we_care_about <- c("Cyanobacteria  -  Cyanobacteriia", "Gemmatimonadota  -  Gemmatimonadetes",
                           "Bacteroidota  -  Bacteroidia", "Verrucomicrobiota  -  Verrucomicrobiae",
                           "Proteobacteria  -  Gammaproteobacteria", "Actinobacteriota  -  Actinobacteria",
                           "Proteobacteria  -  Alphaproteobacteria")

#### HIGHLIGHTED FROM TEMPORAL ANALYSIS ####
data <- rbind(nt, tac, tm) %>% 
  select(field_date, sample_type, site_reach, site, phylum_class, qiime2_relative_abundance) %>% 
  filter(phylum_class %in% classes_we_care_about)

# summarize
data <- data %>% 
  group_by(site, field_date, sample_type, phylum_class) %>% 
  dplyr::summarize(rel_abundance = mean(qiime2_relative_abundance))

# split data based on site
data_site <- split(data, data$site)

# add in grouping factor

# make plots
ggplot(data_site$`RUS`, aes(x = field_date, y = sample_type, fill = grouping2_factor)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  facet_wrap(~taxa, ncol = 1) +
  geom_text(aes(label = round(relative_abundance, 1))) +      # add text on top of tile
  coord_fixed() +
  scale_fill_discrete(palette = c("#206349", "#368f6c", "#5aa185", "#97ccb7", "#bae6d4", "#ccf0e1", "#b0b0b0"))