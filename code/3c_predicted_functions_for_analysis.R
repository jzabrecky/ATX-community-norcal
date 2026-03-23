#### Further selecting functional groups for samples
### Jordan Zabrecky
## last edited: 03.23.2026

# This script processes the further processes the predicted KEGG orthologs 
# from the previous script and only keeps orthologs in processes we care about
# Then, it matches those orthologs to KEGG pathways

#### (1) Loading data & libraries ####

# loading libraries
lapply(c("tidyverse", "ggpicrust2"), require, character.only = T)

# loading data
raw_data <- read.csv("./data/molecular/PICRUSt2_predicted_KO_all.csv")

#### (2) Filtering for Functions ####

# Functional groups we care about that may or may not be important to anatoxin production"
# nitrogen fixation N2 -> NH3  (many papers) 
nfixers <- c("K02588", "K02586", "K02591", "K00531", "K22896", "K22897", "K22898", "K22899")
# nitrification NH3 -> NO2/NO3 (many papers)
nitrification <- c("K28504", "K10945", "K10946", "K10535")
# Phosphorus Pathways- TBD w/ Raina but will start with 00440
# Thiamine synthesis [pathway 00730] (Tee et al. 2021; Junier et al. 2024)
# specifically TenA [orthologs K03707], thiD [ortholog K00941], thiK [ortholog K07251], 
# (Tee et al. 2021, Junier et al. 2024)
thiamine <- c("K03707", "K00941", "K07251")
# Cobalamin, cbiG (B12) synthesis [orthologs K13541,K02189] (Junier et al. 2024)
cobalamin <- c("K13541", "K02189")
# Pyridoxal, pdxH (B6) synthesis [orthologs K00275] (Junier et al. 2024)
pyridoxal <- c("K00275")

# put all in a list
important_orthologs <- c(nfixers, nitrification, thiamine, cobalimin, pyridoxal)

# filter out from reference provided by "ggpicrust2"
ref <- ko_to_kegg_reference %>% 
  filter(ko_id %in% important_orthologs)

# determine if all ko's in ref are unique
length(unique(ref$ko_id)) #  7 not unique

# some double matches with nitrogen pathways here need to remove
ref <- ref %>% 
  filter(!pathway_number %in% c("00625", "00680", "03000"))

# filter data
data <- raw_data %>% 
  filter(ko_id %in% ref$ko_id) %>%  # this reduces the original df to 0.15%
  left_join(ref, by = "ko_id")

#### (3) Add in More Detail & Save ####

# add in more detailed groups
data <- data %>% 
  mutate(my_grouping = case_when(ko_id %in% cobalimin ~ "cobalamin (B12)",
                                 ko_id %in% nfixers ~ "nitrogen fixation",
                                 ko_id %in% nitrification ~ "nitrification", 
                                 ko_id %in% pyridoxal ~ "pyridoxal",
                                 ko_id %in% thiamine ~ "thiamine"))

# fix Russian date (sampling was supposed to be on 7/6 but issues got it pushed to 7/7)
data$field_date[which(data$field_date == "7/7/2022")] <- "7/6/2022"

# save
write.csv(data, "./data/molecular/PICRUSt2_predicted_KO_select.csv")

#### (4) Preliminary Investigation

# summary 
summary <- data %>% 
  dplyr::group_by(site, sample_type, my_grouping) %>% 
  dplyr::summarize(mean_predicted_gene_abundance = mean(predicted_gene_abundance)) %>% 
  ungroup()

# prelim comparisons by site across sample type
ggplot(summary, aes(y = mean_predicted_gene_abundance, x = my_grouping, fill = site)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~sample_type) + 
  scale_x_discrete(guide = guide_axis(angle = 60))
 
# temporal summary
temporal_summary <- data %>% 
  dplyr::group_by(site, sample_type, field_date, my_grouping) %>% 
  dplyr::summarize(mean_predicted_gene_abundance = mean(predicted_gene_abundance)) %>% 
  ungroup()

# among time
ggplot(temporal_summary, aes(y = mean_predicted_gene_abundance, x = field_date, fill = site)) +
  geom_bar(stat = "identity") +
  facet_grid(my_grouping~sample_type, scale = "free_y") + 
  scale_x_discrete(guide = guide_axis(angle = 60))
