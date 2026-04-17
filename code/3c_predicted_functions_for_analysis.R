#### Further selecting functional groups for samples
### Jordan Zabrecky
## last edited: 04.08.2026

# This script processes the further processes the predicted KEGG orthologs 
# from the previous script and only keeps orthologs in processes we care about
# Then, it matches those orthologs to KEGG pathways

#### (1) Loading data & libraries ####

# loading libraries
lapply(c("tidyverse", "ggpicrust2"), require, character.only = T)

# loading data (all and then target taxa w/o respective taxa ASV reads!)
raw_data <- read.csv("./data/molecular/PICRUSt2_predicted_KO_all_all.csv")
raw_data_tm_nomicro <-read.csv("./data/molecular/PICRUSt2_predicted_KO_all_tm_nomicro.csv")
raw_data_tac_noanacyl <-read.csv("./data/molecular/PICRUSt2_predicted_KO_all_tac_noanacyl.csv")

# put them in list
data_list <- list(raw_data, raw_data_tm_nomicro, raw_data_tac_noanacyl)
names(data_list) <- c("all", "tm_nomicro", "tac_noanacyl")

#### (2) Filtering for Functions ####

# Functional groups we care about that may or may not be important to anatoxin production"
# nitrogen fixation N2 -> NH3  (many papers) 
nfixers <- c("K02588", "K02586", "K02591", "K00531", "K22896", "K22897", "K22898", "K22899")
# nitrification NH3 -> NO2/NO3 (many papers)
nitrification <- c("K28504", "K10945", "K10946", "K10535")
# Phosphorus Pathways- phoD and phoA (decided w/ Raina)
phosphatase_transporters <- c("K01113", "K01077")
# Thiamine synthesis [pathway 00730] (Tee et al. 2021; Junier et al. 2024)
# specifically TenA [orthologs K03707], thiD [ortholog K00941], thiK [ortholog K07251], 
# (Tee et al. 2021, Junier et al. 2024)
thiamine <- c("K03707", "K00941", "K07251")
# Cobalamin, cbiG (B12) synthesis [orthologs K13541,K02189] (Junier et al. 2024)
cobalamin <- c("K13541", "K02189")
# Pyridoxal, pdxH (B6) synthesis [orthologs K00275] (Junier et al. 2024)
pyridoxal <- c("K00275")

# put all in a list
important_orthologs <- c(nfixers, nitrification, phosphatase_transporters, thiamine, cobalamin, pyridoxal)

# filter out from reference provided by "ggpicrust2"
ref <- ko_to_kegg_reference %>% 
  filter(ko_id %in% important_orthologs)

# determine if all ko's in ref are unique
length(unique(ref$ko_id)) #  20 not unique

# some double matches with nitrogen pathways here need to remove
# and for phosphatase transporters
ref <- ref %>% 
  filter(!pathway_number %in% c("00625", "00680", "03000", "00790", "00537", "04147")) %>%
  filter(!(pathway_number == "00730" & ko_id == "K01077"))

# filter data
data_list <- lapply(data_list, function(x) {
  data = x %>% 
    filter(ko_id %in% ref$ko_id) %>%  # this reduces the original df to 0.15%
    left_join(ref, by = "ko_id")
  })

#### (3) Add in More Detail & Save ####

# add in more detailed groups & sum genes in each group
data_list <- lapply(data_list, function(x) {
  data = x %>% 
    mutate(my_grouping = case_when(ko_id %in% cobalamin ~ "cobalamin_B12",
                                   ko_id %in% nfixers ~ "nitrogen_fixation",
                                   ko_id %in% nitrification ~ "nitrification", 
                                   ko_id %in% pyridoxal ~ "pyridoxal",
                                   ko_id %in% thiamine ~ "thiamine",
                                   ko_id %in% phosphatase_transporters ~ "phosphatase_transporters")) %>% 
    dplyr::group_by(site, site_reach, field_date, sample_type, my_grouping) %>% 
    dplyr::summarize(predicted_gene_abundance = sum(predicted_gene_abundance))

  # also fix Russian date (sampling was supposed to be on 7/6 but issues got it pushed to 7/7)
  data$field_date[which(data$field_date == "7/7/2022")] <- "7/6/2022"
  
  return(data)
})

# save
write.csv(data_list$all, "./data/molecular/PICRUSt2_predicted_KO_select.csv", row.names = FALSE)
write.csv(data_list$tm_nomicro, "./data/molecular/PICRUSt2_predicted_KO_select_tm_nomicro.csv", row.names = FALSE)
write.csv(data_list$tac_noanacyl, "./data/molecular/PICRUSt2_predicted_KO_select_tac_noanacyl.csv", row.names = FALSE)

#### (4) Preliminary Investigation

# summary 
summary <- lapply(data_list, function(x) {
  data = x %>% 
    dplyr::group_by(site, sample_type, my_grouping) %>% 
    dplyr::summarize(mean_predicted_gene_abundance = mean(predicted_gene_abundance)) %>% 
    ungroup()
  return(data)
})

# prelim comparisons by site across sample type
lapply(summary, function(x) {
  plot = ggplot(x, aes(y = mean_predicted_gene_abundance, x = my_grouping, fill = site)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~sample_type, scale = "free") + 
    scale_x_discrete(guide = guide_axis(angle = 60))
  print(plot)
})

# temporal summary
temporal_summary <- lapply(data_list, function(x) {
  data = x %>% 
    dplyr::group_by(site, sample_type, field_date, my_grouping) %>% 
    dplyr::summarize(mean_predicted_gene_abundance = mean(predicted_gene_abundance)) %>% 
    ungroup()
  return(data)
})
 
# among time
lapply(temporal_summary, function(x) {
  plot = ggplot(x, aes(y = mean_predicted_gene_abundance, x = field_date, fill = site)) +
    geom_bar(stat = "identity") +
    facet_grid(my_grouping~sample_type, scale = "free_y") + 
    scale_x_discrete(guide = guide_axis(angle = 60))
  print(plot)
})
