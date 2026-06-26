#### Supplemental figure comparing nif gene predictions to observed known n-fixing cyanobacteria
### Jordan Zabrecky
## 06.24.2026

# This script creates a supplemental figure comparing predicted nif gene abundances
# and relative abundance of nitrogen fixers in samples and assesses the correlation

#### (1) Loading libraries & data ####

# loading libraries
lapply(c("tidyverse", "ggtext"), require, character.only = T)

## (a) molecular data

# read in processed molecular data in long format
nt <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_FINAL.csv") %>% 
  filter(sample_type == "NT")
tm <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_TM_nomicro.csv")
tac <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_TAC_noanacyl.csv")

# add into list
molec_data <- list(nt, tm, tac)
names(molec_data) <- c("NT", "TM", "TAC")

# filter for known nitrogen fixers & pivot_wider
molec_data <- lapply(molec_data, function(x) x %>% 
                       filter(order == "Nostocales" | order == "Endosymbiotic Diazoplast") %>% 
                       select(site_reach, site, field_date, sample_type, order, picrust2_relative_abundance) %>% 
                       # need to gather from different observations before pivoting wider
                       dplyr::group_by(site_reach, site, field_date, sample_type, order) %>% 
                       dplyr::summarize(picrust2_relative_abundance = sum(picrust2_relative_abundance)) %>% 
                       pivot_wider(names_from = "order", values_from = "picrust2_relative_abundance"))

## (b) predicted abundance data

# load in nitrogen fixation data
nt <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select.csv") %>% 
  filter(sample_type == "NT")
tm <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tm_nomicro.csv")
tac <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tac_noanacyl.csv")

# put all dataframes into a list
func_select <- list(nt, tac, tm)
names(func_select) <- c("NT", "TAC", "TM")

# join different KEGG orthologs
func_select <- lapply(func_select, function(x) x %>% 
                        dplyr::group_by(site_reach, site, field_date, sample_type, functional_grouping) %>% 
                        dplyr::summarize(predicted_gene_abundance = sum(predicted_gene_abundance)))

# join in with molecular data
molec_data_combined <- lapply(names(data_select), function(x) {
  all = left_join(molec_data[[x]], func_select[[x]] %>% 
                    filter(functional_grouping == "nitrogen_fixation"), 
                  by = c("field_date", "sample_type", "site", "site_reach"))
})
names(molec_data_combined) <- c("NT", "TAC", "TM")

# joining all data together
final_data <- rbind(molec_data_combined$NT, molec_data_combined$TAC, molec_data_combined$TM) %>% 
  mutate(nfixers = Nostocales + `Endosymbiotic Diazoplast`)

#### (2) Making plots ####

# set universal plot theme
theme_set(theme_bw() + theme(strip.background = element_blank(),
                             plot.title = element_text(hjust = 0.5), legend.text = element_markdown(),
                             text = element_text(size = 8), strip.text = element_text(size = 8),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# make plots
ggplot(data = final_data, aes(x = nfixers, y = predicted_gene_abundance)) +
  geom_point(size = 3) +
  facet_wrap(sample_type~site, scale = "free")

cor((final_data %>% filter(site == "SAL" & sample_type == "NT"))$nfixers, 
    (final_data %>% filter(site == "SAL" & sample_type == "NT"))$predicted_gene_abundance,
    use = "na.or.complete")

?cor()
