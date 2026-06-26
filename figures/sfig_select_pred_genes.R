#### Supplemental figure of boxplots of predicted gene abundances for select functional groups
### Jordan Zabrecky
## last edited: 06.25.2026

# This script creates a supplemental figure showing the range of predicted
# gene abundances for select genes among rivers

#### (1) Loading Libraries & Data ####

# loading libraries
lapply(c("tidyverse", "ggtext"), require, character.only = T)

# load data (for selected functional groups)
nt <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select.csv") %>% 
  filter(sample_type == "NT")
tm <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tm_nomicro.csv")
tac <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tac_noanacyl.csv")

# put all dataframes into a list
data_select <- list(nt, tm, tac)
names(data_select) <- c("nt", "tm", "tac")

# need to group different orthologs in same functional group for selected genes
# and customize order of functional groupings
data_select <- lapply(data_select, function(x) {
  y = x %>% 
    dplyr::group_by(site, site_reach, field_date, sample_type, functional_grouping) %>% 
    dplyr::summarize(predicted_gene_abundance = sum(predicted_gene_abundance)) %>% 
    ungroup() %>% 
    mutate(functional_group_factor = factor(functional_grouping, levels = c("nitrogen_fixation",
                                                                             "nitrification",
                                                                             "phosphatase_transporters",
                                                                             "thiamine",
                                                                             "pyridoxal",
                                                                             "cobalamin_B12")))
})

#### (2) Making Box Plots ####

# set universal plot theme
theme_set(theme_bw() + theme(strip.background = element_blank(),
                             plot.title = element_markdown(hjust = 0.5), 
                             text = element_text(size = 10), strip.text = element_markdown(size = 10),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# plot titles
plot_titles <- list("Non-Target Periphyton Assemblages", "*Microcoleus* Assemblages",
                    "*Anabaena* Assemblages")
names(plot_titles) <- c("nt", "tm", "tac")

# make boxplots
# set seed for consistent jitter :)
set.seed(6)
boxplots <- lapply(names(plot_titles), function(x) {
  plot <- ggplot(data_select[[x]], aes(x = site, y = predicted_gene_abundance, fill = site)) +
              geom_boxplot(alpha = 0.8) +
              geom_jitter(alpha = 0.5, aes(fill = site, shape = sample_type), color = "black", size = 2, shape = 21) +
              scale_fill_manual(values = c("SAL" = "#81bbfc",
                                           "SFE-M" = "#416f16",
                                           "RUS" = "#ab9f00")) +
              scale_x_discrete(labels = c("Russian<br>River", "Salmon<br>River", "South Fork<br>Eel River")) + 
              theme(axis.text.x = element_markdown(size = 9, color = "#333333")) +
              facet_wrap(~functional_group_factor, scale = "free",
                         labeller = as_labeller(c(`nitrogen_fixation` = "Nitrogen Fixation", 
                                                  `nitrification`= "Nitrification",
                                                  `phosphatase_transporters` = "Phosphatase Transporters",
                                                  `thiamine` = "Thiamine (B<sub>1</sub>)",
                                                  `pyridoxal` = "Pyrodixal (B<sub>6</sub>)",
                                                  `cobalamin_B12` = "Cobalamin (B<sub>12</sub>)"))) +
              ggtitle(paste(plot_titles[[x]], sep ="")) +
              theme(legend.position = "none") +
              labs(x = NULL, y = "Predicted Gene Abundances")
  return(plot)
})

lapply(boxplots, print)

## TO-DO ADD SIGNIFICANCE DIFFERENCE
