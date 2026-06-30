#### Supplemental figure of boxplots of predicted gene abundances for select functional groups
### Jordan Zabrecky
## last edited: 06.27.2026

# This script creates a supplemental figure showing the range of predicted
# gene abundances for select genes among rivers

#### (1) Loading Libraries & Data ####

# loading libraries
lapply(c("tidyverse", "ggtext", "plyr", "FSA", "cowplot"), require, character.only = T)

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
                                                                             "cobalamin_B12")),
           # also, add in site factor so for TAC we can move Salmon to the right as it's excluded
           # from the statistics
           site_factor = factor(site, levels = c("RUS", "SAL", "SFE-M")))
})

# alter for TAC
data_select$tac <- data_select$tac %>% 
  mutate(site_factor = factor(site, levels = c("RUS", "SFE-M", "SAL")))

# load kruskal-wallis test results from previous script
kruskal_results <- ldply(list.files(path = "./data/kruskal_wallis_results/", pattern = "_selectorthologs.csv"), 
                         function(filename) {
  d <- read.csv(paste("./data/kruskal_wallis_results/", filename, sep = ""))
  d$sample_type = filename %>% stringr::str_remove("_selectorthologs.csv")
  d$sample_type = str_remove(d$sample_type, "Q1_")
  return(d)
})

# add in significance labels to kruskal results
kruskal_results <- kruskal_results %>% 
  mutate(label = case_when(kruskal_test >= 0.05 ~ "n.s.",
                           kruskal_test < 0.05 & kruskal_test >= 0.01 ~ "*",
                           kruskal_test < 0.01 & kruskal_test >= 0.001 ~ "**",
                           kruskal_test < 0.001 ~ "***"))

# follow-up tests results (will put these in manually)
dunnTest(predicted_gene_abundance ~ site,
         data=data_select$nt %>% filter(functional_grouping == "cobalamin_B12"),
         method="bonferroni") # n.s. between SAL & SFE, all sig different from RUS
dunnTest(predicted_gene_abundance ~ site,
         data=data_select$nt %>% filter(functional_grouping == "nitrification"),
         method="bonferroni") # n.s. between RUS & SAL, all sig different from SFE
dunnTest(predicted_gene_abundance ~ site,
         data=data_select$nt %>% filter(functional_grouping == "phosphatase_transporters"),
         method="bonferroni") # n.s. between SAL & SFE, all sig different from RUS
dunnTest(predicted_gene_abundance ~ site,
         data=data_select$nt %>% filter(functional_grouping == "pyridoxal"),
         method="bonferroni") # n.s. between SAL & SFE, all sig different from RUS
dunnTest(predicted_gene_abundance ~ site,
         data=data_select$nt %>% filter(functional_grouping == "thiamine"),
         method="bonferroni") # n.s. between all (will not add second bar then)

# OUTLINE TO DO:
# factor sites so that Salmon River sample is shown but is on right for TAC
# and labels based on SAL, SFE, RUS in plot
# follow-up NT tests
# save plot and check line-ups
# TAC labels & followups
# TM labels & followups

#### (2) Making Box Plots ####

# set universal plot theme
theme_set(theme_bw() + theme(strip.background = element_blank(),
                             plot.title = element_markdown(hjust = 0.5), 
                             text = element_text(size = 9), strip.text = element_markdown(size = 9),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# plot titles
plot_titles <- list("Non-Target Periphyton Assemblages", "*Microcoleus* Assemblages",
                    "*Anabaena* Assemblages")
names(plot_titles) <- c("nt", "tm", "tac")

# significance bars (referencing kruskal results manually :) )
sig_bars <- list()
sig_bars[["nt"]] <- data.frame(functional_group_factor = c(kruskal_results$function_groups[which(kruskal_results$sample_type == "nt")],
                                # dunn's test results
                                c("cobalamin_B12", "nitrification", "phosphatase_transporters", "pyridoxal")),
                               x = c(1, 1, 1, 1, 1, 1, 
                                     2, 1, 2, 2),
                               xend = c(3, 3, 3, 3, 3, 3, 
                                        3, 2, 3, 3),
                               y = c(217000, 2100, 82000, 200000, 145000, 125000,
                                     190000, 1800, 170000, 125000))
sig_bars[["tac"]] <- data.frame(functional_group_factor = kruskal_results$function_groups[which(kruskal_results$sample_type == "tac_nosal")],
                               x = c(1, 1, 1, 1, 1, 1),
                               xend = c(2, 2, 2, 2, 2, 2),
                               y = c(115000, 1900, 160000, 160000, 110000, 85000))
sig_bars[["tm"]] <- data.frame(functional_group_factor = kruskal_results$function_groups[which(kruskal_results$sample_type == "tm")],
                                x = c(1, 1, 1, 1, 1, 1),
                                xend = c(2, 2, 2, 2, 2, 2),
                                y = c(30000, 1100, 50000, 70000, 57000, 35000))

# stars & "n.s."
stars <- list()
stars[["nt"]] <- data.frame(functional_group_factor = c(kruskal_results$function_groups[which(kruskal_results$sample_type == "nt")],
                              # dunn's test
                              c("cobalamin_B12", "nitrification", "phosphatase_transporters", "pyridoxal")),
                            x = c(2, 2, 2, 2, 2, 2,
                                  2.5, 1.5, 2.5, 2.5),
                            y = c(222000, 2170, 90000, 205000, 150000, 129000,
                                  206000, 1970, 185000, 135000),
                            label = c(kruskal_results$label[which(kruskal_results$sample_type == "nt")],
                                      "n.s.", "n.s.", "n.s.", "n.s."))
stars[["tac"]] <- data.frame(functional_group_factor = kruskal_results$function_groups[which(kruskal_results$sample_type == "tac_nosal")],
                                x = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5),
                                y = c(125000, 1950, 166000, 164000, 120000, 91500),
                             label = kruskal_results$label[which(kruskal_results$sample_type == "tac_nosal")])
stars[["tm"]] <- data.frame(functional_group_factor = kruskal_results$function_groups[which(kruskal_results$sample_type == "tm")],
                               x = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5),
                               y = c(32500, 1200, 54000, 76000, 61500, 38000),
                               label = kruskal_results$label[which(kruskal_results$sample_type == "tm")])

# and lastly, scale y axis customized so that we don't have to have the asterisk bump up against the top
y_lims <- data.frame(sample_type = c(rep("nt", 6), rep("tac", 6), rep("tm", 6)),
                     # in alphabetical order
                     ylim = c(230000, 2270, 90000, 210000, 155000, 135000,
                              128000, 2050, 169000, 166000, 122000, 94000,
                              33800, 1250, 55000, 76000, 62000, 39000))


# split into list based on sample type
y_lims <- split(y_lims, y_lims$sample_type)

# make boxplots
# set seed for consistent jitter :)
set.seed(6)
boxplots <- lapply(names(plot_titles), function(x) {
  
  # if anabeana, remove showing the stat summary bar from the single Salmon River sample
  # and make single salmon point not a jitter so it lines up
  if(x == "tac") {
    jitter_data = data_select[[x]] %>% filter(site != "SAL")
  } else {
    jitter_data = data_select[[x]]
  }
  
  # make plot
  plot <- ggplot(data_select[[x]], aes(x = site_factor, y = predicted_gene_abundance, fill = site)) +
              stat_summary(data = jitter_data, fun = median, show.legend = FALSE, geom = "crossbar", color = "#424242") +
              geom_jitter(data = jitter_data, alpha = 0.5,
                          aes(fill = site, shape = sample_type), color = "black", size = 2, shape = 21) +
              scale_fill_manual(values = c("SAL" = "#81bbfc",
                                           "SFE-M" = "#416f16",
                                           "RUS" = "#ab9f00")) +
              scale_x_discrete(labels=c("SAL" = "Salmon<br>River", "SFE-M" = "South Fork<br>Eel River",
                                                         "RUS" = "Russian<br>River")) + 
              theme(axis.text.x = element_markdown(size = 9, color = "#333333")) +
              facet_wrap(~functional_group_factor, scale = "free",
                         labeller = as_labeller(c(`nitrogen_fixation` = "Nitrogen Fixation", 
                                                  `nitrification`= "Nitrification",
                                                  `phosphatase_transporters` = "Phosphatase Transporters",
                                                  `thiamine` = "Thiamine (B<sub>1</sub>)",
                                                  `pyridoxal` = "Pyrodixal (B<sub>6</sub>)",
                                                  `cobalamin_B12` = "Cobalamin (B<sub>12</sub>)"))) +
              ggtitle(paste(plot_titles[[x]], sep ="")) +
              geom_segment(data = sig_bars[[x]], aes(x = x, xend = xend, y = y, yend = y),
                           inherit.aes = FALSE) +
              geom_text(data = stars[[x]], aes(x = x, y = y, label = label), inherit.aes = FALSE) +
              theme(legend.position = "none") +
              labs(x = NULL, y = "Predicted Gene Abundances")
  
  # add in single salmon point (not as jitter) if anabaena plot
  if(x == "tac") {
    plot = plot + geom_point(data = data_select[[x]] %>% filter(site == "SAL"), alpha = 0.5,
                             color = "black", fill = "#81bbfc", size = 2, shape = 21)
  }
  
  # lastly add customized scales so asterisks are not hitting the top too closely
  # note that this makes us loose our factor ordering but ... whatever
  plot = plot + ggh4x::facetted_pos_scales(
    y = list(scale_y_continuous(limits = c(0, y_lims[[x]]$ylim[1])),
             scale_y_continuous(limits = c(0, y_lims[[x]]$ylim[2])),
             scale_y_continuous(limits = c(0, y_lims[[x]]$ylim[3])),
             scale_y_continuous(limits = c(0, y_lims[[x]]$ylim[4])),
             scale_y_continuous(limits = c(0, y_lims[[x]]$ylim[5])),
             scale_y_continuous(limits = c(0, y_lims[[x]]$ylim[6])))
  )
    
  # print plot & return
  print(plot)
  return(plot)
})

#### (3) Putting Together Final Figure ####

# nt only (not necessary but will just throw it into plot grid for consistency w/ target...)
nt <- plot_grid(boxplots[[1]], ncol = 1)
nt

# save
ggsave("./figures/tiff_files/sfig_select_pred_genes_nt.tiff", dpi = 500,
       width = 18, height = 12, unit="cm")

# tac & tm
target <- plot_grid(boxplots[[2]], boxplots[[3]], ncol = 1)
target

# save
ggsave("./figures/tiff_files/sfig_select_pred_genes_target.tiff", dpi = 500,
       width = 18, height = 21, unit = "cm")
