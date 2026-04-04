#### Comparing molecular 16s data among rivers
### Jordan Zabrecky
## last edited: 03.23.2026

# This code compares normalized select orthologs/functions predicted via PICRUSt2-SC,
# from NT, TM, and TAC samples across rivers to answer Q1.
# Data is analyzed using Kruskal-Wallis Tests and visualizations

# TBD on keeping ISA results

#### (1) Loading libraries & data ####

# load libraries
lapply(c("tidyverse", "plyr", "cowplot", "indicspecies"), require, character.only = T)

# load data 
data <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select.csv")

#### (2) Plotting #####

# create summary table with means and plus or minuses
summary <- data %>% 
  dplyr::group_by(site, my_grouping, sample_type) %>% 
  dplyr::summarize(mean = mean(predicted_gene_abundance),
                   sd = sd(predicted_gene_abundance),
                   total = n(),
                   se = sd / sqrt(total))

function_plot <- ggplot(summary, aes(x = site, y = mean, fill = site)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(x=site, ymin=mean-se, ymax=mean+se, color = site), width = 0.3) +
  facet_grid(my_grouping~sample_type, scale = "free") + 
  scale_x_discrete(guide = guide_axis(angle = 75))
function_plot # will make nicer later

#### (3) Kruskal Wallis Tests ####

# split into list for sample types
data_list <- split(data, data$sample_type)

kruskal_test_results <- lapply(data_list, function(x) {
  function_groups = unique(x$my_grouping)
  results = data.frame(function_groups = NA,
                       kruskal_test = NA)
  for(i in function_groups) {
    results = rbind(results, data.frame(function_groups = i,
               kruskal_test = (kruskal.test(site~predicted_gene_abundance, data = (x %>% filter(my_grouping == i))))$p.value))
  }
  
  return(results[-1,])
})

lapply(kruskal_test_results, function(x) x[which(x$kruskal_test < 0.1),])
# nitrification is the only one significantly different among groups!

# post-hoc Dunns test for nitrification for TM and NT