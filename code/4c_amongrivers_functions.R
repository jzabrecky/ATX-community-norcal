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
nt <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select.csv") %>% 
  filter(sample_type == "NT")
tm <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tm_nomicro.csv")
tac <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tac_noanacyl.csv")

# put all dataframes into a list
data <- list(nt, tm, tac)
names(data) <- c("nt", "tm", "tac")

#### (2) Plotting #####

# make box plots!!!
boxplots <- lapply(data, function(x) {
  plot = ggplot(x, aes(x = my_grouping, y = predicted_gene_abundance, fill = site)) +
    geom_boxplot()
  print(plot)
  return(plot)
})

#### (3) Kruskal Wallis Tests ####

# run tests
kruskal_test_results <- lapply(data, function(x) {
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
# none significantly different among groups
