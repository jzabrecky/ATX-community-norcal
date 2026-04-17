
# loading libraries
library(tidyverse)

# load data 
nt <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select.csv") %>% 
  filter(sample_type == "NT") %>% 
  mutate(field_date = mdy(field_date))
tm <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tm_nomicro.csv") %>% 
  mutate(field_date = mdy(field_date))
tac <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tac_noanacyl.csv") %>% 
  mutate(field_date = mdy(field_date))

nt_list <- split(nt, nt$site)
tm_list <- split(tm, tm$site)
tac_list <- split(tac, tac$site)

# plot across time

lapply(nt_list, function(x) {
  title = x$site[1]
  ggplot(x, aes(x = field_date, y = predicted_gene_abundance)) +
    geom_point(aes(color = site_reach, shape = site_reach), size = 3) +
    facet_wrap(~my_grouping, scales = "free") +
    labs(title = title)
})

lapply(tm_list, function(x) {
  title = x$site[1]
  ggplot(x, aes(x = field_date, y = predicted_gene_abundance)) +
    geom_point(aes(color = site_reach, shape = site_reach), size = 3) +
    facet_wrap(~my_grouping, scales = "free") +
    labs(title = title)
})

lapply(tac_list, function(x) {
  title = x$site[1]
  ggplot(x, aes(x = field_date, y = predicted_gene_abundance)) +
    geom_point(aes(color = site_reach, shape = site_reach), size = 3) +
    facet_wrap(~my_grouping, scales = "free") +
    labs(title = title)
})
