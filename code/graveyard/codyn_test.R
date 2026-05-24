install.packages("codyn")


library(codyn)

data("collins08")

head(collins08[1:3,])
rate_change(df = collins08, 
         time.var = "year", species.var = "species",
         abundance.var = "abundance", replicate.var = "replicate")
rate_change_interval(df = collins08, 
                     time.var = "year", species.var = "species",
                     abundance.var = "abundance", replicate.var = "replicate")

test <- read.csv("./data/molecular/transformed/16s_nochimera_rarefied_95_TM_nomicro_sqrttransformed.csv") %>% 
  select(site == "SFE-M")

data_longest <- pivot_longer(data = data$tac, cols = c(8:ncol(data$tac)), names_to = "ASV", values_to = "transformed_rel_abun")

rate_change_interval(df = data_longest, time.var = "event_no",
                     species.var = "ASV", abundance.var = "transformed_rel_abun", replicate.var = "site_reach")
