#### (1) Microscopy ####


# also load un-transformed relative abundances and make it longer for bar plots through time
nt <- read.csv("./data/morphological/nt_algalonly.csv") %>% 
tm <- read.csv("./data/morphological/tm_algalonly_nomicro.csv")
tac <- read.csv("./data/morphological/tac_algalonly_noanacylgreenalgae.csv")

microscopy <- list(nt, tm, tac)
names(microscopy) <- c("nt", "tm", "tac")

microscopy <- lapply(microscopy, function(x) {
  x$field_date[which(x$field_date == "2022-07-07")] <- "2022-07-06"
  
  y <- x %>% 
    mutate(field_date = ymd(field_date)) %>% 
    filter(year(field_date) == 2022) %>% 
    pivot_longer(cols = c(5:ncol(.)), names_to = "taxa", values_to = "relative_abundance") %>% 
    group_by(field_date, site, taxa) %>% 
    dplyr::summarize(mean = mean(relative_abundance),
                     sd = sd(relative_abundance))
  return(y) 
})

source("./code/supplemental_code/S4c_grouping_func.R")

# add in broader groupings
microscopy$nt <- nontarget_broader(microscopy$nt)
microscopy$tac <- target_broader(microscopy$tac)
microscopy$tm <- target_broader(microscopy$tm)

source("./code/supplemental_code/S4b_community_analyses_func.R")
# barplots through time
lapply(microscopy, function(x) {
  plot <- barplot(data = x,
        x = "field_date", y  = "mean", fill = "broader") +
  facet_wrap(~site, scales = "free", ncol = 1)
  
  print(plot)
})
