#### Comparing plots of PCoA and NMDS
### Jordan Zabrecky
## 05.27.2026

# This script compares using NMDS versus PCoA for plot ordinations

#### (1) Loading libraries & data ####

# load libraries

# microscopy data
microscopy <- lapply(list.files(path = "./data/morphological/transformed/", pattern = ".csv"),
                     function(x) read.csv(paste("./data/morphological/transformed/", x, sep = "")))
names(microscopy) <- c("NT", "TAC", "TM")

# molecular data 
molecular <- lapply(list.files(path = "./data/molecular/transformed/", pattern = "16s_nochimera"),
                    function(x) read.csv(paste("./data/molecular/transformed/", x, sep = "")))
names(molecular) <- c("NT", "TAC", "TM")

# functional data
functional <- lapply(list.files(path = "./data/molecular/transformed/", pattern = "KO_all"),
                     function(x) read.csv(paste("./data/molecular/transformed/", x, sep = "")))
names(functional) <- c("NT", "TAC", "TM")

#### (2) Q1: Among Rivers ####

# get community analyses functions
source("./code/supplemental_code/S4b_community_analyses_func.R")

# set start_col
micro_start <- 5
molec_start <- 6
func_start <- 5

# function to make plots to compare the two
compare <- function(data, start_col, shape, color, end_col = NA) {
  
  # set up NMDS (won't run loadings)
  NMDS_data = getNMDSdata(data, start_col, end_col = end_col, ASV = TRUE)
  NMDS_plot = makeNMDSplot(NMDS_data, FALSE, FALSE, shape = shape, color = color) +
    ggtitle("NMDS") +
    theme(legend.position = "bottom")
  
  # set up PCoA
  PCoA_data = getPCoAdata(data, start_col, end_col)
  PCoA_plot = makePCoAplot(PCoA_data, shape = shape, color = color) +
    ggtitle("PCoA") + theme(legend.position = "bottom")
  
  # put two plots together
  final = plot_grid(NMDS_plot, PCoA_plot, nrow = 1)
  print(final)
}

# compare plots!
set.seed(1)
lapply(microscopy, function(x) compare(x, micro_start, shape = "month", color = "site"))
# no major changes!
lapply(molecular, function(x) compare(x, molec_start, shape = "month", color = "site"))
# TAC appear more different; however we know on a class level, they are still similar 
lapply(functional, function(x) compare(x, func_start, shape = "month", color = "site"))
# outlier still annoying with PCoA, not appealing

# removing outlier
compare(functional$NT %>% filter(!c(site_reach == "RUS-1S" & field_date == "2022-07-20")),
        func_start, shape = "month", color = "site")

#### (3) Q2: W/ Varying ATX Conc. ####

# load anatoxin data
atx <- read.csv("./data/field_and_lab/atx_w_categorical_groupings.csv")

# join in microscopy data
microscopy_atx <- lapply(microscopy, function(x) 
  left_join(x, atx, by = c("field_date", "site_reach", "site", "sample_type")))
molecular_atx <- lapply(molecular, function(x) 
  left_join(x, atx, by = c("field_date", "site_reach", "site", "sample_type")))
functional_atx <- lapply(functional, function(x) 
  left_join(x, atx, by = c("field_date", "site_reach", "site", "sample_type")))

# south fork eel samples first
set.seed(1)
lapply(c("NT", "TAC", "TM"), function(x){
  end_col = ncol(microscopy[[x]]) # need to set an end_col as we left_joined in data
  data = microscopy_atx[[x]] %>% filter(site == "SFE-M") # filter for SFE
  compare(data, micro_start, end_col = end_col, shape = "atx_group", color = "atx_detected")
})
# all good
lapply(c("NT", "TAC", "TM"), function(x){
  end_col = ncol(molecular[[x]]) # need to set an end_col as we left_joined in data
  data = molecular_atx[[x]] %>% filter(site == "SFE-M") %>% na.omit() # filter for SFE
  compare(data, molec_start, end_col = end_col, shape = "atx_group", color = "atx_detected")
})
# able to plot
lapply(c("NT", "TAC", "TM"), function(x){
  end_col = ncol(functional[[x]]) # need to set an end_col as we left_joined in data
  data = functional_atx[[x]] %>% filter(site == "SFE-M") %>% na.omit() # filter for SFE
  compare(data, func_start, end_col = end_col, shape = "atx_group", color = "atx_detected")
})

# now, russian samples
set.seed(1)
lapply(c("NT", "TAC"), function(x){
  end_col = ncol(microscopy[[x]]) # need to set an end_col as we left_joined in data
  data = microscopy_atx[[x]] %>% filter(site == "RUS") %>% na.omit() # filter for SFE
  compare(data, micro_start, end_col = end_col, shape = "atx_group", color = "atx_detected")
})
lapply(c("NT", "TAC"), function(x){
  end_col = ncol(molecular[[x]]) # need to set an end_col as we left_joined in data
  data = molecular_atx[[x]] %>% filter(site == "RUS") %>% na.omit() # filter for SFE
  compare(data, molec_start, end_col = end_col, shape = "atx_group", color = "atx_detected")
})
# still have one outlier here that would be nice to remove
lapply(c("NT", "TAC"), function(x){
  end_col = ncol(functional[[x]]) # need to set an end_col as we left_joined in data
  data = functional_atx[[x]] %>% filter(site == "RUS") %>% na.omit() # filter for SFE
  compare(data, func_start, end_col = end_col, shape = "site_reach", color = "month")
})
# same outlier sample is annoying here

# removing outlier in molecular sample
compare(molecular_atx$NT %>% filter(site == "RUS") %>% filter(!c(site_reach == "RUS-1S" & field_date == "2022-09-15")) %>% na.omit(),
        molec_start, end_col = ncol(molecular$NT), shape = "site_reach", color = "month")
