#### Identifying other anatoxin-associated taxa in microscopy target samples
### Jordan Zabrecky
## last edited: 10.15.2026

# This code identifies other anatoxin-associated taxa as listed by
# Wood et al. (2020), Chapman & Foss (2020), and Bernard et al. (2017) in microscopy samples

#### (1) Load data and libraries ####

# load libraries
lapply(c("tidyverse"), require, character.only = T)

# load list of toxic taxa from papers
producers <- read.csv("./data/potential_toxin_producers.csv")

# load microscopy data & pivot longer
tm <- read.csv("./data/morphological/tm_algalonly_nomicro.csv") %>% 
  pivot_longer(c(5:ncol(.)), names_to = "taxa", values_to = "percent") %>% 
  filter(year(ymd(field_date)) == 2022)
tac <- read.csv("./data/morphological/tac_algalonly_noanacylgreenalgae.csv")%>% 
  pivot_longer(c(5:ncol(.)), names_to = "taxa", values_to = "percent") %>% 
  filter(year(ymd(field_date)) == 2022)

# put data into list
data_longer <- list(tac, tm)
names(data_longer) <- c("tac", "tm")

#### (2) Identify other Anatoxin-Associated Taxa ####

# see our taxa groups
unique(tac$taxa)
unique(tm$taxa)

# subset for only anatoxin producing taxa, calculate number of samples with that taxa,
# and mean and standard deviation of taxa in samples
atx_taxa_only <- lapply(data_longer, function(x) {
  
  # make dataframe with only taxa in list above 
  # (only writing what taxa we recorded from that list)
  df = x %>% 
    # filter for in potential producers list (plus ones that have multiple taxa listed)
    filter(taxa %in% c(tolower(producers$Genus),
                       "anabaena_and_cylindrospermum", "miscellaneous_oscillatoriales",
                       "phormidium_unknown", "anabaena_and_cylindrospermum", 
                       "leptolyngbya_and_geitlerinema")) %>% 
    mutate(sample_name = paste("(", month(field_date), "-", day(field_date), ") ", site_reach, sep = ""))
  
  
  # get number of samples
  n_samples = length(unique(df$sample_name))
  
  # determine percent of samples with each taxa as well as mean percent for each taxa
  df2 = df %>% 
    filter(percent != 0) %>% 
    dplyr::group_by(taxa) %>% 
    dplyr::summarize(count = length(taxa)) %>% 
    ungroup() %>% 
    mutate(percent_samples = count / n_samples * 100)
  
  # get mean and sd
  df3 = df %>% 
    # care about means when present?
    filter(percent != 0) %>% 
    dplyr::group_by(taxa) %>% 
    dplyr::summarize(mean = mean(percent),
                     sd = sd(percent)) %>% 
    ungroup()
  
  # return a list of dataframes
  final <- list(df, df2, df3)
  names(final) <- c("raw", "percent_samples", "mean_sd")
  return(final)
})

# view results
view(atx_taxa_only$tm$percent_samples)
