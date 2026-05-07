#### Identifying other anatoxin-associated taxa in microscopy target samples
### Jordan Zabrecky
## last edited: 05.06.2026

# This code identifies other anatoxin-associated taxa as listed by
# Christensen & Khan (2019) and Wood et al. (2020) in microscopy samples

#### (1) Load data and libraries ####

# load libraries
lapply(c("tidyverse"), require, character.only = T)

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

# Let's look at only other anatoxin associated taxa in all samples
# using list from Christensen & Khan et al. (2019): Anabaena, Aphanizomenon,
# Aphanothece, Arthospira, Cylindrospermopsis, Cylindrospermum, Gomphosphaeria,
# Limnothrix, Lyngbya, Microcystis, Nostoc, Oscillatoria, Phormidium/Microcoleus,
# Planktothrix, Planktolyngbia, Synechocystis, Psuedoanabaena,
# Raphidopsis, Tychonema
# noting that it doesn't have Geilerinema, so we should also include list
# of ATX producers from Wood et al. (2020) which adds: Fisherella, 
# Geitlerinema, Leptolyngbya, Microseira (formerly Lyngbya), planktothrix

# see our taxa groups
unique(tac$taxa)
unique(tm$taxa)

# subset for only anatoxin producing taxa, calculate number of samples with that taxa,
# and mean and standard deviation of taxa in samples
atx_taxa_only <- lapply(data_longer, function(x) {
  
  # make dataframe with only taxa in list above 
  # (only writing what taxa we recorded from that list)
  df = x %>% 
    filter(taxa %in% c("aphanothece", "anabaena_and_cylindrospermum",
                       "other_coccoids", "leptolyngbya_and_geitlerinema", 
                       "miscellaneous_oscillatoriales", "other_coccoids", "nostoc",
                       "oscillatoria", "phormidium_unknown", "microcoleus")) %>% 
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
view(atx_taxa_only$tac$percent_samples)
