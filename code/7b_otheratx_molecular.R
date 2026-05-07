#### Identifying other anatoxin-associated taxa in microscopy target samples
### Jordan Zabrecky
## last edited: 05.06.2026

# This code identifies other anatoxin-associated taxa as listed by
# Christensen & Khan (2019) and Wood et al. (2020) in molecular samples

# MAYBE REWRITE WITHOUT TARGET TAXA TAKEN OUT - think about supplemental table organization

#### (1) Load data and libraries ####

# load libraries
lapply(c("tidyverse"), require, character.only = T)

# load microscopy data & pivot longer
tm <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_TM_nomicro.csv")
tac <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_TAC_noanacyl.csv")

# put into list
data_longer <- list(tac, tm)
names(data_longer) <- c("tac", "tm")

#### (8) Misc. Questions ####

## How present are other anatoxin associated taxa?
# using ATX taxa as identified in Christensen & Khan (2019) and Wood et al. (2020)
# using list from Christensen & Khan et al. (2019): Anabaena, Aphanizomenon,
# Aphanothece, Arthospira, Cylindrospermopsis, Cylindrospermum, Gomphosphaeria,
# Limnothrix, Lyngbya, Microcystis, Nostoc, Oscillatoria, Phormidium/Microcoleus,
# Planktothrix, Planktolyngbia, Synechocystis, Psuedoanabaena,
# Raphidopsis, Tychonema
# noting that it doesn't have Geilerinema, so we should also include list
# of ATX producers from Wood et al. (2020) which adds: Fisherella, 
# Geitlerinema, Leptolyngbya, Microseira (formerly Lyngbya)

# note: searched within dataframes for genuses of the above and made the list below

atx_taxa_only <- lapply(data_longer, function(x) {
  
  # make dataframe with only taxa in list above 
  # (only writing what taxa we recorded from that list)
  df = x %>% 
    filter(genus %in% c("Anabaena", "Aphanizomenon", "Cylindrospermum", "Cylindrospermopsis",
                        "Limnothrix", "Leptolyngbya", "Microcystis", "Nostoc", "Oscillatoria",
                        "Phormidium", "Microcoleus", "Planktothrix", "Synechocystis", 
                        "Geitlerinema", "Pseudanabaena"
    )) %>% 
    # searching the TAC and TM dataframes do not yield: aphanothece, arthospira, gomphosphaeria,
    # lyngbya, planktolyngbya, microseira
    mutate(field_date = mdy(field_date),
           sample_name = paste("(", month(field_date), "-", day(field_date), ") ", site_reach, sep = "")) %>% 
    select(sample_name, field_date, sample_type, site, site_reach, picrust2_relative_abundance, order, genus) %>% 
    mutate(order_genus = paste(order, " - ", genus, sep = ""))
  
  
  # get number of samples
  n_samples = length(unique(df$sample_name))
  
  # determine percent of samples with each taxa as well as mean percent for each taxa
  df2 = df %>% 
    # need to group ASVs of same genus and calculate total abundance for each genus
    group_by(sample_name, site, genus) %>% 
    dplyr::summarize(relative_abundance = sum(picrust2_relative_abundance)) %>% 
    ungroup() %>% 
    # now, can caluate percent of samples with each genus present
    filter(relative_abundance != 0) %>% 
    dplyr::group_by(genus) %>% 
    dplyr::summarize(count = length(genus)) %>% 
    ungroup() %>% 
    mutate(percent_samples = count / n_samples * 100)
  
  # get mean and sd
  df3 = df %>% 
    # need to group ASVs of same genus and calculate total abundance for each genus
    group_by(sample_name, site, genus) %>% 
    dplyr::summarize(relative_abundance = sum(picrust2_relative_abundance)) %>% 
    ungroup() %>% 
    # now can calculate means and sds
    # care about means when present?
    filter(relative_abundance != 0) %>% 
    dplyr::group_by(genus) %>% 
    dplyr::summarize(mean = mean(relative_abundance),
                     sd = sd(relative_abundance)) %>% 
    ungroup()
  
  # return a list of dataframes
  final <- list(df, df2, df3)
  names(final) <- c("raw", "percent_samples", "mean_sd")
  return(final)
})

# view dataframes (adjust as needed)
view(atx_taxa_only$tm$percent_samples)

# some discrepancies here between and microscopy data
# could be limits of 16s with genus level or the poor resolution of database
# for more, see Dvorak et al. (2025) or, of course, poor resolution with microscopy!