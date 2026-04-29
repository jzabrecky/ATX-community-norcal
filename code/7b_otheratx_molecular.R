
### old code from script 4b

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
atx_taxa_only <- lapply(data_long, function(x) {
  
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
    mutate(sample_name = paste("(", month(field_date), "-", day(field_date), ") ", site_reach, sep = "")) %>% 
    select(sample_name, field_date, sample_type, site, site_reach, relative_abundance, order, genus) %>% 
    mutate(order_genus = paste(order, " - ", genus, sep = ""))
  
  # make bar plot (show each sample individually)
  plot <- ggplot(data = df, aes(x = sample_name, y = relative_abundance / 100, fill = genus)) +
    geom_bar(stat = "identity") +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(x = NULL, y = "Relative Abundance") +
    facet_wrap(~site, scales = "free_x")
  print(plot) # view plot
  
  # return a list including dataframe, then plot
  return(list(df, plot))
})

# some discrepancies here between and microscopy data
# could be limits of 16s with genus level or the poor resolution of database
# for more, see Dvorak et al. (2025)