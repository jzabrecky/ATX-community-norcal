

## OLD CODE SNIPPETS BELOW TO REORGANIZE :)

#### could move the below to Q4 ####

## Is Leptolyngbya/Geitlerinema present in all TM samples?
data$tm$leptolyngbya_geitlerinema
count(data$tm$leptolyngbya_geitlerinema > 0)
# present in 17 out of 23

## What about the presence of Leptolyngbya/Geitlerinema in TAC samples?
data$tac$leptolyngbya_geitlerinema
count(data$tac$leptolyngbya_geitlerinema > 0)
# all of them which is crazy at >1%!

## How about Microcoleus in TAC samples? 
## (particularly interested in Russian River where we did not obsere M. macroscopically)
data$tac$microcoleus
count(data$tac$microcoleus > 0 & data$tac$site == "RUS")
# true for 22/28
# 10 of those are russian river
count(data$tac$site == "RUS") # of 15 samples

# Let's look at only other anatoxin associated taxa in all samples
# using list from Christensen & Khan et al. (2019): Anabaena, Aphanizomenon,
# Aphanothece, Arthospira, Cylindrospermopsis, Cylindrospermum, Gomphosphaeria,
# Limnothrix, Lyngbya, Microcystis, Nostoc, Oscillatoria, Phormidium/Microcoleus,
# Planktothrix, Planktolyngbia, Synechocystis, Psuedoanabaena,
# Raphidopsis, Tychonema
# noting that it doesn't have Geilerinema, so we should also include list
# of ATX producers from Wood et al. (2020) which adds: Fisherella, 
# Geitlerinema, Leptolyngbya, Microseira (formerly Lyngbya), planktothrix
lapply(data, function(x) colnames(x[,5:ncol(x)]))
atx_taxa_only <- lapply(data_longer, function(x) {
  
  # make dataframe with only taxa in list above 
  # (only writing what taxa we recorded from that list)
  df = x %>% 
    filter(taxa %in% c("aphanothece", "anabaena_and_cylindrospermum",
                       "other_coccoids", "geitlerinema", "leptolyngbya", 
                       "lyngbya", "other_coccoids", "nostoc",
                       "oscillatoria", "phormidium_unknown", "microcoleus")) %>% 
    mutate(sample_name = paste("(", month(field_date), "-", day(field_date), ") ", site_reach, sep = ""))
  
  # make bar plot (show each sample individually)
  plot <- ggplot(data = df, aes(x = sample_name, y = percent / 100, fill = taxa)) +
    geom_bar(stat = "identity") +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(x = NULL, y = "Relative Abundance") +
    facet_wrap(~site, scales = "free_x")
  print(plot) # view plot
  
  # return a list including dataframe, then plot
  return(list(df, plot))
})

# NOTE: may shove this to a later script