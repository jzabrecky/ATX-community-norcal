##### Comparing PERMANOVA (& visualizations) for our question across phylogenetic ranks
### 05.28.2026
## Jordan Zabrecky

# comparing results of PERMANOVA and PERMDISP and PCoA visualizations
# across phylogenetic ranks for each question and sample type

#### (1) Loading libraries & data ####

# load libraries
lapply(c("tidyverse", "vegan", "cowplot"), require, character.only = T)

# load molecular data (untransformed because we want more than ASV)
nt <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_FINAL.csv") %>% 
  filter(sample_type == "NT")
tm <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_TM_nomicro.csv")
tac <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_TAC_noanacyl.csv")

# join together in list
data <- list(nt, tm, tac)

# square-root transform abundances (and change field_date)
data <- lapply(data, function(x) x %>% mutate(relative_abun = sqrt(picrust2_relative_abundance),
                                              field_date = mdy(field_date)))
names(data) = c("NT", "TM", "TAC")

# list of ranks we care about
ranks <- c("phylum", "class", "order", "feature_ID")

# make wide dataframes for phylums, classes, order, & ASVs 
# (skipping genera because there will be a lot more NA there)
wide_data <- lapply(ranks, function(x) {
  # create empty list
  final <- list()
  
  # iterate through sample types and create wide dataframe
  for(i in c("NT", "TM", "TAC")) {
    data_wide = data[[i]] %>% 
      select(all_of(c("site_reach", "site", "field_date", "sample_type", "relative_abun", x))) %>% 
      # remove rows with NA (aka no phylum given)
      na.omit() %>% 
      # now, group different observations of phylogenetic rank and total
      dplyr::group_by_at(all_of(c("site_reach", "site", "field_date", "sample_type", x))) %>% 
      dplyr::summarize(relative_abun = sum(relative_abun)) %>% 
      ungroup() %>% 
      # finally, pivot wider
      pivot_wider(names_from = x, values_from = "relative_abun", values_fill = 0)
    
    # add to list
    final[[i]] = data_wide
  }
  
  # return list
  return(final)
})
names(wide_data) = c("phylum", "class", "order", "feature_ID")

#### (2) Q1: Among Rivers? ####

# load community analyses functions
source("./code/supplemental_code/S4b_community_analyses_func.R")

# set start col
rivers_start = 5

## (a) NT

# PERMANOVAs
set.seed(1)
lapply(ranks, function(x) runPERMANOVA(wide_data[[x]]$NT, rivers_start,
                                       group = wide_data[[x]]$NT$site))
# phylum: p = 0.003
# class: p = 0.001
# order: p = 0.001
# feature_ID: p = 0.001
# all significant *** !

# PERMDISP
set.seed(1)
lapply(ranks, function(x) anova(betadisper(vegdist(wide_data[[x]]$NT[,rivers_start:ncol(wide_data[[x]]$NT)], 
                                                   method = "bray"), wide_data[[x]]$NT$site)))
# phylum: p = 0.76
# class: p = 0.77
# order: p = 0.85
# feature_ID: p = 0.63
# none are signficant, so PERMANOVA result not influenced

# PCoA
plot <- list()
for(i in ranks) {
  PCoA_data <- getPCoAdata(wide_data[[i]]$NT, start_col = rivers_start)
  plot[[i]] <- makePCoAplot(PCoA_data, color = "site", shape = "site") +
    ggtitle(i)
}
plot_grid(plot[[1]], plot[[2]], plot[[3]], plot[[4]], nrow = 2)

## (b) TAC

# PERMANOVAs
set.seed(1)
lapply(ranks, function(x) runPERMANOVA(wide_data[[x]]$TAC, rivers_start,
                                       group = wide_data[[x]]$TAC$site))
# phylum: p = 0.33
# class: p = 0.221
# order: p = 0.082
# feature_ID: p = 0.001
# only signficant at feature_ID level

# PERMDISP
set.seed(1)
lapply(ranks, function(x) anova(betadisper(vegdist(wide_data[[x]]$TAC[,rivers_start:ncol(wide_data[[x]]$TAC)], 
                                                   method = "bray"), wide_data[[x]]$TAC$site)))
# phylum: p = 0.01
# class: p = 0.009
# order: p = 0.004
# feature_ID: p = 0.001
# significant, so does influence PERMANOVA results but only signficant PERMANOVA
# we saw was feature_ID

# PCoA
plot <- list()
for(i in ranks) {
  PCoA_data <- getPCoAdata(wide_data[[i]]$TAC, start_col = rivers_start)
  plot[[i]] <- makePCoAplot(PCoA_data, color = "site", shape = "site") +
    ggtitle(i)
}
plot_grid(plot[[1]], plot[[2]], plot[[3]], plot[[4]], nrow = 2)
# much more similar at ranks higher than feature ID

## (c) TM

# PERMANOVAs
set.seed(1)
lapply(ranks, function(x) runPERMANOVA(wide_data[[x]]$TM, rivers_start,
                                       group = wide_data[[x]]$TM$site))
# phylum: p = 0.001
# class: p = 0.001
# order: p = 0.001
# feature_ID: p = 0.001
# all a very significant differences

# PERMDISP
set.seed(1)
lapply(ranks, function(x) anova(betadisper(vegdist(wide_data[[x]]$TM[,rivers_start:ncol(wide_data[[x]]$TM)], 
                                                   method = "bray"), wide_data[[x]]$TM$site)))
# phylum: p = 0.38
# class: p = 0.27
# order: p = 0.06
# feature_ID: p = 0.007
# significant, so does influence PERMANOVA results but the groups are plotted very far apart
# on the ones that have significant dispersion differences!

# PCoA
plot <- list()
for(i in ranks) {
  PCoA_data <- getPCoAdata(wide_data[[i]]$TM, start_col = rivers_start)
  plot[[i]] <- makePCoAplot(PCoA_data, color = "site", shape = "site") +
    ggtitle(i)
}
plot_grid(plot[[1]], plot[[2]], plot[[3]], plot[[4]], nrow = 2)
# much more similar at ranks higher than feature ID

#### (3) Q2: W/ Varying Anatoxin concentrations ####

# load in anatoxin data
atx <- read.csv("./data/field_and_lab/atx_w_categorical_groupings.csv") %>% 
  mutate(field_date = ymd(field_date))

# join dataframes with anatoxin data
wide_data_atx <- lapply(wide_data, function(x) {
  for(i in c("NT", "TM", "TAC")) {
   x[[i]] <- left_join(atx, x[[i]], by = c("site", "site_reach", "field_date", "sample_type")) %>% na.omit()
  }
  return(x)
})

# set start column for dataframe with atx
atx_start <- 8

## (a) TM in South Fork Eel

# PERMANOVAs - detection
set.seed(1)
lapply(ranks, function(x) runPERMANOVA(wide_data_atx[[x]]$TM %>% filter(site == "SFE-M"), atx_start,
                                       group = (wide_data_atx[[x]]$TM %>% filter(site == "SFE-M"))$atx_detected))
# phylum: p = 0.582
# class: p = 0.392
# order: p = 0.177
# feature_ID: p = 0.002
# significant only at ASV level

# PERMDISP - detection
set.seed(1)
lapply(ranks, function(x) anova(betadisper(vegdist((wide_data_atx[[x]]$TM %>% filter(site == "SFE-M"))[,atx_start:ncol(wide_data[[x]]$TM)], 
                                                   method = "bray"), (wide_data_atx[[x]]$TM %>% filter(site == "SFE-M"))$atx_detected)))
# phylum: p = 0.582
# class: p = 0.392
# order: p = 0.177
# feature_ID: p = 0.002
# true difference at ASV level

# PERMANOVAs - high versus low
set.seed(1)
lapply(ranks, function(x) runPERMANOVA(wide_data_atx[[x]]$TM %>% filter(site == "SFE-M" & atx_group != "none"), atx_start,
                                       group = (wide_data_atx[[x]]$TM %>% filter(site == "SFE-M" & atx_group != "none"))$atx_group))
# phylum: p = 0.77
# class: p = 0.75
# order: p = 0.85
# feature_ID: p = 0.78
# no difference with varying anatoxin concentrations

# PCoA
plot <- list()
for(i in ranks) {
  PCoA_data <- getPCoAdata(wide_data_atx[[i]]$TM %>% filter(site == "SFE-M"), start_col = atx_start)
  plot[[i]] <- makePCoAplot(PCoA_data, color = "atx_detected", shape = "atx_group") +
    ggtitle(i)
}
plot_grid(plot[[1]], plot[[2]], plot[[3]], plot[[4]], nrow = 2)

## (b) TAC in South Fork Eel

# PERMANOVAs - detection
set.seed(1)
lapply(ranks, function(x) runPERMANOVA(wide_data_atx[[x]]$TAC %>% filter(site == "SFE-M"), atx_start,
                                       group = (wide_data_atx[[x]]$TAC %>% filter(site == "SFE-M"))$atx_detected))
# phylum: p = 0.225
# class: p = 0.224
# order: p = 0.019*
# feature_ID: p = 0.013*
# significant only at ASV level

# PERMDISP - detection
set.seed(1)
lapply(ranks, function(x) anova(betadisper(vegdist((wide_data_atx[[x]]$TAC %>% filter(site == "SFE-M"))[,atx_start:ncol(wide_data[[x]]$TM)], 
                                                   method = "bray"), (wide_data_atx[[x]]$TAC %>% filter(site == "SFE-M"))$atx_detected)))
# phylum: p = 0.02*
# class: p = 0.04*
# order: p = 0.887
# feature_ID: p = 0.3797
# dispersion not an issue for significant result at ASV level!

# PERMANOVAs - high versus low
set.seed(1)
lapply(ranks, function(x) runPERMANOVA(wide_data_atx[[x]]$TAC %>% filter(site == "SFE-M" & atx_group != "none"), atx_start,
                                       group = (wide_data_atx[[x]]$TAC %>% filter(site == "SFE-M" & atx_group != "none"))$atx_group))
# phylum: p = 0.26
# class: p = 0.22
# order: p = 0.186
# feature_ID: p = 0.034*
# difference in composition at ASV level between low and high!
# interesting given we suspect that Anabaena is not an anatoxin producer!

# PCoA
plot <- list()
for(i in ranks) {
  PCoA_data <- getPCoAdata(wide_data_atx[[i]]$TAC %>% filter(site == "SFE-M"), start_col = atx_start)
  plot[[i]] <- makePCoAplot(PCoA_data, color = "atx_detected", shape = "atx_group") +
    ggtitle(i)
}
plot_grid(plot[[1]], plot[[2]], plot[[3]], plot[[4]], nrow = 2)
# supports above

## (c) SFE-M NT

# PERMANOVAs - detection
set.seed(1)
lapply(ranks, function(x) runPERMANOVA(wide_data_atx[[x]]$NT %>% filter(site == "SFE-M"), atx_start,
                                       group = (wide_data_atx[[x]]$NT %>% filter(site == "SFE-M"))$atx_detected))
# phylum: p = 0.004**
# class: p = 0.003**
# order: p = 0.007**
# feature_ID: p = 0.004**
# significant difference at all levels

# PERMDISP - detection
set.seed(1)
lapply(ranks, function(x) anova(betadisper(vegdist((wide_data_atx[[x]]$NT %>% filter(site == "SFE-M"))[,atx_start:ncol(wide_data[[x]]$TM)], 
                                                   method = "bray"), (wide_data_atx[[x]]$NT %>% filter(site == "SFE-M"))$atx_detected)))
# phylum: p = 0.92
# class: p = 0.85
# order: p = 0.70
# feature_ID: p = 0.16
# dispersion not an issue for any results above

# PERMANOVAs - high versus low
set.seed(1)
lapply(ranks, function(x) runPERMANOVA(wide_data_atx[[x]]$NT %>% filter(site == "SFE-M" & atx_group != "none"), atx_start,
                                       group = (wide_data_atx[[x]]$NT %>% filter(site == "SFE-M" & atx_group != "none"))$atx_group))
# phylum: p = 0.20
# class: p = 0.16
# order: p = 0.13
# feature_ID: p = 0.01**
# difference in composition at ASV level between low and high!

# PCoA
plot <- list()
for(i in ranks) {
  PCoA_data <- getPCoAdata(wide_data_atx[[i]]$NT %>% filter(site == "SFE-M"), start_col = atx_start)
  plot[[i]] <- makePCoAplot(PCoA_data, color = "atx_detected", shape = "atx_group") +
    ggtitle(i)
}
plot_grid(plot[[1]], plot[[2]], plot[[3]], plot[[4]], nrow = 2)
# supports above

## (d) RUS TAC

# PERMANOVAs - detection
set.seed(1)
lapply(ranks, function(x) runPERMANOVA(wide_data_atx[[x]]$TAC %>% filter(site == "RUS"), atx_start,
                                       group = (wide_data_atx[[x]]$TAC %>% filter(site == "RUS"))$atx_detected))
# phylum: p = 0.27
# class: p = 0.288
# order: p = 0.262
# feature_ID: p = 0.111
# not significant

# PERMANOVAs - high versus low
set.seed(1)
lapply(ranks, function(x) runPERMANOVA(wide_data_atx[[x]]$TAC %>% filter(site == "RUS" & atx_group != "none"), atx_start,
                                       group = (wide_data_atx[[x]]$TAC %>% filter(site == "RUS" & atx_group != "none"))$atx_group))
# phylum: p = 0.1
# class: p = 0.1
# order: p = 0.1
# feature_ID: p = 0.1
# no differences, but also, seems like not enough points to compare

# PCoA
plot <- list()
for(i in ranks) {
  PCoA_data <- getPCoAdata(wide_data_atx[[i]]$TAC %>% filter(site == "RUS"), start_col = atx_start)
  plot[[i]] <- makePCoAplot(PCoA_data, color = "atx_detected", shape = "atx_group") +
    ggtitle(i)
}
plot_grid(plot[[1]], plot[[2]], plot[[3]], plot[[4]], nrow = 2)
# supports above

## (e) RUS NT

# PERMANOVAs - detection
set.seed(1)
lapply(ranks, function(x) runPERMANOVA(wide_data_atx[[x]]$NT %>% filter(site == "RUS"), atx_start,
                                       group = (wide_data_atx[[x]]$NT %>% filter(site == "RUS"))$atx_detected))
# phylum: p = 0.348
# class: p = 0.353
# order: p = 0.621
# feature_ID: p = 0.487
# not significant

# PERMANOVAs - high versus low
set.seed(1)
lapply(ranks, function(x) runPERMANOVA(wide_data_atx[[x]]$NT %>% filter(site == "RUS" & atx_group != "none"), atx_start,
                                       group = (wide_data_atx[[x]]$NT %>% filter(site == "RUS" & atx_group != "none"))$atx_group))
# phylum: p = 0.2667
# class: p = 0.333
# order: p = 0.2667
# feature_ID: p = 0.4
# no differences, but also, seems like not enough points to compare

# PCoA
plot <- list()
for(i in ranks) {
  PCoA_data <- getPCoAdata(wide_data_atx[[i]]$NT %>% filter(site == "RUS"), start_col = atx_start)
  plot[[i]] <- makePCoAplot(PCoA_data, color = "atx_detected", shape = "atx_group") +
    ggtitle(i)
}
plot_grid(plot[[1]], plot[[2]], plot[[3]], plot[[4]], nrow = 2)
# supports above (one outlier is particularly annoying me)
