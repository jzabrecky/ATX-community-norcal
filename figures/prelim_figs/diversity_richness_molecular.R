


environmental <- read.csv("./data/field_and_lab/environmental_covariates.csv") %>% 
  mutate(field_date = ymd(field_date)) %>% 
  filter(year(field_date) == 2022)
  

molecular <- read.csv("./data/molecular/16s_nochimera_rarefied_90_FINAL.csv") %>% 
  select(site_reach, site, field_date, sample_type, feature_ID, relative_abundance) %>% 
  group_by(site_reach, site, field_date, sample_type, feature_ID) %>% 
  dplyr::summarize(relative_abundance = sum(relative_abundance)) %>% 
  pivot_wider(names_from = feature_ID, values_from = relative_abundance, values_fill = 0)

# separate based on sample type
molecular_type <- split(molecular, molecular$sample_type)

## richness
molecular_species_richness <- lapply(molecular_type, function(x)
  x %>% mutate(richness = rowSums(x[c(5:ncol(x))] > 0)))

# plot
richness_plot <- list()
for(i in 1:length(molecular_species_richness)) {
  richness_plot[[i]] = ggplot(data = molecular_species_richness[[i]], aes(x = site, y = richness)) +
    geom_boxplot(aes(color = site)) +
    labs(title = paste(names(molecular_species_richness[i])," [richness via ASVs]", sep = ""))
  print(richness_plot[[i]])
}
richness <- plot_grid(richness_plot[[1]], richness_plot[[2]], richness_plot[[3]])
richness


# plot over time
richness_plot_time <- list()
for(i in 1:length(molecular_species_richness)) {
  richness_plot_time[[i]] = ggplot(data = molecular_species_richness[[i]], aes(x = field_date, y = richness)) +
    #facet_wrap(~year, scales = "free", ncol = 1) +
    geom_point(aes(color = site)) +
    labs(title = names(molecular_species_richness[i]))
  print(richness_plot_time[[i]])
}
richness_nt <- microscopy_species_richness$NT %>% 
  select(field_date, site_reach, site, richness)
richness_nt <- left_join(richness_nt, environmental, by = c("field_date", "site_reach", "site"))


## diversity
molecular_species_diversity <- lapply(molecular_type, function(x) {
  x$diversity = diversity(x[,5:ncol(x)])
  return(x)})

# plot
diversity_plot <- list()
for(i in 1:length(molecular_species_diversity)) {
  diversity_plot[[i]] = ggplot(data = molecular_species_diversity[[i]], aes(x = site, y = diversity)) +
    geom_boxplot(aes(color = site)) +
    labs(title = names(molecular_species_diversity[i]))
  print(diversity_plot[[i]])
}
diversity <- plot_grid(diversity_plot[[1]], diversity_plot[[2]], diversity_plot[[3]])
diversity


diversity_plot_time <- list()
for(i in 1:length(molecular_species_diversity)) {
  diversity_plot_time[[i]] = ggplot(data = molecular_species_diversity[[i]], aes(x = field_date, y = diversity)) +
    #facet_wrap(~year, scales = "free", ncol = 1) +
    geom_point(aes(color = site)) +
    labs(title = names(molecular_species_diversity[i]))
  print(diversity_plot_time[[i]])
}

# split based on river- does composition change over time? via adonis
molecular_sites <- lapply(molecular_type, function(x) split(x, x$site))
adonis2(molecular_sites$NT$`SFE-M`[5:ncol(molecular_sites$NT$`SFE-M`)] ~ field_date, 
        data = molecular_sites$NT$`SFE-M`)
# significant change in NT composition SFE-M over time (p = 0.018)
adonis2(molecular_sites$TM$`SFE-M`[5:ncol(molecular_sites$TM$`SFE-M`)] ~ field_date, 
        data = molecular_sites$TM$`SFE-M`)
# no significant change in TM composition SFE-M over time (p = 0.249)
adonis2(molecular_sites$TAC$`SFE-M`[5:ncol(molecular_sites$TAC$`SFE-M`)] ~ field_date, 
        data = molecular_sites$TAC$`SFE-M`)
# significant change in TAC composition over time in SFE-M (p = 0.004)

adonis2(molecular_sites$NT$`SAL`[5:ncol(molecular_sites$NT$`SAL`)] ~ field_date, 
        data = molecular_sites$NT$`SAL`)
# significant change in NT composition SAL over time (p = 0.009)
adonis2(molecular_sites$TM$`SAL`[5:ncol(molecular_sites$TM$`SAL`)] ~ field_date, 
        data = molecular_sites$TM$`SAL`)
# no significant change in TM composition SFE-M over time (p = 0.2)


adonis2(molecular_sites$TAC$`RUS`[5:ncol(molecular_sites$TAC$`RUS`)] ~ field_date, 
        data = molecular_sites$TAC$`RUS`)
# no significant change in NT composition RUS over time (p = 0.582)
adonis2(molecular_sites$TAC$`RUS`[5:ncol(molecular_sites$TAC$`RUS`)] ~ field_date, 
        data = molecular_sites$TAC$`RUS`)
# no significant change in TAC composition RUS over time (p = 0.564)
