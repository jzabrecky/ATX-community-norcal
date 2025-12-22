library(vegan)
library(tidyverse)
library(plyr)
library(cowplot)

environmental <- read.csv("./data/field_and_lab/environmental_covariates.csv") %>% 
  mutate(field_date = ymd(field_date))

microscopy <- ldply(list.files(path = "./data/morphological/", pattern = "algalonly.csv"), 
                    function(filename) {
  d <- read.csv(paste("data/morphological/", filename, sep = ""))
  return(d)
})

microscopy <- microscopy %>% 
  mutate(field_date = ymd(field_date),
         year = year(field_date)) %>% 
  relocate(year, .before = "field_date")

# separate based on sample type
microscopy_list <- split(microscopy, microscopy$sample_type)

# remove columns with NA (different categories for NT and T samples)
microscopy_list <- lapply(microscopy_list, function(x) 
  x[, apply(x, 2, function(y) !all(is.na(y)))]
)

#### species number at each reach on each day ####
species_rich_nt <- specnumber(microscopy_list$NT)
# alt
microscopy_species_richness <- lapply(microscopy_list, function(x)
  x %>% mutate(richness = rowSums(x[c(6:ncol(x))] > 0)))

# plot
richness_plot <- list()
for(i in 1:length(microscopy_species_richness)) {
  richness_plot[[i]] = ggplot(data = microscopy_species_richness[[i]], aes(x = site, y = richness)) +
                          geom_boxplot(aes(color = site)) +
                          labs(title = names(microscopy_species_richness[i]))
  print(richness_plot[[i]])
}
richness <- plot_grid(richness_plot[[1]], richness_plot[[2]], richness_plot[[3]])
richness

# plot over time
richness_plot_time <- list()
for(i in 1:length(microscopy_species_richness)) {
  richness_plot_time[[i]] = ggplot(data = microscopy_species_richness[[i]], aes(x = field_date, y = richness)) +
    facet_wrap(~year, scales = "free", ncol = 1) +
    geom_point(aes(color = site)) +
    labs(title = names(microscopy_species_richness[i]))
  print(richness_plot_time[[i]])
}
richness_nt <- microscopy_species_richness$NT %>% 
  select(field_date, site_reach, site, richness)
richness_nt <- left_join(richness_nt, environmental, by = c("field_date", "site_reach", "site"))

summary(aov(richness ~ site, data = richness_nt)) # significant difference based on river

##### diversity #####
microscopy_species_diversity <- lapply(microscopy_list, function(x) {
  x$diversity = diversity(x[,6:ncol(x)])
  return(x)})
diversity(microscopy_list$NT[,6:ncol(microscopy_list$NT)])

# plot
diversity_plot <- list()
for(i in 1:length(microscopy_species_diversity)) {
  diversity_plot[[i]] = ggplot(data = microscopy_species_diversity[[i]], aes(x = site, y = diversity)) +
    geom_boxplot(aes(color = site)) +
    labs(title = names(microscopy_species_diversity[i]))
  print(diversity_plot[[i]])
}
diversity <- plot_grid(diversity_plot[[1]], diversity_plot[[2]], diversity_plot[[3]])
diversity


diversity_plot_time <- list()
for(i in 1:length(microscopy_species_diversity)) {
  diversity_plot_time[[i]] = ggplot(data = microscopy_species_diversity[[i]], aes(x = field_date, y = diversity)) +
    facet_wrap(~year, scales = "free", ncol = 1) +
    geom_point(aes(color = site)) +
    labs(title = names(microscopy_species_diversity[i]))
  print(diversity_plot_time[[i]])
}
diversity_nt <- microscopy_species_diversity$NT %>% 
  select(field_date, site_reach, site, diversity)
diversity_nt <- left_join(diversity_nt, environmental, by = c("field_date", "site_reach", "site"))

summary(aov(diversity ~ site, data = diversity_nt)) # significant difference based on river

#### composition ####
nt_2022 <- microscopy_list$NT %>% 
  filter(year == 2022)
ac_2022 <- microscopy_list$TAC %>% 
  filter(year == 2022)

nt_2022_site <- split(nt_2022, nt_2022$site)
ac_2022_site <- split(ac_2022, ac_2022$site)
env_site <- split(environmental, environmental$site)

adonis2(nt_2022[,6:ncol(microscopy_list$NT)] ~ site, data = environmental %>% 
                            filter(year(field_date) == 2022))

# what about composition over sampling dates 
adonis2(nt_2022_site$RUS[,6:ncol(microscopy_list$NT)] ~ field_date, data = env_site$RUS) 
adonis2(ac_2022_site$RUS[,6:ncol(microscopy_list$TAC)] ~ field_date, data = env_site$RUS %>% 
          filter(field_date %in% ac_2022_site$RUS$field_date & site_reach %in% ac_2022_site$RUS$site_reach))
# minorly significant, more significant at month level
# whereas TAC is consistent

## Russian 
# field date NT: 0.014
# month NT: 0.055
# field date TAC: 0.604
# month TAC: 0.491

adonis2(nt_2022_site$`SFE-M`[,6:ncol(microscopy_list$NT)] ~ month(field_date), data = env_site$`SFE-M` %>% filter(year(field_date) == 2022)) 
adonis2(ac_2022_site$`SFE-M`[,6:ncol(microscopy_list$TAC)] ~ field_date, data = env_site$`SFE-M` %>%
          filter(field_date %in% ac_2022_site$`SFE-M`$field_date & site_reach %in% ac_2022_site$`SFE-M`$site_reach))

# not a significant composition difference
# field date NT: 
# month NT: 0.291