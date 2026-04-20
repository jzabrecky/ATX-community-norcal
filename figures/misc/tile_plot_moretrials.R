# significant time taxa from script 5a
# nostoc, anabaena, epithemia, microcoleus, geitlerinema, nodularia (may drop this one)
# dropping nodularia bc it's not recorded in NT samples

# read in algal data
tac <- read.csv("./data/morphological/tac_algalonly_noanacyl.csv") %>% 
  select(site, site_reach, field_date, sample_type, nostoc, e_diatoms, microcoleus, leptolyngbya_geitlerinema) %>% 
  filter(year(ymd(field_date)) == 2022)
tm <- read.csv("./data/morphological/tm_algalonly_nomicro.csv") %>% 
  select(site, site_reach, field_date, sample_type, nostoc, e_diatoms, anabaena_and_cylindrospermum, leptolyngbya_geitlerinema) %>% 
  filter(year(ymd(field_date)) == 2022)
nt <- read.csv("./data/morphological/nt_algalonly.csv") %>% 
  select(site, site_reach, field_date, sample_type, nostoc, epithemia, anabaena_and_cylindrospermum, microcoleus, leptolyngbya_geitlerinema) %>% 
  dplyr::rename(e_diatoms = epithemia) %>% 
  filter(year(ymd(field_date)) == 2022)

# list
data <- list(tac, tm, nt)

# pivot longer
data <- lapply(data, function(x) x %>% pivot_longer(cols = c(5:ncol(x)), names_to = "taxa", values_to = "relative_abundance"))

# average across site reach
data <- lapply(data, function(x) x %>% dplyr::group_by(site, field_date, sample_type, taxa) %>% 
                 dplyr::summarize(relative_abundance = mean(relative_abundance)))

# join into one
all_data <- rbind(data[[1]], data[[2]], data[[3]]) %>% 
  mutate(grouping = case_when(relative_abundance >= 40 ~ ">= 40%",
                              relative_abundance >= 30 & relative_abundance < 40 ~ "30-39%",
                              relative_abundance >= 20 & relative_abundance < 30  ~ "20-29%",
                              relative_abundance >= 10 & relative_abundance < 20 ~ "10-19%",
                              relative_abundance < 10 & relative_abundance > 0 ~ "<10%",
                              relative_abundance == 0 ~ "absent"))

# split into list based on site
all_data_split <- split(all_data, all_data$site)

ggplot(all_data_split$`SFE-M`, aes(x = field_date, y = sample_type, fill = grouping)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  facet_wrap(~taxa, ncol = 1) +
  coord_fixed() +
  scale_fill_discrete(palette = c("#97ccb7", "#5aa185", "#368f6c", "#124732", "#b0b0b0"))

# play with groupings more??? look at quartile ranges

# continuous groupings
ggplot(all_data_split$`SFE-M`, aes(x = field_date, y = sample_type, fill = relative_abundance)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  facet_wrap(~taxa, ncol = 1) +
  coord_fixed() +
  scale_fill_continuous()
