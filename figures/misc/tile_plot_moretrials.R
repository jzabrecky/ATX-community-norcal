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
  mutate(grouping = case_when(relative_abundance >= 40 ~ ">=40%",
                              relative_abundance >= 30 & relative_abundance < 40 ~ "30-39%",
                              relative_abundance >= 20 & relative_abundance < 30  ~ "20-29%",
                              relative_abundance >= 10 & relative_abundance < 20 ~ "10-19%",
                              relative_abundance >= 5 & relative_abundance < 10 ~ "5-10%",
                              relative_abundance < 5 & relative_abundance > 0 ~ "<5%",
                              relative_abundance == 0 ~ "absent")) %>% 
  mutate(grouping_factor = factor(grouping, levels = c(">=40%", "30-39%", "20-29%", "10-19%",
                                                       "5-10%", "<5%", "absent")))

# split into list based on site
all_data_split <- split(all_data, all_data$site)

ggplot(all_data_split$`SFE-M`, aes(x = field_date, y = sample_type, fill = grouping_factor)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  facet_wrap(~taxa, ncol = 1) +
  geom_text(aes(label = round(relative_abundance, 1))) +      # add text on top of tile
  coord_fixed() +
  scale_fill_discrete(palette = c("#206349", "#368f6c", "#5aa185", "#97ccb7", "#ccf0e1", "#b0b0b0"))

# play with groupings more??? look at quartile ranges
summary(all_data$relative_abundance) # 3rd quatile is 9.38
# maybe do >25 and then groups of 5%

all_data <- all_data %>% 
  mutate(grouping2 = case_when(relative_abundance >= 25 ~ ">25%",
                              relative_abundance >= 20 & relative_abundance < 25 ~ "20-25%",
                              relative_abundance >= 15 & relative_abundance < 20  ~ "15-20%",
                              relative_abundance >= 10 & relative_abundance < 15 ~ "10-15%",
                              relative_abundance >= 5 & relative_abundance < 10 ~ "5-10%",
                              relative_abundance < 5 & relative_abundance > 0 ~ "<5%",
                              relative_abundance == 0 ~ "absent")) %>% 
  mutate(grouping2_factor = factor(grouping2, levels = c(">25%", "20-25%", "15-20%", "10-15%",
                                                       "5-10%", "<5%", "absent")))

# split into list based on site
all_data_split <- split(all_data, all_data$site)

ggplot(all_data_split$`SFE-M`, aes(x = field_date, y = sample_type, fill = grouping2_factor)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  facet_wrap(~taxa, ncol = 1) +
  geom_text(aes(label = round(relative_abundance, 1))) +      # add text on top of tile
  coord_fixed() +
  scale_fill_discrete(palette = c("#206349", "#368f6c", "#5aa185", "#97ccb7", "#bae6d4", "#ccf0e1", "#b0b0b0"))

ggplot(all_data_split$`RUS`, aes(x = field_date, y = sample_type, fill = grouping2_factor)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  facet_wrap(~taxa, ncol = 1) +
  geom_text(aes(label = round(relative_abundance, 1))) +      # add text on top of tile
  coord_fixed() +
  scale_fill_discrete(palette = c("#206349", "#368f6c", "#5aa185", "#97ccb7", "#bae6d4", "#ccf0e1", "#b0b0b0"))

ggplot(all_data_split$`SAL`, aes(x = field_date, y = sample_type, fill = grouping2_factor)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  facet_wrap(~taxa, ncol = 1) +
  geom_text(aes(label = round(relative_abundance, 1))) +      # add text on top of tile
  coord_fixed() +
  scale_fill_discrete(palette = c("#206349", "#368f6c", "#5aa185", "#97ccb7", "#bae6d4", "#ccf0e1", "#b0b0b0"))
