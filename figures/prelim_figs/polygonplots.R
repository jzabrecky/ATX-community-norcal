#### preliminary figure
# polygon plots through time

#### (1) Loading data ####

library(tidyverse)

# molecular
molecular <- read.csv("./data/molecular/16s_nochimera_rarefied_95_processed.csv") %>% 
  mutate(field_date = mdy(field_date),
         year = year(field_date))


# morphological
tm <- read.csv("./data/morphological/tm_algalonly_with_covar.csv")
tac <- read.csv("./data/morphological/tac_algalonly_with_covar.csv")
target <- rbind(tm, tac) %>% 
  mutate(field_date = ymd(field_date),
         year = year(field_date))
nt <- read.csv("./data/morphological/nt_algalonly_with_covar.csv") %>% 
  mutate(field_date = ymd(field_date),
         year = year(field_date))

# pivot longer
target_long <- target[,c(1:26, ncol(target))]
target_long <- pivot_longer(target_long, cols = c(5:26), values_to = "percent", names_to = "group")
nt_long <- nt[,c(1:47, ncol(nt))]
nt_long <- pivot_longer(nt_long, cols = c(5:47), values_to = "percent", names_to = "group")

# group across year and site
target_group <- target_long %>% 
  dplyr::group_by(field_date, site, year, sample_type, group) %>% 
  dplyr::summarize(mean = mean(percent))
nt_group <- nt_long %>% 
  dplyr::group_by(field_date, site, year, sample_type, group) %>% 
  dplyr::summarize(mean = mean(percent))
# make into list
target_long_list <- split(target_long, target_long$sample_type)
target_group_list <- split(target_group, target_group$sample_type)
# remove smallest categories and put into other???


#### polygon plots 

# to do: plot vertical line for highest atx, consider putting minimized taxa into an other bin

theme_set(theme_bw())

tm_morpho_reach <- ggplot(data = target_long_list$TM, aes(x = field_date, y = percent, fill = group)) +
  geom_area() + 
  facet_wrap(~ site_reach + year, scales = "free")
tm_morpho_reach

tm_morpho_site <- ggplot(data = target_group_list$TM, aes(x = field_date, y = mean, fill = group)) +
  geom_area() + 
  facet_wrap(~ site + year, scales = "free")
tm_morpho_site

# tac plots
tac_morpho_reach <- ggplot(data = target_long_list$TAC, aes(x = field_date, y = percent, fill = group)) +
  geom_area() + 
  facet_wrap(~ site_reach + year, scales = "free")
tac_morpho_reach

tac_morpho_site <- ggplot(data = target_group_list$TAC, aes(x = field_date, y = mean, fill = group)) +
  geom_area() + 
  facet_wrap(~ site + year, scales = "free")
tac_morpho_site

# nt plots
nt_morpho_site <- ggplot(data = nt_long, aes(x = field_date, y = percent, fill = group)) +
  geom_area() + 
  facet_wrap(~ site_reach + year, scales = "free")
nt_morpho_site

nt_morpho_site <- ggplot(data = nt_group, aes(x = field_date, y = mean, fill = group)) +
  geom_area() + 
  facet_wrap(~ site + year, scales = "free")
nt_morpho_site

#### molecular polygon plots

molecular_phylums <- molecular %>% 
  tidyr::complete(phylum, nesting(field_date, sample_type, site_reach, year),
                  fill = list(relative_abundance = 0))

# tm 
tm_molec_reach <- ggplot(data = molecular_phylums,
                        aes(x = field_date, y = relative_abundance, fill = phylum)) +
  geom_area() + 
  facet_wrap(~ site_reach + year, scales = "free")
tm_molec_reach

tm_bars_reach <- ggplot(data = molecular_phylums %>% filter(sample_type == "TM"),
                         aes(x = field_date, y = relative_abundance, fill = phylum)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~ site_reach + year, scales = "free")
tm_bars_reach

tac_bars_reach <- ggplot(data = molecular_phylums %>% filter(sample_type == "TAC"),
                        aes(x = field_date, y = relative_abundance, fill = phylum)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~ site_reach + year, scales = "free")
tac_bars_reach

# will still try to figure out polygon plots ugh not sure what is up


# curious about raw read for RUSSIAN TAC 9-15