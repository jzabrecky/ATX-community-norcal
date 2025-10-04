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

tm_morpho_reach_sal <- ggplot(data = target_long_list$TM %>% 
                          filter(site == "SAL"), aes(x = field_date, y = percent, fill = group)) +
  geom_area() + 
  facet_wrap(~ site_reach + year, scales = "free")
tm_morpho_reach_sal

tm_morpho_reach_sfk <- ggplot(data = target_long_list$TM %>% 
                                filter(site == "SFE-M"), aes(x = field_date, y = percent, fill = group)) +
  geom_area() + 
  facet_wrap(~ site_reach + year, scales = "free")
tm_morpho_reach_sfk

tm_morpho_site <- ggplot(data = target_group_list$TM, aes(x = field_date, y = mean, fill = group)) +
  geom_area() + 
  facet_wrap(~ site + year, scales = "free")
tm_morpho_site

# tac plots
tac_morpho_reach_rus <- ggplot(data = target_long_list$TAC %>% 
                             filter(site == "RUS"), aes(x = field_date, y = percent, fill = group)) +
  geom_area() + 
  facet_wrap(~ site_reach + year, scales = "free")
tac_morpho_reach_rus

tac_morpho_reach_sfe <- ggplot(data = target_long_list$TAC %>% 
                                 filter(site == "SFE-M"), aes(x = field_date, y = percent, fill = group)) +
  geom_area() + 
  facet_wrap(~ site_reach + year, scales = "free")
tac_morpho_reach_sfe

tac_morpho_site <- ggplot(data = target_group_list$TAC, aes(x = field_date, y = mean, fill = group)) +
  geom_area() + 
  facet_wrap(~ site + year, scales = "free")
tac_morpho_site

# nt plots
nt_morpho_reach_rus <- ggplot(data = nt_long %>% 
                                filter(site == "RUS"), aes(x = field_date, y = percent, fill = group)) +
  geom_area() + 
  facet_wrap(~ site_reach + year, scales = "free")
nt_morpho_reach_rus

nt_morpho_reach_sfe <- ggplot(data = nt_long %>% 
                                filter(site == "SFE-M"), aes(x = field_date, y = percent, fill = group)) +
  geom_area() + 
  facet_wrap(~ site_reach + year, scales = "free")
nt_morpho_reach_sfe

nt_morpho_reach_sal <- ggplot(data = nt_long %>% 
                                filter(site == "SAL"), aes(x = field_date, y = percent, fill = group)) +
  geom_area() + 
  facet_wrap(~ site_reach + year, scales = "free")
nt_morpho_reach_sal

nt_morpho_site <- ggplot(data = nt_group, aes(x = field_date, y = mean, fill = group)) +
  geom_area() + 
  facet_wrap(~ site + year, scales = "free")
nt_morpho_site

#### molecular polygon plots

## have to readjust 

molecular_phylums <- molecular %>% 
  select(site_reach, site, field_date, year, sample_type, phylum, relative_abundance) %>% 
  group_by(site_reach, field_date, year, sample_type, phylum) %>% 
  dplyr::summarize(relative_abundance_total = sum(relative_abundance))
# need to figure this out

# tm 
tm_molec_reach <- ggplot(data = molecular_phylums %>% filter(sample_type == "TM"),
                        aes(x = field_date, y = relative_abundance_total, fill = phylum)) +
  geom_area() + 
  facet_wrap(~ site_reach + year, scales = "free")
tm_molec_reach

# tac
tac_molec_reach <- ggplot(data = molecular_phylums %>% filter(sample_type == "TAC"),
                          aes(x = field_date, y = relative_abundance_total, fill = phylum)) +
  geom_area() + 
  facet_wrap(~ site_reach + year, scales = "free")
tac_molec_reach

# nt
nt_molec_reach <- ggplot(data = molecular_phylums %>% filter(sample_type == "NT"),
                          aes(x = field_date, y = relative_abundance_total, fill = phylum)) +
  geom_area() + 
  facet_wrap(~ site_reach + year, scales = "free")
nt_molec_reach

## class

molecular_class <- molecular %>% 
  select(site_reach, site, field_date, year, sample_type, class, relative_abundance) %>% 
  group_by(site_reach, field_date, year, sample_type, class) %>% 
  dplyr::summarize(relative_abundance_total = sum(relative_abundance))
  
# tm 
tm_molec_reach_class <- ggplot(data = molecular_class %>% filter(sample_type == "TM"),
                           aes(x = field_date, y = relative_abundance_total, fill = class)) +
  geom_area() + 
  facet_wrap(~ site_reach + year, scales = "free")
tm_molec_reach_class

# tac
tac_molec_reach_class <- ggplot(data = molecular_class %>% filter(sample_type == "TAC"),
                          aes(x = field_date, y = relative_abundance_total, fill = class)) +
  geom_area() + 
  facet_wrap(~ site_reach + year, scales = "free")
tac_molec_reach_class

# nt
nt_molec_reach_class <- ggplot(data = molecular_class %>% filter(sample_type == "NT"),
                         aes(x = field_date, y = relative_abundance_total, fill = class)) +
  geom_area() + 
  facet_wrap(~ site_reach + year, scales = "free")
nt_molec_reach_class

## cyanobacteria orders

## class

molecular_orders <- molecular %>% 
  filter(phylum == "Cyanobacteria") %>% 
  select(site_reach, site, field_date, year, sample_type, order, relative_abundance) %>% 
  group_by(site_reach, field_date, year, sample_type, order) %>% 
  dplyr::summarize(relative_abundance_total = sum(relative_abundance))

# recalculate relative abundance of orders within cyanobacteria phylum
molecular_order_total <- molecular_orders %>% 
  group_by(site_reach, field_date, year, sample_type) %>% 
  dplyr::summarize(total_cyano_phylum = sum(relative_abundance_total)) %>%
  ungroup()

# calculate relative abundance within cyanobacteria
molecular_orders <- left_join(molecular_orders, molecular_order_total, 
                              by = c("site_reach", "field_date", "year", "sample_type")) %>% 
  mutate(relative_abun_within_cyano = relative_abundance_total / total_cyano_phylum)

# tm 
tm_molec_reach_order <- ggplot(data = molecular_orders %>% filter(sample_type == "TM"),
                               aes(x = field_date, y = relative_abun_within_cyano, fill = order)) +
  geom_area() + 
  facet_wrap(~ site_reach + year, scales = "free")
tm_molec_reach_order

# tac
tac_molec_reach_class <- ggplot(data = molecular_orders %>% filter(sample_type == "TAC"),
                                aes(x = field_date, y = relative_abun_within_cyano, fill = order)) +
  geom_area() + 
  facet_wrap(~ site_reach + year, scales = "free")
tac_molec_reach_class

# nt
nt_molec_reach_order <- ggplot(data = molecular_orders %>% filter(sample_type == "NT"),
                               aes(x = field_date, y = relative_abun_within_cyano, fill = order)) +
  geom_area() + 
  facet_wrap(~ site_reach + year, scales = "free")
nt_molec_reach_order

# will still try to figure out polygon plots ugh not sure what is up


# curious about raw read for RUSSIAN TAC 9-15
