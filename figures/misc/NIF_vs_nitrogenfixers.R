# nif versus observed nitrogen fixers (both microscopically & molecularly)

library(tidyverse)

#### TM ####

# algal
tm_algal <- read.csv("./data/morphological/tm_algalonly.csv ") %>% 
  mutate(nfixers = anabaena_and_cylindrospermum + e_diatoms + calothrix + gloeotrichia + nodularia + nostoc + rivularia) %>% 
  mutate(field_date = ymd(field_date)) %>% 
  filter(year(field_date) == 2022)

# 16s
tm_16s <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_TM_nomicro.csv") #%>% 
  filter(order == "Nostocales" | order == "Endosymbiotic Diazoplast") %>% 
  select(site_reach, site, field_date, order, picrust2_relative_abundance) %>% 
  group_by(site_reach, site, field_date, order) %>% 
  dplyr::summarize(total = sum(picrust2_relative_abundance)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "order", values_from = "total") %>% 
  mutate(field_date = mdy(field_date),
         total_nfixers = case_when(is.na(`Endosymbiotic Diazoplast`) ~ Nostocales,
                                   TRUE ~ `Endosymbiotic Diazoplast` + Nostocales))

# functional
tm_picrust <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tm_nomicro.csv") %>% 
  pivot_wider(names_from = my_grouping, values_from = predicted_gene_abundance) %>% 
  mutate(field_date = mdy(field_date)) 
# pivot wider

# join
tm <- left_join(tm_algal, tm_picrust, by = c("site", "site_reach", "field_date"))
tm <- left_join(tm, tm_16s, by = c("site", "site_reach", "field_date"))


# plots
ggplot(data = tm, aes(x = nfixers, y = nitrogen_fixation)) +
  geom_point() +
  facet_wrap(~site, scale = "free")

ggplot(data = tm, aes(x = total_nfixers, y = nitrogen_fixation)) +
  geom_point() +
  facet_wrap(~site, scale = "free")

#### TAC ####


# algal
tac_algal <- read.csv("./data/morphological/tac_algalonly_noanacyl.csv ") %>% 
  mutate(nfixers = e_diatoms + calothrix + gloeotrichia + nodularia + nostoc + rivularia) %>% 
  mutate(field_date = ymd(field_date)) %>% 
  filter(year(field_date) == 2022)

# 16s
tac_16s <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_tac_noanacyl.csv") %>% 
  filter(order == "Nostocales" | order == "Endosymbiotic Diazoplast") %>% 
  select(site_reach, site, field_date, order, picrust2_relative_abundance) %>% 
  group_by(site_reach, site, field_date, order) %>% 
  dplyr::summarize(total = sum(picrust2_relative_abundance)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "order", values_from = "total") %>% 
  mutate(field_date = mdy(field_date),
         total_nfixers = case_when(is.na(`Endosymbiotic Diazoplast`) ~ Nostocales,
                                   TRUE ~ `Endosymbiotic Diazoplast` + Nostocales))

# functional
tac_picrust <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tac_noanacyl.csv") %>% 
  pivot_wider(names_from = my_grouping, values_from = predicted_gene_abundance) %>% 
  mutate(field_date = mdy(field_date)) 
# pivot wider

# join
tac <- left_join(tac_algal, tac_picrust, by = c("site", "site_reach", "field_date"))
tac <- left_join(tac, tac_16s, by = c("site", "site_reach", "field_date"))


# plots
ggplot(data = tac, aes(x = nfixers, y = nitrogen_fixation)) +
  geom_point(size = 3) +
  facet_wrap(~site, scale = "free")

ggplot(data = tac, aes(x = total_nfixers, y = nitrogen_fixation)) +
  geom_point(size = 3) +
  facet_wrap(~site, scale = "free")

#### NT ####

# algal
nt_algal <- read.csv("./data/morphological/nt_algalonly.csv ") %>% 
  mutate(nfixers = anabaena_and_cylindrospermum + epithemia + rhopalodia + calothrix + gloeotrichia + nostoc + rivularia) %>% 
  mutate(field_date = ymd(field_date)) %>% 
  filter(year(field_date) == 2022)

# 16s
nt_16s <- read.csv("./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_FINAL.csv") %>% 
  filter(sample_type == "NT") %>% 
  filter(order == "Nostocales" | order == "Endosymbiotic Diazoplast") %>% 
  select(site_reach, site, field_date, order, picrust2_relative_abundance) %>% 
  group_by(site_reach, site, field_date, order) %>% 
  dplyr::summarize(total = sum(picrust2_relative_abundance)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "order", values_from = "total") %>% 
  mutate(field_date = mdy(field_date),
         total_nfixers = case_when(is.na(`Endosymbiotic Diazoplast`) ~ Nostocales,
                                   TRUE ~ `Endosymbiotic Diazoplast` + Nostocales))

# functional
nt_picrust <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select.csv") %>% 
  filter(sample_type == "NT") %>% 
  pivot_wider(names_from = my_grouping, values_from = predicted_gene_abundance) %>% 
  mutate(field_date = mdy(field_date)) 
# pivot wider

# join
nt <- left_join(nt_algal, nt_picrust, by = c("site", "site_reach", "field_date"))
nt <- left_join(nt, nt_16s, by = c("site", "site_reach", "field_date"))


# plots
ggplot(data = nt, aes(x = nfixers, y = nitrogen_fixation)) +
  geom_point(size = 3) +
  facet_wrap(~site, scale = "free")

ggplot(data = nt, aes(x = total_nfixers, y = nitrogen_fixation)) +
  geom_point(size = 3) +
  facet_wrap(~site, scale = "free")
