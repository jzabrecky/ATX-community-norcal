#### Nitrogen fixers


#### microscopy ####
NT <- read.csv("./data/morphological/nt_algalonly.csv") %>% 
  mutate(field_date = ymd(field_date),
         year = year(field_date)) %>% 
  select(site, site_reach, sample_type, field_date, year, anabaena_and_cylindrospermum, epithemia,
         nostoc, calothrix, gloeotrichia, rhopalodia, lyngbya, rivularia, scytonema)

NT_long <- pivot_longer(NT, cols = (6:ncol(NT)), names_to = "taxon", values_to = "percent")

ggplot(NT_long, aes(x=field_date, y=percent, fill=taxon)) + 
  geom_area() +
  facet_wrap(~site_reach + year, scales = "free_x")

TM <- read.csv("./data/morphological/tm_algalonly.csv") %>% 
  mutate(field_date = ymd(field_date),
         year = year(field_date)) %>% 
  select(site, site_reach, sample_type, field_date, year, anabaena_and_cylindrospermum, e_diatoms,
         nostoc, calothrix, tolypothrix, gloeotrichia, lyngbya, rivularia, nodularia, scytonema)

TM_long <- pivot_longer(TM, cols = (6:ncol(TM)), names_to = "taxon", values_to = "percent")

ggplot(TM_long, aes(x=field_date, y=percent, fill=taxon)) + 
  geom_area() +
  facet_wrap(~site_reach + year, scales = "free_x")

TAC <- read.csv("./data/morphological/tac_algalonly.csv") %>% 
  mutate(field_date = ymd(field_date),
         year = year(field_date)) %>% 
  select(site, site_reach, sample_type, field_date, year, anabaena_and_cylindrospermum, e_diatoms,
         nostoc, calothrix, tolypothrix, gloeotrichia, lyngbya, rivularia, nodularia, scytonema)

TAC_long <- pivot_longer(TAC, cols = (6:ncol(TAC)), names_to = "taxon", values_to = "percent")

ggplot(TAC_long, aes(x=field_date, y=percent, fill=taxon)) + 
  geom_area() +
  facet_wrap(~site_reach + year, scales = "free_x")

response <- read.csv("./data/field_and_lab/environmental_covariates_and_toxins.csv") %>% 
  mutate(field_date = ymd(field_date),
         year = year(field_date))

NT_summed <- NT_long %>% 
  dplyr::group_by(field_date, site_reach) %>% 
  dplyr::summarize(total_n_fixers = sum(percent)) %>% 
  left_join(response, by = c("site_reach", "field_date"))

ggplot(NT_summed, aes(x = total_n_fixers, y = TM_ATX_all_ug_orgmat_g)) +
  geom_point(aes(color = site)) +
  facet_wrap(~site_reach + year, scales = "free")

ggplot(NT_summed, aes(x = total_n_fixers, y = TAC_ATX_all_ug_orgmat_g)) +
  geom_point(aes(color = site)) +
  facet_wrap(~site_reach + year, scales = "free")

TAC_summed <- TAC_long %>% 
  dplyr::group_by(field_date, site_reach) %>% 
  dplyr::summarize(total_n_fixers = sum(percent)) %>% 
  left_join(response, by = c("site_reach", "field_date"))

ggplot(TAC_summed, aes(x = total_n_fixers, y = TAC_ATX_all_ug_orgmat_g)) +
  geom_point(aes(color = site)) +
  facet_wrap(~site_reach + year, scales = "free")

TM_summed <- TM_long %>% 
  dplyr::group_by(field_date, site_reach) %>% 
  dplyr::summarize(total_n_fixers = sum(percent)) %>% 
  left_join(response, by = c("site_reach", "field_date"))

ggplot(TM_summed, aes(x = total_n_fixers, y = TM_ATX_all_ug_orgmat_g)) +
  geom_point(aes(color = site)) +
  facet_wrap(~site_reach + year, scales = "free")

 
#### molecular