#### how does presence of other toxin producers relate to atx concentrations?

## defining anatoxin producer from Christensen and Khan list
## anabaena, aphanozemonon, aphanothece, arthrospira, cylindrospermum, limnothrix,
## lynbya, nostoc, oscillatoria, tychonema, planktothrix, phormidium, microcoleus

# note HCB-2 list is cylindrospermum, geitlerinema, kamptonema, microcoleus,
# nostoc, oscillatoria, tychonema

#### microscopy ####
TM <- read.csv("./data/morphological/tm_algalonly.csv") %>% 
  mutate(field_date = ymd(field_date),
         year = year(field_date)) %>% 
  select(site, site_reach, sample_type, field_date, year, anabaena_and_cylindrospermum,
         nostoc, geitlerinema, microcoleus, phormidium_unknown, 
         oscillatoria, lyngbya, aphanothece)

TM_long <- pivot_longer(TM, cols = (6:ncol(TM)), names_to = "taxon", values_to = "percent")

ggplot(TM_long %>% 
         filter(taxon != "microcoleus"), aes(x=field_date, y=percent, fill=taxon)) + 
  geom_area() +
  facet_wrap(~site_reach + year, scales = "free_x")

TAC <- read.csv("./data/morphological/tac_algalonly.csv") %>% 
  mutate(field_date = ymd(field_date),
         year = year(field_date)) %>% 
  select(site, site_reach, sample_type, field_date, year, anabaena_and_cylindrospermum,
         nostoc, geitlerinema, microcoleus, phormidium_unknown, 
         oscillatoria, lyngbya, aphanothece)

TAC_long <- pivot_longer(TAC, cols = (6:ncol(TM)), names_to = "taxon", values_to = "percent")

ggplot(TAC_long %>% 
         filter(taxon != "anabaena_and_cylindrospermum"), aes(x=field_date, y=percent, fill=taxon)) + 
  geom_area() +
  facet_wrap(~site_reach + year, scales = "free_x")

response <- read.csv("./data/field_and_lab/environmental_covariates_and_toxins.csv") %>% 
  mutate(field_date = ymd(field_date),
         year = year(field_date))

TAC_summed <- TAC_long %>% 
  dplyr::group_by(field_date, site_reach) %>% 
  dplyr::summarize(total_toxin_producers = sum(percent)) %>% 
  left_join(response, by = c("site_reach", "field_date"))

TAC_joined <- TAC_long %>% 
  left_join(response, by = c("site_reach", "field_date"))

ggplot(TAC_summed, aes(x = total_toxin_producers, y = TAC_ATX_all_ug_orgmat_g)) +
  geom_point(aes(color = site)) +
  facet_wrap(~site_reach + year, scales = "free")

ggplot(TAC_joined, 
       aes(x = percent, y = TAC_ATX_all_ug_orgmat_g)) +
  geom_point(aes(color = taxon)) +
  facet_wrap(~site_reach + year.x, scales = "free")

TM_summed <- TM_long %>% 
  filter(taxon != "microcoleus") %>% 
  dplyr::group_by(field_date, site_reach) %>% 
  dplyr::summarize(total_toxin_producers = sum(percent)) %>% 
  left_join(response, by = c("site_reach", "field_date"))

TM_joined <- TM_long %>% 
  left_join(response, by = c("site_reach", "field_date"))

ggplot(TM_summed, aes(x = total_toxin_producers, y = TM_ATX_all_ug_orgmat_g)) +
  geom_point(aes(color = site)) +
  facet_wrap(~site_reach + year, scales = "free")

ggplot(TM_joined, 
       aes(x = percent, y = TM_ATX_all_ug_orgmat_g)) +
  geom_point(aes(color = taxon)) +
  facet_wrap(~site_reach + year.x, scales = "free")
