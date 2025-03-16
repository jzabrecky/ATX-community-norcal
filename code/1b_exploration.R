## libraries

library(tidyverse)
library(lubridate)


## looking at data
tm <- read.csv("./data/morphological/tm_microscopy_with_covar.csv")
tac <- read.csv("./data/morphological/tac_microscopy_with_covar.csv")

eel_tm <- tm %>% 
  mutate(field_date = ymd(field_date),
         year = year(field_date),
         oPhos_mg_P_L = oPhos_ug_P_L / 1000,
         amm_plus_nitrate = ammonium_mg_N_L + nitrate_mg_N_L,
         n_to_p = amm_plus_nitrate / oPhos_mg_P_L) %>% # hypothetical?
  filter(site == "SFE-M" | site == "SFE-SH")

eel_tac <- tac %>% 
  mutate(field_date = ymd(field_date),
         year = year(field_date),
         oPhos_mg_P_L = oPhos_ug_P_L / 1000,
         amm_plus_nitrate = ammonium_mg_N_L + nitrate_mg_N_L,
         n_to_p = amm_plus_nitrate / oPhos_mg_P_L) %>% # hypothetical?
  filter(site == "SFE-M" | site == "SFE-SH")


plot <- ggplot(data = eel, aes(x = field_date, y = amm_plus_nitrate)) +
  geom_point() +
  facet_wrap(~year, scale = "free")
plot

plot <- ggplot(data = eel, aes(x = field_date, y = oPhos_ug_P_L)) +
  geom_point() +
  facet_wrap(~year, scale = "free")
plot

plot <- ggplot(data = eel, aes(x = field_date, y = n_to_p)) +
  geom_point() +
  facet_wrap(~year, scale = "free")
plot

plot <- ggplot(data = eel_tac, aes(x = field_date)) +
  geom_point(aes(y = amm_plus_nitrate), color = "blue") +
  geom_point(aes(y = e_diatoms / 500), color = "brown") +
  scale_y_continuous(sec.axis = (~ . * 500)) +
  ggtitle("Anabaena Samples") +
  facet_wrap(~year, scale = "free")
plot

plot <- ggplot(data = eel_tm, aes(x = field_date)) +
  geom_point(aes(y = amm_plus_nitrate), color = "blue") +
  geom_point(aes(y = e_diatoms / 500), color = "brown") +
  scale_y_continuous(sec.axis = (~ . * 500)) +
  ggtitle("Microcoleus Samples") +
  facet_wrap(~year, scale = "free")
plot

plot(eel$amm_plus_nitrate, eel$ATX_all_ug_g)
cor(eel$e_diatoms, eel$ATX_all_ug_g)