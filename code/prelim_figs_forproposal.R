### preliminary chapter 2 figures ####

library(tidyverse)
library(ggplot2)

TM_microscopy <- read.csv("./data/EDI_data_package/microscopy_target_samples.csv") %>% 
  mutate(field_date = ymd(field_date)) %>% 
  mutate(year = year(field_date))

TM_summary <- TM_microscopy %>% 
  filter(sample_type == "TM") %>% 
  dplyr::group_by(field_date, site_reach, year) %>% 
  dplyr::summarize(microcoleus = mean(microcoleus),
            anabaena = mean(anabaena_and_cylindrospermum),
            n_fixers = mean(anabaena_and_cylindrospermum) + mean(gloeotrichia) + mean(rivularia) + 
              mean(nostoc) + mean(calothrix) + mean(tolypothrix) + mean(e_diatoms),
            e_diatoms = mean())

anatoxins <- read.csv("./data/cyano_atx.csv") %>% 
  mutate(field_date = ymd(field_date)) %>% 
  filter(sample_type == "TM")

TM_all <- left_join(TM_summary, anatoxins, by = c("site_reach", "field_date"))

n_fix_time <- ggplot(data = TM_summary, aes(x = field_date, y = n_fixers)) +
  geom_point() +
  facet_wrap(~year)
n_fix_time

# fill na
TM_all$ATX_all_ug_orgmat_g <- replace_na(TM_all$ATX_all_ug_orgmat_g, 0)

n_fix_atx <- ggplot(data = TM_all, aes(x = n_fixers, y = ATX_all_ug_orgmat_g)) +
  geom_point() +
  scale_y_log10()
n_fix_atx

cor.test(TM_all$n_fixers, y = TM_all$ATX_all_ug_orgmat_g)

daily_avg <- TM_all %>%
  filter(site_reach == "SFE-M-1S" | site_reach == "SFE-M-2" | site_reach == "SFE-M-3"
         | site_reach == "SFE-M-4" | site_reach == "SFE-SH-1S") %>% 
  group_by(year, field_date) %>% 
  summarize(n_fixers = mean(n_fixers),
            e_diatoms = mean(e_diatoms),
            atx = mean(ATX_all_ug_orgmat_g))

ggplot(data = daily_avg, aes(x = field_date, y = atx)) +
  geom_point(aes(y = e_diatoms), color = "red") +
  scale_y_log10() +
  facet_wrap(~year) +
  geom_point()
