# waffle plots as discussed with Joanna

# remotes::install_github("hrbrmstr/waffle")

#### Algal ####
 
# read in algal data
tac <- read.csv("./data/morphological/tac_algalonly_noanacyl.csv")
tm <- read.csv("./data/morphological/tm_algalonly_nomicro.csv")
nt <- read.csv("./data/morphological/nt_algalonly.csv")

# species we care about

# 
library(ggplot2)
install.packages("plotly")
library(plotly)

# put together for nostoc
nostoc <- rbind(tac %>% select(site, site_reach, field_date, sample_type, nostoc),
                tm %>% select(site, site_reach, field_date, sample_type, nostoc),
                nt %>% select(site, site_reach, field_date, sample_type, nostoc)) %>% 
  mutate(present = case_when(nostoc > 0 ~ "y",
                             TRUE ~ "n")) %>% 
  filter(site == "SFE-M") %>% 
  filter(year(ymd(field_date)) == 2022)

ggplot(nostoc, aes(x = field_date, y = sample_type, fill = present)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  facet_wrap(~site_reach, ncol = 1) +
  scale_fill_discrete(palette = c("#E8E8E8", "#3DD17F")) + 
  coord_fixed()


epithemia <- rbind(tac %>% select(site, site_reach, field_date, sample_type, e_diatoms),
                tm %>% select(site, site_reach, field_date, sample_type, e_diatoms),
                nt %>% dplyr::rename(e_diatoms = epithemia) %>% 
                  select(site, site_reach, field_date, sample_type, e_diatoms)) %>% 
  mutate(present = case_when(e_diatoms > 20 ~ ">20%",
                             e_diatoms < 20 & e_diatoms > 10 ~ ">10%",
                             e_diatoms > 0 & e_diatoms < 10 ~ "<10%",
                             TRUE ~ "n")) %>% 
  filter(site == "SFE-M") %>% 
  filter(year(ymd(field_date)) == 2022)

ggplot(epithemia, aes(x = field_date, y = sample_type, fill = present)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  facet_wrap(~site_reach, ncol = 1) +
  scale_fill_discrete(palette = c("#E8E0BA", "#C4B676", "#8F7917")) + 
  coord_fixed()
