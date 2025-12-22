library(tidyverse)
library(vegan)

# REDO WITH STANDARDIZED COVARIATES!?!?!

# load data- only care about community and not environmental here
tm <- read.csv("./data/morphological/tm_algalonly.csv") %>% 
  mutate(month = month(field_date),
         year = year(field_date))
tac <- read.csv("./data/morphological/tac_algalonly.csv") %>% 
  mutate(month = month(field_date),
         year = year(field_date))
nt <- read.csv("./data/morphological/nt_algalonly.csv") %>% 
  mutate(month = month(field_date),
         year = year(field_date))
env <- read.csv("./data/field_and_lab/environmental_covariates.csv")

# function to see which columns have nothing other than zero
check_zeros <- function(df) {
  for(i in 1:ncol(df)) {
    check <- any(which(df[i] > 0))
    if(check == FALSE) {
      return(colnames(df)[i])
    }
  }
}
## microcoleus in SFE
check_zeros(tm %>% filter(site == "SFE-M" | site == "SFE-SH"))
tm <- tm %>% 
  select(!aphanothece)

# join in microcoleus with environmental parameters
tm_env <- left_join(tm, env, by = c("field_date", "site_reach", "site"))

tm_sub_sfk <- tm_env %>% 
  filter(site == "SFE-M" | site == "SFE-SH") %>% 
  select(c(anabaena_and_cylindrospermum:unknown))
tm_matrix <- as.matrix(tm_sub_sfk)
env_sub_sfk <- tm_env %>% 
  filter(site == "SFE-M" | site == "SFE-SH") %>% 
  select(assumed_pH, temp_C, DO_mg_L, cond_uS_cm, nitrate_mg_N_L, ammonium_mg_N_L,
         TDC_mg_L, DOC_mg_L, TM_ATXa_ug_g, TM_dhATXa_ug_g, TM_ATX_all_ug_orgmat_g)
# curious if environmental parameters should be log-scaled or standardized

nmds_tm <- metaMDS(tm_matrix,
                   distance = "bray",
                   trymax = 500,
                   autotransform = TRUE)
tm_env <- envfit(nmds_tm, env_sub_sfk, permutations = 999, na.rm = TRUE)
tm_env

env_coord <- as.data.frame(scores(tm_env, "vectors")) * ordiArrowMul(tm_env)

nmds_scores <- as.data.frame(scores(nmds_tm, "sites"))
nmds_data <- cbind(nmds_scores, tm %>% filter(site == "SFE-M" | site == "SFE-SH") %>% 
                     select(site, site_reach, field_date, month, year)) %>% 
  mutate(month = as.factor(month))

tm_sfk <- ggplot(nmds_data, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = month, shape = site_reach), size = 4) +
  stat_ellipse(aes(color = month), type = "t", linetype = 2, size = 0.5) +
  theme_minimal() +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = env_coord, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = env_coord, aes(x = NMDS1, y = NMDS2), colour = "grey30",
            fontface = "bold", label = row.names(env_coord))
  facet_wrap(~year, ncol = 1)
tm_sfk

# do communities change based on anatoxins?
adonis2(tm_sub_sfk ~ TM_ATX_all_ug_orgmat_g, data = env_sub_sfk)
# yes..... p = 0.01 but what characterizes this change?
# also probably need to do this on standardized data