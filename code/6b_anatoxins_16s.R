

# set seed for reproducibility
set.seed(2025)

# libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot"), require, character.only = T)

# read in files (data transformed in previous script, "4a_amongrivers_microscopy.R")
data <- lapply(list.files(path = "./data/molecular/transformed/", pattern = ".csv"),
               function(x) read.csv(paste("./data/molecular/transformed/", x, sep = "")))
names(data) <- c("nt", "tac", "tm")

# fix stupid date
data <- lapply(data, function(x) x = x %>% 
                 mutate(field_date = mdy(field_date)))

env_tog <- read.csv("./data/field_and_lab/env_tox_standardized_together.csv") %>% 
  mutate(field_date = ymd(field_date))
env_sep <- read.csv("./data/field_and_lab/env_tox_standardized_byriver.csv") %>% 
  mutate(field_date = ymd(field_date))

# lastly, merge with our community data to get rows formatted the same
data_tog <- lapply(data, function(x) left_join(x, env_tog, by = c("site", "site_reach", "field_date")))
data_sep <- lapply(data, function(x) left_join(x, env_sep, by = c("site", "site_reach", "field_date")))

# for microcoleus & anabaena remove rows where there are NAs for anatoxins
which(is.na(data_tog$tm$TM_ATX_all_ug_orgmat_g)) #22
which(is.na(data_tog$tac$TAC_ATX_all_ug_orgmat_g)) #none!

# 22 was a "not sure if this is microcoleus sample" and there was not enough material
# to analyze for anatoxins
data_tog$tm <- data_tog$tm[-7,]
data_sep$tm <- data_sep$tm[-7,]

# split into lists by river
for(i in 1:length(data_sep)) {
  data_sep[[i]] <- split(data_sep[[i]], data_sep[[i]]$site)
}

#### Attempting TITAN2

library("TITAN2")
## (a) TM

# first establish that anatoxin is associated with significant changes in composition
adonis2(data = data_tog$tm, data_tog$tm[,6:ncol(data$tm)] ~ 
          TM_ATX_all_ug_orgmat_g) # yes**

# rivers separately
adonis2(data = data_sep$tm$`SFE-M`, data_sep$tm$`SFE-M`[,6:ncol(data$tm)] ~ 
          TM_ATX_all_ug_orgmat_g) # yes**
adonis2(data = data_sep$tm$`SAL`, data_sep$tm$`SAL`[,6:ncol(data$tm)] ~ 
          TM_ATX_all_ug_orgmat_g)
# no, but only one sample had small amounts of atx so dropping this site

## (b) TAC

# first establish that anatoxin is associated with significant changes in composition
adonis2(data = data_tog$tac, data_tog$tac[,6:ncol(data$tac)] ~ 
          TAC_ATX_all_ug_orgmat_g) # yes*

# rivers separately
adonis2(data = data_sep$tac$`SFE-M`, data_sep$tac$`SFE-M`[,6:ncol(data$tac)] ~ 
          TAC_ATX_all_ug_orgmat_g) # yes**
adonis2(data = data_sep$tac$`RUS`, data_sep$tac$`RUS`[,6:ncol(data$tac)] ~ 
          TAC_ATX_all_ug_orgmat_g) # no
# no, but only one sample had small amounts of atx so dropping this site

test <- data_tog$tm

titan.PDIR <- titan(env = test$TM_ATX_all_ug_orgmat_g,
                    txa =  test[,5:500])
plot_taxa_ridges(titan.PDIR,
                 xlabel = "PDIR")

# okay not sure if this is the best

### attempting differential abundance analysis

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ALDEx2")

# left off here

# Q from Grant Johnson: 

# how does bacterial diversity scale with anatoxins?

# calculate diversity for each dataframe
diversity <- lapply(data, function(x) {
  x = x %>% 
    mutate(shannon_diversity = calc_diversity(x, start_col)) %>% 
    select(field_date, site_reach, site, shannon_diversity)})

# plot diversity as boxplots
for(i in 1:length(diversity)) {
  boxplot = ggplot(data = diversity[[i]], aes(x = site, y = shannon_diversity, fill = site)) +
    geom_boxplot()
  print(boxplot)
}

#### Grant Q: Does diversity differ across rivers?



# calculate diversity for each dataframe
diversity <- lapply(list.files(path = "./data/molecular/shannon_diversity/", pattern = ".csv"),
                     function(x) read.csv(paste("./data/molecular/shannon_diversity/", x, sep = "")))
names(diversity) <- c("nt", "tac", "tm")

# plot diversity as boxplots
for(i in 1:length(diversity)) {
  boxplot = ggplot(data = diversity[[i]], aes(x = site, y = shannon_diversity, fill = site)) +
    geom_boxplot()
  print(boxplot)
}

# plot diversity as time plot
for(i in 1:length(diversity)) {
  diversity[[i]] = diversity[[i]] %>% 
    mutate(field_date = ymd(field_date))
  time = ggplot(data = diversity[[i]], aes(x = field_date, y = shannon_diversity, color = site)) +
    geom_point()
  print(time)
}

# Does diversity differ across rivers?
lapply(diversity, function(x) kruskal.test(shannon_diversity~site, data = x))
# not significantly different for any group but close for TM (p = 0.06)

lapply(diversity, function(x) kruskal.test(shannon_diversity~site, data = x))
# not significantly different for any group but close for TM (p = 0.06)

# how about diversity with anatoxins
diversity$tm <- left_join(diversity$tm, env_tog, by = c("field_date", "site_reach", "site"))
ggplot(data = diversity$tm, aes(x = TM_ATX_all_ug_orgmat_g.y, y = shannon_diversity, color = site)) +
  geom_point()
cor.test(diversity$tm$shannon_diversity[-7], diversity$tm$TM_ATX_all_ug_orgmat_g.y[-7])
# p = 0.2667
