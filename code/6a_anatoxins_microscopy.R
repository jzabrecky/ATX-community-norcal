#### Comparing microscopy data with regard to anatoxin concentrations
### Jordan Zabrecky
## last edited: 12.20.2025

# This script examines how communities as identified by microscopy
# change with increasing anatoxin concentrations by <INSERT>

#### (1) Loading libraries & data ####

# set seed for reproducibility
set.seed(2025)

# libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot"), require, character.only = T)

# read in relative abundance files (data transformed in previous script, "4a_amongrivers_microscopy.R")
data <- lapply(list.files(path = "./data/morphological/transformed/", pattern = ".csv"),
               function(x) read.csv(paste("./data/morphological/transformed/", x, sep = "")))
names(data) <- c("nt", "tac", "tm")

# read in environmental covariates & toxin data
env <- read.csv("data/field_and_lab/environmental_covariates_and_toxins.csv") %>% 
  # analyze congeners all together, so removed separted congeners
  select(!c("pH", "TM_ATXa_ug_g", "TAC_ATXa_ug_g", "TM_dhATXa_ug_g", 
            "TAC_dhATXa_ug_g", "TM_HTXa_ug_g", "TAC_HTXa_ug_g", "TM_Chla_ug_g",
            "TAC_Chla_ug_g", "TM_Pheo_ug_g", "TAC_Pheo_ug_g", "TM_percent_organic_matter",
            "TAC_percent_organic_matter", "TM_ATX_all_ug_chla_ug", "TAC_ATX_all_ug_chla_ug",
            "TM_ATX_all_ug_g", "TAC_ATX_all_ug_g")) %>% 
  filter(year(ymd(field_date)) == 2022)

#### (2) Predictor Data Transformation ####

# replace 0's with small value for log-transformation (to avoid outweighing outliers)
# DO NOT replace NAs because in those circumstances, no mat was present for us to sample
# first look at histogram
hist(env$TM_ATX_all_ug_orgmat_g, breaks = 30)
hist(env$TAC_ATX_all_ug_orgmat_g, breaks = 30)

# lowest non-zero value is 0.03, how about 0.005 as zero replacement?
log(0.03) # -3.506558
log(0.01) # -4.60517

# replace zeros
env$TM_ATX_all_ug_orgmat_g <- replace(env$TM_ATX_all_ug_orgmat_g , 
                                      env$TM_ATX_all_ug_orgmat_g == 0, 0.01)
env$TAC_ATX_all_ug_orgmat_g <- replace(env$TAC_ATX_all_ug_orgmat_g , 
                                       env$TAC_ATX_all_ug_orgmat_g == 0, 0.01)

# log-transform variables 
env[,4:ncol(env)] <- log(env[,4:ncol(env)])
# center & scale
# data all-together
env_tog <- env
env_tog[,4:ncol(env_tog)] <- apply(env_tog[,4:ncol(env_tog)], 2, scale)
# rivers separately (group-level variable)
env_sep <- env
env_sep[,4:ncol(env_sep)] <- apply(env_sep[,4:ncol(env_sep)], 2,
                                   function(x) ave(x, env_sep$site, FUN = scale))

# compare histograms of the two
hist(env_tog$TM_ATX_all_ug_orgmat_g)
hist(env_sep$TM_ATX_all_ug_orgmat_g)
hist(env_tog$TAC_ATX_all_ug_orgmat_g)
hist(env_sep$TAC_ATX_all_ug_orgmat_g)

# save to avoid doing this in future scripts (run once)
#write.csv(env_tog, "./data/field_and_lab/env_tox_standardized_together.csv", 
#          row.names = FALSE)
#write.csv(env_sep, "./data/field_and_lab/env_tox_standardized_byriver.csv", 
#          row.names = FALSE)

# lastly, merge with our community data to get rows formatted the same
data_tog <- lapply(data, function(x) left_join(x, env_tog, by = c("site", "site_reach", "field_date")))
data_sep <- lapply(data, function(x) left_join(x, env_sep, by = c("site", "site_reach", "field_date")))

#### (3) Redundancy Analysis ####

# for microcoleus & anabaena remove rows where there are NAs for anatoxins
which(is.na(data_tog$tm$TM_ATX_all_ug_orgmat_g))
which(is.na(data_tog$tac$TAC_ATX_all_ug_orgmat_g))
# we need to remove one for microcoleus!

# also remove TM for 

# what is our community here is it transformed I can't even remember

test.rda <- rda(data_tog$tac[,5:ncol(data$tac)] ~ 
                  data_tog$tac[,c("TAC_ATX_all_ug_orgmat_g", "nitrate_mg_N_L")])
# bruh idk
summary(test.rda)

model <- ordiplot(test.rda, type = "none", scaling = 2, cex=10, xlab = "RDA1 (26.2%)", ylab = "RDA2 (9.8%)", cex.lab=1.25)
points(test.rda, col="darkgrey", cex=1)
points(test.rda, dis="sp", col="blue")
text(test.rda, dis="sp", col="blue")
text(test.rda, dis="bp", col="black")

# maybe RDA is more about just exploring potential explanatory varaibles- 
# would want to confrim that we see differences with toxins though

# how do I make this look nice:
# https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html
# basically, Anabaena samples that had higher toxins had a higher amount of Microcoleus???
# consider bioindicators of toxin blooms as done in Jansen et al. 2025 using the ALDEx2 package?

# may want to transform predictor variables
