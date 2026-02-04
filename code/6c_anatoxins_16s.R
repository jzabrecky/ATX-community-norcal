#### Comparing 16s molecular data with regard to anatoxin concentrations
### Jordan Zabrecky
## last edited: 02.02.2026

## insert description later

#### (1) Loading libraries & data ####

# libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot", "indicspecies"), require, character.only = T)

# read in files (data transformed in previous script, "4b_amongrivers_16s.R")
data <- lapply(list.files(path = "./data/molecular/transformed/", pattern = ".csv"),
               function(x) read.csv(paste("./data/molecular/transformed/", x, sep = "")))
names(data) <- c("nt", "tac", "tm")

# load environmental data
env_tog <- read.csv("./data/field_and_lab/env_tox_standardized_together.csv")
env_sep <- read.csv("./data/field_and_lab/env_tox_standardized_byriver.csv")

# lastly, merge with our community data to get rows formatted the same
data_tog <- lapply(data, function(x) left_join(x, env_tog, by = c("site", "site_reach", "field_date")))
data_sep <- lapply(data, function(x) left_join(x, env_sep, by = c("site", "site_reach", "field_date")))

# for microcoleus & anabaena remove rows where there are NAs for anatoxins
which(is.na(data_tog$tm$TM_ATX_all_ug_orgmat_g)) #7
which(is.na(data_tog$tac$TAC_ATX_all_ug_orgmat_g)) #none!

# 7 was a "not sure if this is microcoleus sample" and there was not enough material
# to analyze for anatoxins
data_tog$tm <- data_tog$tm[-7,]
data_sep$tm <- data_sep$tm[-7,]

data_sep_list <- list()
for(i in 1:length(data_sep)) {
  data_sep_list[[i]] <- split(data_sep[[i]], data_sep[[i]]$site)
}
names(data_sep_list) <- names(data_sep)

# set starting column for community data
start_col <- 6

#### (2) Functions for Analyses ####

# load community analyses functions from other script
source("./code/supplemental_code/S4a_community_analyses_func.R")

# set start column of community data
start_col <- 6

#### (3) Looking at Diversity Changes #### 

#### (4) PERMANOVA and BETADISP ####

## Lastly, what about with a binary?

# (a) just doing single within rivers
# (omitting Salmon because only one sample was toxic and very minorly)
set.seed(1)
runPERMANOVA(data_sep_list$tm$`SFE-M`, start_col, end_col = ncol(data$tm),
             group = data_sep_list$tm$`SFE-M`$TM_atx_detected)
anova(betadisper(vegdist(data_sep_list$tm$`SFE-M`[,start_col:ncol(data$tm)], method = "bray"), 
                 data_sep_list$tm$`SFE-M`$TM_atx_detected))
# SFE TM: significant: **, dispersion is not!
set.seed(1)
runPERMANOVA(data_sep_list$tac$`SFE-M`, start_col, end_col = ncol(data$tac),
             group = data_sep_list$tac$`SFE-M`$TAC_atx_detected)
anova(betadisper(vegdist(data_sep_list$tac$`SFE-M`[,start_col:ncol(data$tac)], method = "bray"), 
                 data_sep_list$tac$`SFE-M`$TAC_atx_detected))
# SFE TAC: significant *, dispersion is not!
set.seed(1)
runPERMANOVA(data_sep_list$tac$`RUS`, start_col, end_col = ncol(data$tac),
             group = data_sep_list$tac$`RUS`$TAC_atx_detected)
anova(betadisper(vegdist(data_sep_list$tac$`RUS`[,start_col:ncol(data$tac)], method = "bray"), 
                 data_sep_list$tac$`RUS`$TAC_atx_detected))
# RUS TAC: not significant, nor is dispersion

#### (5) dbRDA ####

## see variable selection in script "6b_anatoxins_microscopy.R"

## (?) apply to each river separately with separate standardization

# russian TAC
russian_TAC = run_dbRDA(data_sep_list$tac$`SFE-M`, 
                        start_col = start_col, 
                        end_col = ncol(data$tac),
                        mat_atx = "tac")

# sfkeel TAC
sfkeel_TAC = run_dbRDA(data_sep_list$tac$`SFE-M`,
                       start_col = start_col,
                       end_col = ncol(data$tac),
                       mat_atx = "tac")

# sfkeel TM
sfkeel_TM = run_dbRDA(data_sep_list$tm$`SFE-M`,
                      start_col = start_col,
                      end_col = ncol(data$tm),
                      mat_atx = "tm")

# add into list
individual_t_rdas <- list(russian_TAC, sfkeel_TAC, sfkeel_TM)
names(individual_t_rdas) <- c("russian_TAC", "sfkeel_TAC", "sfkeel_TM")

# plot
lapply(individual_t_rdas, function(x) plot(x$object))

# are models significant?
lapply(individual_t_rdas, function(x) print(x$model_sig))
# russian tac no, sfkeel tac no, sfkeel tm yes
# this seems to be a reverse of the situation for microscopy data!

# what is the rsquared?
lapply(individual_t_rdas, function(x) print(x$rsquared))
# adjusted: russian tac 0.17, sfkeel tac 0.17, sfkeel tm 0.22

# significant variables?
lapply(individual_t_rdas, function(x) print(x$sig_variables))
# ATX for all * plus oPhos for TM sfkeel

#### (5) NMDS Visualization ####

# South Fork Eel Microcoleus
makeNMDSplot(getNMDSdata(data_sep_list$tm$`SFE-M`, start_col, end_col = ncol(data$tm), 
                         ASV = TRUE),
             color = "TM_atx_detected", shape = "TM_atx_category", loading = FALSE)

# Russian Anabaena
makeNMDSplot(getNMDSdata(data_sep_list$tac$`RUS`, start_col, end_col = ncol(data$tac),
                         ASV = TRUE),
             color = "TAC_atx_detected", shape = "TAC_atx_category", loading = FALSE)

# South Fork Eel Anabaena
makeNMDSplot(getNMDSdata(data_sep_list$tac$`SFE-M`, start_col, end_col = ncol(data$tac),
                         ASV = TRUE),
             color = "TAC_atx_detected", shape = "TAC_atx_category", loading = FALSE)
# stress is weirdly 0- maybe because of an outlier sample?
# one none detected is away from the others, let's try
which(data_sep_list$tac$`SFE-M`$TAC_atx_detected == "none") # 6 or 8

# it is apparently sample 6 which is 2022-09-06 SFE-M-1S
makeNMDSplot(getNMDSdata(data_sep_list$tac$`SFE-M`[-8,], start_col, end_col = ncol(data$tac),
                         ASV = TRUE),
             color = "TAC_atx_detected", shape = "TAC_atx_category", loading = FALSE)

#### (6) 

