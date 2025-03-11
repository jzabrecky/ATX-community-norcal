#### Processing microscopy data and adding in environmental & atx data
### Jordan Zabrecky
## last edited 03.10.2025

## This code processes microscopy data from the EDI data release 
## (averaging across slides and evaluting rsd) and adds in water 
## chemistry data and anatoxin concentrations

## For code processing of water chemistry and target sample anatoxin 
## data from original EDI csv's, see
## https://github.com/jzabrecky/ATX-synchrony-norcal/blob/main/code/2b_water_chemistry.R
## https://github.com/jzabrecky/ATX-synchrony-norcal/blob/main/code/2e_target_sample_anatoxins.R

## need to make another script to process NT anatoxins :)
## maybe do it in S1a?

#### (1) Loading in libraries & data ####

# loading libraries
lapply(c("tidyverse", "lubridate"), require, character.only = T)

# read in microscopy data
target <- read.csv("./data/EDI_data_package/microscopy_target_samples.csv")
nontarget <- read.csv("./data/EDI_data_package/microscopy_non_target_samples.csv")

# read in processed anatoxin & water chemistry data
atx_target <- read.csv("./data/field_and_lab/cyano_atx.csv") # processed version
water_chemistry <- read.csv("./data/field_and_lab/water_chemistry.csv") %>% 
  select(!c(time, reach))

#### (2) Analyzing microscopy data (w/ multiple slides) ####

## (a) target samples

# look at average, sd, & rsd for each taxa across slides
target_multislides <- target %>% 
  filter(slide_rep != "Final") %>% 
  pivot_longer(microcoleus:aphanothece, names_to = "taxa", values_to = "percent") %>% 
  dplyr::group_by(taxa, site_reach, field_date, sample_type) %>% 
  dplyr::summarize(mean = mean(percent),
                   sd = sd(percent),
                   rsd = (sd * 100) / mean)
view(target_multislides) # the highest are all samples where one slide had 0.1% cover
mean(target_multislides$rsd, na.rm = TRUE) # 44.3% high but let's look at taxa 
                                           # that have >5% mean cover
plot(target_multislides$mean, target_multislides$rsd)
# yup, samples with a higher mean have a lower rsd

# looking at taxa that have >5% mean cover
target_multislides_sub <- target_multislides %>% 
  filter(mean > 5)
view(target_multislides_sub)
mean(target_multislides_sub$rsd) # 21.7% high but better

# pivot back wider to join back into processed
target_multislides_wider <- target_multislides %>% 
  select(!(sd:rsd)) %>% 
  pivot_wider(names_from = "taxa", values_from = "mean")

# for data that were already averaged across three slides,
# we need to match it up to the processed data frame
target_processed <- target %>% 
  filter(slide_rep == "Final") %>% 
  select(!(slide_rep:method)) # remove columns we don't need

# match up column order
target_processed <- target_processed[,colnames(target_multislides_wider)]

# rbind to join dataframes
target_processed <- rbind(target_processed, target_multislides_wider)

## (b) non-target samples

# look at average, sd, & rsd for each taxa across slides
# no samples here are averaged over as indicated by "Final"
nontarget_multislides <- nontarget %>% 
  pivot_longer(non_algal:unknown, names_to = "taxa", values_to = "percent") %>% 
  dplyr::group_by(taxa, site_reach, field_date, sample_type) %>% 
  dplyr::summarize(mean = mean(percent),
                   sd = sd(percent),
                   rsd = (sd * 100) / mean)
view(nontarget_multislides) # the highest are all samples where one slide had 0.1% cover
mean(nontarget_multislides$rsd, na.rm = TRUE) # 97.6% high but let's look at taxa 
                                              # that have >5% mean cover
plot(nontarget_multislides$mean, nontarget_multislides$rsd)
# less strong of a relationship with non-target samples

# looking at taxa that have >5% mean cover
nontarget_multislides_sub <- nontarget_multislides %>% 
  filter(mean > 5)
view(nontarget_multislides_sub)
mean(nontarget_multislides_sub$rsd) # 42.8% high but better

# pivot back wider to join with anatoxin and environmental data
nt_processed <- nontarget_multislides %>% 
  select(!(sd:rsd)) %>% 
  pivot_wider(names_from = "taxa", values_from = "mean") 

#### (3) Joining in Anatoxin & Water Chemistry Data ####

# separate out TAC and TM 
tm_processed <- target_processed %>% 
 filter(sample_type == "TM")
tac_processed <- target_processed %>% 
  filter(sample_type == "TAC")

# replace 2022-09-06 microcoleus sample with one from 2022-09-08
# (substitute technician accidentally took 2022-09-06 with different methods)
tm_processed$field_date[which(tm_processed$field_date == "2022-09-08")] <- "2022-09-06"
atx_target$field_date[which(atx_target$field_date == "2022-09-08")] <- "2022-09-06"

# put all processed data into a list
processed <- list(tm_processed, tac_processed, nt_processed)
names(processed) <- c("tm", "tac", "nt")

# left join in water data and change column organization
for(i in 1:length(processed)) {
  processed[[i]] <- left_join(processed[[i]], water_chemistry, 
                              by = c("site_reach", "field_date")) %>%  
    relocate(site, .after = site_reach)
}

# left join in anatoxin data (TM gets TM, TAC gets TAC, and NT gets both)
processed$tm <- left_join(processed$tm, atx_target %>% filter(sample_type == "TM"), 
                          by = c("site_reach", "field_date", "sample_type"))
processed$tac <- left_join(processed$tac, atx_target %>% filter(sample_type == "TAC"), 
                           by = c("site_reach", "field_date", "sample_type"))
processed$nt <- left_join(processed$nt, atx_target %>% filter (sample_type == "TM") 
                          %>% select(!sample_type),
                          by = c("site_reach", "field_date")) %>% 
  rename(TM_ATX_all_ug_chla_ug = ATX_all_ug_chla_ug,
         TM_ATX_all_ug_orgmat_g = ATX_all_ug_orgmat_g) %>% 
  select(!(ATXa_ug_g:percent_organic_matter))
processed$nt <- left_join(processed$nt, atx_target %>% filter (sample_type == "TAC") 
                       %>% select(!sample_type),
                       by = c("site_reach", "field_date")) %>% 
  rename(TAC_ATX_all_ug_chla_ug = ATX_all_ug_chla_ug,
         TAC_ATX_all_ug_orgmat_g = ATX_all_ug_orgmat_g) %>% 
  select(!(ATXa_ug_g:percent_organic_matter))

# turns out we don't have ATX data for TM SFE-M-4 08-23-2022
# I think we weren't clear if it would be Microcoleus and did not have enough
# so let's just remove that row
processed$tm <- processed$tm[-22,]

# saving csv's
path <- paste(getwd(), "/data/morphological/", sep = "")
lapply(names(processed), function(x) write.csv(processed[[x]], file = paste(path, x, "_microscopy_with_covar.csv", 
                                                                              sep = ""), row.names = FALSE))
