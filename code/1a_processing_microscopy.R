#### Processing microscopy data and adding in environmental & atx data
### Jordan Zabrecky
## last edited 03.06.2025

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

# read in anatoxin & water chemistry data
atx_target <- read.csv("./data/field_and_lab/cyano_atx.csv")
water_chemistry <- read.csv("./data/field_and_lab/water_chemistry.csv") %>% 
  select(!time)

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

# put all processed data into a list
processed <- list(tm_processed, tac_processed, nt_processed)

# left join in water data and change column organization
for(i in 1:length(processed)) {
  processed[[i]] <- left_join(processed[[i]], water_chemistry, by = c("site_reach", "field_date")) %>%  
    dplyr::select(!reach) # don't need a column that just has reach
}

# left join in water chemistry data
tm_processed <- left_join(tm_processed, water_chemistry, by = c("site_reach", "field_date")) %>% 
  dplyr::relocate(site, .after = site_reach)
tac_processed <- left_join(tac_processed, water_chemistry, by = c("site_reach", "field_date"))
nt_processed <- left_join(nt_processed, water_chemistry, by = c("site_reach", "field_date"))
# also want to add in TM and TAC ATX to NT dataframe