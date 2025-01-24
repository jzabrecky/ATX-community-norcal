#### Processing miscellaneous Non-target data
### Jordan Zabrecky
## last edited 01.25.2024

# This code analyzes the triplicates and blanks of non-target data
# to create a more processed final csv that may or may not be used

#### (1) Loading libraries and data ####

# loading libraries
lapply(c("tidyverse", "lubridate"), require, character.only = T)

# read in data
nt_misc <- read.csv("./data/EDI_data_package/non_target_sample_miscellaneous.csv")

# convert field_date to date object
nt_misc$field_date <- ymd(nt_misc$field_date)

# separate out into single sample info (area scraped, etc.), afdm, and chlorophyll-a data
samples <- nt_misc %>% 
  select(site_reach, field_date, area_scraped_cm_2, volume_scraped_mL, macroalgal_displacement_mL) %>% 
  na.omit()
afdm <- nt_misc %>% 
  select(site_reach, field_date, triplicate, dry_sample_filter_tin_g, ash_filter_tin_g, afdm_g,
         afdm_sample_filtered_mL)
chla <- nt_misc %>% 
  select(site_reach, field_date, triplicate, Chla_ug_L, Pheo_ug_L, Chla_Pheo_flag, neg_Pheo_flag,
         chla_sample_filtered_mL)

#### (2) Processing AFDM data ####

## (a) process blanks

# gather blanks of samples and do some editing to the dataframe
afdm_blanks <- afdm %>% 
  filter(site_reach == "blank") %>% 
  dplyr::rename(lab_date = field_date,
                blank_afdm_g = afdm_g) %>% 
  select(lab_date, blank_afdm_g) %>% 
  mutate(field_date = lab_date - 1)

# creating data frame for more processing of afdm data
afdm_processing <- afdm %>% 
  filter(site_reach != "blank")

# left_join processing and blanks data
afdm_processing <- left_join(afdm_processing, afdm_blanks, by = "field_date")

# we didn't do a blank for sampling that occurred on 6/27/2024 so 
# let's just assume it's similar to the blank from the next sampling on the Salmon
afdm_processing$blank_afdm_g[2:3] <- afdm_processing$blank_afdm_g[12]

# also sample from 7/7/2024 was processed the same day (7/7/2024) as those on 7/6/2024
afdm_processing$blank_afdm_g[11] <- afdm_processing$blank_afdm_g[10]

# subtract blank from afdm measurement
afdm_processing$ash_free_dry_mass_g <- afdm_processing$afdm_g - 
  afdm_processing$blank_afdm_g

# calculate afdm per milliters of sample
afdm_processing$afdm_g_mL <- afdm_processing$ash_free_dry_mass_g / afdm_processing$afdm_sample_filtered_mL

## (b) process triplicates

# calculate summary statistics for triplicates
triplicates_afdm <- afdm_processing %>% 
  filter(triplicate == "y") %>% 
  dplyr::group_by(field_date, site_reach) %>% 
  dplyr::summarize(mean = mean(afdm_g_mL),
                   sd = sd(afdm_g_mL),
                   rsd = sd * 100 / mean) %>% 
  ungroup() %>% 
  distinct()

# look at triplicate results
view(triplicates_afdm) # RSD max is 38.3% (the first one); we got better (<16% afterwards)
mean(triplicates_afdm$rsd) # mean is 11.3%

# take dataset and select for columns we care about before merging
triplicates_afdm_final <- triplicates_afdm %>% 
  dplyr::rename(afdm_g_mL = mean) %>% 
  mutate(afdm_g_mL = round(afdm_g_mL, 7)) %>% # round to nearest 2 decimals
  # (briefly looking at the data that is roughly our sig figs...)
  select(field_date, site_reach, afdm_g_mL)

# remove triplicates from original dataset before joining the two
afdm_final <- afdm_processing %>% 
  filter(triplicate == "n") %>% 
  mutate(afdm_g_mL = round(afdm_g_mL, 7)) %>% # also round to the nearest 2 decimals
  select(field_date, site_reach, afdm_g_mL)

# add processed/averaged triplicates back in
afdm_final <- rbind(afdm_final, triplicates_afdm_final)

#### (3) Processing Chlorophyll-a data ####

# convert chlorophyll-a & pheophytin in ug/L to ug/mL
chla$Chla_ug_mL <- chla$Chla_ug_L / 1000
chla$Pheo_ug_mL <- chla$Pheo_ug_L / 1000

## (a) process blanks

# check blanks to see if they had detected chlorophyll-a
chla$Chla_Pheo_flag[which(chla$site_reach == "blank")]

# two blanks are above detection limit?
chla_blanks <- chla[(which(chla$site_reach == "blank")),] %>% 
  # rename field_date as lab_date (which is what it truly is) 
  dplyr::rename(lab_date = field_date) %>% 
  # if measured Chla & Pheo below detection, replace with 0
  mutate(Chla_blank_ug_mL = case_when(Chla_Pheo_flag == "y" ~ 0,
                                      TRUE ~ Chla_ug_mL),
         Pheo_blank_ug_mL = case_when(Chla_Pheo_flag == "y" ~ 0,
                                      TRUE ~ Pheo_ug_mL),
         field_date = ymd(lab_date) - 1) %>% 
  select(field_date, Chla_blank_ug_mL, Pheo_blank_ug_mL)

# left join blanks
chla <- left_join(chla %>% filter(site_reach != "blank"), chla_blanks, by = "field_date")

# did not run blanks after 6/27/2022 Salmon samples but no other blanks after salmon samples
# were detectable; also 7/7/2022 Russian sample was ran on same day with 7/6/2022 samples
# which did not have a detectable blank-- fill NA of column with 0
chla$Chla_blank_ug_mL <- replace(chla$Chla_blank_ug_mL, which(is.na(chla$Chla_blank_ug_mL)), 0)
chla$Pheo_blank_ug_mL <- replace(chla$Pheo_blank_ug_mL, which(is.na(chla$Pheo_blank_ug_mL)), 0)

# calculate final chlorophyll-a & pheophytin minus detectable blanks
chla$Chla_ug_mL_final <- chla$Chla_ug_mL - chla$Chla_blank_ug_mL
chla$Pheo_ug_mL_final <- chla$Pheo_ug_mL - chla$Pheo_blank_ug_mL

## (b) process triplicates

# filter out triplicates and calculate summary statistics
triplicates_chla <- chla %>% 
  filter(triplicate == "y") %>% 
  dplyr::group_by(field_date, site_reach) %>% 
  dplyr::summarize(mean_chla = mean(Chla_ug_mL_final),
                   sd_chla = sd(Chla_ug_mL_final),
                   rsd_chla = sd_chla * 100 / mean_chla,
                   mean_pheo = mean(Pheo_ug_mL_final)) %>% 
  ungroup() %>% 
  distinct()

# look at triplicate results
view(triplicates_chla) # RSD max is 38.8% which is pretty high, the next is 34.7%
mean(triplicates_chla$rsd_chla) # average is 19.6%

# take dataset and select for columns we care about before merging
triplicates_chla_final <- triplicates_chla %>% 
  dplyr::rename(Chla_ug_mL_final = mean_chla,
                Pheo_ug_mL_final = mean_pheo) %>% 
  select(field_date, site_reach, Chla_ug_mL_final, Pheo_ug_mL_final)

# remove triplicates from original dataset before joining the two
chla_final <- chla %>% 
  filter(triplicate == "n") %>% 
  select(field_date, site_reach, Chla_ug_mL_final, Pheo_ug_mL_final)

# add processed/averaged triplicates back in & do some conversions
chla_final <- rbind(chla_final, triplicates_chla_final) %>% 
  select(field_date, site_reach, Chla_ug_mL_final, Pheo_ug_mL_final) %>%
  dplyr::rename(Chla_ug_mL = Chla_ug_mL_final, # renaming for final df
                Pheo_ug_mL = Pheo_ug_mL_final) %>% 
  # add back in a flag for data with negative Pheophytin (means chl-a is overestimated)
  mutate(flag_Pheo = case_when(Pheo_ug_mL <  0 ~ 1,
                               TRUE ~ 0))

#### (4) Converting Chl-a and afdm per mL to per area scraped ####

# calculate ratio of mL of sample to area sampled
samples$vol_mL_to_area_cm_2 <- samples$volume_scraped_mL / samples$area_scraped_cm_2

# merge afdm and chla data into area_volume
all <- left_join(samples, afdm_final, by = c("field_date", "site_reach"))
all <- left_join(all, chla_final, by = c("field_date", "site_reach")) 

# calculate afdm, chla, & pheophytin per area
all$afdm_g_cm_2 <- all$afdm_g_mL * all$vol_mL_to_area_cm_2
all$Chla_ug_cm_2 <- all$Chla_ug_mL * all$vol_mL_to_area_cm_2
all$Pheo_ug_cm_2 <- all$Pheo_ug_mL * all$vol_mL_to_area_cm_2

# final csv only keeping columns we care about
all_final <- all %>%
  mutate(afdm_mg_cm_2 = afdm_g_cm_2 * 1000) %>%  # convert g/cm^2 to mg/cm^2
  select(field_date, site_reach, macroalgal_displacement_mL, area_scraped_cm_2, volume_scraped_mL,
         vol_mL_to_area_cm_2, afdm_mg_cm_2, afdm_g_cm_2, Chla_ug_cm_2, Pheo_ug_cm_2, flag_Pheo)

# save final dataset
write.csv(all_final, "./data/NT_misc/NT_misc_data_processed.csv", row.names = FALSE)
