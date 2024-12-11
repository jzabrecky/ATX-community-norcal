#### Processing miscellaneous Non-target data
### Jordan Zabrecky
## last edited 12.10.2024

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
  select(site_reach, field_date, triplicate, dry_sample_filter_tin_g, ash_filter_tin_g, afdm_filter_g,
         afdm_sample_filtered_mL)
chla <- nt_misc %>% 
  select(site_reach, field_date, triplicate, Chla_ug_L, Pheo_ug_L, Chla_Pheo_flag, neg_Pheo_flag,
         chla_sample_filtered_mL)

#### (2) Processing AFDM data ####

# (a) process blanks

# gather blanks of samples and do some editing to the dataframe
afdm_blanks <- afdm %>% 
  filter(site_reach == "blank") %>% 
  dplyr::rename(lab_date = field_date,
                blank_afdm_g = afdm_filter_g) %>% 
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
afdm_processing$ash_free_dry_mass_g <- afdm_processing$afdm_filter_g - 
  afdm_processing$blank_afdm_g

# calculate afdm per milliters of sample
afdm_processing$afdm_g_mL <- afdm_processing$ash_free_dry_mass_g / afdm_processing$afdm_sample_filtered_mL

# (b) process triplicates

# calculate summary statistics for triplicates
triplicates <- afdm_processing %>% 
  filter(triplicate == "y") %>% 
  dplyr::group_by(field_date, site_reach) %>% 
  dplyr::summarize(mean = mean(afdm_g_mL),
                   sd = sd(afdm_g_mL),
                   rsd = sd * 100 / mean) %>% 
  ungroup() %>% 
  distinct()

# look at triplicate results
view(triplicates) # RSD max is 38.3% (the first one); we got better (<16% afterwards)
mean(triplicates$rsd) # mean is 11.3%

# take dataset and select for columns we care about before merging
triplicates_final <- triplicates %>% 
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
afdm_final <- rbind(afdm_final, triplicates_final)

#### (3) Processing Chlorophyll-a data ####

# convert chlorophyll-a & pheophytin in ug/L to ug/mL
chla$Chla_ug_mL <- chla$Chla_ug_L / 1000
chla$Pheo_ug_mL <- chla$Pheo_ug_L / 1000

# (a) process blanks

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
chla_blanks$Chla_blank_ug_mL <- replace(chla_blanks$Chla_blank_ug_mL, which(is.na(chla_blanks$Chla_blank_ug_mL)), 0)
chla_blanks$Pheo_blank_ug_mL <- replace(chla_blanks$Pheo_blank_ug_mL, which(is.na(chla_blanks$Pheo_blank_ug_mL)), 0)

# calculate final chlorophyll-a & pheophytin minus detectable blanks
chla_pheo_calculations$Chla_ug_mL_final <- chla_pheo_calculations$Chla_ug_mL - chla_pheo_calculations$Chla_blank_ug_mL
chla_pheo_calculations$Pheo_ug_mL_final <- chla_pheo_calculations$Pheo_ug_mL - chla_pheo_calculations$Pheo_blank_ug_mL

### NEED TO FIGURE OUT NEGATIVE PHEO FLAG WHAT TO DO

#### PROCESSING CHL-A CAN MOVE TO ANOTHER SCRIPT

# convert chlorophyll-a & pheophytin in ug/L to ug/mL
chla_pheo_calculations$Chla_ug_mL <- chla_pheo_calculations$Chla_ug_L / 1000
chla_pheo_calculations$Pheo_ug_mL <- chla_pheo_calculations$Pheo_ug_L / 1000

# check blanks to see if they had detected chlorophyll-a
chla_pheo_calculations$RFU_bdl[which(chla_pheo_calculations$site_reach == "blank")]

# two blanks are above detection limit?
pos_blanks <- chla_pheo_calculations[(which(chla_pheo_calculations$RFU_bdl == 0 &
                                              chla_pheo_calculations$site_reach == "blank")),] %>% 
  dplyr::rename(lab_date = field_date,
                Chla_blank_ug_mL = Chla_ug_mL,
                Pheo_blank_ug_mL = Pheo_ug_mL) %>% 
  mutate(field_date = mdy(lab_date) - 1) %>% # for blanks field date is actually the lab date
  select(field_date, Chla_blank_ug_mL, Pheo_blank_ug_mL)

# convert field date to Date object to join
chla_pheo_calculations$field_date <- mdy(chla_pheo_calculations$field_date)

# left join blanks
chla_pheo_calculations <- left_join(chla_pheo_calculations, pos_blanks, by = "field_date")

# fill non-detectable blanks with zero
chla_pheo_calculations$Chla_blank_ug_mL <- replace(chla_pheo_calculations$Chla_blank_ug_mL,
                                                   is.na(chla_pheo_calculations$Chla_blank_ug_mL), 0)
chla_pheo_calculations$Pheo_blank_ug_mL <- replace(chla_pheo_calculations$Pheo_blank_ug_mL,
                                                   is.na(chla_pheo_calculations$Pheo_blank_ug_mL), 0)

# calculate final chlorophyll-a & pheophytin minus detectable blanks
chla_pheo_calculations$Chla_ug_mL_final <- chla_pheo_calculations$Chla_ug_mL - chla_pheo_calculations$Chla_blank_ug_mL
chla_pheo_calculations$Pheo_ug_mL_final <- chla_pheo_calculations$Pheo_ug_mL - chla_pheo_calculations$Pheo_blank_ug_mL

# now we can remove the blanks!
chla_pheo_calculations <- chla_pheo_calculations[-which(chla_pheo_calculations$site_reach == "blank"),]

# cannot rerun samples with negative pheophytin values (sample gone), 
# but will add a flag at the end to indicate unreliable data

## (b) analyze triplicates

# filter out triplicates and calculate summary statistics
triplicates <- chla_pheo_calculations %>% 
  filter(triplicate == "y") %>% 
  dplyr::group_by(field_date, site_reach) %>% 
  dplyr::summarize(mean_chla = mean(Chla_ug_mL_final),
                   sd_chla = sd(Chla_ug_mL_final),
                   rsd_chla = sd_chla * 100 / mean_chla,
                   mean_pheo = mean(Pheo_ug_mL_final)) %>% 
  ungroup() %>% 
  distinct()

# look at triplicate results
view(triplicates) # RSD max is 38.8% which is pretty high, the next is 34.7%
mean(triplicates$rsd_chla) # average is 19.6%

## (c) put together final dataset and save

# take dataset and select for columns we care about before merging
triplicates_final <- triplicates %>% 
  dplyr::rename(Chla_ug_mL_final = mean_chla,
                Pheo_ug_mL_final = mean_pheo) %>% 
  select(field_date, site_reach, Chla_ug_mL_final, Pheo_ug_mL_final)

# remove triplicates from original dataset before joining the two
final <- chla_pheo_calculations %>% 
  filter(triplicate == "n") %>% 
  select(field_date, site_reach, Chla_ug_mL_final, Pheo_ug_mL_final)

# add processed/averaged triplicates back in & do some conversions
final <- rbind(final, triplicates_final) %>% 
  select(field_date, site_reach, Chla_ug_mL_final, Pheo_ug_mL_final) %>%
  dplyr::rename(Chla_ug_mL = Chla_ug_mL_final, # renaming for final df
                Pheo_ug_mL = Pheo_ug_mL_final) %>% 
  mutate(flag_Pheo = case_when(Pheo_ug_mL <  0 ~ 1,
                               TRUE ~ 0))

# save csv
write.csv(final, "../rock_scrape_chla.csv", row.names = FALSE)