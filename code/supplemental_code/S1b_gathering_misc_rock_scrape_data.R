#### Code to process and organize miscellaneous rock scrape/NT for EDI data package
### Jordan Zabrecky
## last edited 10.07.2024

# This code gathers miscellaneous rock scrape data that is not currently used for 
# either ATX projects. Data includes total area scraped, total volume of sample, 
# displacement of preserved sample, chlorophyll-a per area scraped, 
# AFDM per area scraped, and, if present, detected anatoxins of sample

#### (1) Loading libraries & data ####

# loading libraries
lapply(c("tidyverse", "lubridate", "plyr"), require, character.only = T)

# area scraped and sample volume data
area_volume <- read.csv("./data/NT_rock_scrapes_misc/raw/volumes_and_area.csv")

# filtered chlorophyll-a data (processed in script "1a_rock_scrape_filter_chla.R")
chla <- read.csv("./data/NT_rock_scrapes_misc/rock_scrape_chla.csv")

# filtered ash-free dry mass data
afdm <- ldply(list.files(path = "./data/NT_rock_scrapes_misc/raw/", 
                         pattern = "angelo"), function(filename) {
                           d <- read.csv(paste("./data/NT_rock_scrapes_misc/raw/", filename, sep = ""))
                           d$field_date <- mdy(d$field_date)
                           return(d)
                         })


# rock scrape anatoxin data
anatoxins_full <- read.csv("./data/NT_rock_scrapes_misc/raw/SUNY_ESF_LCMS_MS.csv")
metadata <- read.csv("./data/NT_rock_scrapes_misc/raw/SUNY_ESF_metadata.csv")

# changing date fields to date objects
area_volume$field_date <- mdy(area_volume$field_date)
afdm$field_date <- ymd(afdm$field_date)
metadata$field_date <- mdy(metadata$field_date)
chla$field_date <- ymd(chla$field_date)

#### (2) Processing ash-free dry mass data ####

# calculate afdm for all samples
afdm$ash_free_dry_mass_g_raw <- afdm$dry_sample_filter_tin_g - afdm$ash_filter_tin_g

# gather blanks of samples and do some editing to the dataframe
afdm_blanks <- afdm %>% 
  filter(site_reach == "blank") %>% 
  dplyr::rename(lab_date = field_date,
                blank_afdm_g = ash_free_dry_mass_g_raw) %>% 
  select(lab_date, blank_afdm_g) %>% 
  mutate(field_date = lab_date - 1)

# creating data frame for more processing of afdm data
afdm_processing <- afdm %>% 
  filter(site_reach != "blank")

# change class of mL filtered from character to integer
afdm_processing$filtered_amount_mL <- as.integer(afdm_processing$filtered_amount_mL)

# left_join processing and blanks data
afdm_processing <- left_join(afdm_processing, afdm_blanks, by = "field_date")

# we didn't do a blank for sampling that occurred on 6/27/2024 so 
# let's just assume it's similar to the blank from the next sampling on the Salmon
afdm_processing$blank_afdm_g[2:3] <- afdm_processing$blank_afdm_g[12]

# also sample from 7/7/2024 was processed the same day (7/7/2024) as those on 7/6/2024
afdm_processing$blank_afdm_g[11] <- afdm_processing$blank_afdm_g[10]

# subtract blank from afdm measurement
afdm_processing$ash_free_dry_mass_g <- afdm_processing$ash_free_dry_mass_g_raw - 
  afdm_processing$blank_afdm_g

# calculate afdm per milliters of sample
afdm_processing$afdm_g_mL <- afdm_processing$ash_free_dry_mass_g / afdm_processing$filtered_amount_mL

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

#### (3) Processing anatoxins data ####

## (a) cleaning ESF csv
# lots of cleaning to be done with SUNY ESF data sheet...
view(anatoxins_full)

# columns 2, 26, 27, 31, 32, 33, 36, 38, 40, 42 are important
anatoxins_reduced <- anatoxins_full[c(2, 26, 27, 34, 37, 39, 41, 43)]
# (I already know all of our samples had none of the 
# three cylindrospermopsin derivatives tested, so not including them here)

# rename these columns to be more computer-legible
labels <- c("ESF_ID", "MCY_ug_g", "MCY_det_limit_full", "ATX_det_limit_full", 
            "ATXa_ug_g", "HTXa_ug_g", "dhATXa_ug_g", "dhHTXa_ug_g")

# adding labels as column names
colnames(anatoxins_reduced) <- labels

# convert columns that are character class to numeric 
# (this will replace the "-" with NA)
anatoxins_reduced[5:8] <- sapply(anatoxins_reduced[5:8], as.numeric)

# parse detection limit column and convert to numeric
anatoxins_reduced <- anatoxins_reduced %>%
  mutate(MCY_det_limit = as.numeric(str_sub(anatoxins_reduced$MCY_det_limit_full,3, 6)),
         ATX_det_limit = as.numeric(str_sub(anatoxins_reduced$ATX_det_limit_full,3, 8))) %>% 
  select(-MCY_det_limit_full, -ATX_det_limit_full)

## (b) processing data values
# values below detection limit were included in our report, but we will not be
# reporting those values, so let's replace them with 0
anatoxins_processed <- anatoxins_reduced %>% 
  mutate(
    MCY_ug_g = case_when(MCY_ug_g < MCY_det_limit ~ 0,
                         TRUE ~ MCY_ug_g),
    ATXa_ug_g = case_when(ATXa_ug_g < ATX_det_limit ~ 0,
                          TRUE ~ ATXa_ug_g),
    HTXa_ug_g = case_when(HTXa_ug_g < ATX_det_limit ~ 0,
                          TRUE ~ HTXa_ug_g),
    dhATXa_ug_g = case_when(dhATXa_ug_g < ATX_det_limit ~ 0,
                            TRUE ~ dhATXa_ug_g),
    dhHTXa_ug_g = case_when(dhHTXa_ug_g < ATX_det_limit ~ 0,
                            TRUE ~ dhHTXa_ug_g)
  )

# fill NA values with 0 
anatoxins_processed <- replace(anatoxins_processed, is.na(anatoxins_processed), 0)

# calculate total anatoxins
anatoxins_processed <- anatoxins_processed %>% 
  mutate(ATX_all_ug_g = ATXa_ug_g + HTXa_ug_g + dhATXa_ug_g + dhHTXa_ug_g)

# match ESF_ID with metadata
combined <- left_join(metadata, anatoxins_processed, by = "ESF_ID") %>% 
  # only calculating NT/rock scrape samples ATX here (or blanks)
  filter(sample_type == "NT" | sample_type == "BLANK")

# check blanks to make sure they had no anatoxin detections
combined$ATX_all_ug_g[which(combined$sample_type == "BLANK")]

# blanks are zero- can remove them!
combined <- combined[-which(combined$sample_type == "BLANK"),]

# no triplicates to process here- only 14 samples had detectable ATX
# all from 2022 and all but one from South Fork Eel (one was from Russian)
# no MCY, HTX-a, or dhHTX-a detected
anatoxins_final <- combined %>% 
  select(site_reach, field_date, MCY_ug_g, MCY_det_limit, ATXa_ug_g, dhATXa_ug_g, ATX_det_limit)

#### (4) Converting chl-a and afdm per mL to per area ####

# calculate ratio of mL of sample to area sampled
area_volume$vol_mL_to_area_cm_2 <- area_volume$volume_scraped_sample_mL / area_volume$area_scraped_cm_2

# merge afdm and chla data into area_volume
all <- left_join(area_volume, afdm_final, by = c("field_date", "site_reach"))
all <- left_join(all, chla, by = c("field_date", "site_reach"))
# no chl-a data for samples RUS-3 7/20/2022 and SAL-3 7/26/2022; forgot to acidify

# calculate afdm, chla, & pheophytin per area
all$afdm_g_cm_2 <- all$afdm_g_mL * all$vol_mL_to_area_cm_2
all$Chla_ug_cm_2 <- all$Chla_ug_mL * all$vol_mL_to_area_cm_2
all$Pheo_ug_cm_2 <- all$Pheo_ug_mL * all$vol_mL_to_area_cm_2

# merge in anatoxin data
all <- left_join(all, anatoxins_final, by = c("field_date", "site_reach"))

# keep only columns we care about
all_final <- all %>%
  mutate(afdm_mg_cm_2 = afdm_g_cm_2 * 1000) %>%  # convert g/cm^2 to mg/cm^2
  select(field_date, site_reach, macroalgal_displacement_mL, area_scraped_cm_2, volume_scraped_sample_mL,
         vol_mL_to_area_cm_2, afdm_mg_cm_2, Chla_ug_cm_2, Pheo_ug_cm_2, flag_Pheo, MCY_ug_g,
         MCY_det_limit, ATXa_ug_g, dhATXa_ug_g, ATX_det_limit)

# save final dataset
write.csv(all_final, "./data/NT_rock_scrapes_misc/NT_misc_data.csv", row.names = FALSE)
