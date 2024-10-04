#### Code to process and organize miscellaneous rock scrape/NT for EDI data package
### Jordan Zabrecky
## last edited 10.04.2024

# This code gathers miscellaneous rock scrape data that is not currently used for 
# either ATX projects. Data includes area scraped, total volume of sample, 
# displacement of preserved sample, chlorophyll-a of filtered sample, 
# AFDM of filtered sample, and, if present, detected anatoxins of sample

#### (1) Loading libraries & data ####

# loading libraries
lapply(c("tidyverse", "lubridate", "plyr"), require, character.only = T)

# area scraped and sample volume data
area_volume <- read.csv("./data/NT_rock_scrapes_misc/raw/volumes_and_area.csv")

# filtered chlorophyll-a data (processed in script "1a_rock_scrape_filter_chla.R")

# filtered ash-free dry mass data
afdm <- ldply(list.files(path = "./data/NT_rock_scrapes_misc/raw/", 
                         pattern = "angelo"), function(filename) {
                           d <- read.csv(paste("./data/NT_rock_scrapes_misc/raw/", filename, sep = ""))
                           d$field_date <- mdy(d$field_date)
                           return(d)
                         })


# rock scrape anatoxin data

# changing date fields to date objects
area_volume$field_date <- mdy(area_volume$field_date)
afdm$field_date <- ymd(afdm$field_date)

#### (2) Processing ash-free dry mass data ####

# calculate afdm for all samples
afdm$ash_free_dry_mass_g_raw <- afdm$dry_sample_filter_tin_g - afdm$ash_filter_tin_g

# gather blanks of samples and do some editing to the dataframe
afdm_blanks <- afdm %>% 
  filter(site_reach == "blank") %>% 
  dplyr::rename(lab_date = field_date,
                blank_afdm_g = ash_free_dry_mass_raw) %>% 
  select(lab_date, )

# thoughts gather blanks and left join via date? as column

#### (3) Processing anatoxins data ####c

## change field_date to year_month_day format before saving!!