#### Code to process and organize miscellaneous rock scrape/NT for EDI data package
### Jordan Zabrecky
## last edited 01.23.2025

# This code gathers miscellaneous rock scrape data that is not currently used for 
# either ATX projects. Data includes total area scraped, total volume of sample, 
# displacement of preserved sample, chlorophyll-a and pheophytin of sample and amount
# filtered for chlorophyll-a processing, and dry mass, ashed mass, ash-free dry mass
# and amount filtered for afdm processing.

#### (1) Loading libraries & data ####

# loading libraries
lapply(c("tidyverse", "lubridate", "plyr"), require, character.only = T)

# area scraped and sample volume data
area_volume <- read.csv("./data/NT_misc/raw/volumes_and_area.csv")

# filtered chlorophyll-a data (processed in script "1a_rock_scrape_filter_chla.R")
chla <- read.csv("./data/NT_misc/NT_chlorophylla_notprocessed.csv")

# filtered ash-free dry mass data
afdm <- ldply(list.files(path = "./data/NT_misc/raw/", 
                         pattern = "angelo"), function(filename) {
                           d <- read.csv(paste("./data/NT_misc/raw/", filename, sep = ""))
                           d$field_date <- mdy(d$field_date)
                           return(d)
                         })

# changing date fields to date objects
area_volume$field_date <- mdy(area_volume$field_date)
afdm$field_date <- ymd(afdm$field_date)
chla$field_date <- mdy(chla$field_date)

# changing label to clarify what analyses the filtered sample amount was for
afdm <- afdm %>% 
  dplyr::rename(afdm_sample_filtered_mL = filtered_amount_mL)
chla <- chla %>% 
  dplyr::rename(chla_sample_filtered_mL = volume_filtered_mL)

#### (2) Joining csv's together ####

# add number tags to triplicates so we done have 3x3 when joining
assign_triplicate_numbers <- function(df) {
  
  # separate triplicates and non_triplicates into two data frames
  list <- split(df, df$triplicate)
  
  # assign triplicate_no column which is 0 is not a triplicate and alternating 1-3 for triplicate
  list$n$triplicate_no <- rep(0, length(nrow(list$n)))
  list$y$triplicate_no <- rep(1:3, length(nrow(list$y)))
  
  # reduce list and reorder by date
  final <- rbind(list$n, list$y) %>% 
    arrange(field_date) %>% 
    relocate(triplicate_no, .after = triplicate)
  
  # return data frame
  return(final)
}

# apply function to afdm and chlorophyll-a data frames
afdm <- assign_triplicate_numbers(afdm)
chla <- assign_triplicate_numbers(chla)

# calculate afdm with dry mass and ash (inorganic) mass
afdm <- afdm %>% 
  mutate(afdm_g = dry_sample_filter_tin_g - ash_filter_tin_g) %>% 
  relocate(afdm_g, .after = ash_filter_tin_g)

# join chlorophyll-a and afdm
all <- join(afdm, chla, by = c("field_date", "site_reach", "triplicate", "triplicate_no"))

# create temporary column to join with area_volume dataframe (and not have duplicates)
all <- all %>% 
  mutate(first_sample = case_when(triplicate_no == 0 ~ "y",
                                  triplicate_no == 1 ~ "y",
                                  TRUE ~ "n"))
area_volume$first_sample <- "y"

# joining together all three
final <- left_join(all, area_volume, by = c("site_reach", "field_date", "first_sample")) %>% 
  select(!first_sample) %>% 
  relocate(macroalgal_displacement_mL, .after = triplicate_no) %>% 
  relocate(volume_scraped_mL, .after = triplicate_no) %>% 
  relocate(area_scraped_cm_2, .after = triplicate_no)

# save final csv
write.csv(final, "./data/EDI_data_package/non_target_sample_miscellaneous.csv", row.names = FALSE)
