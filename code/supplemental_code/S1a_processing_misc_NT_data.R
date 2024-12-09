#### Processing miscellaneous Non-target data
### Jordan Zabrecky
## last edited 12.09.2024

# This code analyzes the triplicates and blanks of non-target data
# to create a more processed final csv that may or may not be used

#### (1) Loading libraries and data



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