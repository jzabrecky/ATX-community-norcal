#### Processing chlorophyll-a and pheophytin concentrations for freeze-dried cyanobacteria mass
### Jordan Zabrecky (modified from Joanna Blaszczak)
## last modified 12.09.2024

# This code calculates chlorophyll-a and pheophytin concentrations for filtered
# non-target/rock scrape from RFUs obtained from the Blaszczak Lab's Trilogy 
# Fluorometer using the 2023-2024 calibration (2022 samples were analyzed in 
# early 2023 and 2023 samples in late 2023). Chl-a is given in mg/L of the sample
# (conversion to mg/cm^2 is in script "pre_EDI_1c_gathering_misc_rock_scrape_data.R")

# link to calculation power point: https://docs.google.com/presentation/d/1qoOF3yhsAzrzcR9l_J0bSWqni8gjoPfq/edit#slide=id.p3

#### (1) Loading libraries & data ####

# import packages
lapply(c("plyr","dplyr", "tidyverse", "lubridate"), require, character.only=T)

# set working directory to lead in data
setwd("./data/NT_misc/raw/chla")

# RFU data
rawRFU <- ldply(list.files(pattern = "_RFUdata.csv"), function(filename) {
  d <- read.csv(filename)
  d$file <- filename
  d$analysis_date <- sapply(strsplit(d$file,"_"), `[`, 1)
  return(d)
})

# metadata
metadata <- ldply(list.files(pattern = "_metadata.csv"), function(filename) {
  d <- read.csv(filename)
  d$file <- filename
  d$analysis_date <- sapply(strsplit(d$file,"_"), `[`, 1)
  return(d)
})

# check data frame formatting
sapply(rawRFU, class)
sapply(metadata, class)

# fill in volume filtered for blanks as 15 mL as we typically filtered that amount of tap water
metadata$volume_filtered_mL <- replace(metadata$volume_filtered_mL, 
                                       is.na(metadata$volume_filtered_mL), 15)

# split into list of dataframes
rawRFU_list <- split(rawRFU, rawRFU$analysis_date)
metadata_list <- split(metadata, metadata$analysis_date)

# merge all by analysis date and Rack_ID
merged_list <- list()
for(i in 1:length(names(rawRFU_list))){
  merged_list[[i]] <- left_join(rawRFU_list[[i]][,c("Fluorometer_ID","RFU","Units","Unacidified_acidified","Rack_ID","Extraction_vol")],
                                metadata_list[[i]][,c("Rack_ID", "site_reach", "field_date", "triplicate", "volume_filtered_mL")],
                                by = "Rack_ID", multiple = "all")
}

# add analysis date as name for each dataframe in merged list
names(merged_list) <- names(rawRFU_list)

# check an individual dataframe from analysis date
View(merged_list$'20230407')

# remove samples where acidified was forgotten (cannot calculate chl-a on those)
merged_list$"20230407" <- merged_list$"20230407"[-which(merged_list$"20230407"$Rack_ID == 25),]
merged_list$"20230407" <- merged_list$"20230407"[-which(merged_list$"20230407"$Rack_ID == 31),]

#### (2) Define the standard curves for low and high ####

# RFU cutoff between low and high standard curve = 1617913.38
# calibration curve for 2/1/2023 - 1/31/2024: https://docs.google.com/spreadsheets/d/1qO1-9TI8dwxcrY6I3o4K0haUapuAcIKJRSivkrpgS8o/edit#gid=400761777
low_Fm <- 1.79; high_Fm <- 1.8
low_Fa <- 49339.41; high_Fa <- 453486.34
low_Fb <- 85555.48; high_Fb <- 807198.44
low_std_C1 <- 25; high_std_C1 <- 250
low_blank <- 3547.67; high_blank <- 8740.07

#### (3) Calculating chlorophyll-a and pheophytin concentrations ####

# function to calculate chlorophyll-a and pheophytin
chla_pheo_calc <- function(ex){
  
  # define volume of the solvent and mass of sample
  vol_solvent <- ex$Extraction_vol[1]
  vol_sample <- ex$volume_filtered_mL[1]
  
  # next, determine which curve type to use
  ex$curve_type <- ifelse(ex$RFU <= 1617913.38,
                          yes = "low", no = "high")
  
  # determine if acidified or unacidified RFU below detection limit 
  # (lowest standard from lowest cal curve)
  ex$RFU_bdl_sep <- (ex$RFU < 85555.48)
  # for reference, highest standard is 18570952 RFU so over detection
  # limit shouldn't be an issue
  
  # define parameters for low and high curves
  Int_B_low <- low_std_C1*((ex[which(ex$Unacidified_acidified == "U"),]$RFU - low_blank)/(low_Fb - low_blank))
  Int_A_low <- low_std_C1*((ex[which(ex$Unacidified_acidified == "A"),]$RFU - low_blank)/(low_Fb - low_blank))
  Int_B_high <- high_std_C1*((ex[which(ex$Unacidified_acidified == "U"),]$RFU - high_blank)/(high_Fb - high_blank))
  Int_A_high <- high_std_C1*((ex[which(ex$Unacidified_acidified == "A"),]$RFU - high_blank)/(high_Fb - high_blank))
  
  # apply interpolation based on whether its a high or low curve
  ex$Chla_ug_L_raw <- ifelse(ex$curve_type == "low",
                              yes = (low_Fm/(low_Fm - 1))*(Int_B_low - Int_A_low)*(vol_solvent/vol_sample), # (vol mL/vol mL is unitless factor)
                              no = (high_Fm/(high_Fm - 1))*(Int_B_high - Int_A_high)*(vol_solvent/vol_sample))
  
  ex$Pheo_ug_L_raw <- ifelse(ex$curve_type == "low",
                              yes = (low_Fm/(low_Fm - 1))*((low_Fm*Int_A_low) - Int_B_low)*(vol_solvent/vol_sample),
                              no = (high_Fm/(high_Fm - 1))*((high_Fm*Int_A_high) - Int_B_high)*(vol_solvent/vol_sample))
  
  # return true if either acidified or unacidified RFU is below detection limit
  ex$RFU_bdl <- ifelse(ex$RFU_bdl_sep[1] == TRUE | ex$RFU_bdl_sep[2] == TRUE, TRUE, FALSE)
  
  # preserve dilution information to potentially adjust for future chlorophyll-a runs
  output <- as.data.frame(cbind(ex$Chla_ug_L[1], ex$Pheo_ug_L[1], ex$RFU_bdl[1]))
  
  return(output)
}

# function to apply the above function across list of dataframes
chla_pheo_processing <- function(dat){
  # split data by Rack_ID
  df_rack_split <- split(dat, dat$Rack_ID)
  # apply chlorophyll-a and pheophytin calculation function
  chla_Pheo_df <- ldply(lapply(df_rack_split, function(x) chla_pheo_calc(x)), data.frame)
  # column names for output
  colnames(chla_Pheo_df) <- c("Rack_ID", "Chla_ug_L","Pheo_ug_L", "RFU_bdl")
  
  return(chla_Pheo_df)
}

# applying function to list of merged data frames
processed_data <- ldply(lapply(merged_list, function(x) chla_pheo_processing(x)))

# adding in name of first column
colnames(processed_data)[1] <- "analysis_date"

# merge again with metadata
processed_data$Rack_ID <- as.character(processed_data$Rack_ID) # convert to character to match
metadata$Rack_ID <- as.character(metadata$Rack_ID)
chla_pheo_calculations <- left_join(na.omit(processed_data),
                                    metadata[,c("Rack_ID","site_reach","field_date", "triplicate", "volume_filtered_mL", "analysis_date")],
                                    by = c("analysis_date","Rack_ID"))

#### (4) Final processing of chlorophyll-a and pheophytin calculations ####

# flag for negative pheophytin values and only keep columns we care about
chla_pheo_final <- chla_pheo_calculations %>% 
  mutate(neg_Pheo_flag = case_when(Pheo_ug_L < 0 ~ "y",
                                   TRUE ~ "n"),
         Chla_Pheo_flag = case_when(RFU_bdl == 1 ~ "y",
                                    TRUE ~ "n")) %>% 
  select(field_date, site_reach, triplicate, Chla_ug_L, Pheo_ug_L, Chla_Pheo_flag, neg_Pheo_flag, volume_filtered_mL)

# note: no chl-a data for samples RUS-3 7/20/2022 and SAL-3 7/26/2022; forgot to acidify

# save final csv
write.csv(chla_pheo_final, "../../../NT_misc/NT_chlorophylla_notprocessed.csv", row.names = FALSE)
