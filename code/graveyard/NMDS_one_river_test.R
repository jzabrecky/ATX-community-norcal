library(tidyverse)

#### morphological

# read in files (data transformed in previous script, "4a_amongrivers_microscopy.R")
data <- lapply(list.files(path = "./data/morphological/transformed/", pattern = ".csv"),
               function(x) read.csv(paste("./data/morphological/transformed/", x, sep = "")))
names(data) <- c("nt", "tac", "tm")

add_event_no <- function(data) {
  data %>% 
    mutate(field_date = ymd(field_date),
           month = month(field_date)) %>% 
    mutate(event_no = case_when((field_date >= ymd("2022-06-24") & field_date <= ymd("2022-06-29")) ~ 1,
                                (field_date >= ymd("2022-07-06") & field_date <= ymd("2022-07-14")) ~ 2,
                                (field_date >= ymd("2022-07-20") & field_date <= ymd("2022-07-28")) ~ 3,
                                (field_date >= ymd("2022-08-02") & field_date <= ymd("2022-08-10")) ~ 4,
                                (field_date >= ymd("2022-08-17") & field_date <= ymd("2022-08-23")) ~ 5,
                                (field_date >= ymd("2022-09-01") & field_date <= ymd("2022-09-06")) ~ 6,
                                (field_date >= ymd("2022-09-15") & field_date <= ymd("2022-09-22")) ~ 7)) %>% 
    relocate(event_no, .before = "field_date") %>% 
    relocate(month, .before = "field_date")
}

data <- lapply(data, function(x) add_event_no(x))

setdiff(colnames(data$nt), colnames(data$tm)) # gives differences in column names
colnames(data$tm) # what to match it into
nt_unified <- data$nt %>% 
  pivot_longer(c(7:ncol(data$nt)), names_to = "nt_taxa", values_to = "hellinger_value") %>% 
  mutate(tac_tm_taxa = case_when(nt_taxa == "ankistrodesmus" ~ "green_algae",
                                 nt_taxa == "chantransia" ~ "unknown", # treating unknown as other
                                 nt_taxa == "cladophora" ~ "green_algae",
                                 nt_taxa == "closterium" ~ "green_algae",
                                 nt_taxa == "coelastrum" ~ "green_algae",
                                 nt_taxa == "cosmarium" ~ "green_algae",
                                 nt_taxa == "desmodesmus_spines" ~ "green_algae",
                                 nt_taxa == "epithemia" ~ "e_diatoms",
                                 nt_taxa == "euglenoid" ~ "unknown",
                                 nt_taxa == "gloeocystis" ~ "green_algae",
                                 nt_taxa == "lacunastrum" ~ "green_algae",
                                 nt_taxa == "mougeotia" ~ "green_algae",
                                 nt_taxa == "non_e_r_diatoms" ~ "non_e_diatoms",
                                 nt_taxa == "oedogonium" ~ "green_algae",
                                 nt_taxa == "oocystis" ~ "green_algae",
                                 nt_taxa == "pediastrum" ~ "green_algae",
                                 nt_taxa == "rhopalodia" ~ "non_e_diatoms",
                                 nt_taxa == "scenedesmus_no_spines" ~ "green_algae",
                                 nt_taxa == "spirogyra" ~ "green_algae",
                                 nt_taxa == "stauridium" ~ "green_algae",
                                 nt_taxa == "stigeoclonium" ~ "green_algae",
                                 nt_taxa == "tetraedron" ~ "green_algae",
                                 nt_taxa == "ulothrix" ~ "green_algae",
                                 nt_taxa == "unknown_green_algae" ~ "green_algae",
                                 nt_taxa == "zygnema" ~ "zygnema",
                                 TRUE ~ nt_taxa)) %>% 
  select(!nt_taxa) %>% 
  pivot_wider(names_from = tac_tm_taxa, values_from = hellinger_value, values_fn = sum)

# separate out by river
tm <- split(data$tm, data$tm$site)
tac <- split(data$tac, data$tac$site)
nt <- split(nt_unified, nt_unified$site)

# put together by river
sfkeel <- bind_rows(tm$`SFE-M`, tac$`SFE-M`, nt$`SFE-M`)
sfkeel <- replace(sfkeel, is.na(sfkeel), 0)

# source functions
source("./code/supplemental_code/S4b_community_analyses_func.R")

eelNMDS <- getNMDSdata(sfkeel, 7, ASV = TRUE) # some issue with vectors, so leaving it out
makeNMDSplot(eelNMDS, FALSE, FALSE, color = "month", shape = "sample_type")
# event no not released with current NMDS function

eelTMNMDS <- getNMDSdata(tm$`SFE-M`, 7, ASV = TRUE)
makeNMDSplot(eelTMNMDS, FALSE, FALSE, color = "month", shape = "month")

rusNMDS <- getNMDSdata(tac$RUS, 7, ASV = TRUE)
makeNMDSplot(rusNMDS, FALSE, FALSE, color = "month", shape = "month")

eeltac_NMDS <- getNMDSdata(tac$`SFE-M`, 7, ASV = TRUE)
makeNMDSplot(eeltac_NMDS, FALSE, FALSE, color = "month", shape = "month")


eelNT_NMDS <- getNMDSdata(nt$`SFE-M`, 7, ASV = TRUE)
makeNMDSplot(eelNT_NMDS, FALSE, FALSE, color = "month", shape = "month")
