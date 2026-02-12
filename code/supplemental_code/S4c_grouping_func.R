#### Script of functions used to add broader classifications to microscopy data
### Jordan Zabrecky
## last edited: 02.05.2026

# This script hosts functions to add broader classifications to microscopy data

#### (1) Loading libraries ####

# just one!
library(tidyverse)

#### (2) Functions ####

## (a) target_broader
# @param data is relative abundance data in long format of target taxa (TM or TAC samples)
target_broader <- function(data) {
  # add in column for broader classification
  new_data <- data %>% 
    mutate(broader = case_when(taxa == "lyngbya" | taxa == "nodularia" |  taxa == "calothrix" |
                                 taxa == "scytonema" | taxa == "gloeotrichia" | taxa == "rivularia" |
                                 taxa == "tolypothrix" ~ "Other N-fixing Cyanobacteria",
                               taxa == "nostoc" ~ "Nostoc",
                               taxa == "chroococcus" | taxa == "other_coccoids" | taxa == "aphanothece"
                               ~ "Unicellullar Cyanobacteria",
                               taxa == "anabaena_and_cylindrospermum" ~ "Anabaena or Cylindrospermum",
                               taxa == "e_diatoms" ~ "Epithemia",
                               taxa == "geitlerinema" ~ "Geitlerinema",
                               taxa == "green_algae" ~ "Green Algae",
                               taxa == "oscillatoria" | taxa == "phormidium_unknown" |
                                 taxa == "leptolyngbya" | taxa == "homoeothrix"
                               ~ "Other Filamentous Cyanobacteria",
                               taxa == "microcoleus" ~ "Microcoleus",
                               taxa == "non_e_diatoms" ~ "Diatoms Other than Epithemia",
                               taxa == "unknown" ~ "Unknown"))
  
  # return new df
  return(new_data)
}

## (b) nontarget_broader
# @param data is relative abundance data in long format of non-target taxa (NT samples)
nontarget_broader <- function(data) {
  # add in column for broader classification
  new_data <- data %>% 
    mutate(broader = case_when(taxa == "lyngbya" | taxa == "nodularia" |  taxa == "calothrix" |
                                 taxa == "scytonema" | taxa == "gloeotrichia" | taxa == "rivularia" |
                                 taxa == "tolypothrix"
                               ~ "Other N-fixing Cyanobacteria",
                               taxa == "nostoc" ~ "Nostoc",
                               taxa == "chroococcus" | taxa == "other_coccoids" | taxa == "aphanothece"
                               ~ "Unicellullar Cyanobacteria",
                               taxa == "anabaena_and_cylindrospermum" ~ "Anabaena or Cylindrospermum",
                               taxa == "oscillatoria" | taxa == "phormidium_unknown" |
                                 taxa == "leptolyngbya" | taxa == "homoeothrix" | taxa == "geitlerinema"
                               ~ "Other Filamentous Cyanobacteria",
                               taxa == "microcoleus" ~ "Microcoleus",
                               taxa == "non_e_r_diatoms" ~ "Diatoms Other than Epithemia or Rhopalodia",
                               taxa == "unknown" | taxa == "chantransia" | taxa == "euglenoid" |
                                 taxa == "unknown_green_algae"
                               ~ "Misc. Other",
                               taxa == "ankistrodesmus" | taxa == "gloeocystis" | taxa == "lacunastrum" | 
                                 taxa == "oocystis" | taxa == "pediastrum" | taxa == "scenedesmus_no_spines" |
                                 taxa == "cosmarium" | taxa == "desmodesmus_spines" | taxa == "closterium" |
                                 taxa == "coelastrum" | taxa == "stauridium" | taxa == "tetraedron"
                               ~ "Other Green Algae",
                               taxa == "cladophora" ~ "Cladophora",
                               # filamentous green algae (now grouping w/ unicellular)
                               taxa == "mougeotia" | taxa == "ulothrix" | taxa == "zygnema" |
                                 taxa == "stigeoclonium" | taxa == "oedogonium"
                               ~ "Other Green Algae",
                               taxa == "rhopalodia" | taxa == "epithemia" ~ "Epithemia or Rhopalodia",
                               taxa == "spirogyra" ~ "Spirogyra"
    ))
  
  # return new df
  return(new_data)
}

## (c) microbial_grouping
# @param data is molecular data in long format
# @rank is the column to group_by
# @param threshold is the mean relative abundance under which the phylum 
# gets categorized under "Other"
# @param cyanobacteria is binary for if we want to filter for cyanobacteria phylum only
microbial_grouping <- function(data, rank, threshold, cyanobacteria = FALSE) {
  
  # optional- filter for cyanobacteria
  if(cyanobacteria) {
      data <- data %>% 
        filter(phylum == "Cyanobacteria")
  }
  
  # calculate average abundance for each phylum
  grouped <- data %>% 
    dplyr::group_by(.data[[rank]]) %>% 
    dplyr::summarize(mean = mean(relative_abundance)) %>% 
    # put into "Other" category if NA or % is less than #%
    mutate(broader = case_when(mean < threshold  ~ "Other"))
  
  # put into "Other" if NA
  grouped[which(apply(grouped[,1], 1, function(x) str_detect(x, "NA"))), 3] <- "Other"
  # keep rank if not NA
  grouped[which(is.na(grouped[,3])), 3] <- grouped[which(is.na(grouped[,3])), 1]
  
  # add into regular data (after removing mean column)
  new <- left_join(data, grouped %>% select(!mean), by = c(rank))
  
  return(new)
}

