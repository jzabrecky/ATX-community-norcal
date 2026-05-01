#### Updating microscopy classification clarified by Rosalina
### Jordan Zabrecky
## 05.01.2026

# This script updates microscopy classifications from the current EDI data package 
# for this work (ver 3) by doing some refinement based on categories with the help 
# of Dr. Rosalina Stancheva to be saved into a new csv for EDI ver 4 for this work.

# In summary the changes:
# (1) Leptolyngbya and Geitlerinema will be grouped together as they are very difficult to distinguish
# (2) Most of our Homoeothrix identifications were Calothrix, so those observations will be moved to Calothrix
# (3) We acknowledge that our Rivularia observations may have been Gloeotrichia before formation of akinetes,
# so we change the name of this category to acknowledge that 
# (4) Many of our Lyngbya observations contained Homoeothrix Juliana, Lyngbya, and Coleofasciculaceae
# so we are changing that category to Misc. Oscillatoriales

#### (1) Read in Old Versions (EDI ver. 3) ####

# reading in dataframes
nontarget <- read.csv("./data/EDI_data_package/ver_3/microscopy_non_target_samples.csv")
target <- read.csv("./data/EDI_data_package/ver_3/microscopy_target_samples.csv")

# put in a list to process both at once
data <- list(nontarget, target)
names(data) <- c("nontarget", "target")

#### (2) Make changes ####

# make changes in lapply
final_data <- lapply(data, function(x) {
  y <- x %>% 
    # Change 1: merge Leptolygnbya and Geitlerinema
    mutate(leptolyngbya_and_geitlerinema = leptolyngbya + geitlerinema) %>% 
    select(!c("leptolyngbya", "geitlerinema")) %>% 
    # Change 2: move Homoeothrix to Calothrix
    mutate(calothrix = case_when(grepl("calothrix", paste(colnames(.), collapse = " ")[1]) ~  calothrix + homoeothrix)) %>% 
    select(!any_of(c("homoeothrix"))) %>% 
    # Change 3: acknowledge Rivularia may be Gloeotrichia pre akinete formation
    rename(rivularia_or_early_stage_gloeotrichia = rivularia) %>% 
    # Change 4: change Lyngbya to Miscellaneous Oscillatoriales
    dplyr::rename(miscellaneous_oscillatoriales = lyngbya)
})

# make sure columns still add to 100
lapply(final_data, function(x) any(rowSums(x[,c(7:ncol(x))]) != 100)) # all add to 100!

#### (3) Save data ####

# save csvs
write.csv(final_data$nontarget, "./data/EDI_data_package/microscopy_non_target_samples.csv", row.names = FALSE)
write.csv(final_data$target, "./data/EDI_data_package/microscopy_target_samples.csv", row.names = FALSE)
