#### draft script for processing PICRUSt2-SC data

# ideal overview: read in data, match to 95% rarefied data, 
# figure out broader functional grouping and add in

# see how abundance matches on each dataframe?
# consider saving outputs from script 2d elsewhere
# and having a script 3c where predicted functional groups are verified

library(ggpicrust2)
library(tidyverse)

# read in files (following ggpicrust2 tutorial)
abundance_file <- "./data/molecular/picrust2_outputs/pred_metagenome_unstrat1.tsv"

# get list of file names
plate_files <- list.files(path = "./data/molecular/metadata/", pattern = "plate")

# create empty list (different plates have same IDs so want to keep them separate!)
plate_data <- list()

# fill in list with dataframes
for(i in 1:length(plate_files)) {
  plate_data[[i]] <- read.table(paste("data/molecular/metadata/", plate_files[i], sep = ""))
  colnames(plate_data[[i]]) <- c("plate_ID", "vial_ID")
}

test = read_tsv("./data/molecular/picrust2_outputs/pred_metagenome_unstrat1.tsv") %>% 
  mutate(`function` = str_remove(`function`, "ko:"))
test2 = read_tsv("./data/molecular/picrust2_outputs/pred_metagenome_unstrat2.tsv") %>% 
  mutate(`function` = str_remove(`function`, "ko:"))
#test = ko2kegg_abundance("./data/molecular/picrust2_outputs/pred_metagenome_unstrat1.tsv")
#test = cbind(rownames(test), test)
#test2 = ko2kegg_abundance("./data/molecular/picrust2_outputs/pred_metagenome_unstrat2.tsv")
#test2 = cbind(rownames(test2), test2)
#colnames(test)[1] = "pathway_id"
#colnames(test2)[1] = "pathway_id"

view(test)

# pivot longer
test_longer <- test %>% 
  pivot_longer(c(2:ncol(test)), values_to = "abundance", names_to = "plate_ID")
test_longer2 <- test2 %>% 
  pivot_longer(c(2:ncol(test2)), values_to = "abundance", names_to = "plate_ID")

together1 <- left_join(test_longer, plate_data[[1]],
                       by = "plate_ID")
together2 <- left_join(test_longer2, plate_data[[2]],
                       by = "plate_ID")

together_all <- rbind(together1, together2)
together_all <- left_join(together_all, metadata,
                          by = "vial_ID")
colnames(together_all)[1] <- "ko_id"

# need to figure out what each kegg pathway is; use ggpicrust2
view(ko_to_kegg_reference)

together_all <- left_join(together_all, ko_to_kegg_reference, by = "ko_id")

# struggling to match the ko and kegg abundances??

# things to do:
# look at NMDS plots

test = pathway_annotation(file = "./data/molecular/picrust2_outputs/pred_metagenome_unstrat1.tsv",
                          pathway = "KO",
                          ko_to_kegg = FALSE)
