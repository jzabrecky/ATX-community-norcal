#### Processing predicted orthologs for all samples
### Jordan Zabrecky
## last edited: 04.09.2026

# This script processes the predicted functional groups (as KEGG orthologs 
# predicted genes) from PICRUSt2-SC standard pipeline. 
# This script saves all KO's which are further processed for ones we care about in
# the script 3c
# Note: This will only be done for rarefied (95% threshold) & filtered .biom data

#### (1) Loading libraries & data ####

# loading libraries
lapply(c("tidyverse", "ggpicrust2"), require, character.only = T)

## (a) getting metadata

# get list of file names
plate_files <- list.files(path = "./data/molecular/metadata/", pattern = "plate")

# create empty list (different plates have same IDs so want to keep them separate!)
plate_data <- list()

# fill in list with dataframes
for(i in 1:length(plate_files)) {
  plate_data[[i]] <- read.table(paste("data/molecular/metadata/", plate_files[i], sep = ""))
  colnames(plate_data[[i]]) <- c("plate_ID", "vial_ID")
}

# preface plate ID's with p2 or p3 if 2nd or 3rd plate respectively
plate_data[[2]]$plate_ID <- paste("p2", plate_data[[2]]$plate_ID, sep = "")
plate_data[[3]]$plate_ID <- paste("p3", plate_data[[3]]$plate_ID, sep = "")

# put into one dataframe
plate_data_final <- plate_data[[1]]
for(i in 2:length(plate_data)) {
  plate_data_final <- rbind(plate_data_final, plate_data[[i]])
}


# read in sample metadata & add into plate data
metadata <- read.csv("./data/molecular/metadata/16s_sample_metadata.csv")

# join in plate_data_final with metadata
metadata <- left_join(metadata, plate_data_final, by = c("vial_ID"))

## (b) reading in PICRUSt2-SC predicted genes
picrust_predfunc_all <- read_tsv("./data/molecular/picrust2_outputs/all/pred_metagenome_unstrat.tsv") %>% 
  mutate(`function` = str_remove(`function`, "ko:")) %>% 
  dplyr::rename(ko_id = `function`) %>% 
  pivot_longer(cols = c(2:ncol(.)) ,  names_to = "plate_ID",
               values_to = "abundance") %>% filter(abundance != 0)
picrust_predfunc_tm <- read_tsv("./data/molecular/picrust2_outputs/tm_nomicro/pred_metagenome_unstrat.tsv") %>% 
  mutate(`function` = str_remove(`function`, "ko:")) %>% 
  dplyr::rename(ko_id = `function`) %>% 
  pivot_longer(cols = c(2:ncol(.)) ,  names_to = "plate_ID",
               values_to = "abundance") %>% filter(abundance != 0)
picrust_predfunc_tac <- read_tsv("./data/molecular/picrust2_outputs/tac_noanacyl/pred_metagenome_unstrat.tsv") %>% 
  mutate(`function` = str_remove(`function`, "ko:")) %>% 
  dplyr::rename(ko_id = `function`) %>% 
  pivot_longer(cols = c(2:ncol(.)) ,  names_to = "plate_ID",
               values_to = "abundance") %>% filter(abundance != 0)

# put into list
picrust_predfunc_list <- list(picrust_predfunc_all, picrust_predfunc_tm, picrust_predfunc_tac)
names(picrust_predfunc_list) <- c("all", "tm", "tac")

#### (2) Joining Data Together ####

# pivot longer and join with plate metadata
picrust_predfunc_list <- lapply(picrust_predfunc_list, function(x) left_join(x, metadata, by = "plate_ID"))

# deal with a couple samples as also dealt with in qiime2 processing:

# (i)  very unsure about TAC at SAL-2 on 9/22/22 (only took a little for 16s)
# in the field, so will probably just remove (not as clearly Anabaena as the sample at SAL-3)
# despite some reads being anabaena! (is vial 230)
picrust_predfunc_list <- lapply(picrust_predfunc_list, function(x) {
  y <- x %>% 
    filter(vial_ID != 230)
  return(y)
})

# (ii) fix field date label issue
# remove TM sample taken on 9/6 at SFE-M-1S (taken inconsistently from other samples; turkey-baster, lots of water)
# & replace with sample taken on 9/8 (taken correctly; mostly mat )
picrust_predfunc_list$all <- picrust_predfunc_list$all[-which(picrust_predfunc_list$all$field_date == "9/6/2022" & picrust_predfunc_list$all$site_reach == "SFE-M-1S" & picrust_predfunc_list$all$sample_type == "TM"),]
picrust_predfunc_list$all[which(picrust_predfunc_list$all$field_date == "9/8/2022" & picrust_predfunc_list$all$site_reach == "SFE-M-1S"& picrust_predfunc_list$all$sample_type == "TM"),]$field_date <- "9/6/2022"
picrust_predfunc_list$tm <- picrust_predfunc_list$tm[-which(picrust_predfunc_list$tm$field_date == "9/6/2022" & picrust_predfunc_list$tm$site_reach == "SFE-M-1S" & picrust_predfunc_list$tm$sample_type == "TM"),]
picrust_predfunc_list$tm[which(picrust_predfunc_list$tm$field_date == "9/8/2022" & picrust_predfunc_list$tm$site_reach == "SFE-M-1S"& picrust_predfunc_list$tm$sample_type == "TM"),]$field_date <- "9/6/2022"

#### (4) Processing Triplicates ####

# filter out for triplicates and not triplicates
triplicates <- lapply(picrust_predfunc_list, function(x) {
  y <- x %>% 
    filter(triplicate == "y") %>% 
    mutate(full_sample_name = paste(site_reach, field_date, sample_type))
  return(y)
})

# remove odd samples as identified in prior processing scripts
IDs_to_remove <- c(110, 225, 50, 123)
triplicates_adjusted <- lapply(triplicates, function(x) {
  y <- x %>% 
    filter(!vial_ID %in% IDs_to_remove)
  return(y)
})

# average across triplicates
triplicates_adjusted <- lapply(triplicates_adjusted, function(x) {
  y <- x %>% 
    dplyr::group_by(site_reach, site, field_date, sample_type, triplicate, full_sample_name, ko_id, abundance) %>% 
    # do mean of predicted genes
    dplyr::summarize(mean_abundance = mean(abundance)) %>%  
    ungroup() %>% 
    select(site_reach, site, field_date, sample_type, ko_id, mean_abundance)
  return(y)
})

# remove excess files
rm(triplicates)

#### (5) Putting Data Together ####

# remove triplicates from original data
final_data <- lapply(picrust_predfunc_list, function(x) {
  y <- x %>% 
    filter(triplicate == "n") %>% 
    # filter out fake_targets and blanks %>% 
    filter(sample_type!= "blank" & fake_target == "n") %>% 
    select(site_reach, site, field_date, sample_type, ko_id, abundance)
})

# merge in averaged & re-relativized
final_data <- lapply(names(final_data), function(x) {
  y <- rbind(final_data[[x]] %>% dplyr::rename(predicted_gene_abundance = abundance),
             triplicates_adjusted[[x]] %>% dplyr::rename(predicted_gene_abundance = mean_abundance))
  return(y)
})
names(final_data) <- c("all", "tm_nomicro", "tac_noanacyl")

# save data
lapply(names(final_data), function(x) write.csv(final_data[[x]],
                                                paste("./data/molecular/PICRUSt2_predicted_KO_all_", x, ".csv", sep = ""),
                                                row.names = FALSE))

#### (6) Quickly Look at Prominent Pathways ####

# for analyses, I will only look at pathways I care about, 
# but am curious what the most abundant pathways are broadly speaking

# note: single KEGG ortholog IDs can match to multiple pathways
# and the function below makes weird choices sometimes, 
# but since we want to get a quick overview on everything first,
# we will just use that function
kegg_pathways <- ko2kegg_abundance("./data/molecular/picrust2_outputs/all/pred_metagenome_unstrat.tsv") %>%
  mutate(pathway_id = rownames(.)) %>% 
  relocate(pathway_id, .before = 1) %>% 
  pivot_longer(cols = c(2:ncol(.)) ,  names_to = "plate_ID",
               values_to = "abundance") %>% filter(abundance != 0) %>% 
  # join in metadata
  left_join(metadata, by = "plate_ID") %>% 
  # join in pathway information
  left_join(ko_to_kegg_reference %>% select(pathway_id, level1, level2, level3) %>% distinct(), by = "pathway_id") 

# Q1: What is most abundant at level 1?
levelone <- kegg_pathways %>% 
  dplyr::group_by(sample_type, site, level1) %>% 
  dplyr::summarize(total = sum(abundance)) %>% 
  slice(1:5)
# for all (1) Metabolism, (2) Genetic Information Processing, (3) Environmental Information Processing,
# (4) Cellular Processes, (5) Organismal Systems
# function for filtering for prokayotes should not have organismal systems so high?
# again, ggpicrust2 makes weird choices, so I am assuming those ko's match other things..

# Q2: What is most abundant at level 2?
leveltwo <- kegg_pathways %>% 
  dplyr::group_by(sample_type, site, level2) %>% 
  dplyr::summarize(total = sum(abundance)) %>% 
  slice(1:5)
# carbohydrate metabolism, energy metabolism, lipid metabolism, nucleotide metabolism,
# amino acid metabolism for all

# Q3: what is most abundant at level 3?
levelthree <- kegg_pathways %>% 
  dplyr::group_by(sample_type, site, level3) %>% 
  dplyr::summarize(total = sum(abundance)) %>% 
  slice(1:5)
