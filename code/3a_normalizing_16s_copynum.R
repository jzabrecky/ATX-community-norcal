#### Normalizing relative abundance with predicted 16s gene copy numbers
### Jordan Zabrecky
## last edited: 03.20.2026

# This script pulls in relative abundance based on the predicted 16s gene 
# copy numbers for each ASV obtained via PICRUSt2-SC and compares it to 
# the relative abundance without doing so (QIIME2 outputs obtained and
# processed in step 2 of code)

# Note: This will only be done for rarefied (95% threshold) & processed data

#### (1) Loading libraries & data ####

# loading libraries
lapply(c("tidyverse", "cowplot"), require, character.only = T)

# functions to quickly plot data
plot_phylum <- function(data, abundance_var, x_var) {
  plot <- ggplot(data = data) +
    geom_bar(aes(x = .data[[x_var]], y = .data[[abundance_var]], fill = phylum), 
             stat = "identity", position = "fill")
  return(plot)
}
plot_cyano <- function(data, abundance_var, x_var) {
  # filter out for phylum cyanobacteria
  cyano_only = data %>%
    filter(phylum == "Cyanobacteria")
  plot <- ggplot(data = cyano_only) +
    geom_bar(aes(x = .data[[x_var]], y = .data[[abundance_var]], fill = genus), 
             stat = "identity", position = "fill")
  return(plot)
  
}

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
metadata <- left_join(read.csv("./data/molecular/metadata/16s_sample_metadata.csv"), plate_data_final,
                      by = "vial_ID")

## (b) reading in PICRUSt2-SC NSTI values (whereby lower values closer to 0 indicate better match)
picrust_nsti <- read_tsv("./data/molecular/picrust2_outputs/weighted_nsti.tsv") %>% 
  dplyr::rename(plate_ID = sample)

## (c) reading in PICRUSt2-SC estimated 16s rRNA copy number values
picrust_predabun <- read_tsv("./data/molecular/picrust2_outputs/seqtab_norm.tsv") %>% 
  pivot_longer(cols = c(2:ncol(.)),  names_to = "plate_ID",
               values_to = "picrust2_abundance") %>% 
  # remove filler zero values from wide format
  filter(picrust2_abundance != 0)

## (d) reading in QIIME2 outputs from previous processing steps (code step 2)
qiime2 <- read.csv("./data/molecular/16s_nochimera_rarefied_95_endcode2.csv")

#### (2) Join Data Together ####

# merge in metadata with nsti values
data <- left_join(picrust_nsti, metadata, by = "plate_ID")

# join in abundances normalized by 16s copy number
data <- left_join(data, picrust_predabun, by = "plate_ID") %>% 
  # when we read in table feature ID got put under normalized!
  dplyr::rename(feature_ID = normalized)

# missing vial check
setdiff(unique(metadata$vial_ID), unique(data$vial_ID)) 
# 4, 18, 34, 40, 43, 58, 99, 102, 104, 115, 135, 137, 142, 208, 215

# this is also assessed in "2b_reading_in_QIIME_outputs_rarefied.R", copy & pasted results below:

# we know from looking at nonrarefied data that:
# 104, 108, and 208 were lost in translation
# 18, 34, and 135 had poor sequencing quality

## lost from 95% sequencing depth rarefaction:
# 4: SAL-3 blank 6-27-2022 (inconsequential!)
# 40: SAL-2 TM 7-26-2022
# 43: SFE-M-3 TAC 7-14-2022
# 58: SFE-M-3 TAC 7-28-2022 (triplicate- inconsequential!)
# 99: RUS-1S TM 8-17-2022 (fake target- inconsequential!)
# 102: RUS-1S blank 8-17-2022 (blank- inconsequential!)
# 115: RUS-2 TAC 9-1-2022
# 137: SFE-M-1S blank 9-6-2022 (blank- inconsequential!)
# 142: RUS-1S TM 9-1-2022 (fake target- inconsequential!)
# 215: RUS-1S TM 9-15-2022 (fake target- inconsequential!)

# lastly, deal with a couple samples as also dealt with in qiime2 processing:

# (i)  very unsure about TAC at SAL-2 on 9/22/22 (only took a little for 16s)
# in the field, so will probably just remove (not as clearly Anabaena as the sample at SAL-3)
# despite some reads being anabaena! (is vial 230)
data <- data %>% 
  filter(vial_ID != 230)

# (ii) fix field date label issue
# remove TM sample taken on 9/6 at SFE-M-1S (taken inconsistently from other samples; turkey-baster, lots of water)
# & replace with sample taken on 9/8 (taken correctly; mostly mat )
data <- data[-which(data$field_date == "9/6/2022" & data$site_reach == "SFE-M-1S" & data$sample_type == "TM"),]
data[which(data$field_date == "9/8/2022" & data$site_reach == "SFE-M-1S"& data$sample_type == "TM"),]$field_date <- "9/6/2022"

#### (3) Relativizing Normalized Abundances ####

# first save normalized abundances that are not relativized and include triplicates
write.csv(data, "./data/molecular/intermediate_csvs/16s_picrustnormalizedabundances_rarefied95_filtered.csv")

# calculate total abundances per vial
total_abundance_per_vial <- data %>% 
  dplyr::group_by(vial_ID) %>% 
  dplyr::summarize(total_abundance = sum(picrust2_abundance)) %>% 
  ungroup()

# left join in this data to full dataframe and calculate relative abundance
data <- left_join(data, total_abundance_per_vial, by = c("vial_ID")) %>% 
  mutate(picrust2_relative_abundance = picrust2_abundance / total_abundance * 100) %>% # doing % to avoid really small decimals
  relocate(total_abundance, .before = feature_ID) %>% 
  relocate(picrust2_relative_abundance, .before = feature_ID)

# check to make sure all add up to 100%
relativized_check <- data %>% 
  dplyr::group_by(vial_ID) %>% 
  dplyr::summarize(total = sum(picrust2_relative_abundance)) %>% 
  ungroup()
relativized_check$total[which(relativized_check$total != 100)] 
# all that don't are very close to 100/reading in as 100, probably a float issue

# remove dataframe
rm(relativized_check)

#### (4) Processing Triplicates ####

# filter out for triplicates and not triplicates
triplicates <- data %>% 
  filter(triplicate == "y") %>% 
  mutate(full_sample_name = paste(site_reach, field_date, sample_type))

# add in taxonomy data
triplicates <- left_join(triplicates, unique(qiime2 %>% select(feature_ID, phylum, genus)), by = "feature_ID")

# split out into a list for plotting purposes
triplicates_list <- split(triplicates, triplicates$full_sample_name)

# visually look at differences
for(i in 1:length(triplicates_list)) {
  test_data = triplicates_list[[i]]
  title_label = paste(test_data$site_reach[i],
                      test_data$sample_type[i],
                      test_data$field_date[i], sep = " ")
  
  print(plot_phylum(test_data %>% mutate(vial_ID = as.character(vial_ID)), "picrust2_relative_abundance", "vial_ID") + 
          labs(title = title_label))
  print(plot_cyano(test_data %>% mutate(vial_ID = as.character(vial_ID)), "picrust2_relative_abundance", "vial_ID") + 
          labs(title = title_label))
}

## same findings from QIIME2 processing here- 
# generally samples look good
# samples that appear to have one-odd-one-out:
# RUS-3 8/17/22 NT, SAL-3 9/22/22 NT, SFE-M-1S 7/28/22 TM, SFE-M-3 9/6/22 NT
# will remove odd-one-out:
IDs_to_remove <- c(110, 225, 50, 123)
triplicates_adjusted <- triplicates %>% 
  filter(!vial_ID %in% IDs_to_remove)

# average across triplicates for a single weighted NSTI
triplicates_nsti <- triplicates_adjusted %>% 
  dplyr::group_by(site_reach, site, field_date, sample_type) %>% 
  dplyr::summarize(weighted_NSTI = mean(weighted_NSTI)) %>% 
  ungroup()

# average across triplicates
triplicates_adjusted <- triplicates_adjusted %>% 
  dplyr::group_by(site_reach, site, field_date, sample_type, triplicate, feature_ID, full_sample_name) %>% 
  # need to do mean of relative abundance as each vial has different number of reads!
  dplyr::summarize(relative_abundance_means = mean(picrust2_relative_abundance)) %>%  # may not sum to 100!
  ungroup()

# to ensure they still sum to 100, take the total of relative abundances
total_relative_abundances <- triplicates_adjusted %>% 
  dplyr::group_by(full_sample_name) %>% 
  dplyr::summarize(relative_abundance_totals = sum(relative_abundance_means)) %>% 
  ungroup()

# then re-relativize triplicates
triplicates_adjusted_final <- left_join(triplicates_adjusted, total_relative_abundances, 
                                        by = c("full_sample_name")) %>% 
  mutate(picrust2_relative_abundance = relative_abundance_means / relative_abundance_totals * 100) %>% 
  # join in NSTI data
  left_join(triplicates_nsti, by = c("site_reach", "site", "field_date", "sample_type")) %>% 
  select(site_reach, site, field_date, sample_type, weighted_NSTI, feature_ID, picrust2_relative_abundance)

# check that they all sum to 100
triplicate_relativize_check <- triplicates_adjusted_final %>% 
  mutate(full_sample_name = paste(site_reach, field_date, sample_type)) %>% 
  dplyr::group_by(full_sample_name) %>% 
  dplyr::summarize(total = sum(picrust2_relative_abundance)) %>%  # all summing to 100!
  ungroup()
triplicate_relativize_check$total[which(triplicate_relativize_check$total != 100)]
# all that don't are very close to 100/reading in as 100, probably a float issue

# remove excess files
rm(triplicate_relativize_check, total_relative_abundances, triplicates_adjusted)

#### (5) Putting Data Together ####

# remove triplicates from original data
final_data <- data %>% 
  filter(triplicate == "n") %>% 
  # filter out fake_targets and blanks %>% 
  filter(sample_type!= "blank" & fake_target == "n") %>% 
  select(site_reach, site, field_date, sample_type, weighted_NSTI, feature_ID, picrust2_relative_abundance)

# merge in averaged & re-relativized
final_data <- rbind(final_data, triplicates_adjusted_final)

# add in processed QIIME2 data
final_data <- left_join(final_data, qiime2, by = c("site_reach", "site", "sample_type", "field_date", "feature_ID")) %>% 
  dplyr::rename(qiime2_relative_abundance = relative_abundance)

# confirm there are no samples missing qiime2 data
any(is.na(final_data$qiime2_relative_abundance)) # none!

#### (6) Comparing QIIME2 Outputs versus PICRUSt2-SC Normalized ####

## (a) first, plotting NSTI values
ggplot(data = final_data %>% select(site_reach, site, field_date, sample_type, weighted_NSTI) %>% unique(), 
       aes(y = weighted_NSTI, x = field_date)) +
 geom_point(aes(color = sample_type, shape = site))
# standard pipeline removes values above 2; our average around 0.05-0.1 with TM being lower than TAC and NT

## (b) take random samples to compare

# grab random 12 samples with metadata (that are not fake targets)
set.seed(6)
test_samples <- (metadata %>% filter(fake_target == "n" & sample_type != "blank"))[sample(1:73, 12),] %>%  
  mutate(site_reach_date_type = paste(site_reach, field_date, sample_type))

## (c) plot those samples with QIIME2 outputs on the left and PICRUSt2-SC normalized on the right
for(i in 1:nrow(test_samples)) {
  
  # get data for a single sample
  data = final_data %>% 
    mutate(site_reach_date_type = paste(site_reach, field_date, sample_type)) %>% 
    filter(site_reach_date_type == test_samples$site_reach_date_type[i]) %>% 
    pivot_longer(cols = c("picrust2_relative_abundance", "qiime2_relative_abundance"),
                 values_to = "relative_abundance", names_to = "processing")
  
  # print plots for phylums
  print(plot_phylum(data, "relative_abundance", "site_reach_date_type") +
          facet_wrap(~processing))
  
  # print plots for cyanos
  print(plot_cyano(data, "relative_abundance", "site_reach_date_type") +
          facet_wrap(~processing))
}
# not any crazy differences!

#### (7) Save csv ####

# save full csv
write.csv(final_data, "./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_FINAL.csv",
          row.names = FALSE)

# save versions of TM and TAC with Microcoleus 
# and Anabaena/Cylindrospermum/Trichormus removed respectively

# for microcoleus
TM_data <- final_data %>% 
  filter(genus != "Microcoleus") %>% 
  filter(sample_type == "TM")
write.csv(TM_data, "./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_TM_nomicro.csv", row.names = FALSE)

# for anabaena/cylindrospermum
TAC_data <- final_data %>% 
# note: cylindrospermopsis is also reading a highly close match to anabaena in BLAST so will
# remove that as well
  filter(! genus %in% c("Anabaena","Cylindrospermum","Trichormus",  "Cylindrospermopsis")) %>% 
  filter(sample_type == "TAC")
write.csv(TAC_data, "./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_TAC_noanacyl.csv", row.names = FALSE)
