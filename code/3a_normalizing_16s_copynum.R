#### Normalizing relative abundance with predicted 16s gene copy numbers
### Jordan Zabrecky
## last edited: 03.11.2026

# This script pulls in relative abundance based on the predicted 16s gene 
# copy numbers for each ASV obtained via PICRUSt2-SC and compares it to 
# the relative abundance without doing so (QIIME2 outputs obtained and
# processed in step 2 of code)

# Note: This will only be done for rarefied (95% threshold) & processed data

# TO DO: rerun QIIME2 scripts and save 95 rarefied instead and change name
# remove anacyl and microcoleus at end, will need to add "p" 

#### (1) Loading libraries & data ####

# loading libraries
lapply(c("tidyverse"), require, character.only = T)

# functions to quickly plot data
plot_phylum <- function(data, abundance_var) {
  plot <- ggplot(data = data) +
    geom_bar(aes(x = vial_ID, y = .data[[abundance_var]], fill = phylum), 
             stat = "identity", position = "fill")
  return(plot)
}
plot_cyano <- function(data, abundance_var) {
  # filter out for phylum cyanobacteria
  cyano_only = data %>%
    filter(phylum == "Cyanobacteria")
  plot <- ggplot(data = cyano_only) +
    geom_bar(aes(x = vial_ID, y = .data[[abundance_var]], fill = genus), 
             stat = "identity", position = "fill")
  return(plot)
  
}

## (a) getting metadata

# get list of plates in order (will read in 1, 2, then 3)
plate_data <- lapply(list.files("./data/molecular/metadata/",
                                  pattern = "plate"), function(x) { 
                                      y = read.table(paste("./data/molecular/metadata/", x, sep = ""))
                                      colnames(y) <- c("plate_ID", "vial_ID")
                                      return(y)
                                  })

# read in sample metadata & add into plate data
metadata <- read.csv("./data/molecular/metadata/16s_sample_metadata.csv")

## (b) reading in PICRUSt2-SC NSTI values (whereby lower values closer to 0 indicate better match)
picrust_nsti <- lapply(list.files("./data/molecular/picrust2_outputs/nochimera_nonrarefied",
                                  pattern = "nsti"), function(x) 
                                  read_tsv(paste("./data/molecular/picrust2_outputs/nochimera_nonrarefied/", x, sep = "")) %>% 
                                             dplyr::rename(plate_ID = sample))

## (c) reading in PICRUSt2-SC estimated 16s rRNA copy number values
picrust_predcopynum <- lapply(list.files("./data/molecular/picrust2_outputs/nochimera_nonrarefied",
                                        pattern = "seqtab"), function(x) 
                                          read_tsv(paste("./data/molecular/picrust2_outputs/nochimera_nonrarefied/", x, sep = "")))

## (d) reading in QIIME2 outputs from previous processing steps (code step 2)
qiime2 <- read.csv("./data/molecular/16s_nochimera_rarefied_95_endcode2.csv")

#### (2) Join Data Together ####

# merge in plate data with nsti values
merged_data <- list()
for(i in 1:length(plate_data)) {
  merged_data[[i]] <- left_join(plate_data[[i]], picrust_nsti[[i]], by = "plate_ID")
}

# join in abundances normalized by 16s copy number
for(i in 1:length(plate_data)) {
  # pivot picrust2 normalized abundances longer
  temp = picrust_predcopynum[[i]] %>% pivot_longer(cols = (2:ncol(.)), names_to = "plate_ID", values_to = "picrust2_abundance") %>% 
    filter(picrust2_abundance != 0)
  merged_data[[i]] <- left_join(merged_data[[i]], temp, by = "plate_ID")
}

# convert into one dataframe
merged_data = rbind(merged_data[[1]], merged_data[[2]], merged_data[[3]])

# join in metadata
merged_data <- left_join(metadata, merged_data, by = c("vial_ID")) %>% 
  # rename feature_IDs (which is under "normalized" after reading in the ASV)
  dplyr::rename(feature_ID = normalized)

#### TO-DO: checks for missing vials, etc. when we have more complete version of data

#### (3) Relativizing Normalized Abundances ####

# first save normalized abundances that are not relativized and include triplicates
write.csv(merged_data, "./data/molecular/intermediate_csvs/16s_picrustnormalizedabundances_rarefied95_filtered.csv")

# calculate total abundances per vial
total_abundance_per_vial <- merged_data %>% 
  dplyr::group_by(vial_ID) %>% 
  dplyr::summarize(total_abundance = sum(picrust2_abundance)) %>% 
  ungroup()

# left join in this data to full dataframe and calculate relative abundance
merged_data <- left_join(merged_data, total_abundance_per_vial, by = c("vial_ID")) %>% 
  mutate(picrust2_relative_abundance = picrust2_abundance / total_abundance * 100) %>% # doing % to avoid really small decimals
  relocate(total_abundance, .before = feature_ID) %>% 
  relocate(picrust2_relative_abundance, .before = feature_ID)

# check to make sure all add up to 100%
relativized_check <- merged_data %>% 
  dplyr::group_by(vial_ID) %>% 
  dplyr::summarize(total = sum(picrust2_relative_abundance)) %>% 
  ungroup()
relativized_check$total[which(relativized_check$total != 100)] 
# all that don't are very close to 100/reading in as 100, probably a float issue

# remove dataframe
rm(relativized_check)

#### (4) Processing Triplicates ####

# filter out for triplicates and not triplicates
triplicates <- merged_data %>% 
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
  
  print(plot_phylum(test_data, "picrust2_relative_abundance") + 
          labs(title = title_label))
  print(plot_cyano(test_data, "picrust2_relative_abundance") + 
          labs(title = title_label))
}

#### TO-DO: resolve this ####
# seeming to be missing a lot of ASV labels?
# could be difference in rarefied versus non rarefied ATM?
# or processing because the one test just showed up as assigned to bacteria domain

## code from QIIME2 processing 
# generally samples look good
# samples that appear to have one-odd-one-out:
# RUS-3 8/17/22 NT, SAL-3 9/22/22 NT, SFE-M-1S 7/28/22 TM, SFE-M-3 9/6/22 NT
# will remove odd-one-out:
IDs_to_remove <- c(110, 225, 50, 123)
triplicates_adjusted <- triplicates %>% 
  filter(!vial_ID %in% IDs_to_remove)
## end code from QIIME2 processing

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
rm(triplicate_relativize_check, total_relative_abundances, temp, triplicates_adjusted)

#### (5) Putting Data Together ####

# remove triplicates from original data
final_data <- merged_data %>% 
  filter(triplicate == "n") %>% 
  # filter out fake_targets and blanks %>% 
  filter(sample_type!= "blank" & fake_target == "n") %>% 
  select(site_reach, site, field_date, sample_type, weighted_NSTI, feature_ID, picrust2_relative_abundance)

# merge in averaged & re-relativized
final_data <- rbind(final_data, triplicates_adjusted_final)

# add in processed QIIME2 data
final_data <- left_join(final_data, qiime2, by = c("site_reach", "site", "sample_type", "field_date", "feature_ID"))

#### TO-DO: check for NAs after we have final processed version! ####

#### (6) Comparing QIIME2 Outputs versus PICRUSt2-SC Normalized ####

## (a) first, plotting NSTI values
ggplot(data = final_data %>% select(site_reach, site, field_date, sample_type, weighted_NSTI) %>% unique(), 
       aes(y = weighted_NSTI, x = field_date)) +
 geom_point(aes(color = sample_type, shape = site))
# standard pipeline removes values above 2; our average around 0.05-0.1 with TM being lower than TAC and NT

## (b) take random samples to compare

## (c) plot those samples with QIIME2 outputs on the left and PICRUSt2-SC normalized on the right

#### (7) Save csv ####

# save full csv
write.csv(final_data, "./data/molecular/16s_nochimera_rarefied_95_copynum_normalized_FINAL.csv")

# save versions of TM and TAC with Microcoleus and Anabaena/Cylindrospermum/Trichormus removed respectively

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
