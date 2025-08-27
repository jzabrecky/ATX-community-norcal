#### Further processing of QIIME2 outputs
### Jordan Zabrecky
## last edited: 08.27.2025

## This code reads in the csv of assembled QIIME2 outputs and metadata
## and further processes it by removing reads that are "Mitochondria"
## or "Chloroplasts", removes low confidence (<0.85) reads, removes "fake"
## field target samples, and analyzes/processes blanks and triplicates
## and saves to a final csv

# left off: making plotting functions, comparing 80-90% confidence filtering for 
# meeting update powerpoint

#### (1) Loading libraries and data ####

# loading libraries
lapply(c("tidyverse"), require, character.only = T)

# load in data
data <- ldply(list.files(path = "./data/molecular/", pattern = "unfiltered"), function(filename) {
  d <- read.csv(paste("data/molecular/", filename, sep = ""))
  d$file <- filename
  return(d)
})

# decide on random 12 samples to compare during processing steps
set.seed(23)
indeces <- sample(1:73, 12)

# get information for these random samples
test_samples <- data %>% filter(sample_type != "blank" & fake_target == "n" 
       & file == "16s_nochimera_rarefied_90_unfiltered.csv", triplicate == "n") %>% 
  select(site_reach, sample_type, field_date) %>%
  unique()
test_samples <- test_samples[indeces,]
test_samples$test_index <- seq(1:12)

# join into data
data <- left_join(data, test_samples, by = c("site_reach", "sample_type", "field_date"))

# lastly, make plotting functions to easily look at samples as we are processing
plot_phylum <- function(data, facet_wrap, group) {
}
plot_cyano <- function(data, facet_wrap, group) {
  
}

#### (2) Filter out chloroplast and mitochondria assignments ####

data_ver2_filtered <- data %>% 
  filter(!grepl("mitochondria", family, ignore.case =TRUE)) %>% 
  filter(!grepl("chloroplast", family, ignore.case = TRUE)) %>% 
  filter(domain != "Eukaryota") # exclude anything matching to Eurkaryota
# seems like family level removes all cases of mitochondria and chloroplast from the dataset!

##### (3) Removing low confidence reads ####

## dataframe size:
# all: 466395
# >0.7: 466395 (100% of original retained)
# >0.8: 405259 (86.8% of original retained)
# lots of uncultured bacteria at .8 level
# >0.9: 335626 (72.0 % of original retained)

# visual test looking at test samples
data_confint_test_80 <- data_ver2_filtered %>% 
  filter(!is.na(test_index)) %>% 
  filter(confidence > 0.80) %>% 
  mutate(filtered = "> 0.80")
data_confint_test_85 <- data_ver2_filtered %>% 
  filter(!is.na(test_index)) %>% 
  filter(confidence > 0.85) %>% 
  mutate(filtered = "> 0.85")
data_confint_test_90 <- data_ver2_filtered %>% 
  filter(!is.na(test_index)) %>% 
  filter(confidence > 0.90) %>% 
  mutate(filtered = "> 0.90")
data_confint_test_95 <- data_ver2_filtered %>% 
  filter(!is.na(test_index)) %>% 
  filter(confidence > 0.95) %>% 
  mutate(filtered = "> 0.95")
data_confint_test <- rbind(data_confint_test_80, data_confint_test_85, data_confint_test_90,
                           data_confint_test_95) %>% 
  filter(file == "16s_nochimera_rarefied_95_unfiltered.csv")

## MAYBE TURN THESE INTO FUNCTIONS

# will go only through test samples; MAYBE WANT TO CALCULATE RELATIVE ABUNDANCE FIRST!?!
for(i in 1:length(indeces)) {
  test_data = data_confint_test %>% filter(test_index == i)
  title_label <- paste(as.character(test_samples[i,]), sep = " ")
  
  plot <- ggplot(data = test_data) +
    geom_bar(aes(x = vial_ID, y = abundance, fill = phylum), 
             stat = "identity", position = "fill") +
    labs(title = title_label) +
    facet_wrap(~filtered)
  
  print(plot)
}

for(i in 1:length(indeces)) {
  cyano_only = data_confint_test %>% filter(test_index == i) %>%  # filter out for cyanobacteria
    filter(phylum == "Cyanobacteria")
  title_label <- paste(as.character(test_samples[i,]), sep = " ")

  plot <- ggplot(data = cyano_only) +
    geom_bar(aes(x = vial_ID, y = abundance, fill = genus), 
           stat = "identity", position = "fill") +
    labs(title = title_label) +
    facet_wrap(~filtered)
  
  print(plot)
}

# will just go with 85% for now because there is a lot of anabaena at 85% :)

## double-check to see if we lost any samples from doing this
# as a reminder, unrarefied should have X, 95% rarefied should have 95, and
# 90% rarefied should have 93
data_ver3_conf %>% filter(sample_type != "blank" & fake_target == "n") %>% 
  select(site_reach, sample_type, field_date, file) %>% 
  dplyr::group_by(file) %>% 
  unique() %>%
  dplyr::summarize(count = n())
# lost 1 for nonrarefied but won't be using that anyways, let's continue!

#### (4) Removing "fake" target samples ####

## These were samples that were sometimes taken in absence of macroscopially identifiable
## Microcoleus or Anabaena/Cylindrospermum

# see if we have any Microcoleus in samples
fakes <- data_ver3_conf %>% 
  filter(fake_target == "y")

view(fakes)
# yes, we had presence in some Russian samples and in Salmon september samples even though
# it was not visually apparent

# however, this is probably too confusing for methods and non-target samples will
# show presence of "macroscopically absent" taxa we care about, so focus only on "true" samples
data_ver4_true <- data_ver3_conf %>% 
  filter(fake_target == "n")

#### (5) Relativizing Abundances #####

# calculate total abundances/reads per vial
total_abundance_per_vial <- data_ver4_true %>% 
  dplyr::group_by(vial_ID) %>% 
  dplyr::summarize(total_reads = sum(abundance))

# summary of number of reads
mean(total_abundance_per_vial$total_reads) # 87131
min(total_abundance_per_vial$total_reads) # 1945
max(total_abundance_per_vial$total_reads) # 270140

# left join in this data to full dataframe and calculate relative abundance
data_ver5_relativized <- left_join(data_ver4_true, total_abundance_per_vial, by = c("vial_ID")) %>% 
  mutate(relative_abundance = abundance / total_reads) %>% 
  relocate(total_reads, .before = feature_ID) %>% 
  relocate(relative_abundance, .before = feature_ID)

# reads per sample
reads_per_sample <- data_ver5_relativized %>% 
  select(site_reach, field_date, total_reads, sample_type) %>% 
  unique() 
# 3 of the 4 lowest are blanks
# 3 of the 4 highest are NT... and one is a blank...?

# plot to look at relationship with reads 
ggplot(data = reads_per_sample) +
  geom_boxplot(aes(y = total_reads, color = sample_type)) +
  theme_bw()
# can also probably be just dependent on how extraction of sample went?

#### (6) Processing blanks ####

# look at blanks
blanks <- data_ver5_relativized %>% 
  filter(sample_type == "blank")
blank_reads_per_sample <- reads_per_sample %>% 
  filter(sample_type == "blank")
# half are below average number of reads in a sample but seemingly no pattern of what makes
# number of total reads low or high (e.g. site or date)

# indicates our while out in the field sample processing was not perfect
# (which involved a quick bleach of the filtering equipment then rinse)
# however this is to be expected; would have been interesting to check
# blanks after in the more extensive in-lab bleaching process

# remove blanks from data
data_ver6_noblanks <- data_ver5_relativized %>% 
  filter(sample_type != "blank")

#### (7) Cleaning names ####

# look at names at different levels
unique(data_ver6_noblanks$domain)
unique(data_ver6_noblanks$phylum)
unique(data_ver6_noblanks$class)
unique(data_ver6_noblanks[which(data_ver6_noblanks$phylum == "Cyanobacteria"),]$family)
unique(data_ver6_noblanks[which(data_ver6_noblanks$phylum == "Cyanobacteria"),]$order)
unique(data_ver6_noblanks[which(data_ver6_noblanks$phylum == "Cyanobacteria"),]$genus)

# clean accordingly
# have some _ to fix for domain

#### (8) Processing Triplicates ####

# filter out for triplicates and not triplicates
triplicates <- data_ver6_noblanks %>% 
  filter(triplicate == "y") %>% 
  mutate(full_sample_name = paste(site_reach, field_date, sample_type)) %>% 
  mutate(vial_ID = as.character(vial_ID))
  
# split out into a list for plotting purposes
triplicates_list <- split(triplicates, triplicates$full_sample_name)

# visually look at differences across phyla
for(i in 1:length(triplicates_list)) {
  title_label <- triplicates_list[[i]]$full_sample_name[1]
  
  plot <- ggplot(data = triplicates_list[[i]]) +
    geom_bar(aes(x = vial_ID, y = relative_abundance, fill = phylum), stat = "identity") +
    labs(title = title_label)
  
  print(plot)
}
# generally samples look good
# samples that appear to have one-odd-one-out:
# RUS-3 8/17/22 NT, SAL-3 9/22/22 NT, SFE-M-1S 7/28/22 TM, SFE-M-3 9/6/22 NT

# visually look at differences within cyanobacteria
for(i in 1:length(triplicates_list)) {
  
  cyano_only = triplicates_list[[i]] %>% # filter out for cyanobacteria
    filter(phylum == "Cyanobacteria")
  title_label <- triplicates_list[[i]]$full_sample_name[1]
  
  plot <- ggplot(data = cyano_only) +
    geom_bar(aes(x = vial_ID, y = abundance, fill = genus), 
             stat = "identity", position = "fill") +
    labs(title = title_label)
  
  print(plot)
}
# more samples here that appear to have one-odd-one-out:
# RUS-1 9/15/2022 TAC, RUS-3 8/17/2022 NT, SAL-3 NT 9/22/22, SAL-3 9/22/22 TAC,
# SFE-M-1S 7/28/2022 TM, SFE-M-3 7/28/22 TAC, SFE-M-3 9/6/22 NT

# save as 16s_nochimera_processed.csv
#### TO-DO anything else I should do???

# average triplicate groups together