#### Further processing of QIIME2 outputs
### Jordan Zabrecky
## last edited: 10.10.2025

## This code reads in the csv of assembled QIIME2 outputs and metadata
## and further processes it by removing reads that are "Mitochondria"
## or "Chloroplasts", removes low confidence (<0.85) reads, removes "fake"
## field target samples, and analyzes/processes blanks and triplicates,
## and saves to a final csv (one with raw reads [triplicates not averaged] 
# and one with relative)

#### (1) Loading libraries and data ####

# loading libraries
lapply(c("tidyverse", "plyr"), require, character.only = T)

# load in data
data <- ldply(list.files(path = "./data/molecular/intermediate_csvs", 
                         pattern = "unfiltered"), function(filename) {
  d <- read.csv(paste("data/molecular/intermediate_csvs/", filename, sep = ""))
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

# functions to quickly plot data
plot_phylum <- function(data) {
  plot <- ggplot(data = data) +
    geom_bar(aes(x = vial_ID, y = abundance, fill = phylum), 
             stat = "identity", position = "fill")
  return(plot)
}
plot_cyano <- function(data) {
  # filter out for phylum cyanobacteria
  cyano_only = data %>%
    filter(phylum == "Cyanobacteria")
  plot <- ggplot(data = cyano_only) +
    geom_bar(aes(x = vial_ID, y = abundance, fill = genus), 
             stat = "identity", position = "fill")
  return(plot)
  
}

# compare non-rarefied versus rarefied data for random samples
for(i in 1:length(indeces)) {
  test_data = data %>% filter(test_index == i)
  title_label = paste(test_samples[i,], collapse = " ")
  
  print(plot_phylum(test_data) + 
          labs(title = title_label) +
          facet_wrap(~file))
  print(plot_cyano(test_data) + 
          labs(title = title_label) +
          facet_wrap(~file))
}
# really no easily visible change in relative abundance!

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
data_confint_test_75 <- data_ver2_filtered %>% 
  filter(confidence > 0.75) %>% 
  mutate(filtered = "> 0.75")
data_confint_test_80 <- data_ver2_filtered %>% 
  filter(confidence > 0.80) %>% 
  mutate(filtered = "> 0.80")
data_confint_test_85 <- data_ver2_filtered %>% 
  filter(confidence > 0.85) %>% 
  mutate(filtered = "> 0.85")
data_confint_test_90 <- data_ver2_filtered %>% 
  filter(confidence > 0.90) %>% 
  mutate(filtered = "> 0.90")
data_confint_test <- rbind(data_confint_test_75,data_confint_test_80, 
                           data_confint_test_85, data_confint_test_90) %>% 
  filter(!is.na(test_index)) %>% 
  filter(file == "16s_nochimera_rarefied_90_unfiltered.csv")

for(i in 1:length(indeces)) {
  test_data = data_confint_test %>% filter(test_index == i)
  title_label <- paste(test_samples[i,], collapse = " ")
  
  print(plot_phylum(test_data) + 
          labs(title = title_label) +
          facet_wrap(~filtered))
  print(plot_cyano(test_data) + 
          labs(title = title_label) +
          facet_wrap(~filtered))
}
# really no easily visible change in relative abundance!

# will just go with 85% because there is a lot of anabaena at 85% 
# & that gets rid of our misc. "uncultured bacterium" reads :)
data_ver3_conf <- data_ver2_filtered %>% 
  filter(confidence > 0.85) %>% 
  mutate(filtered = "> 0.85")

## double-check to see if we lost any samples from doing this
# as a reminder, unrarefied should have 102, 95% rarefied should have 99, and
# 90% rarefied should have 93
data_ver3_conf %>% filter(sample_type != "blank" & fake_target == "n") %>% 
  select(site_reach, sample_type, field_date, file) %>% 
  dplyr::group_by(file) %>% 
  unique() %>%
  dplyr::summarize(count = n())
# lost 1 for non-rarefied but won't be using that anyways, let's continue!

#### (4) Removing "fake" target samples ####

## These were samples that were sometimes taken in absence of macroscopially identifiable
## Microcoleus or Anabaena/Cylindrospermum

# see if we have any Microcoleus in samples
fakes <- data_ver3_conf %>% 
  filter(fake_target == "y")

view(fakes)
# yes, we had presence in some Russian samples and in Salmon september samples even though
# it was not visually apparent (note searched "tychonema")

# we had one sample on 7-14-2022 for TM at SFE-M-4 that was a "maybe"
maybe <- data_ver3_conf %>% 
  filter(fake_target == "maybe")

view(maybe)
# single read of microcoleus/tychonema!

# however, these "fake targets" are probably too confusing for methods and non-target samples will
# show presence of "macroscopically absent" taxa we care about, so focus only on "true" samples
data_ver4_true <- data_ver3_conf %>% 
  filter(fake_target == "n")

# was also very unsure about TAC at SAL-2 on 9/22/22 (only took a little for 16s)
# in the field, so will probably just remove (not as clearly Anabaena as the sample at SAL-3)
# despite some reads being anabaena! (is vial 230)
data_ver4_true <- data_ver4_true %>% 
  filter(vial_ID != 230)

#### (5) Processing blanks ####

# look at blanks
blanks <- data_ver4_true %>% 
  filter(sample_type == "blank") %>% 
  mutate(full_sample_name = paste(site_reach, field_date, sample_type))

view(blanks)
# not a ton of microcoleus or anabaena which is good

# remove blanks
data_ver5_noblanks <- data_ver4_true %>% 
  filter(sample_type != "blank")

#### (5) Relativizing Abundances #####

# first, let's save the raw reads so we can compare ASV's across samples before merging
# triplicates (which requires relativizing)
raw_reads <- split(data_ver5_noblanks, data_ver5_noblanks$file)
names(raw_reads) <- lapply(names(raw_reads), function(x) gsub("_unfiltered.csv*.","", x))
lapply(names(raw_reads), function(x) write.csv(raw_reads[[x]] %>% select(!file), 
                                           paste("./data/molecular/intermediate_csvs/", x, 
                                                 "_filtered_rawreads.csv", sep = ""), 
                                           row.names = FALSE))


# calculate total abundances/reads per vial for each file
total_abundance_per_vial <- data_ver5_noblanks %>% 
  dplyr::group_by(vial_ID, file) %>% 
  dplyr::summarize(total_reads = sum(abundance))

# left join in this data to full dataframe and calculate relative abundance
data_ver6_relativized <- left_join(data_ver5_noblanks, total_abundance_per_vial, by = c("vial_ID", "file")) %>% 
  mutate(relative_abundance = abundance / total_reads * 100) %>% # doing % to avoid really small decimals
  relocate(total_reads, .before = feature_ID) %>% 
  relocate(relative_abundance, .before = feature_ID)

# check to make sure all add up to 100%
relativized_check <- data_ver6_relativized %>% 
  dplyr::group_by(vial_ID, file) %>% 
  dplyr::summarize(total = sum(relative_abundance))

#### (7) Processing Triplicates ####

# filter out for triplicates and not triplicates
triplicates <- data_ver6_relativized %>% 
  filter(triplicate == "y") %>% 
  mutate(full_sample_name = paste(site_reach, field_date, sample_type))

# split out into a list for plotting purposes
triplicates_list <- split(triplicates, triplicates$full_sample_name)

# visually look at differences
for(i in 1:length(triplicates_list)) {
  test_data = triplicates_list[[i]] %>% filter(file == "16s_nochimera_rarefied_90_unfiltered.csv")
  title_label = paste(test_data$site_reach[i],
                      test_data$sample_type[i],
                      test_data$field_date[i], sep = " ")
  
  print(plot_phylum(test_data) + 
          labs(title = title_label))
  print(plot_cyano(test_data) + 
          labs(title = title_label))
}
# generally samples look good
# samples that appear to have one-odd-one-out:
# RUS-3 8/17/22 NT, SAL-3 9/22/22 NT, SFE-M-1S 7/28/22 TM, SFE-M-3 9/6/22 NT
# will remove odd-one-out:
IDs_to_remove <- c(110, 225, 50, 123)
triplicates_adjusted <- triplicates %>% 
  filter(!vial_ID %in% IDs_to_remove)

# average across triplcates
triplicates_adjusted <- triplicates_adjusted %>% 
  dplyr::group_by(site_reach, site, field_date, sample_type, triplicate, feature_ID, taxon_full, 
                  domain, phylum, class, order, family, genus, species, file, full_sample_name) %>% 
  # need to do mean of relative abundance as each vial has different number of reads!
  dplyr::summarize(relative_abundance_means = mean(relative_abundance), # may not sum to 100!
                   confidence = mean(confidence)) %>% 
  mutate(vial_ID = NA)

# to ensure they still sum to 100, take the total of relative abundances
total_relative_abundances <- triplicates_adjusted %>% 
  dplyr::group_by(full_sample_name, file) %>% 
  dplyr::summarize(relative_abundance_totals = sum(relative_abundance_means))

# then re-relativize triplicates
triplicates_adjusted_final <- left_join(triplicates_adjusted, total_relative_abundances, 
                                        by = c("full_sample_name", "file")) %>% 
  mutate(relative_abundance = relative_abundance_means / relative_abundance_totals * 100,
         abundance = relative_abundance) # for plotting function purposes!

# check that they all sum to 100
triplicate_relativize_check <- triplicates_adjusted_final %>% 
  dplyr::group_by(full_sample_name, file) %>% 
  dplyr::summarize(total = sum(relative_abundance)) # all summing to 100!

# plot merged triplicates
triplicates_adjusted_list <- split(triplicates_adjusted_final, triplicates_adjusted_final$full_sample_name) 
for(i in 1:length(triplicates_adjusted_list)) {
  test_data = triplicates_adjusted_list[[i]] %>%  filter(file == "16s_nochimera_rarefied_90_unfiltered.csv")
  title_label = paste(test_data$site_reach[i],
                      test_data$sample_type[i],
                      test_data$field_date[i], sep = " ")
  
  print(plot_phylum(test_data) + 
          labs(title = title_label))
  print(plot_cyano(test_data) + 
          labs(title = title_label))
}

#### (8) Putting final data together

# get data without triplicates and trim columns
data_ver7_final <- data_ver6_relativized %>% 
  filter(triplicate == "n") %>% 
  select(site_reach, site, field_date, sample_type, triplicate, relative_abundance, feature_ID,
         taxon_full, domain, phylum, class, order, family, genus, species, confidence, file)

# adjust column names for triplicates
triplicates_final <- triplicates_adjusted_final %>% 
  relocate(relative_abundance, .before = "feature_ID") %>% 
  select(!c("vial_ID", "full_sample_name", "relative_abundance_means", 
            "relative_abundance_totals", "abundance")) %>% 
  relocate(confidence, .before = file)

# double-check that column names match those of triplicate
eval(colnames(data_ver7_final) == colnames(triplicates_final))

# rbind non-triplicates and triplicates together
data_ver7_final <- rbind(data_ver7_final, triplicates_final) %>% 
  mutate(file = str_remove(file, "_unfiltered.csv"))

# instead of saving several decimals, let's just save up to the minimum required for non-zero values
min((data_ver7_final %>% filter(relative_abundance != 0))[,"relative_abundance"]) # 0.0004859393

# save to the 5th decimal place
data_ver7_final <- data_ver7_final %>% 
  mutate(relative_abundance = round(relative_abundance, 5))

# split into list based on file name
final <- split(data_ver7_final, data_ver7_final$file)

# save all as individual csv's
lapply(names(final), function(x) write.csv(final[[x]] %>% select(!file), 
                                           paste("./data/molecular/intermediate_csvs/", x, 
                                                 "_filtered_relativized.csv", sep = ""), 
                                           row.names = FALSE))
