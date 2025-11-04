#### Comparing microscopy data among rivers
### Jordan Zabrecky
## last edited: 11.03.2025

#### (1) Loading libraries & data ####

# libraries
lapply(c("tidyverse", "plyr"), require, character.only = T)

# read in files (note: two target taxa csvs- one with target taxa included, other with it excluded)
# doing as a list rather than one csv as each has a different # of columns
files <- list.files(path = "./data/morphological/", pattern = ".csv")
data_wide <- lapply(files, function(x) read.csv(paste("./data/morphological/", x, sep = "")))
names(data_wide) <- files

#### (2) Data Transformation

# pivot longer
data <- lapply(data_wide, function(x) x %>% pivot_longer(cols = c(5:ncol(x)), names_to = "taxa",
                                                    values_to = "percent"))

# focus on 2022 data for the three river comparison
data <- lapply(data, function(x) x %>% mutate(year = year(ymd(field_date))) %>% filter(year == 2022) %>% 
                 relocate(year, .before = "field_date"))

# histogram of raw percent (highly right-skew)
lapply(data, function(x) hist(x$percent))

# log-transforming for multivariate analyses since data is highly right-skewed
# minimum nonzero value across all dataframes is...
lapply(data, function(x) min(x$percent[x$percent != 0])) # 0.02079

# beforehand, curious if any taxa has 0 for all samples
lapply(data_wide, function(x) which(colSums(x[,5:ncol(x)]) == 0))

# no aphanothece in microcoleus samples
data$tm_algalonly.csv <- data$tm_algalonly.csv %>% filter(taxa != "aphanothece")
data$tm_algalonly_nomicro.csv <- data$tm_algalonly_nomicro.csv %>% filter(taxa != "aphanothece")
 
# log-transform percent values adding a small amount for zero
data_filtered <- lapply(data_filtered, function(x) x %>% 
                          mutate(percent = case_when(percent == 0 ~ 0.01,
                                                     TRUE ~ percent),
                          log_percent = log(percent)))

# updated histogram
lapply(data, function(x) hist(x$log_percent))

## CONSIDER REMOVING TAXA THAT COMPOSE < than a certain percent!

#### (3) TM Samples ####

## (a) average across site bar plot
TM_bar <- ggplot(data = data$tm_algalonly.csv, aes(x = site, y = percent, color = taxa)) +
  geom_bar()
TM_bar
 
## (b) NMDS

## (c) misc. questions

#### (4) TAC Samples ####

#### (5) NT Samples ####