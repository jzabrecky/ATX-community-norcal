#### TM data 2022

library(tidyverse)
# can full join to join everything

# function to make reach a character
change_df <- function(df) {
  df$reach <- as.character(df$reach)
  return(df)
}

# salmon 2022
sal_1s_2022 <- read.csv("data/microscopy/raw/2022_TM_sal-1s.csv")
sal_2_2022 <- change_df(read.csv("data/microscopy/raw/2022_TM_sal-2.csv"))
sal_3_2022 <- change_df(read.csv("data/microscopy/raw/2022_TM_sal-3.csv"))

tm_2022 <- full_join(sal_1s_2022, sal_2_2022)
tm_2022 <- full_join(tm_2022, sal_3_2022)

# south fork eel 2022
sfe_1s_2022 <- read.csv("data/microscopy/raw/2022_TM_sfe-m-1s.csv")
sfe_3_2022 <- change_df(read.csv("data/microscopy/raw/2022_TM_sfe-m-3.csv"))
sfe_4_2022 <- change_df(read.csv("data/microscopy/raw/2022_TM_sfe-m-4.csv"))

tm_2022 <- full_join(tm_2022, sfe_1s_2022)
tm_2022 <- full_join(tm_2022, sfe_3_2022)
tm_2022 <- full_join(tm_2022, sfe_4_2022)


## checking 2023 TM data
sal_1s_2023 <-read.csv("data/microscopy/raw/2023_TM_sal-1s.csv")
sal_2_2023 <- change_df(read.csv("data/microscopy/raw/2023_TM_sal-3.csv"))
sal_3_2023 <- change_df(read.csv("data/microscopy/raw/2023_TM_sal-3.csv"))

tm_2023 <- full_join(sal_1s_2023, sal_2_2023)
tm_2023 <- full_join(tm_2023, sal_3_2023)


library(plyr)
TM_samples <- ldply(list.files(path = "./data/microscopy/raw/", pattern = "TM"), function(filename) {
  d <- read.csv(paste("data/microscopy/raw/", filename, sep = ""))
  return(d)
})

# NEED TO CHECK THAT THEY SUM TO 100

TM_samples <- replace(TM_samples, is.na(TM_samples), 0)

totals <- rowSums(TM_samples[8:29])

TM_samples$total <- totals

which(TM_samples$total != 100)

TA_samples <- ldply(list.files(path = "./data/microscopy/raw/", pattern = "TA"), function(filename) {
  d <- read.csv(paste("data/microscopy/raw/", filename, sep = ""))
  return(d)
})
# leaving in phormidium unknown instead of moving it into unknown because amounts are >1

TA_samples <- replace(TA_samples, is.na(TA_samples), 0)

colnames(TA_samples[8:28])

totals_TA <- rowSums(TA_samples[8:28])

TA_samples$total <- totals_TA

which(TA_samples$total != 100)