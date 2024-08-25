#### TM data 2022

library(tidyverse)
# can full join to join everything

# function to make reach a character
change_df <- function(df) {
  df$reach <- as.character(df$reach)
  return(df)
}

##### OLD CODE #####

# salmon 2022
sal_1s_2022 <- read.csv("data/morphological/raw/2022_TM_sal-1s.csv")
sal_2_2022 <- change_df(read.csv("data/morphological/raw/2022_TM_sal-2.csv"))
sal_3_2022 <- change_df(read.csv("data/morphological/raw/2022_TM_sal-3.csv"))

tm_2022 <- full_join(sal_1s_2022, sal_2_2022)
tm_2022 <- full_join(tm_2022, sal_3_2022)

# south fork eel 2022
sfe_1s_2022 <- read.csv("data/morphological/raw/2022_TM_sfe-m-1s.csv")
sfe_3_2022 <- change_df(read.csv("data/morphological/raw/2022_TM_sfe-m-3.csv"))
sfe_4_2022 <- change_df(read.csv("data/morphological/raw/2022_TM_sfe-m-4.csv"))

tm_2022 <- full_join(tm_2022, sfe_1s_2022)
tm_2022 <- full_join(tm_2022, sfe_3_2022)
tm_2022 <- full_join(tm_2022, sfe_4_2022)


## checking 2023 TM data
sal_1s_2023 <-read.csv("data/morphological/raw/2023_TM_sal-1s.csv")
sal_2_2023 <- change_df(read.csv("data/morphological/raw/2023_TM_sal-3.csv"))
sal_3_2023 <- change_df(read.csv("data/morphological/raw/2023_TM_sal-3.csv"))

tm_2023 <- full_join(sal_1s_2023, sal_2_2023)
tm_2023 <- full_join(tm_2023, sal_3_2023)

### USE CODE STARTING HERE #####

library(plyr)
TM_samples <- ldply(list.files(path = "./data/morphological/raw/", pattern = "TM"), function(filename) {
  d <- read.csv(paste("data/morphological/raw/", filename, sep = ""))
  return(d)
})
# need to decide on chroococuss and aphanothene toxin producing coccoids

TM_samples <- replace(TM_samples, is.na(TM_samples), 0)

colnames(TM_samples[8:29])

totals <- rowSums(TM_samples[8:29])

TM_samples$total <- totals

which(TM_samples$total != 100)

TAC_samples <- ldply(list.files(path = "./data/morphological/raw/", pattern = "TAC"), function(filename) {
  d <- read.csv(paste("data/morphological/raw/", filename, sep = ""))
  return(d)
})

TAC_samples <- replace(TAC_samples, is.na(TAC_samples), 0)

colnames(TAC_samples[8:32]) # need to fix chroococcus label

totals_TAC <- rowSums(TAC_samples[8:32])

TAC_samples$total <- totals_TAC

which(TAC_samples$total != 100)

## NT samples

NT_samples <- ldply(list.files(path = "./data/morphological/raw/", pattern = "NT"), function(filename) {
  d <- read.csv(paste("data/morphological/raw/", filename, sep = ""))
  return(d)
})