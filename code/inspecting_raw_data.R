#### microscopy data check and read in

library(tidyverse)
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
# q for taryn on desmodesmus_spines and scenedesmus_no_spines
# also change coccoid label to exclude microcystis, chrocconous, aphanthene, etc.
# other_coccoid_cyanos??? label name?

NT_samples <- replace(NT_samples, is.na(NT_samples), 0)

colnames(NT_samples[8:57]) # some issue with aphanothece but will fix later after talking w/ taryn

totals_NT <- rowSums(NT_samples[8:57])

NT_samples$total <- totals_NT

which(NT_samples$total != 100)
