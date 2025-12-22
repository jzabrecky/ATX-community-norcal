#### Venn Diagrams for ASVs shared across sample types
### Jordan Zabrecky
## 12.16.2025

# This code determines which ASVs are shared across sample types and rivers
# and illustrates those findings with Venn Diagrams

# will make VennDiagrams prettier later and maybe move this script to figures folder?
# may also want to include the number of samples per each

#### (1) Loading libraries & data ####

# load libraries
lapply(c("VennDiagram", "tidyverse"), require, character.only = T)

# load data
data <- read.csv("./data/molecular/16s_nochimera_rarefied_95_FINAL.csv")

# split by sample type and river
TM <- split(data %>% filter(sample_type == "TM"), (data %>% filter(sample_type == "TM"))$site)
TAC <- split(data %>% filter(sample_type == "TAC"), (data %>% filter(sample_type == "TAC"))$site)
NT <- split(data %>% filter(sample_type == "NT"), (data %>% filter(sample_type == "NT"))$site)
eel <- split(data %>% filter(site == "SFE-M"), (data %>% filter(site == "SFE-M"))$sample_type)
rus <- split(data %>% filter(site == "RUS"), (data %>% filter(site == "RUS"))$sample_type)
sal <- split(data %>% filter(site == "SAL"), (data %>% filter(site == "SAL"))$sample_type)

#### (2) Making Venn Diagrams ####

# figures save instead of popping up in viewer window, so set working directory
setwd("./figures/venn_diagrams")

## (a) comparing rivers within a sample type

# Microcoleus 
venn.diagram(
  x = list(TM$SAL$feature_ID, TM$`SFE-M`$feature_ID),
  filename = "tm_venn_diagram.png",
  category.names = c("Salmon Microcoleus" , "South Fork Eel Microcoleus"),
  print.mode = c("raw", "percent"),
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#62a7f8", "#416f16"),
)

# Anabaena/Cylindrospermum
venn.diagram(
  x = list(TAC$RUS$feature_ID, TAC$SAL$feature_ID, TAC$`SFE-M`$feature_ID),
  filename = "tac_venn_diagram.png",
  category.names = c("Russian Anabaena", "Salmon Anabaena" , "South Fork Eel Anabaena"),
  print.mode = c("raw", "percent"),
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#ebdf38", "#62a7f8", "#416f16"),
)
# maybe omit salmon single sample?

# Non-Target
venn.diagram(
  x = list(NT$RUS$feature_ID, NT$SAL$feature_ID, NT$`SFE-M`$feature_ID),
  filename = "nt_venn_diagram.png",
  category.names = c("Russian Non-Target", "Salmon Non-Target" , "South Fork Eel Non-Target"),
  print.mode = c("raw", "percent"),
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#ebdf38", "#62a7f8", "#416f16"),
)

## (b) comparing sample types within a river

# South Fork Eel River
venn.diagram(
  x = list(eel$NT$feature_ID, eel$TM$feature_ID, eel$TAC$feature_ID),
  filename = "sfkeel_venn_diagram.png",
  category.names = c("Non-target", "Microcoleus" , "Anabaena"),
  print.mode = c("raw", "percent"),
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#416f16", "#6eab35", "#9ddb63"),
)

# Salmon River
venn.diagram(
  x = list(sal$NT$feature_ID, sal$TM$feature_ID, sal$TAC$feature_ID),
  filename = "salmon_venn_diagram.png",
  category.names = c("Non-target", "Microcoleus" , "Anabaena"),
  print.mode = c("raw", "percent"),
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#62a7f8", "#7eb8fc", "#bfddff"),
)

# Russian River
venn.diagram(
  x = list(rus$NT$feature_ID, rus$TAC$feature_ID),
  filename = "russian_venn_diagram.png",
  category.names = c("Non-target", "Anabaena"),
  print.mode = c("raw", "percent"),
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#ebdf38", "#f7ef7c"),
)


# end of script, set working directory back to original
setwd("../..")
