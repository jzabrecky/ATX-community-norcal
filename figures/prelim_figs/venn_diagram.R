
data <- read.csv("./data/molecular/16s_nochimera_rarefied_90_filtered_relativized.csv")

sites <- split(data, data$site)
sample_type <- lapply(sites, function(x) split(x, x$sample_type))
TM <- lapply(sites, function(x) x %>% filter(sample_type == "TM"))
TAC <- lapply(sites, function(x) x %>% filter(sample_type == "TAC"))

library(VennDiagram)

# save to figures
setwd("./figures/prelim_figs/venn_diagrams")

venn.diagram(
  x = list(sites$RUS$feature_ID, sites$SAL$feature_ID, sites$`SFE-M`$feature_ID),
  filename = "rivers_venn_diagram.png",
  category.names = c("Russian" , "Salmon " , "South Fork Eel"),
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#ebdf38", "#62a7f8", "#416f16"),
)

venn.diagram(
  x = list(TM$SAL$feature_ID, TM$`SFE-M`$feature_ID),
  filename = "tm_venn_diagram.png",
  category.names = c("Salmon Microcoleus" , "South Fork Eel Microcoleus"),
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#62a7f8", "#416f16"),
)

venn.diagram(
  x = list(TAC$RUS$feature_ID, TAC$SAL$feature_ID, TAC$`SFE-M`$feature_ID),
  filename = "tac_venn_diagram.png",
  category.names = c("Russian Anabaena", "Salmon Anabaena" , "South Fork Eel Anabaena"),
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#ebdf38", "#62a7f8", "#416f16"),
)


venn.diagram(
  x = list(sample_type$`SFE-M`$NT$feature_ID, sample_type$`SFE-M`$TM$feature_ID, 
           sample_type$`SFE-M`$TAC$feature_ID),
  filename = "sfkeel_venn_diagram.png",
  category.names = c("Non-target", "Microcoleus" , "Anabaena"),
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#416f16", "#6eab35", "#9ddb63"),
)

venn.diagram(
  x = list(sample_type$SAL$NT$feature_ID, sample_type$SAL$TM$feature_ID, 
           sample_type$SAL$TAC$feature_ID),
  filename = "salmon_venn_diagram.png",
  category.names = c("Non-target", "Microcoleus" , "Anabaena"),
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#62a7f8", "#7eb8fc", "#bfddff"),
)

venn.diagram(
  x = list(sample_type$RUS$NT$feature_ID, 
           sample_type$RUS$TAC$feature_ID),
  filename = "russian_venn_diagram.png",
  category.names = c("Non-target", "Anabaena"),
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#ebdf38", "#f7ef7c"),
)
