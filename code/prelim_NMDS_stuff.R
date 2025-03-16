####
# loading packages
library(tidyverse)
library(vegan)

# load data- only care about community and not environmental here
tm <- read.csv("./data/morphological/tm_algalonly_with_covar.csv") %>% 
  dplyr::select(!c(pH:Pheo_ug_g)) %>% 
  dplyr::select(!c(percent_organic_matter, ATX_all_ug_chla_ug, ATX_all_ug_orgmat_g)) %>% 
  mutate(month = month(field_date),
         year = year(field_date))
tac <- read.csv("./data/morphological/tac_algalonly_with_covar.csv") %>% 
  dplyr::select(!c(pH:Pheo_ug_g)) %>% 
  dplyr::select(!c(percent_organic_matter, ATX_all_ug_chla_ug, ATX_all_ug_orgmat_g)) %>% 
  mutate(month = month(field_date),
         year = year(field_date))
nt <- read.csv("./data/morphological/nt_algalonly_with_covar.csv") %>% 
  dplyr::select(!c(pH:TAC_ATX_all_ug_orgmat_g)) %>% 
  mutate(month = month(field_date),
         year = year(field_date))

# function to see which columns have nothing other than zero
check_zeros <- function(df) {
  for(i in 1:ncol(df)) {
    check <- any(which(df[i] > 0))
    if(check == FALSE) {
      return(colnames(df)[i])
    }
  }
}
# should probably move this code block to processing code
check_zeros(tm)
check_zeros(tac)
check_zeros(nt)
# remove aphanothece to avoid errors
tm <- tm %>% 
  select(!aphanothece)

# create matrix
tm_sub<- tm[,c(4:(ncol(tm) - 2))] %>% 
  dplyr::select(!sample_type)
tm_matrix <- as.matrix(tm_sub)

tac_sub <- tac[,c(4:(ncol(tac) - 2))] %>% 
  dplyr::select(!sample_type)
nt_sub <- nt[,c(4:(ncol(nt) - 2))] %>% 
  dplyr::select(!sample_type)
tac_matrix <- as.matrix(tac_sub)
nt_matrix <- as.matrix(nt_sub)

# nmds

nmds_tm <- metaMDS(tm_matrix,
                   distance = "bray",
                   trymax = 500,
                   autotransform = TRUE)

nmds_scores <- as.data.frame(scores(nmds_tm, "sites"))
nmds_data <- cbind(nmds_scores, tm %>% select(site_reach, site, field_date, month, year))
nmds_data$month <- as.character(nmds_data$month)

nmds_tac <- metaMDS(tac_matrix,
                    distance = "bray",
                    trymax = 500,
                    autotransform = TRUE)
nmds_scores_tac <- as.data.frame(scores(nmds_tac, "sites"))
nmds_data_tac <- cbind(nmds_scores_tac, tac %>% select(site_reach, site, field_date, month, year))
nmds_data_tac$month <- as.character(nmds_data_tac$month)

nmds_nt <- metaMDS(nt_matrix,
                    distance = "bray",
                    trymax = 500,
                    autotransform = TRUE)
nmds_scores_nt <- as.data.frame(scores(nmds_nt, "sites"))
nmds_data_nt <- cbind(nmds_scores_nt, nt %>% select(site_reach, site, field_date, month, year))
nmds_data_nt$month <- as.character(nmds_data_nt$month)


# loadings attempt
vf <- envfit(nmds_tm, tm_matrix, perm = 999)
vf_tac <- envfit(nmds_tac, tac_matrix, perm = 999)
vf_nt <- envfit(nmds_nt, nt_matrix, perm = 999)

vf_coord <- as.data.frame(scores(vf, "vectors")) * ordiArrowMul(vf)
vf_coord_tac <- as.data.frame(scores(vf_tac, "vectors")) * ordiArrowMul(vf_tac)
vf_coord_nt <- as.data.frame(scores(vf_nt, "vectors")) * ordiArrowMul(vf_nt)

# only keep those > 0.7
vf_coord <- vf_coord %>% 
  filter(abs(NMDS1) >= 0.75 | abs(NMDS2) >= 0.75)
vf_coord_tac <- vf_coord_tac %>% 
  filter(abs(NMDS1) >= 0.8 | abs(NMDS2) >= 0.8)
vf_coord_nt <- vf_coord_nt %>% 
  filter(abs(NMDS1) >= 0.8 | abs(NMDS2) >= 0.8)
# remove unknown taxa
vf_coord_tac <- vf_coord_tac[-c(10,12),]

# change anabaena & cyl to just anabaena for plotting purposees
rownames(vf_coord)[1] <- "anabaena"

# text specific locations change for graphs
vf_coord_text <- vf_coord
rownames(vf_coord_text) <- c("Anabaena", "Calothrix", "Epithemia", "Geitlerinema",
                             "Green Algae", "Microcoleus", "Other Diatoms", "Coccoid Cyanobacteria")
vf_coord_text["Microcoleus",][1,1] <- (-1.1)
vf_coord_text["Epithemia",][1,1] <- (1.2)
vf_coord_text["Anabaena",][1,1] <- (1.2)

vf_coord_text_tac <- vf_coord_tac
rownames(vf_coord_text_tac) <- c("Calothrix", "Chroococcus", "Epithemia", "Gloeotrichia",
                                 "Homoeothrix", "Leptolyngbya", "Microcoleus", "Nodularia",
                                 "Oscillatoria")
vf_coord_text_tac["Chroococcus",][1,1] <- (1.22)
vf_coord_text_tac["Leptolyngbya",][1,2] <- (-.09)

vf_coord_text_nt <- vf_coord_nt
rownames(vf_coord_text_nt) <- c("Epithemia", "Leptolynbya", "Microcoleus", "Mougeotia",
                                "Other Diatoms", "Oedogonium", "Spirogyra")
vf_coord_text_nt["Spirogyra",][1,2] <- (-0.7)
vf_coord_text_nt["Spirogyra",][1,1] <- (-0.75)
vf_coord_text_nt["Epithemia",][1,1] <- (-0.9)
vf_coord_text_nt["Epithemia",][1,2] <- (-0.75)

# graphs

ggplot(nmds_data, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = site, shape = month), size = 4) +
  stat_ellipse(aes(color = site), type = "t", linetype = 2, size = 0.5) +
  scale_color_manual(values = c("SAL" = "#62a7f8",
                                "SFE-M" = "#416f16",
                                "SFE-SH" = "#a8ff82")) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = vf_coord, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = vf_coord_text, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(vf_coord_text)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "right") +
  labs(title = "NMDS of Morphologically-Identified Microcoleus Mat Community",
       subtitle = paste("Stress:", round(nmds_tm$stress, 3)),
       x = "NMDS Axis 1",
       y = "NMDS Axis 2")

ggplot(nmds_data_tac, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = site, shape = month), size = 4) +
  stat_ellipse(aes(color = site), type = "t", linetype = 2, size = 0.5) +
  scale_color_manual(values = c("SAL" = "#62a7f8",
                                "SFE-M" = "#416f16",
                                "SFE-SH" = "#a8ff82",
                                "RUS" = "#bdb000")) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = vf_coord_tac, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = vf_coord_text_tac, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(vf_coord_text_tac)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "right") +
  labs(title = "NMDS of Morphologically-Identified Anabaena Mat Community",
       subtitle = paste("Stress:", round(nmds_tm$stress, 3)),
       x = "NMDS Axis 1",
       y = "NMDS Axis 2")

ggplot(nmds_data_nt, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = site, shape = month), size = 4) +
  stat_ellipse(aes(color = site), type = "t", linetype = 2, size = 0.5) +
  scale_color_manual(values = c("SAL" = "#62a7f8",
                                "SFE-M" = "#416f16",
                                "RUS" = "#bdb000")) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = vf_coord_nt, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = vf_coord_text_nt, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(vf_coord_text_nt)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "right") +
  labs(title = "NMDS of Morphologically-Identified Periphyton Community",
       subtitle = paste("Stress:", round(nmds_tm$stress, 3)),
       x = "NMDS Axis 1",
       y = "NMDS Axis 2")
