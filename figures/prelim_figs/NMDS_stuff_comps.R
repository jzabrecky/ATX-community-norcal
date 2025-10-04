####
# loading packages
library(tidyverse)
library(vegan)

set.seed(666)

# load data- only care about community and not environmental here
tm <- read.csv("./data/morphological/tm_algalonly.csv") %>% 
  mutate(month = month(field_date),
         year = year(field_date))
tac <- read.csv("./data/morphological/tac_algalonly.csv") %>% 
  mutate(month = month(field_date),
         year = year(field_date))
nt <- read.csv("./data/morphological/nt_algalonly.csv") %>% 
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
green_algae <- vf_coord_tac["green_algae",]
geitlerinema <- vf_coord_nt["geitlerinema",]
ulothrix <- vf_coord_nt["ulothrix",]
homoeothrix <- vf_coord_nt["homoeothrix",]

# remove unknowns
vf_coord_tac <- vf_coord_tac[-c(18,22),]
vf_coord <- vf_coord[-c(17, 21),]
vf_coord_nt <- vf_coord_nt[-c(41,42),]

# remove ones I don't care about or find redundant
vf_coord <- vf_coord %>% 
  filter(abs(NMDS1) >= 0.3 | abs(NMDS2) >= 0.3)
vf_coord_tac <- vf_coord_tac %>% 
  filter(abs(NMDS1) >= 0.3 | abs(NMDS2) >= 0.3)
vf_coord_nt <- vf_coord_nt %>% 
  filter(abs(NMDS1) >= 0.35 | abs(NMDS2) >= 0.35)

# add back in green algae to TAC and stuff to NT
vf_coord_tac <- rbind(vf_coord_tac, green_algae)
vf_coord_nt <- rbind(vf_coord_nt, geitlerinema)
vf_coord_nt <- rbind(vf_coord_nt, ulothrix)
vf_coord_nt <- rbind(vf_coord_nt, homoeothrix)

# change anabaena & cyl to just anabaena for plotting purposees

# text specific locations change for graphs
vf_coord_text <- vf_coord
rownames(vf_coord_text) <- c("Anabaena", "Calothrix", "Chroococcus", "Epithemia", 
                             "Geitlerinema", "Green Algae", "Leptolyngbya", "Lyngbya",
                             "Microcoleus", "Other Diatoms", "Other Coccoid Cyanobacteria")
vf_coord_text["Chroococcus",][1,2] <- (-.26)

vf_coord_text_tac <- vf_coord_tac
rownames(vf_coord_text_tac) <- c("Aphanothece", "Calothrix", "Chroococcus", "Epithemia", 
                                 "Gloeotrichia", "Homoeothrix", "Leptolyngbya", "Lyngbya", 
                                 "Microcoleus", "Nodularia", "Oscillatoria", "Scytonema",
                                 "Green Algae")
vf_coord_text_tac["Chroococcus",][1,1] <- (0.7)
vf_coord_text_tac["Leptolyngbya",][1,2] <- (-.05)
vf_coord_text_tac["Gloeotrichia",][1,2] <- (.48)
vf_coord_text_tac["Scytonema",][1,2] <- (.52)
vf_coord_text_tac["Scytonema",][1,1] <- (.19)

vf_coord_text_nt <- vf_coord_nt
rownames(vf_coord_text_nt) <- c("Anabaena", "Cladophora", "Coelastrum", "Epithemia", "Leptolynbya", 
                                "Microcoleus", "Mougeotia", "Other Diatoms", "Nostoc", "Oedogonium", 
                                "Rhopalodia", "Spirogyra", "Stauridium", "Geitlerinema", "Ulothrix",
                                "Homoeothrix")
vf_coord_text_nt["Spirogyra",][1,2] <- (-0.52)
vf_coord_text_nt["Epithemia",][1,1] <- (-0.83)
vf_coord_text_nt["Cladophora",][1,1] <- (-0.53)
vf_coord_text_nt["Coelastrum",][1,1] <- (-0.04)
vf_coord_text_nt["Ulothrix",][1,2] <- (0.25)

# graphs

tm <- ggplot(nmds_data, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = site, shape = month), size = 4) +
  stat_ellipse(aes(color = site), type = "t", linetype = 2, size = 0.5) +
  scale_color_manual(values = c("SAL" = "#62a7f8",
                                "SFE-M" = "#416f16",
                                "SFE-SH" = "#a8ff82")) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = vf_coord_text, size =1, alpha = 0.5, colour = "grey30") +
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

tac <- ggplot(nmds_data_tac, aes(x = NMDS1, y = NMDS2)) +
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

nt <- ggplot(nmds_data_nt, aes(x = NMDS1, y = NMDS2)) +
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


library(cowplot)
plot_grid(tm, tac, labels = c('A', 'B'), align = "v")
