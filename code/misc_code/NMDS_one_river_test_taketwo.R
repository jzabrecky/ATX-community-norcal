## take two, plotting all sample types of one river, together
# cannot omit target taxa!!!!

library(tidyverse)

#### morphology ####
tac <- read.csv("./data/morphological/tac_algalonly.csv") %>% 
  filter(year(field_date) == 2022)
tm <- read.csv("./data/morphological/tm_algalonly.csv") %>% 
  filter(year(field_date) == 2022)
nt <- read.csv("./data/morphological/nt_algalonly.csv") %>% 
  filter(year(field_date) == 2022)

# need to match column IDs for NT
colnames(tac)
colnames(tm)
colnames(nt)
setdiff(colnames(nt), colnames(tm))
setdiff(colnames(tm), colnames(nt))
setdiff(colnames(tac), colnames(tm))
nt_modified <- nt %>% 
  mutate(green_algae = ankistrodesmus + cladophora + closterium + coelastrum + cosmarium + 
           desmodesmus_spines + gloeocystis + lacunastrum + mougeotia + oedogonium + 
           oocystis + pediastrum + scenedesmus_no_spines + spirogyra + stauridium + 
           stigeoclonium + tetraedron + ulothrix + unknown_green_algae + zygnema,
         e_diatoms = epithemia,
         non_e_diatoms = non_e_r_diatoms + rhopalodia,
         nodularia = 0, # nodularia in some TAC and TM
         unknown = euglenoid + chantransia) %>%  # treating unkown like other
  select(colnames(tm))

# okay now join all together
morpho <- rbind(tac, tm, nt_modified)
morpho <- split(morpho, morpho$site)

# NMDS plots

# function to add event_no
source("./code/supplemental_code/S4b_community_analyses_func.R")

eelNMDS <- getNMDSdata(morpho$`SFE-M`, 5, ASV = FALSE) # some issue with vectors, so leaving it out

# copying code from Q2 figure script
processing <- function(x) {
  # add even number to NMDS dataframe and calculate mean & sd of samples from three reaches
  final = x[["nmds"]] %>% 
    add_event_no() %>% 
    group_by(site, month, event_no, sample_type) %>% 
    dplyr::summarize(mean1 = mean(NMDS1),
                     mean2 = mean(NMDS2),
                     sd1 = sd(NMDS1),
                     sd2 = sd(NMDS2)) %>% 
    ungroup() %>% 
    arrange(sample_type)
  
  # add in arrows
  final$start_mean1 = c(NA, final$mean1[-length(final$mean1)])
  final$start_mean2 = c(NA, final$mean2[-length(final$mean2)])
  
  # no arrow pointing to first point of each site
  final$start_mean1[which(final$event_no == 1)] = NA
  final$start_mean2[which(final$event_no == 1)] = NA
  
  return(final)
}

NMDSaverageddata <- processing(eelNMDS)
view(NMDSaverageddata)
# no TAC first date, so day 2 needs to have NA
NMDSaverageddata[which(NMDSaverageddata$sample_type == "TAC" & 
                         NMDSaverageddata$event_no == 2), 
                 c("start_mean1", "start_mean2")] <- NA

# copied from plotting code in supplemental code
loadings = eelNMDS$coord
pvalues = as.data.frame(eelNMDS$vs$vectors$pvals)
colnames(pvalues) = "pvalue"
loadings = cbind(loadings, pvalues) %>% 
  filter(pvalue < 0.05)

ggplot(NMDSaverageddata, aes(x = mean1, y = mean2)) +
  geom_errorbar(aes(ymin = mean2 - sd2, ymax = mean2 + sd2, color = sample_type), 
                alpha = 0.2, linewidth = 1.2) +
  geom_errorbarh(aes(xmin = mean1 - sd1, xmax = mean1 + sd1, color = sample_type), 
                 alpha = 0.2, linewidth = 1.2) +
  geom_point(aes(color = sample_type), size = 4) +
  geom_segment(aes(x = start_mean1, y = start_mean2, xend = mean1, yend = mean2),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"), # Adds an arrow head
               color = "grey") +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
                 data = loadings, linewidth =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = loadings, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = rownames(loadings)) +
  ggtitle("South fork Eel")

# see significant loadings

## (b) RUssian

rusNMDS <- getNMDSdata(morpho$`RUS`, 5, ASV = FALSE) # some issue with vectors, so leaving it out
NMDSprocessedrus <- processing(rusNMDS)

view(NMDSprocessedrus)
# no TAC first date, so day 2 needs to have NA
NMDSprocessedrus[which(NMDSprocessedrus$sample_type == "TAC" & 
                         NMDSprocessedrus$event_no == 3), 
                 c("start_mean1", "start_mean2")] <- NA

# copied from plotting code in supplemental code
loadings = rusNMDS$coord
pvalues = as.data.frame(rusNMDS$vs$vectors$pvals)
colnames(pvalues) = "pvalue"
loadings = cbind(loadings, pvalues) %>% 
  filter(pvalue < 0.05)

ggplot(NMDSprocessedrus, aes(x = mean1, y = mean2)) +
  geom_errorbar(aes(ymin = mean2 - sd2, ymax = mean2 + sd2, color = sample_type), 
                alpha = 0.2, linewidth = 1.2) +
  geom_errorbarh(aes(xmin = mean1 - sd1, xmax = mean1 + sd1, color = sample_type), 
                 alpha = 0.2, linewidth = 1.2) +
  geom_point(aes(color = sample_type), size = 4) +
  geom_segment(aes(x = start_mean1, y = start_mean2, xend = mean1, yend = mean2),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"), # Adds an arrow head
               color = "grey") +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = loadings, linewidth =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = loadings, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = rownames(loadings)) +
  ggtitle("Russian")

## (b) Salmon

salNMDS <- getNMDSdata(morpho$`SAL`, 5, ASV = FALSE) # some issue with vectors, so leaving it out
NMDSprocessedsal <- processing(salNMDS)

view(NMDSprocessedsal)
# no TAC first date, so day 2 needs to have NA
NMDSprocessedsal[which(NMDSprocessedsal$sample_type == "TM" & 
                         NMDSprocessedsal$event_no == 2), 
                 c("start_mean1", "start_mean2")] <- NA
NMDSprocessedsal[which(NMDSprocessedsal$sample_type == "TAC"), 
                 c("start_mean1", "start_mean2")] <- NA

# copied from plotting code in supplemental code
loadings = salNMDS$coord
pvalues = as.data.frame(salNMDS$vs$vectors$pvals)
colnames(pvalues) = "pvalue"
loadings = cbind(loadings, pvalues) %>% 
  filter(pvalue < 0.05)

ggplot(NMDSprocessedsal, aes(x = mean1, y = mean2)) +
  geom_errorbar(aes(ymin = mean2 - sd2, ymax = mean2 + sd2, color = sample_type), 
                alpha = 0.2, linewidth = 1.2) +
  geom_errorbarh(aes(xmin = mean1 - sd1, xmax = mean1 + sd1, color = sample_type), 
                 alpha = 0.2, linewidth = 1.2) +
  geom_point(aes(color = sample_type), size = 4) +
  geom_segment(aes(x = start_mean1, y = start_mean2, xend = mean1, yend = mean2),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"), # Adds an arrow head
               color = "grey") +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = loadings, linewidth =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = loadings, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = rownames(loadings)) +
  ggtitle("Salmon")

#### Molecular plots ####

# all sample types here
molec <- read.csv("./data/molecular/16s_nochimera_rarefied_95_FINAL.csv") %>% 
  mutate(field_date = mdy(field_date))

# split by site
molec <- split(molec, molec$site)

# pivot wider and remove columns that equal zero
molec_wide <- lapply(molec, function(x) {
  y = x %>% 
    select(site, site_reach, sample_type, field_date, feature_ID, relative_abundance) %>% 
    pivot_wider(values_from = relative_abundance, names_from = feature_ID, values_fill = 0)
  return(y)
})

# remove columns where ASV is not observed
test = colSums(molec_wide$RUS[,5:ncol(molec_wide$RUS)])
view(test) # not seeming to be any?
which(colSums(molec_wide$SAL[,5:ncol(molec_wide$SAL)]) == 0)
which(colSums(molec_wide$`SFE-M`[,5:ncol(molec_wide$`SFE-M`)]) == 0)

eelNMDS_molec <- getNMDSdata(molec_wide$`SFE-M`, 5, ASV = TRUE) 
NMDSaverageddata <- processing(eelNMDS_molec)
view(NMDSaverageddata)

# no TAC first date, so day 2 needs to have NA
NMDSaverageddata[which(NMDSaverageddata$sample_type == "TAC" & 
                         NMDSaverageddata$event_no == 2), 
                 c("start_mean1", "start_mean2")] <- NA

ggplot(NMDSaverageddata, aes(x = mean1, y = mean2)) +
  geom_errorbar(aes(ymin = mean2 - sd2, ymax = mean2 + sd2, color = sample_type), 
                alpha = 0.2, linewidth = 1.2) +
  geom_errorbarh(aes(xmin = mean1 - sd1, xmax = mean1 + sd1, color = sample_type), 
                 alpha = 0.2, linewidth = 1.2) +
  geom_point(aes(color = sample_type), size = 4) +
  geom_segment(aes(x = start_mean1, y = start_mean2, xend = mean1, yend = mean2),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"), # Adds an arrow head
               color = "grey") +
  ggtitle("South fork Eel")


eelNMDS_molec <- getNMDSdata(molec_wide$`SFE-M`, 5, ASV = TRUE) 
NMDSaverageddata <- processing(eelNMDS_molec)
view(NMDSaverageddata)

# no TAC first date, so day 2 needs to have NA
NMDSaverageddata[which(NMDSaverageddata$sample_type == "TAC" & 
                         NMDSaverageddata$event_no == 2), 
                 c("start_mean1", "start_mean2")] <- NA

ggplot(NMDSaverageddata, aes(x = mean1, y = mean2)) +
  geom_errorbar(aes(ymin = mean2 - sd2, ymax = mean2 + sd2, color = sample_type), 
                alpha = 0.2, linewidth = 1.2) +
  geom_errorbarh(aes(xmin = mean1 - sd1, xmax = mean1 + sd1, color = sample_type), 
                 alpha = 0.2, linewidth = 1.2) +
  geom_point(aes(color = sample_type), size = 4) +
  geom_segment(aes(x = start_mean1, y = start_mean2, xend = mean1, yend = mean2),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"), # Adds an arrow head
               color = "grey") +
  ggtitle("South fork Eel")
# interesting that despite algal communities being quite different, that bacterial has
# a more consistent change throughout the river

rusNMDS_molec <- getNMDSdata(molec_wide$`RUS`, 5, ASV = TRUE) 
NMDSaverageddata <- processing(rusNMDS_molec)
view(NMDSaverageddata)

# no TAC first date, so day 2 needs to have NA
NMDSaverageddata[which(NMDSaverageddata$sample_type == "TAC" & 
                         NMDSaverageddata$event_no == 3), 
                 c("start_mean1", "start_mean2")] <- NA

ggplot(NMDSaverageddata, aes(x = mean1, y = mean2)) +
  geom_errorbar(aes(ymin = mean2 - sd2, ymax = mean2 + sd2, color = sample_type), 
                alpha = 0.2, linewidth = 1.2) +
  geom_errorbarh(aes(xmin = mean1 - sd1, xmax = mean1 + sd1, color = sample_type), 
                 alpha = 0.2, linewidth = 1.2) +
  geom_point(aes(color = sample_type), size = 4) +
  geom_segment(aes(x = start_mean1, y = start_mean2, xend = mean1, yend = mean2),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"), # Adds an arrow head
               color = "grey") +
  ggtitle("Russian")

# maybe the question is how do sample_types differ within a river?

# salmon
salNMDS_molec <- getNMDSdata(molec_wide$`SAL`, 5, ASV = TRUE) 
NMDSaverageddata <- processing(salNMDS_molec)
view(NMDSaverageddata)

# no TAC first date, so day 2 needs to have NA
NMDSaverageddata[which(NMDSaverageddata$sample_type == "TAC"), 
                 c("start_mean1", "start_mean2")] <- NA
NMDSaverageddata[which(NMDSaverageddata$sample_type == "TM" & 
                         NMDSaverageddata$event_no == 2), 
                 c("start_mean1", "start_mean2")] <- NA

ggplot(NMDSaverageddata, aes(x = mean1, y = mean2)) +
  geom_errorbar(aes(ymin = mean2 - sd2, ymax = mean2 + sd2, color = sample_type), 
                alpha = 0.2, linewidth = 1.2) +
  geom_errorbarh(aes(xmin = mean1 - sd1, xmax = mean1 + sd1, color = sample_type), 
                 alpha = 0.2, linewidth = 1.2) +
  geom_point(aes(color = sample_type), size = 4) +
  geom_segment(aes(x = start_mean1, y = start_mean2, xend = mean1, yend = mean2),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"), # Adds an arrow head
               color = "grey") +
  ggtitle("Salmon")

# PERMANOVAs

# obviously the algal will be distinct but what about south fork eel river
runPERMANOVA(molec_wide$`SFE-M`, start_col = 5, group = molec_wide$`SFE-M`$sample_type)
# apparently they are significantly different by sample_type
print(anova(betadisper(vegdist(molec_wide$`SFE-M`[,5:ncol(molec_wide$`SFE-M`)], method = "bray"), 
                       molec_wide$`SFE-M`$sample_type)))
# not significantly different dispersion
runPERMANOVA(molec_wide$`SFE-M`, start_col = 5, group = molec_wide$`SFE-M`$sample_type,
             strata = molec_wide$`SFE-M`$field_date)
# still with field_date starata

runPERMANOVA(molec_wide$`RUS`, start_col = 5, group = molec_wide$`RUS`$sample_type)
runPERMANOVA(molec_wide$`SAL`, start_col = 5, group = molec_wide$`SAL`$sample_type)

# we would expect different periphyton/hosts to have distinct bacterial microbiomes
# curious if functional groups are distinct as well