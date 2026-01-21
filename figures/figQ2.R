#### NMDS plots

# take from NMDS
source("./code/4a_amongrivers_microscopy.R")

# function to add event no.

NMDS_morphological <- lapply(NMDS_list, function(x) {
  # add even number to NMDS dataframe and calculate mean & sd of samples from three reaches
  final = x[["nmds"]] %>% 
    add_event_no() %>% 
    group_by(site, month, event_no) %>% 
    dplyr::summarize(mean1 = mean(NMDS1),
                     mean2 = mean(NMDS2),
                     sd1 = sd(NMDS1),
                     sd2 = sd(NMDS2)) %>% 
    ungroup() %>% 
    arrange(site)
  
  # add in arrows
  final$start_mean1 = c(NA, final$mean1[-length(final$mean1)])
  final$start_mean2 = c(NA, final$mean2[-length(final$mean2)])
  
  # no arrow pointing to first point of each site
  final$start_mean1[which(final$event_no == 1)] = NA
  final$start_mean2[which(final$event_no == 1)] = NA
  
  return(final)
})

# plot!
plot <- list()
for(i in 1:length(NMDS_morphological)) {
  plot[[i]] <- ggplot(NMDS_morphological[[i]], aes(x = mean1, y = mean2)) +
                  geom_errorbar(aes(ymin = mean2 - sd2, ymax = mean2 + sd2, color = site), alpha = 0.2) +
                  geom_errorbarh(aes(xmin = mean1 - sd1, xmax = mean1 + sd1, color = site), alpha = 0.2) +
                  geom_point(aes(color = site, shape = as.factor(month)), size = 4) +
                  geom_segment(aes(x = start_mean1, y = start_mean2, xend = mean1, yend = mean2),
                    arrow = arrow(length = unit(0.2, "cm"), type = "closed"), # Adds an arrow head
                    color = "grey") +
                  scale_shape_manual(values = c(16, 17, 15, 18)) +
                  scale_color_manual(values = c("SAL" = "#62a7f8",
                                                "SFE-M" = "#416f16",
                                                "RUS" = "#bdb000"))
}
lapply(plot, print)

## bacterial

source("./code/4b_amongrivers_16s.R")

NMDS_molecular <- lapply(NMDS_list, function(x) {
  # add even number to NMDS dataframe and calculate mean & sd of samples from three reaches
  final = x[["nmds"]] %>% 
    add_event_no() %>% 
    group_by(site, month, event_no) %>% 
    dplyr::summarize(mean1 = mean(NMDS1),
                     mean2 = mean(NMDS2),
                     sd1 = sd(NMDS1),
                     sd2 = sd(NMDS2)) %>% 
    ungroup() %>% 
    arrange(site)
  
  # add in arrows
  final$start_mean1 = c(NA, final$mean1[-length(final$mean1)])
  final$start_mean2 = c(NA, final$mean2[-length(final$mean2)])
  
  # no arrow pointing to first point of each site
  final$start_mean1[which(final$event_no == 1)] = NA
  final$start_mean2[which(final$event_no == 1)] = NA
  
  return(final)
})

# plot!
plot2 <- list()
for(i in 1:length(NMDS_molecular)) {
  plot2[[i]] <- ggplot(NMDS_molecular[[i]], aes(x = mean1, y = mean2)) +
    geom_errorbar(aes(ymin = mean2 - sd2, ymax = mean2 + sd2, color = site), alpha = 0.2) +
    geom_errorbarh(aes(xmin = mean1 - sd1, xmax = mean1 + sd1, color = site), alpha = 0.2) +
    geom_point(aes(color = site, shape = as.factor(month)), size = 4) +
    geom_segment(aes(x = start_mean1, y = start_mean2, xend = mean1, yend = mean2),
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed"), # Adds an arrow head
                 color = "grey") +
    scale_shape_manual(values = c(16, 17, 15, 18)) +
    scale_color_manual(values = c("SAL" = "#62a7f8",
                                  "SFE-M" = "#416f16",
                                  "RUS" = "#bdb000"))
}
lapply(plot2, print)
