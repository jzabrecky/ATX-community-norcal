#### NMDS plots

# take from NMDS
source("./code/4a_amongrivers_microscopy.R")

# function to add event no.

NMDS_morphological <- lapply(NMDS_list, function(x) {
  final = x[["nmds"]] %>% 
    add_event_no() %>% 
    group_by(site, event_no) %>% 
    dplyr::summarize(mean1 = mean(NMDS1),
                     mean2 = mean(NMDS2),
                     sd1 = sd(NMDS1),
                     sd2 = sd(NMDS2)) %>% 
    ungroup()
  return(final)
})

# plot!

ggplot(NMDS_morphological[[2]], aes(x = mean1, y = mean2)) +
  geom_errorbar(aes(ymin = mean2 - sd2, ymax = mean2 + sd2, color = site), alpha = 0.2) +
  geom_errorbarh(aes(xmin = mean1 - sd1, xmax = mean1 + sd1, color = site), alpha = 0.2) +
  geom_point(aes(color = site), size = 4) +
  # adding lines in between points
  stat_ellipse(aes(color = .data[[color]]), type = "t", linetype = 2, linewidth = 0.5) +
  labs(subtitle = paste("Stress:", round(stress, 3)),
       x = "NMDS Axis 1",
       y = "NMDS Axis 2")
