#### Supplemental script with functions to create barplot figures
### Jordan Zabrecky
## last edited: 12.17.2025

## This script has functions to plot relative abundance bar plots

#### (1) Loading libraries & setting plot themes ####

# libraries
lapply(c("tidyverse"), require, character.only = T)

# set universal plot theme
theme_set(theme_bw() + theme(panel.grid = element_blank(),
                             panel.border = element_rect(fill = NA, color = "black"),
                             legend.position = "right"))

#### (2) Functions for plotting ####

# bar plot function
barplot <- function(data, x, y, fill) {
  plot = ggplot(data = data, aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]])) +
    geom_bar(position = "fill", stat = "identity")
  
  return(plot)
}
