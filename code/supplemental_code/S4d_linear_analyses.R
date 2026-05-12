#### Functions for doing lm's and plotting!
### Jordan Zabrecky
## 05.11.2026

# This script creates functions to make build and plot linear models
# such as comparing diversity metrics with anatoxins and predicted
# functional gene counts with anatoxin concentration

#### (1) Loading libraries ####

# library loading
lapply(c("tidyverse", "ggtext"), 
        require, character.only = T)

#### (2) Functions ####

## (a) linear models
# this function makes linear models and gives output of (1) the linear model
# and (2) a plot
# @param data is data frame
# @param y is string title of response variable vector
# @param x is string title of explantory variable vector
# @rsquared is optional to show r-squared value in upper left corner
linear_model <- function(data, x, y, rsquared = FALSE) {
  # make linear model
  model = lm(data[[x]] ~ data[[y]], data = data)
  
  # make plot
  plot = ggplot(data = data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_smooth(method = "lm") +
    geom_point(aes(shape = site_reach), size = 2)
   
  if(rsquared) {
    plot = plot + annotate(geom = 'text', label = paste("r^2 = ", round(summary(model)$r.squared, 2), 
                                                        ", p = ", round(summary(model)$coefficients[2,4], 2), sep = ""), 
                           x = -Inf, y = Inf, hjust = 0, vjust = 1) 
  }

  # return model object and plot in a list
  final <- list(model, plot)
  names(final) <- c("model", "plot")
  return(final)
}
