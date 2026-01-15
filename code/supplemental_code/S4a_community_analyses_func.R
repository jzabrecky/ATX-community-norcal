#### Script of functions used in NMDS (and related) analyses
### Jordan Zabrecky
## last edited: 12.17.2025

# This script hosts functions used to create NMDS plots

#### (1) Loading libraries & set plot theme ####

# libraries
lapply(c("tidyverse", "vegan"), require, character.only = T)

# set universal plot theme
theme_set(theme_bw() + theme(panel.grid = element_blank(),
                             panel.border = element_rect(fill = NA, color = "black"),
                             legend.position = "right"))

#### (2) Functions ####

## (a) getNMDSdata
# creates NMDS data point coordinates and loadings
# @param data is relative abundance data in wide format with environmental/sampling data on left
# @param start_col is index of column for which the abundance data starts 
# @param molecular is if the data is molecular (Aka has ASV's or not); this is included because if 
# we ask for loading from it, it will take forever since there are 100+ ASVs even after trimming
getNMDSdata <- function(data, start_col, ASV = FALSE) {
  # use vegan to calculate NMDS distances
  nmds = metaMDS(as.matrix(data[,start_col:ncol(data)]),
                 distance = "bray",
                 trymax = 500,
                 autotransform = TRUE)
  # bind x & y positions to site information
  nmds_final = cbind(as.data.frame(scores(nmds, "sites")), 
                     data %>% select(site_reach, site, field_date, sample_type)) %>% 
    mutate(field_date = ymd(field_date),
           year = year(field_date),
           month = as.character(month(field_date)))
  
  # get loadings for taxa (if not ASV-based!)
  if(ASV == FALSE) {
    vs = envfit(nmds, as.matrix(data[,6:ncol(data)]), perm = 999)
    coord = as.data.frame(scores(vs, "vectors"))
    stress = nmds$stress
    
    # return a named list with both dataframes
    list <- list(nmds_final, vs, coord, stress)
    names(list) = c("nmds", "vs", "coord", "stress")
    return(list)
  } else {
    stress = nmds$stress 
    
    # return a named list with nmds and stress dataframes only
    list <- list(nmds_final, stress)
    names(list) <- c("nmds", "stress")
    return(list)
  }
}

## (b) makeNMDSplot
# makes NMDS plot
# @param data is list output from function "getNMDSdata"
# @param loading is TRUE/FALSE argument for placing loadings on plot
# @param argument is TRUE/FALSE argument for only including significantly loadings (p < 0.05) 
# (previous loading argument must be TRUE)
# @shape graph aesthetic for column that defines shape of point
# @color graph aesthetic for column that defines point color and ellipses drawn
makeNMDSplot <- function(data, loading, significant, color, shape) {
  
  # separating out data to be able to easily call each
  nmds_data = data$nmds
  stress = data$stress
  
  # only specify loadings if those are requested
  if(loading) {
    loadings = data$coord
    pvalues = as.data.frame(data$vs$vectors$pvals)
    colnames(pvalues) = "pvalue"
  }
  
  # make plot
  plot = ggplot(nmds_data, aes(x = NMDS1, y = NMDS2)) +
    geom_point(aes(color = .data[[color]], shape = .data[[shape]]), size = 4) +
    stat_ellipse(aes(color = .data[[color]]), type = "t", linetype = 2, linewidth = 0.5) +
    labs(subtitle = paste("Stress:", round(stress, 3)),
         x = "NMDS Axis 1",
         y = "NMDS Axis 2")
  
  # add in site color if color argument is "site"
  if(color == "site") {
    plot = plot + scale_color_manual(values = c("SAL" = "#62a7f8",
                                                "SFE-M" = "#416f16",
                                                "RUS" = "#bdb000"))
  }
  
  # add in loadings
  if(loading) {
    
    if(significant) {
      loadings = cbind(loadings, pvalues) %>% 
        filter(pvalue < 0.05)
      
    }
    
    plot = plot + geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
                               data = loadings, size =1, alpha = 0.5, colour = "grey30") +
      geom_text(data = loadings, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
                fontface = "bold", label = rownames(loadings))
  }
  
  return(plot)
}

## (c) runPERMANOVA
# runs PERMANOVA test on inputted data
# @param data is relative abundance data in wide format with environmental/sampling data on left
# @param start_col is index of column for which the abundance data starts 
runPERMANOVA <- function(data, start_col, group, strata = NA) {
  # create distance matrix based on Bray-Curtis distances
  dist_matrix = vegdist(data[,start_col:ncol(data)], method = "bray")
  
  # return PERMANOVA test results
  if(is.na(strata[1])) {
    results = adonis2(dist_matrix ~ group)
  } else {
    results = adonis2(dist_matrix ~ group, strata = strata)
  }
  
  return(results)
}

## (d) add_event_no
# add event number for 2022 data
# @param data is wide dataframe with field_date as a column
add_event_no <- function(data) {
  data %>% 
    mutate(field_date = ymd(field_date),
           month = month(field_date)) %>% 
    mutate(event_no = case_when((field_date >= ymd("2022-06-24") & field_date <= ymd("2022-06-29")) ~ 1,
                                (field_date >= ymd("2022-07-06") & field_date <= ymd("2022-07-14")) ~ 2,
                                (field_date >= ymd("2022-07-20") & field_date <= ymd("2022-07-28")) ~ 3,
                                (field_date >= ymd("2022-08-02") & field_date <= ymd("2022-08-10")) ~ 4,
                                (field_date >= ymd("2022-08-17") & field_date <= ymd("2022-08-23")) ~ 5,
                                (field_date >= ymd("2022-09-01") & field_date <= ymd("2022-09-06")) ~ 6,
                                (field_date >= ymd("2022-09-15") & field_date <= ymd("2022-09-22")) ~ 7)) %>% 
    relocate(event_no, .before = "field_date") %>% 
    relocate(month, .before = "field_date")
}
