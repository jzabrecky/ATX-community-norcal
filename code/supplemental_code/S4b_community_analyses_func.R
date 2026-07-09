#### Script of functions used in NMDS (and related) analyses
### Jordan Zabrecky
## last edited: 07.02.2026

# This script hosts functions used to create relative abundance bar plots &
# NMDS and PCoA plots, run PERMANOVAs, and add "event_no" to sampling date which refers 
# to the visit number to the site

#### (1) Loading libraries & set plot theme ####

# libraries
lapply(c("tidyverse", "vegan"), require, character.only = T)

# set universal plot theme
theme_set(theme_bw() + theme(panel.grid = element_blank(),
                             panel.border = element_rect(fill = NA, color = "black"),
                             legend.position = "right"))

#### (2) Functions ####

## (a) barplot
# function to create barplots
# @param data is data in long format
# @param x is x-axis given in "string"
# @param y is y-axis given in "string"
# @param fill is aesthetic grouping for fill of barplots given in "string"
# @param facet_wrap is aesthetic grouping for facet wrap (if called)
barplot <- function(data, x, y, fill, facet_wrap = NA) {
  plot = ggplot(data = data, aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]])) +
    geom_bar(position = "fill", stat = "identity")
  
  if(!is.na(facet_wrap)) {
    plot = plot + facet_wrap(~.data[[facet_wrap]])
  }
  
  return(plot)
}

## (b) getNMDSdata
# creates NMDS data point coordinates and loadings
# @param data is relative abundance data in wide format with environmental/sampling data on left
# @param start_col is index of column for which the abundance data starts
# @param end_col is index of column for which the abundance data ends (default is end of dataframe)
# @param molecular is if the data is molecular (Aka has ASV's or not); this is included because if 
# we ask for loading from it, it will take forever since there are 100+ ASVs even after trimming
getNMDSdata <- function(data, start_col, end_col = NA, ASV = FALSE) {
  
  # if end_col not given, assume end of the dataframe
  if(is.na(end_col)) {
    end_col = ncol(data)
  }
  
  # use vegan to calculate NMDS distances
  nmds = metaMDS(as.matrix(data[,start_col:end_col]),
                 distance = "bray",
                 trymax = 500,
                 autotransform = FALSE)
  # bind x & y positions to site information
  nmds_final = cbind(as.data.frame(scores(nmds, "sites")), 
                     data %>% select(any_of(c("site_reach", "site", "field_date", 
                                              "sample_type", "atx_detected", "atx_group")))) %>% 
    mutate(field_date = ymd(field_date),
           year = year(field_date),
           month = as.character(month(field_date)))
  
  # get loadings for taxa (if not ASV-based!)
  if(ASV == FALSE) {
    vs = envfit(nmds, as.matrix(data[,start_col:end_col]), perm = 999)
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

## (c) makeNMDSplot
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
    geom_point(aes(color = .data[[color]], shape = .data[[shape]]), size = 2) +
    stat_ellipse(aes(color = .data[[color]]), type = "t", linetype = 2, linewidth = 0.5) +
    labs(subtitle = paste("Stress:", round(stress, 3)),
         x = "NMDS Axis 1",
         y = "NMDS Axis 2") +
    theme(plot.subtitle = element_text(hjust=1, vjust=0.5))
  
  # add in site color if color argument is "site"
  if(color == "site") {
    plot = plot + scale_color_manual(values = c("SAL" = "#81bbfc",
                                                "SFE-M" = "#416f16",
                                                "RUS" = "#ab9f00"))
  }
  # add in same shape if shape argument is "site"
  if(shape == "site") {
    plot = plot + scale_shape_manual(values = c("SAL" = 17,
                                                "SFE-M" = 15,
                                                "RUS" = 16))
  }
  
  # add in loadings
  if(loading) {
    
    if(significant) {
      loadings = cbind(loadings, pvalues) %>% 
        filter(pvalue < 0.05)
      
    }
    
    plot = plot + geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
                               data = loadings, linewidth =1, alpha = 0.5, colour = "grey30") +
      geom_text(data = loadings, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
                fontface = "bold", label = rownames(loadings))
  }
  
  return(plot)
}

## (d) getPCoAdata
# creates PCoA data point coordinates and loadings
# @param data is relative abundance data in wide format with environmental/sampling data on left
# @param start_col is index of column for which the abundance data starts
# @param end_col is index of column for which the abundance data ends (default is end of dataframe)
getPCoAdata <- function(data, start_col, end_col = NA) {
  
  # if end_col not given, assume end of the dataframe
  if(is.na(end_col)) {
    end_col = ncol(data)
  }
  
  # get bray-curtis dissimilarity
  bray_dis = vegdist(data[,start_col:end_col], method = "bray")
  
  # run PCoA (use eigenvalues to calculate % explained, and add makes sure our values are all >0)
  pcoa = cmdscale(bray_dis, k = 2, eig = TRUE, add = TRUE)
  
  # add colnames to points for easier plotting
  colnames(pcoa$points) = c("pcoa1", "pcoa2")
  
  # get percent variation explained for each axis
  per_expl = (100 * pcoa$eig / sum(pcoa$eig))[1:2]
  
  # bind x & y positions to site information
  pcoa_final = cbind(as.data.frame(pcoa$points, "sites"), 
                     data %>% select(any_of(c("pcoa1", "pcoa2", "site_reach", "site", "field_date", 
                                              "sample_type", "atx_detected", "atx_group", 
                                              "ATX_all_ug_org_mat")))) %>% 
    mutate(field_date = ymd(field_date),
           year = year(field_date),
           month = as.character(month(field_date))) 
    
  # add ATX interval groupings, if available
  if(any(colnames(pcoa_final) == "ATX_all_ug_org_mat"))  {
    pcoa_final = pcoa_final %>% 
      mutate(atx_group_interval = case_when(ATX_all_ug_org_mat > 10 ~ "greater than 10",
                                            ATX_all_ug_org_mat >= 1 & ATX_all_ug_org_mat <= 10 ~ "between 1 and 10",
                                            ATX_all_ug_org_mat > 0  & ATX_all_ug_org_mat < 1 ~ "less than 1",
                                            ATX_all_ug_org_mat == 0 ~ "none"))
  }
  
  # return list with (1) PCoA plotting data (2) % explained for each axis
  list <- list(pcoa_final, per_expl)
  names(list) <- c("pcoa", "per_expl")
  return(list)
}

## (e) makePCoAplot
# makes PCoA plot
# @param data is list output from function "getPCoAdata"
# @shape graph aesthetic for column that defines shape of point
# @color graph aesthetic for column that defines point color and ellipses drawn
# @stat_ellipse is TRUE/FALSE which dictates if ellipses are drawn or not
makePCoAplot <- function(data, color, shape, stat_ellipse = TRUE) {
  
  # separating out data to be able to easily call each
  pcoa_data = data$pcoa
  pcoa_per_expl = data$per_expl
  
  # make plot
  plot = ggplot(pcoa_data, aes(x = pcoa1, y = pcoa2)) +
    geom_point(aes(color = .data[[color]], shape = .data[[shape]]), size = 2) +
    labs(x = paste("PCoA1 (", round(pcoa_per_expl[1], 2), "%)", sep = ""),
         y =  paste("PCoA2 (", round(pcoa_per_expl[2], 2), "%)", sep = ""),)
  
  if(stat_ellipse == TRUE) {
    plot = plot + 
      stat_ellipse(aes(color = .data[[color]]), type = "t", linetype = 2, linewidth = 0.5)
  }
  
  # add in site color if color argument is "site"
  if(color == "site") {
    plot = plot + scale_color_manual(values = c("SAL" = "#81bbfc",
                                                "SFE-M" = "#416f16",
                                                "RUS" = "#ab9f00"))
  }
  # add in same shape if shape argument is "site"
  if(shape == "site") {
    plot = plot + scale_shape_manual(values = c("SAL" = 17,
                                                "SFE-M" = 15,
                                                "RUS" = 16))
  }
  
  # add in atx color if color argument is "atx_group"
  if(color == "atx_group") {
    plot = plot + scale_color_manual(values = c("none" = "#949494", 
                                                "low" = "#93d152",
                                                "high" = "#3e700a",
                                                "n" = "#949494",
                                                "y" = "#74B030"))
  }
  
  # add in atx color if color argument is "atx_group_interval"
  if(color == "atx_group_interval") {
    plot = plot + scale_color_manual(values = c("none" = "#949494", 
                                                "less than 1" = "#BAE094",
                                                "between 1 and 10" = "#72A644",
                                                "greater than 10" = "#315E07",
                                                "n" = "#949494",
                                                "y" = "#74B030"))
  }
  
  # add in shape if argument is "atx_group"
  if(shape == "atx_group") {
    plot = plot + scale_shape_manual(values = c("none" = 4,
                                                "low" = 17,
                                                "high" = 16))
  }
      
  return(plot)
}

## (f) runPERMANOVA
# runs PERMANOVA test on inputted data
# @param data is relative abundance data in wide format with environmental/sampling data on left
# @param start_col is index of column for which the abundance data starts 
# @param end_col is index of column for which abundance data ends 
# (automatically set to end of givn dataframe unless stated otherwise)
# @param group is explanatory variable for PERMANOVA test given as data$`col_name`
# @strata is an optional argument to include a group-level effect
# @na.action allows option for to remove NAs (automatically set to fail if NAs)
runPERMANOVA <- function(data, start_col, end_col = NA, group, strata = NA,
                         na.action = "na.fail") {
  
  # if "end_col" not given, use the number of columns in dataframe
  if(is.na(end_col)) {
    end_col = ncol(data)
  }
  
  # create distance matrix based on Bray-Curtis distances
  dist_matrix = vegdist(data[,start_col:end_col], method = "bray")
  
  # return PERMANOVA test results
  if(is.na(strata[1])) {
    results = adonis2(dist_matrix ~ group, na.action = na.action)
  } else {
    results = adonis2(dist_matrix ~ group, strata = strata, na.action = na.action)
  }
  
  return(results)
}

## (g) add_event_no
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
