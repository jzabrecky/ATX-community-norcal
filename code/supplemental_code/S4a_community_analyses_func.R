#### Script of functions used in NMDS (and related) analyses
### Jordan Zabrecky
## last edited: 01.20.2025

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
                 autotransform = TRUE)
  # bind x & y positions to site information
  nmds_final = cbind(as.data.frame(scores(nmds, "sites")), 
                     data %>% select(any_of(c("site_reach", "site", "field_date", 
                                              "sample_type", "TM_atx_category",
                                              "TAC_atx_category", "NT_atx_category",
                                              "TM_atx_detected", "TAC_atx_detected",
                                              "NT_atx_detected")))) %>% 
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
                               data = loadings, linewidth =1, alpha = 0.5, colour = "grey30") +
      geom_text(data = loadings, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
                fontface = "bold", label = rownames(loadings))
  }
  
  return(plot)
}

## (c) runPERMANOVA
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

## (e) run_dbRDA
# run distance-based redundancy analysis
# @param data is wide dataframe of hellinger abundances with environmental covariates on end
# @param start_col is first column of abundance data
# @param end_col is last column of abundance data
# @param mat_atx is the sample type of the community being analyze
# where "tac" matches to TAC atx concentrations, "tm" to TM, and "nt"
# to an average or both or, if only one available, the one available
# TBD about the above though
# note: covariates are hard-coded in ATM and TM vs. TAC is accessible via if/else statements
run_dbRDA <- function(data, start_col, end_col, mat_atx, na.action = "na.fail") {
  
  # run dbRDA
  if(mat_atx == "tac") {
    results = dbrda(data = data, data[,start_col:end_col] ~ 
                      TAC_ATX_all_ug_orgmat_g + DIN_mg_N_L + oPhos_ug_P_L +
                      temp_C + TDC_mg_L + DOC_mg_L + Mg_mg_L, na.action = na.action)
  } else if (mat_atx == "tm") {
    results = dbrda(data = data, data[,start_col:end_col] ~ 
                      TM_ATX_all_ug_orgmat_g + DIN_mg_N_L + oPhos_ug_P_L +
                      temp_C + TDC_mg_L + DOC_mg_L + Mg_mg_L, na.action = na.action)
  } else if (mat_atx == "nt") {
    results = dbrda(data = data, data[,start_col:end_col] ~ 
                      mean_ATX_all_ug_orgmat_g + DIN_mg_N_L + oPhos_ug_P_L +
                      temp_C + TDC_mg_L + DOC_mg_L + Mg_mg_L, 
                      na.action = na.action)
  }
  
  # confirm variance inflation is not high
  test = as.vector(vif.cca(results))
  print(paste("VIF above 30:", any(test > 30), sep = " "))
  
  # is RDA model significant? # MAY MOVE THESE THINGS
  model_sig =  anova.cca(results, parallel=getOption("mc.cores"))
  
  # how much variation explained by RDA? (adjusted r-squared)
  model_rsquared = RsquareAdj(results)
  
  # which variables are significant?
  model_sigvariables = anova.cca(results, step = 1000, by = "term")
  
  # return list of dbRDA object, model significance test, model rsquared, and significant variable test
  final = list(results, model_sig, model_rsquared, model_sigvariables)
  names(final) <- c("object", "model_sig", "rsquared", "sig_variables")
  return(final)
}

plotRDA <- function(dbrda_list, original_data) {
  
  # get (db)RDA object from "run_dbRDA" function
  rda <- rda$object
  
  ggvegan::autoplot(rda) 
  # get % explained by each acces
  perc = round(100*(summary(rda)$cont$importance[2, 1:2]), 2)
  
  # get sample x and y axes
  samples = scores(rda, display="sites", choices=c(1,2), scaling=1)
  
  # get site info and anatoxin groupings from data
  cat_df = original_data %>% 
    select(any_of(c("site", "site_reach", "field_date", "sample_type",
                    "TM_atx_category", "TAC_atx_category", "NT_atx_category",
                    "TM_atx_detected", "TAC_atx_detected", "NT_atx_detected")))
  
  # add in environmental data (in same order!)
  df = cbind(samples, env_df)
  
  # get environmental endpoints from dbRDA
  env <- scores(rda, display="bp", choices=c(1, 2), scaling=1)
  
  plot = ggplot(data = df, aes(x = dbRDA1, y = dbRDA2)) +
            geom_vline(aes(xintercept = 0), color = "gray", linetype = "dotted") +
            geom_hline(aes(yintercept = 0), color = "gray", linetype = "dotted") +
            geom_point(aes(color = site, shape = TM_atx_category)) +
            labs(x = paste("dbRDA1 (", perc[1], "%)", sep = ""), 
                 y = paste("dbRDA2 (", perc[2], "%)", sep = "")) +
            geom_segment(data = env, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), 
                         linewidth = 1, alpha = 0.5, colour = "grey30") +
            geom_text(data = env, aes(x = dbRDA1, y = dbRDA2), colour = "grey30", 
                      fontface = "bold", label = rownames(env))
    plot
    
    # want to fix loadings to be arrows and move labels?
  # add in site color if color argument is "site"
  if(color == "site") {
    plot = plot + scale_color_manual(values = c("SAL" = "#62a7f8",
                                                "SFE-M" = "#416f16",
                                                "RUS" = "#bdb000"))
  }
  plot
  
  # from above- NMDS reference
  plot = ggplot(nmds_data, aes(x = NMDS1, y = NMDS2)) +
    geom_point(aes(color = .data[[color]], shape = .data[[shape]]), size = 4) +
    stat_ellipse(aes(color = .data[[color]]), type = "t", linetype = 2, linewidth = 0.5) +
    labs(subtitle = paste("Stress:", round(stress, 3)),
         x = "NMDS Axis 1",
         y = "NMDS Axis 2")
}
