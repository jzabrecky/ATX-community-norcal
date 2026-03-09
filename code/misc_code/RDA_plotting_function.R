
# OLD RDA code: to get rid of 

## TO-DO finish and fill in description
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
