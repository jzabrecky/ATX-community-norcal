#### Comparing molecular 16s data across time
### Jordan Zabrecky
## last edited: 01.05.2026

# This code compares 16s rRNA (bacterial assemblage) data from NT, TM, and TAC samples
# across rivers to answer Q2

#### (1) Loading libraries & data ####

# set seed for reproducibility
set.seed(2025)

# libraries
lapply(c("tidyverse", "plyr", "vegan", "cowplot"), require, character.only = T)

# read in files (data transformed in previous script, "4a_amongrivers_microscopy.R")
data <- lapply(list.files(path = "./data/molecular/transformed/", pattern = ".csv"),
               function(x) read.csv(paste("./data/molecular/transformed/", x, sep = "")))
names(data) <- c("nt", "tac", "tm")

# also load un-transformed relative abundances and make it longer for bar plots through time
nt <- read.csv("./data/molecular/16s_nochimera_rarefied_95_FINAL.csv") %>% 
  filter(sample_type == "NT")
tm <- read.csv("./data/molecular/16s_nochimera_rarefied_95_TM_nomicro.csv")
tac <- read.csv("./data/molecular/16s_nochimera_rarefied_95_TAC_noanacyl.csv")

# add into list
unaltered_data <- list(nt, tm, tac)
names(unaltered_data) <- c("nt", "tm", "tac")
data_longer <- lapply(unaltered_data, function(x) x <- x %>% # change date format for next step
                        mutate(field_date = mdy(field_date)))
data <- lapply(data, function(x) x <- x %>% # change date format for next step
                        mutate(field_date = mdy(field_date)))

#### (2) Add  Columns for Sampling Event and Broader Groupings ####

# we have field dates for sampling distinguishing samples, however, let's change it
# to be the the number of sampling event (where x is data)
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

# apply to all dataframes
data <- lapply(data, add_event_no)
data_longer <- lapply(data_longer, add_event_no)

# lastly, add in broader categories to data_longer so we don't have several with barplots
# could do this quickly with forcats package, but since we are lumping, may make sense
# to have more customized categories
phylums <- lapply(data_longer, function(x) x %>%
                    # calculate average abundance for each phylum
                    group_by(phylum) %>% 
                    dplyr::summarize(mean = mean(relative_abundance)) %>% 
                    # put into "Other" category if NA or % is less than #%
                    mutate(phylum_cat = case_when(is.na(phylum) ~ "Other",
                                                  mean < 0.03  ~ "Other",
                                                  TRUE ~ phylum)) %>% 
                    select(phylum, phylum_cat))
classes <- lapply(data_longer, function(x) x %>%
                    # calculate average abundance for each class
                    group_by(class) %>% 
                    dplyr::summarize(mean = mean(relative_abundance)) %>% 
                    # put into "Other" category if NA or % is less than #%
                    mutate(classes_cat = case_when(is.na(class) ~ "Other",
                                                   mean < 0.04  ~ "Other",
                                                   TRUE ~ class)) %>% 
                    select(class, classes_cat))
cyano_order <- lapply(data_longer, function(x) x %>%
                        # filter for cyanobacteria phylum
                        filter(phylum == "Cyanobacteria") %>% 
                        # calculate average abundance for each order
                        group_by(order) %>% 
                        dplyr::summarize(mean = mean(relative_abundance)) %>% 
                        # put into "Other" category if NA or % is less than #%
                        mutate(cyano_order_cat = case_when(is.na(order) ~ "Other",
                                                           mean < 0.035  ~ "Other",
                                                           TRUE ~ order)) %>% 
                        select(order, cyano_order_cat))
cyano_genus <- lapply(data_longer, function(x) x %>%
                        # filter for cyanobacteria phylum
                        filter(phylum == "Cyanobacteria") %>% 
                        # calculate average abundance for each genus
                        group_by(genus) %>% 
                        dplyr::summarize(mean = mean(relative_abundance)) %>% 
                        # put into "Other" category if NA or % is less than #%
                        mutate(cyano_genus_cat = case_when(is.na(genus) ~ "Other",
                                                           mean < 0.11  ~ "Other",
                                                           TRUE ~ genus)) %>% 
                        select(genus, cyano_genus_cat))
# add these categories into our long dataframes
data_long_broader <- data_longer
for(i in 1:length(data_longer)) {
  data_long_broader[[i]] <- left_join(data_long_broader[[i]], phylums[[i]], by = "phylum")
  data_long_broader[[i]] <- left_join(data_long_broader[[i]], classes[[i]], by = "class")
  data_long_broader[[i]] <- left_join(data_long_broader[[i]], cyano_order[[i]], by = "order")
  data_long_broader[[i]] <- left_join(data_long_broader[[i]], cyano_genus[[i]], by = "genus")
}


##### (3) Function for Analyses ####

# load from supplemental script
source("./code/supplemental_code/S4a_community_analyses_func.R")
source("./code/supplemental_code/S4c_barplot_func.R")

#### (4) Barplots through Time ####

# phylum & class first
barplot_phylum_plots <- lapply(data_long_broader, function(x) barplot(x, x = "event_no", y = "relative_abundance", 
                                                              fill = "phylum_cat", facet_wrap = "site"))
barplot_class_plots <- lapply(data_long_broader, function(x) barplot(x, x = "event_no", y = "relative_abundance", 
                                                                 fill = "classes_cat", facet_wrap = "site"))

# titles for plots
titles <- c("Non-Target Samples", 
            "Microcoleus Samples (excluding M)",
            "Anabaena/Cylindrospermum Samples (excluding AC & GA)")

# view plots
for(i in 1:length(barplot_phylum_plots)) {
  print(barplot_phylum_plots[[i]] + labs(title = titles[i]))
  print(barplot_class_plots[[i]] + labs(title = titles[i]))
}

# within cyanobacteria, order and genera
barplot_cyanoorder_plots <- lapply(data_long_broader, 
                                   function(x) barplot(x %>% filter(phylum == "Cyanobacteria"), x = "event_no", y = "relative_abundance", 
                                                       fill = "cyano_order_cat", facet_wrap = "site"))
barplot_cyanogenus_plots <- lapply(data_long_broader, 
                                   function(x) barplot(x %>% filter(phylum == "Cyanobacteria"), x = "event_no", y = "relative_abundance", 
                                                       fill = "cyano_genus_cat", facet_wrap = "site"))

# view plots
for(i in 1:length(barplot_phylum_plots)) {
  print(barplot_cyanoorder_plots[[i]] + labs(title = titles[i]))
  print(barplot_cyanogenus_plots[[i]] + labs(title = titles[i]))
}


#### (5) Q: What are the dominant five taxa at each event for each river? ####

# just looking at phylums
temporal_overview <- lapply(data_longer, function(x) {
  summary = x %>% 
    group_by(site, event_no, class, sample_type) %>% 
    dplyr::summarize(mean_abun = mean(relative_abundance)) %>% 
    # pivot_wider(names_from = event_no, values_from = mean_abun) %>% 
    ungroup() %>% 
    dplyr::group_by(site, event_no) %>% 
    arrange(desc(mean_abun)) %>% 
    slice_head(n = 5) %>% 
    ungroup()
  
  #summary = summary %>% # get table
  #  mutate(ranking = rep(seq(1:5), nrow(summary) / 5)) %>% 
  #  select(!c(mean_abun, broader)) %>% 
  #  pivot_wider(names_from = event_no, values_from = taxa)
  #return(summary)
  
  # split to have separate colors for each plot
  summaries = split(summary, summary$site)
  for(i in 1:length(summaries)) {
    plot = ggplot(summaries[[i]], aes(x = event_no, y = mean_abun, color = class)) +
      geom_line() +
      geom_point() +
      ggtitle(paste(summary$sample_type[1], ": ", names(summaries)[i], sep = ""))
    print(plot)
  }
})
# come back to this

#### (6) PERMANOVA ####
# set column where abundance data starts
start_col <- 8

# separate data out by site
data_river <- lapply(data, function(x) split(x, x$`site`))

## (a) with the strata argument for EVENT NO.
permanovas_event <- lapply(data, function(x) runPERMANOVA(x, start_col, x$`event_no`, strata = x$'site'))
lapply(permanovas_event, print)
# TM ***, TAC *, NT ***, so yes for all

## (b) separated out by river

permanovas_event_NT_sep <- lapply(data_river$nt, function(x) runPERMANOVA(x, start_col, x$`event_no`))
lapply(permanovas_event_NT_sep, print)
# NT significantly different for SAL***, RUS*, and SFE-M**

permanovas_event_TM_sep <- lapply(data_river$tm, function(x) runPERMANOVA(x, start_col, x$`event_no`))
lapply(permanovas_event_TM_sep, print)
# TM significant for SFE-M** but not SAL but only 2 days that are not significantly different (p = 0.1)

permanovas_event_TAC_sep <- lapply(data_river$tac, function(x) runPERMANOVA(x, start_col, x$`event_no`))
lapply(permanovas_event_TAC_sep, print)
# TAC significant for SFE-M* but not for RUS (0.353)

## (c) check dispersion in groups for event no.

for(i in 1:length(data_river$nt)) {
  print(names(data_river$nt)[i])
  print(anova(betadisper(vegdist(data_river$nt[[i]][,start_col:ncol(data_river$nt[[i]])], method = "bray"), 
                         data_river$nt[[i]]$event_no)))
}
# dispersion not significant for RUS or SFE but has * for SAL

for(i in 1:length(data_river$tm)) {
  print(names(data_river$tm)[i])
  print(anova(betadisper(vegdist(data_river$tm[[i]][,start_col:ncol(data_river$tm[[i]])], method = "bray"), 
                         data_river$tm[[i]]$event_no)))
}
# neither have signficant differences in dispersion for TM

# skip the salmon sample which is index 2
for(i in c(1,3)) {
  print(names(data_river$tac)[i])
  print(anova(betadisper(vegdist(data_river$tac[[i]][,start_col:ncol(data_river$tac[[i]])], method = "bray"), 
                         data_river$tac[[i]]$event_no)))
}
# SFE-M has significant differences in dispersion ***, but not RUS
# NEED TO VISUALLY ASSESS THIS