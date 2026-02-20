# testing barplot change over time

# also load un-transformed relative abundances and make it longer for bar plots through time
nt <- read.csv("./data/morphological/nt_algalonly.csv")
tm <- read.csv("./data/morphological/tm_algalonly_nomicro.csv")
tac <- read.csv("./data/morphological/tac_algalonly_noanacylgreenalgae.csv")

# add into list
unaltered_data <- list(nt, tm, tac)
names(unaltered_data) <- c("nt", "tm", "tac")

# finally, pivot longer for those bar plots and select only 2022 data
data_longer <- lapply(unaltered_data, 
                      function(x) x %>% pivot_longer(cols = c(5:ncol(x)), values_to = "percent",
                                                     names_to = "taxa") %>% 
                        filter(year(ymd(field_date)) == 2022))

# add event no
data_longer <- lapply(data_longer, function(x) {
  y = add_event_no(x) %>% 
    mutate(event_no = paste("no_", event_no, sep = ""))
})

# group by site_reach and average
summarized <- lapply(data_longer, function(x) {
  y = x %>% 
    dplyr::group_by(site, event_no, taxa) %>% 
    dplyr::summarize(mean = mean(percent),
                     sd = sd(percent),
                     min = mean - sd,
                     max = mean + sd)
})

# pivot wider? consider removing Salmon also
wider_mean <- lapply(summarized, function(x) {
  # get unique event numbers (in case one is missing)
  event_nos = unique(x$event_no)
  
  y = x %>% 
    select(site, event_no, taxa, mean) %>% 
    pivot_wider(names_from = "event_no", values_from = "mean", values_fill = 0) %>%
    mutate(no_1 = ifelse("no_1" %in% names(.), no_1, 0),) %>% 
    select(site, taxa, no_1, no_2, no_3, no_4, no_5, no_6, no_7)
})

calc_change <- function(data) {
  end = ncol(data)
  for(i in 1:(end - 3)) {
    rowname = paste(colnames(data)[i+3], "_minus_", colnames(data)[i+2], sep = "")
    data[,ncol(data) + 1] = data[,i+3] - data[,i+2] 
    colnames(data)[ncol(data)] = rowname
  }
}

change_data <- lapply(wider_mean, calc_change)

all_change = bind_rows(change_data$nt %>% mutate(sample_type = "NT"), 
                       change_data$tm %>% mutate(sample_type = "TM"),
                       change_data$tac %>% mutate(sample_type = "TAC"))

plotting_data <- all_change %>% 
  select(!c( no_1, no_2, no_3, no_4, no_5, no_6, no_7)) %>% 
  relocate(sample_type, .before = "site") %>% 
  pivot_longer(4:ncol(.), names_to = "time", values_to = "change")

ggplot(data = plotting_data %>% filter(taxa == "e_diatoms" | taxa == "epithemia"), 
       aes(x = time, y = change, fill = sample_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~site)
ggplot(data = plotting_data %>% filter(taxa == "nostoc"), 
       aes(x = time, y = change, fill = sample_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~site)
ggplot(data = plotting_data %>% filter(taxa == "anabaena_and_cylindrospermum"), 
       aes(x = time, y = change, fill = sample_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~site)
ggplot(data = plotting_data %>% filter(taxa == "microcoleus"), 
       aes(x = time, y = change, fill = sample_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~site)

# how about percent change?

calc_perc_change <- function(data) {
  end = ncol(data)
  for(i in 1:(end - 3)) {
    rowname = paste(colnames(data)[i+3], "_minus_", colnames(data)[i+2], sep = "")
    data[,ncol(data) + 1] = ((data[,i+3] - data[,i+2]) / data[,i+2])
    colnames(data)[ncol(data)] = rowname
  }
  return(data)
}

per_change_data <- lapply(wider_mean, calc_perc_change)

per_all_change = bind_rows(per_change_data$nt %>% mutate(sample_type = "NT"), 
                           per_change_data$tm %>% mutate(sample_type = "TM"),
                           per_change_data$tac %>% mutate(sample_type = "TAC"))
# true, that doesn't work because we have a lot of zeros lol
per_plotting_data <- all_change %>% 
  select(!c( no_1, no_2, no_3, no_4, no_5, no_6, no_7)) %>% 
  relocate(sample_type, .before = "site") %>% 
  pivot_longer(4:ncol(.), names_to = "time", values_to = "change")


ggplot(data = per_plotting_data %>% filter(taxa == "e_diatoms" | taxa == "epithemia"), 
       aes(x = time, y = change, fill = sample_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~site)
ggplot(data = per_plotting_data %>% filter(taxa == "nostoc"), 
       aes(x = time, y = change, fill = sample_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~site)
