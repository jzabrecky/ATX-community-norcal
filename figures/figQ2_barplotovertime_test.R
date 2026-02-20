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
data_longer <- lapply(data_longer, add_event_no)

# group by site_reach and average
summarized <- 