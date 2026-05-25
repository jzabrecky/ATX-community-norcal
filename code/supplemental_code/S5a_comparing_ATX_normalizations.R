#### Comparing ATX versus predicted gene counts with raw ATX versus OM normalized
### Jordan Zabrecky
## last edited: 05.24.2026

# This code compares results from linear models of ATX versus predicted gene counts
# with both raw ATX values (not normalized) and normalized ATX values (OM normalized)

#### (1) Loading libraries & data ####

# read in ATX data
atx <- read.csv("./data/field_and_lab/atx_w_categorical_groupings.csv")

# also need to get back non-normalized ATX!
raw_atx <- read.csv("./data/field_and_lab/cyano_atx.csv") %>% 
  select(field_date, site_reach, sample_type, ATX_all_ug_g)

# create NT ATX with same method
nt_raw_atx <- raw_atx %>% 
  pivot_wider(names_from = sample_type, values_from = ATX_all_ug_g) %>% 
  mutate(NT_ATX_all_ug_orgmat_g = case_when(is.na(TM) ~ TAC,
                                            is.na(TAC) ~ TM,
                                            TRUE ~ (TM + TAC) / 2)) %>% 
  select(field_date, site_reach, NT_ATX_all_ug_orgmat_g) %>% 
  mutate(sample_type = "NT") %>% 
  relocate(sample_type, .before = "NT_ATX_all_ug_orgmat_g") %>% 
  dplyr::rename(ATX_all_ug_g = NT_ATX_all_ug_orgmat_g)

# join in with raw atx data
raw_atx <- rbind(raw_atx, nt_raw_atx) %>% 
  mutate(log_ATX_all_ug_g = case_when(ATX_all_ug_g == 0 ~ log(0.02),
                                            TRUE ~ log(ATX_all_ug_g)))

# join atx data together
atx <- left_join(atx, raw_atx, by = c("field_date", "site_reach", "sample_type"))

# load data (for selected functional groups)
nt <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select.csv") %>% 
  filter(sample_type == "NT")
tm <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tm_nomicro.csv")
tac <- read.csv("./data/molecular/PICRUSt2_predicted_KO_select_tac_noanacyl.csv")

# put all dataframes into a list
data_select <- list(nt, tm, tac)
names(data_select) <- c("NT", "TM", "TAC")

# need to log predicted genes!
data_select <- lapply(data_select, function(x) x %>% 
                        group_by(field_date, site, site_reach, sample_type, functional_grouping) %>% 
                        # need to also group as multiple KOs are listed for each functional group
                        dplyr::summarize(predicted_gene_abundance = sum(predicted_gene_abundance)) %>% 
                        ungroup() %>% 
                        mutate(log_predicted_gene_abundance = log(predicted_gene_abundance)))

# split
atx <- split(atx, atx$sample_type)

# join w/ atx data
data_select <- lapply(names(data_select), function(x) left_join(data_select[[x]] %>% 
                                                                  mutate(field_date = mdy(field_date)),
                                                                atx[[x]] %>% 
                                                                  mutate(field_date = ymd(field_date)),
                                                                by = c("sample_type", "site_reach", "site",
                                                                       "field_date")))
names(data_select) <- c("NT", "TM", "TAC")

#### (2) Making plots ####

# pivot longer
data_select <- lapply(data_select, function(x) x %>% pivot_longer(cols = c("log_ATX_all_ug_org_mat", "log_ATX_all_ug_g"),
                                                                  values_to = "log_ug_g",
                                                                  names_to = "method"))

eel_select <- rbind(data_select$TM %>% filter(site == "SFE-M"), data_select$TAC %>% filter(site == "SFE-M"))

tm_eel_plot = ggplot(data_select$TM %>% filter(site == "SFE-M"), 
              aes(x = log_predicted_gene_abundance, y = log_ug_g)) +
  geom_smooth(aes(color = method, fill = method), method = "lm") +
  facet_wrap(~functional_grouping, scales = "free") + 
  scale_color_discrete(palette = c("#3f9633", "#6d4275")) +
  scale_fill_discrete(palette = c("#acf0a3", "#c79bcf")) +
  geom_point(aes(shape = method, color = method), size = 2, alpha = 0.7) +
  theme(strip.background = element_blank()) +
  ggtitle("TM South Fork Eel")
tm_eel_plot 

tac_eel_plot = ggplot(data_select$TAC %>% filter(site == "SFE-M"), 
                     aes(x = log_predicted_gene_abundance, y = log_ug_g)) +
  geom_smooth(aes(color = method, fill = method), method = "lm") +
  facet_wrap(~functional_grouping, scales = "free") + 
  scale_color_discrete(palette = c("#3f9633", "#6d4275")) +
  scale_fill_discrete(palette = c("#acf0a3", "#c79bcf")) +
  geom_point(aes(shape = method, color = method), size = 2, alpha = 0.7) +
  theme(strip.background = element_blank()) +
  ggtitle("TAC South Fork Eel")
tac_eel_plot 

nt_eel_plot = ggplot(data_select$NT %>% filter(site == "SFE-M"), 
                      aes(x = log_predicted_gene_abundance, y = log_ug_g)) +
  geom_smooth(aes(color = method, fill = method), method = "lm") +
  facet_wrap(~functional_grouping, scales = "free") + 
  scale_color_discrete(palette = c("#3f9633", "#6d4275")) +
  scale_fill_discrete(palette = c("#acf0a3", "#c79bcf")) +
  geom_point(aes(shape = method, color = method), size = 2, alpha = 0.7) +
  theme(strip.background = element_blank()) +
  ggtitle("NT South Fork Eel")
nt_eel_plot 

tac_rus_plot = ggplot(data_select$TAC %>% filter(site == "RUS"), 
                     aes(x = log_predicted_gene_abundance, y = log_ug_g)) +
  geom_smooth(aes(color = method, fill = method), method = "lm") +
  facet_wrap(~functional_grouping, scales = "free") + 
  scale_color_discrete(palette = c("#3f9633", "#6d4275")) +
  scale_fill_discrete(palette = c("#acf0a3", "#c79bcf")) +
  geom_point(aes(shape = method, color = method), size = 2, alpha = 0.7) +
  theme(strip.background = element_blank()) +
  ggtitle("TAC Russian")
tac_rus_plot 

nt_rus_plot = ggplot(data_select$NT %>% filter(site == "RUS"), 
                     aes(x = log_predicted_gene_abundance, y = log_ug_g)) +
  geom_smooth(aes(color = method, fill = method), method = "lm") +
  facet_wrap(~functional_grouping, scales = "free") + 
  scale_color_discrete(palette = c("#3f9633", "#6d4275")) +
  scale_fill_discrete(palette = c("#acf0a3", "#c79bcf")) +
  geom_point(aes(shape = method, color = method), size = 2, alpha = 0.7) +
  theme(strip.background = element_blank()) +
  ggtitle("NT Russian")
nt_rus_plot 

# Summary: not a huge change; some relationships get stronger while others get weaker
