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
nt_raw_atx <- 

# read in selected orthologs data