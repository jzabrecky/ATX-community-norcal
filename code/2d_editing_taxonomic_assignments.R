#### Editing taxonomic assignments
### Jordan Zabrecky
## last edited: 09.05.2025

## This code takes the processed QIIME outputs (presently, just the 90 rarefied)
## and adjusts taxonomic assignments to make sure they are all clean and correct

#### (1) Loading libraries & data ####

# loading libraries
lapply(c("tidyverse"), require, character.only = T)

# reading in data
original <- read.csv("./data/molecular/16s_nochimera_rarefied_90_filtered.csv") # keep original just in case
data <- original # data dataframe is for altering

#### (2) Changing Tychonema to Microcoleus ####

# our most common sequences that are matching as Tychonema also have a 100%
# match with Microcoleus and past work in northern California uses the Microcoleus
# taxonomy, so let's change that
data[which(grepl("Tychonema", data$genus)),]$order <-  "Oscillatoriales"
data[which(grepl("Tychonema", data$genus)),]$family <-  "Microcoleaceae"
data[which(grepl("Tychonema", data$genus)),]$genus <- "Microcoleus"

#### (3) Checking domain-level ####

# unique values
unique(data$domain) # good

# just curious of ratio of bacteria to archaea in our samples
ggplot(data = data) +
  geom_bar(aes(x = sample_type, y = relative_abundance, fill = domain),
           stat = "identity", position = "fill")
# very much dominated by bacteria!

#### (4) Checking phylum-level ####

# unique values
unique(data$phylum)

# BLASTing the weird names shows that they are mostly random cultured bacteriums
# let's just add them to NA
odd_phylum_names <- c("WPS-2", "SAR324_clade(Marine_group_B)", "MBNT15", "NB1-j", 
                      "RCP2-54", "GAL15", "PAUC34f", "FCPU426", "AncK6", "WS1", 
                      "Rs-K70_termite_group", "Sva0485", "LCP-89")

for(i in 1:length(odd_phylum_names)) {
  data[which(data$phylum == odd_phylum_names[i]),
       c("phylum", "class", "order", "family", "genus", "species")] <- NA
}
unique(data$phylum) # all look good now
# NOTE: can check if these are all real later?

#### (4) Checking orders w/in Cyanobacteria ####

# unique names
unique(data[which(data$phylum == "Cyanobacteria"),]$order)

# names that BLAST to again random uncultured bacteriums move to NA
odd_order_names <- c("SepB-3", "uncultured", "RD011")
for(i in 1:length(odd_order_names)) {
  data[which(data$order == odd_order_names[i]),
       c("order", "family", "genus", "species")] <- NA
}
unique(data[which(data$phylum == "Cyanobacteria"),]$order) # all look good now

# Fr127 has a lot of matches to Leptolyngbyales (order) Leptolyngbyaceae (family)
# (not to genus only to family and order)
data[which(data$order == "Fr127"), "family"] <- "Leptolyngbyaceae"
data[which(data$order == "Fr127"), "order"] <- "Leptolyngbyales"

# "Oxyphotobacteria_Incertae_Sedis" (oxyphotobacteria of uncertain placement)
# weirdly within our dataframe has a lot of cyanobacteria genera such as Calothrix,
# Leptolyngbya, etc. so we will fix the order, family for those when we get to genus

# checking again to make sure all of these are true cyanobacteria orders
unique(data[which(data$phylum == "Cyanobacteria"),]$order)

# Cyanobacteriales is not an order- a lot of what is under it is the family Nostocaceae,
# Microcystaceae, etc. will adjust this when I get to family

# Limnotrichales is not an order in NCBI- all of this is the genus Limnothrix, so change
data[which(data$order == "Limnotrichales"), "family"] <- "Pseudanabaenaceae"
data[which(data$order == "Limnotrichales"), "order"] <- "Pseudanabaenales"

# Thermosynechococcales also not an order, what is listed here as the family is
# actually under the order Acaryochloridales
data[which(data$order == "Thermosynechococcales"), "order"] <- "Acaryochloridales"

# Eurycoccales also not an order, pulling up uncultured bacteriums in BLAST
# so just changing this to NA
data[which(data$order == "Eurycoccales"), c("order", "family", "genus", "species")] <- NA

# lastly, Cyanobacteriia is also not an order
# BLASTing this one has a lot of 100% matches for different cyanobacteria
# since there isn't much of it, just change it to NA
data[which(data$order == "Cyanobacteriia"), c("order", "family", "genus", "species")] <- NA

# final check for now
unique(data[which(data$phylum == "Cyanobacteria"),]$order)
# again, will deal with oxyphotobacteria thingy later

#### (5) Checking families w/in Cyanobacteria ####

# unique names
unique(original[which(data$phylum == "Cyanobacteria"),]$family)

odd_family_names <- c("Unknown_family")

# Cyanobacteriaceae is not a family; genus for these is either Geitlerinema, Geminocystis
# or a bacterium- those can be changed to NA
data[which(data$family == "Cyanobacteriaceae" & grepl("Geitlerinema", data$genus)),
     "order"] <- "Geitlerinematales"
data[which(data$family == "Cyanobacteriaceae" & grepl("Geitlerinema", data$genus)),
     "family"] <- "Geitlerinemataceae"
data[which(data$family == "Cyanobacteriaceae" & grepl("Geitlerinema", data$genus)),
     "order"] <- "Chroococcales"
data[which(data$family == "Cyanobacteriaceae" & grepl("Geminocystis", data$genus)),
     "family"] <- "Geminocystaceae"
data[which(data$family == "Cyanobacteriaceae"),c("order", "family", "genus", "species")] <- NA

# Synechococcales_Incertae_Sedis- all of schizothrix genus which in the Leptolyngbya family
data[which(data$family == "Synechococcales_Incertae_Sedis"), "order"] <- "Leptolyngbyales"
data[which(data$family == "Synechococcales_Incertae_Sedis"), "family"] <- "Schizotrichaceae"

# note Phormidiaceae seems to have been changed to Oscillatoriaceae or Microcoleaceae in NCBI
# depending on genus
data[which(data$family == "Phormidiaceae" & grepl("Trichodesmium", data$genus)), "family"] <- "Microcoleaceae"
data[which(data$family == "Phormidiaceae" & grepl("Kamptonema", data$genus)), "family"] <- "Microcoleaceae"
data[which(data$family == "Phormidiaceae" & grepl("Planktothrix", data$genus)), "family"] <- "Microcoleaceae"
# unknowns can go to Osillatoriaceae since that is what the Phormidiaceae search defaults to?
data[which(data$family == "Phormidiaceae"), "family"] <- "Osillatoriaceae"

# Phormidesmiaceae is also not a family, genus is Phormidesmis which is a new Leptolynbya...
# insert adjustment here


## old script


## process here: look at unique names, readjust naming
# data frame for cleaning names
data_ver7_cleannames <- data_ver6_noblanks

## (a) domain
unique(data_ver6_noblanks$domain) # good

# just curious of ratio of bacteria to archaea in our samples
ggplot(data = data_ver6_noblanks) +
  geom_bar(aes(x = sample_type, y = abundance, fill = domain),
           stat = "identity", position = "fill")
# very much dominated by bacteria!

## (b) phylum
unique(data_ver6_noblanks$phylum)
odd_phylum_names <- c("WPS-2", "SAR324_clade(Marine_group_B)", "MBNT15", "NB1-j", 
                      "RCP2-54", "GAL15", "PAUC34f", "FCPU426", "AncK6", "WS1", 
                      "Rs-K70_termite_group", "Sva0485", "LCP-89")

for(i in 1:length(odd_phylum_names)) {
  data_ver7_cleannames[which(data_ver7_cleannames$phylum == odd_phylum_names[i]),
                       c("phylum", "class", "order", "family", 
                         "genus", "species")] <- NA
}
unique(data_ver7_cleannames$phylum)

# BLASTing most abundant sequence for "weird" labels
# most are coming up as uncultured bacteriums and what not, so let's just shove them to NA
# lapply to go through lists?

# our most common reads matching to Tychnonema also 100% match Microcoleus
# considering the poor resolution here & past work in northern California,
# we will change these to Microcoleus

# also want to test certain taxa to make sure they have correct order

## (c) Order within Cyanobacteria
unique(data_ver7_cleannames[which(data_ver7_cleannames$phylum == "Cyanobacteria"),]$order)
odd_order_names <- c("SepB-3", "uncultured", "RD011")

# convert odd ones to NA
for(i in 1:length(odd_phylum_names)) {
  data_ver7_cleannames[which(data_ver7_cleannames$family == odd_family_names[i]),
                       c("order", "family", "genus", "species")] <- NA
}

# Oxyphotobacteria_Incertae_Sedis refers to oxyphotobacteria of uncertain placement
# a lot of these actually have a true genus (e.g. Calothrix or Phormidium) so adjust
unique(data_ver7_cleannames[which(data_ver7_cleannames$order == "Oxyphotobacteria_Incertae_Sedis"),]$genus)
# will edit those for correct family later when checking genus

# Fr127 has a lot of matches to Leptolyngbyales (order) Leptolyngbyaceae (family)
data_ver7_cleannames[which(data_ver7_cleannames$family == "Fr127"), 
                     "order"] <- "Leptolyngbyales"
data_ver7_cleannames[which(data_ver7_cleannames$family == "Fr127"),
                     "family"] <- "Leptolynbyaceae"

## (c) Order within Cyanobacteria
unique(data_ver7_cleannames[which(data_ver7_cleannames$phylum == "Cyanobacteria"),]$order)


unique(data_ver7_cleannames[which(data_ver7_cleannames$phylum == "Cyanobacteria"),]$family)
odd_family_names <- c("Unknown_Family", "SepB-3", "uncultured", "RD011")

for(i in 1:length(odd_phylum_names)) {
  data_ver7_cleannames[which(data_ver7_cleannames$family == odd_family_names[i]),
                       c("family", "genus", "species")] <- NA
}

# Fr127 has a lot of matches to Leptolyngbyales (order) Leptolyngbyaceae (family)
data_ver7_cleannames[which(data_ver7_cleannames$family == "Fr127"), 
                     "order"] <- "Leptolyngbyales"
data_ver7_cleannames[which(data_ver7_cleannames$family == "Fr127"),
                     "family"] <- "Leptolynbyaceae"

# do all the remaining names make sense?
unique(data_ver7_cleannames[which(data_ver7_cleannames$phylum == "Cyanobacteria"),]$family)
# Cyanobiaceaea- Synechococcus-like species, under a reclassification (see Komarek et al. 2020)
# left as is...
# Limnotrichaceae- also not on NCBI taxonomy, all are genus limnothrix, adjust as such

# clean accordingly
# have some _ to fix for domain