#### Editing taxonomic assignments
### Jordan Zabrecky
## last edited: 02.20.2026

# This code takes the processed QIIME outputs (presently, just the 95% rarefied)
# and adjusts taxonomic assignments to make sure they are all clean and correct
# (e.g., had a lot of Cyanobacteria phylum assignments with a specified, real genus 
# but order "Cyanobacteriales" which does not exist)

# Referencing NCBI Taxonomy database and Strunecky et al. (2023) to do this

#### (1) Loading libraries & data ####

# loading libraries
lapply(c("tidyverse"), require, character.only = T)

# reading in data
original <- read.csv("./data/molecular/intermediate_csvs/16s_nochimera_rarefied_95_filtered_relativized.csv") # keep original just in case
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

# checking all of these in NCBI taxonomy, some notes:

# Chloroflexi (should have ota on the end)
data[which(data$phylum == "Chloroflexi"), "phylum"] <- "Chloroflexiota"
# Campylobacterota instead of Campilobacterota
data[which(data$phylum == "Campilobacterota"), "phylum"] <- "Campylobacterota"
# Deferrisomatota is actually a class in the phylum Desulfobacterota
data[which(data$phylum == "Deferrisomatota"), "phylum"] <- "Desulfobacterota"

# synonyms of some of these groups
# Firmicutes is also known as Bacillota, Modulibacteria as Moduliflexiota, 
# Crenarchaeota as Thermoproteota, Actinobacteriota as Actinomycetota, Desulfobacterota as
# Thermodesulfobacterota, Euryarchaeota as Methanobacteriota, Nanoarchaeota as Nanobdellota,
# 
# Patescibacteria is also known as candidate phyla radiation, highly unresolved group
# Dependentiae is also highly unresolved, matches to Candidatus Babelota


#### (5) Checking class-level ####

# unique names
unique(data$class)

# there is a lot- let's just see if there are any that are common that would
# not get filtered into "other"
classes <- data %>% 
  dplyr::group_by(class) %>% 
  dplyr::summarize(total = sum(relative_abundance)) %>% 
  filter(total > 1) %>%
  ungroup()

# BLASTing the weird names shows that they are mostly random cultured bacteriums
# unassigned passed phylum, let's just add them to NA after class level
odd_class_names <- c("vadinHA49", "OM190", "KD4-96", "uncultured", "Subgroup_22", 
                     "Pla4_lineage", "P9X2b3D02", "BD7-11", "Lineage_IIc", "OLB14", 
                     "bacteriap25", "Pla3_lineage", "TK10", "AT-s3-28", "Subgroup_5",
                     "Subgroup_11", "Subgroup_25", "028H05-P-BN-P5", "ABY1", 
                     "BD2-11_terrestrial_group", "WWE3", "V2072-189E03", "MVP-15", 
                     "JG30-KF-CM66", "4-29-1", "Lineage_IIb", "SJA-28", "MB-A2-108", 
                     "Lineage_IIa", "OC31", "Rs-M47", "CPR2", "Subgroup_18")
for(i in 1:length(odd_class_names)) {
  data[which(data$class == odd_class_names[i]),
       c("class", "order", "family", "genus", "species")] <- NA
}
unique(data$class)

# will wait to look more into confirming the correct, up-to-date spellings if 
# I make a figure using the classes and see what are the main categories showing up :)

#### (6) Checking orders w/in Cyanobacteria ####

# unique names
unique(data[which(data$phylum == "Cyanobacteria"),]$order)

# names that BLAST to again random uncultured bacteriums move to NA
odd_order_names <- c("SepB-3", "RD011")
for(i in 1:length(odd_order_names)) {
  data[which(data$order == odd_order_names[i]),
       c("order", "family", "genus", "species")] <- NA
}

# Fr127 has a lot of matches to Leptolyngbyales (order) Leptolyngbyaceae (family)
# (not to genus only to family and order)
data[which(data$order == "Fr127"), "family"] <- "Leptolyngbyaceae"
data[which(data$order == "Fr127"), "order"] <- "Leptolyngbyales"
data[which(data$genus == "Fr127"), c("genus", "species")] <- NA

# "Oxyphotobacteria_Incertae_Sedis" (oxyphotobacteria of uncertain placement)
# weirdly within our dataframe has a lot of cyanobacteria genera such as Calothrix,
# Leptolyngbya, etc. so we will fix the order, family for those later

# checking again to make sure all of these are true cyanobacteria orders
unique(data[which(data$phylum == "Cyanobacteria"),]$order)

# Cyanobacteriales is not an order- a lot of what is under it is the family Nostocaceae,
# Microcystaceae, etc. will adjust this when I get to family

# Limnotrichales is not an order- all of this is the genus Limnothrix, so change
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
# again, will deal with oxyphotobacteria thingy & cyanobacteriales later

#### (7) Checking families w/in Cyanobacteria ####

# unique names
unique(data[which(data$phylum == "Cyanobacteria"),]$family)

# a lot under "Unknown Family" have defined genus
unique(data[which(data$phylum == "Cyanobacteria" & data$family == "Unknown_Family"),"genus"])
# many are also under the oxyphotobacteria thingy, so will deal with this later

# Cyanobacteriaceae is not a family; genus for these is either Geitlerinema, Geminocystis
# or a bacterium- those can be changed to NA
data[which(data$family == "Cyanobacteriaceae" & grepl("Geitlerinema", data$genus)),
     "order"] <- "Geitlerinematales"
data[which(data$family == "Cyanobacteriaceae" & grepl("Geitlerinema", data$genus)),
     "family"] <- "Geitlerinemataceae"
data[which(data$family == "Cyanobacteriaceae" & grepl("Geminocystis", data$genus)),
     "order"] <- "Chroococcales"
data[which(data$family == "Cyanobacteriaceae" & grepl("Geminocystis", data$genus)),
     "family"] <- "Geminocystaceae"
data[which(data$family == "Cyanobacteriaceae"), c("order", "family", "genus", "species")] <- NA

# Synechococcales_Incertae_Sedis- all of schizothrix genus which in the Leptolyngbya family
data[which(data$family == "Synechococcales_Incertae_Sedis"), "order"] <- "Leptolyngbyales"
data[which(data$family == "Synechococcales_Incertae_Sedis"), "family"] <- "Schizotrichaceae"
# NOTE: ambiguous phylogenetic placement as indicated by Strunecky et al. (2023)

# note Phormidiaceae seems to have been changed to Oscillatoriaceae or Microcoleaceae
# depending on genus
data[which(data$family == "Phormidiaceae" & grepl("Trichodesmium", data$genus)), "family"] <- "Microcoleaceae"
data[which(data$family == "Phormidiaceae" & grepl("Kamptonema", data$genus)), "family"] <- "Microcoleaceae"
data[which(data$family == "Phormidiaceae" & grepl("Planktothrix", data$genus)), "family"] <- "Microcoleaceae"
# unknowns can go to Osillatoriaceae since that is what the Phormidiaceae search defaults to?
data[which(data$family == "Phormidiaceae"), "family"] <- "Osillatoriaceae"

# Phormidesmiaceae is also not a family, genus is Phormidesmis which is a new Leptolynbyales (order)
# NOTE: ambiguous phylogenetic placement as indicated by Strunecky et al. (2023)
data[which(data$family == "Phormidesmiaceae"), "order"] <- "Leptolyngbyales"
data[which(data$family == "Phormidesmiaceae"), "family"] <- "Leptolyngbyaceae"

# Gloeocapsaceae is also not a family, genus Gloeocapsa is in family Chroococcaceae and order Chroococcidiopsidaceae
# NOTE: ambiguous phylogenetic placement as indicated by Strunecky et al. (2023)
# there are also some species listed here that are not in the Gloecapsa genus or different family
data[which(data$family == "Gloeocapsaceae"), "order"] <- "Chroococcidiopsidaceae"
data[which(data$family == "Gloeocapsaceae"), "family"] <- "Chroococcidiopsidales"
data[which(data$genus == "Gleocapsa" & grepl("Limnococcus", data$species)), "genus"] <- "Limnococcus"
data[which(data$genus == "Gleocapsa" & grepl("Limnococcus", data$species)), "family"] <- "Cyanothrichaceae"
data[which(data$genus == "Gleocapsa" & grepl("Chroococcus", data$species)), "genus"] <- "Chroococcus"
data[which(data$genus == "Gleocapsa" & grepl("Chroococcus", data$species)), "family"] <- "Chroococcaceae"

# lastly, the Chroococcales_cyanobacterium should not have identification beyond order according to BLAST/GenBank
data[which(data$genus == "Gleocapsa" & data$species == "Chroococcales_cyanobacterium"), c("family", "genus")] <- NA

# Cyanobiaceaea- Synechococcus-like species, under a reclassification (see Komarek et al. 2020)
# but the genus is Cyanobium which is listed under family Prochlorococcaceae in Strunecky et al. 2023
data[which(data$genu == "Cyanobium"), "family"] <- "Prochlorococcaceae"

# NOTE: Vampirovibriophyceae (under Vampirivibrionia here) 
# is the only other class within cyanobacteria phylum (Strunecky et al. 2023)
# naming conventions are unclear, but likely low abundance, so we will just leave as is

# final check
unique(data[which(data$phylum == "Cyanobacteria"),]$family)
# again, will finish dealing with unknown family later

#### (8) Checking genera w/in Cyanobacteria ####

# unique names
unique(data[which(data$phylum == "Cyanobacteria"),]$genus)

# will assign SU2 symbiont group as N-fixing Diazoplasts (according to Marks et al. 2025)
data[which(data$genus == "SU2_symbiont_group"), c("order", "family", "genus")] <- "Endosymbiotic Diazoplast"

# only two CENA518 reads are Leptolyngbya species
data[which(data$genus == "CENA518"), "order"] <- "Leptolyngbyales"
data[which(data$genus == "CENA518"), "family"] <- "Leptolyngbyaceae"
data[which(data$genus == "CENA518"), "genus"] <- "Leptolyngbya"

# remove uncultured genus to NA sans one "uncultured_Oscillatoriales" that should be in Oscillatoriales
data[which(data$genus == "uncultured" & grepl("uncultured_Oscillatoria", data$species)), "order"] <- "Oscillatoriales"
data[which(data$genus == "uncultured" & grepl("uncultured_Oscillatoria", data$species)), "family"] <- "Oscillatoriaceae"
data[which(data$genus == "uncultured" & grepl("uncultured_Oscillatoria", data$species)), "genus"] <- "Oscillatoria"
data[which(data$genus == "uncultured"), "genus"] <- NA

# LB3-76 does not really match to anything in BLAST
data[which(data$genus == "LB3-76"), c("order", "family", "genus")] <- NA

# Unknown Family?
unique(data[which(data$phylum == "Cyanobacteria" & data$genus == "Unknown_Family"),])
# not seeming to match to anything clearly
data[which(data$phylum == "Cyanobacteria" & data$genus == "Unknown_Family"), c("order", "family", "genus")] <- NA

# lots of the genus names have "_[numbers]" following that can be removed
data[which(data$phylum == "Cyanobacteria"), "genus"] <- str_replace(data[which(data$phylum == "Cyanobacteria"), "genus"], "_.*", "")

# look again at names
unique(data[which(data$phylum == "Cyanobacteria"),]$genus)
# checking to make sure all of these are legit in Strunecky et al. (2023)

# Sericytochromatia not seeming to be anything

# Nostocaceae is a family
unique(data[which(data$phylum == "Cyanobacteria" & data$genus == "Nostocaceae"),])
# they are Calothrix species that are matching to a Calothrix genus, etc. on BLAST
data[which(data$genus == "Nostocaceae" & grepl("Calothrix", data$species)), "order"] <- "Nostocales"
data[which(data$genus == "Nostocaceae" & grepl("Calothrix", data$species)), "family"] <- "Rivulariaceae"
data[which(data$genus == "Nostocaceae" & grepl("Calothrix", data$species)), "genus"] <- "Calothrix"

# candidatus is not a  genus
data[which(data$genus == "Candidatus"), "genus"] <- NA
       
# Leptolyngbyaceae is also a family; only a single read, just removing genus
data[which(data$genus == "Leptolyngbyaceae"), "genus"] <- NA

# Caenarcaniphilales considered as an order, not genus (highly unresolved)
data[which(data$genus == "Caenarcaniphilales"), c("family", "genus")] <- NA

# Obscuribacteraceae considered only at family
data[which(data$genus == "Obscuribacteraceae"), "genus"] <- NA

# Vampirovibrionaceae is also just a family
data[which(data$genus == "Vampirovibrionaceae"), "genus"] <- NA

# also, Vampirovibrionales is not a genus but an order
data[which(data$genus == "Vampirovibrionales"), c("genus")] <- NA

# other notes:
# Sericytochromatia is an unranked clade, going to leave as is but not in Strunecky et al. (2023)
# Gleocapsa not in NCBI but in Phycokey so keeping it

#### (9) Checking consistency across genera ####

## lastly, need to check that the same genus have the same family, order, etc.
## this will also deal with the order "Cyanobacteriales" and the oxyphotobacteria thingy
names <- unique(data[which(data$phylum == "Cyanobacteria"), c("order", "family", "genus")])
view(names) # going through genus in this dataframe

# Aliterella family is wrong and order is the "Cyanobacteriales"
data[which(data$genus == "Aliterella"), "family"] <- "Chroococcidiopsidaceae"
data[which(data$genus == "Aliterella"), "order"] <- "Chroococcidiopsidales"

# Anabaena & Aphanizomenon order is wrong ("Cyanobacteriales")
data[which(data$genus == "Anabaena" | data$genus == "Aphanizomenon"), "order"] <- "Nostocales"

# Aphanizomenon also in wrong family
data[which(data$genus == "Aphanizomenon"), "family"] <- "Aphanizomenonaceae"

# two incorrect labels for Calothrix
data[which(data$genus == "Calothrix"), "family"] <- "Rivulariaceae"
data[which(data$genus == "Calothrix"), "order"] <- "Nostocales"

# Chamaesiphon also has incorrect "Cyanobacteriales" order
data[which(data$genus == "Chamaesiphon"), "order"] <- "Gomontiellales"
# also sometimes listed under wrong family
data[which(data$genus == "Chamaesiphon"), "family"] <- "Chamaesiphonaceae"

# Chroococcidiopsis and Chroococcopsis also have wrong family and incorrect "Cyanobacteriales" order
data[which(data$genus == "Chroococcidiopsis"), "family"] <- "Chroococcidiopsidaceae"
data[which(data$genus == "Chroococcidiopsis"), "order"] <- "Chroococcidiopsidales"
data[which(data$genus == "Chroococcopsis"), "family"] <- "Hyellaceae"
data[which(data$genus == "Chroococcopsis"), "order"] <- "Pleurocapsales"
# Note: seems poorly resolved/documents (Strunecky et al. 2023)

# Cuspidothrix has incorrect "Cyanobacteriales" order
data[which(data$genus == "Cuspidothrix"), "order"] <- "Nostocales"
data[which(data$genus == "Cuspidothrix"), "family"] <- "Aphanizomenonaceae"
# NOTE: ambiguous placement from Strunecky et al. (2023)

# Cyanothece incorrect Cyanobacteriales order
data[which(data$genus == "Cyanothece"), "order"] <- "Gomontiellales"
data[which(data$genus == "Cyanothece"), "family"] <- "Cyanothecaceae"
# NOTE: poorly resolved (Strunecky et al. 2023)

# Cylindrospermum & Cylindrospermopsis have been reassigned to Aphanizomenonaceae family
# also both have incorrect Cyanobacteriales order
data[which(data$genus == "Cylindrospermum" | data$genus == "Cylindrospermopsis"), "family"] <- "Aphanizomenonaceae"
data[which(data$genus == "Cylindrospermum" | data$genus == "Cylindrospermopsis"), "order"] <- "Nostocales"

# Desmonostoc same order issue
data[which(data$genus == "Desmonostoc"), "order"] <- "Nostocales"

# Two reads for Geitlerinema have same issue & oxyphoto issue and a mispelling in Family
data[which(data$genus == "Geitlerinema"), "family"] <- "Geitlerinemataceae"
data[which(data$genus == "Geitlerinema"), "order"] <- "Geitlerinematales"

# "Cyanobacteriales" order issue for Geminocystis as well
data[which(data$genus == "Geminocystis"), "order"] <- "Chroococcales"

# Same order issue and a family issue for Gloeothece
data[which(data$genus == "Gloeothece"), "family"] <- "Aphanothecaceae"
data[which(data$genus == "Gloeothece"), "order"] <- "Chroococcales"

# Same order issue and a family issue for Gloeotrichia
data[which(data$genus == "Gloeotrichia"), "family"] <- "Gloeotrichiaceae"
data[which(data$genus == "Gloeotrichia"), "order"] <- "Nostocales"

# Same order issue for Kamptonema
data[which(data$genus == "Kamptonema"), "order"] <- "Oscillatoriales"

# One read for Leptolyngbya with Oxyphoto order issue
data[which(data$genus == "Leptolyngbya"), "family"] <- "Leptolyngbyales"
data[which(data$genus == "Leptolyngbya"), "order"] <- "Leptolyngbyaceae"

# Merismopedia with family and order ("Cyanobacteriales) issues
# noting that a lot of the Cyanobacteriales also have incorrect families that are Microcystaceae
data[which(data$genus == "Merismopedia"), "family"] <- "Merismopediaceae"
data[which(data$genus == "Merismopedia"), "order"] <- "Synechococcales"

# Microcoleus incorrect family and order ("Cyanobacteriales")
# (I guess this is the original Microcoleus reads that weren't Tychonema)
data[which(data$genus == "Microcoleus"), "family"] <- "Microcoleaceae"
data[which(data$genus == "Microcoleus"), "order"] <- "Oscillatoriales"

# Microcystis wrong order ("Cyanobacteriales")
data[which(data$genus == "Microcystis"), "order"] <- "Chroococcales"

# Incorrect order for Nodosilinea
data[which(data$genus == "Nodosilinea"), "order"] <- "Nodosilineales"

# Incorrect "Cyanobacteriales" order for both Nodularia and Nostoc plus incorrect family for Nodularia
data[which(data$genus == "Nostoc" | data$genus == "Nodularia"), "order"] <- "Nostocales"
data[which(data$genus == "Nodularia"), "family"] <- "Nodulariaceae"

# Incorrect oxyphotobacteria & unknown family for Oscillatoria & Phormidium
data[which(data$genus == "Oscillatoria" | data$genus == "Phormidium"), "family"] <- "Oscillatoriaceae"
data[which(data$genus == "Oscillatoria" | data$genus == "Phormidium"), "order"] <- "Oscillatoriales"

# Incorrect "Cyanobacteriales" order for Planktothricoides & Planktothrix
data[which(data$genus == "Planktothricoides" | data$genus == "Planktothrix"), "order"] <- "Oscillatoriales"

# Incorrect family & "Cyanobacteriales" order for Pleurocapsa
data[which(data$genus == "Pleurocapsa"), "order"] <- "Pleurocapsales"
data[which(data$genus == "Pleurocapsa"), "family"] <- "Chroococcales"

# Incorrect family & "Cyanobacteriales" order for Potamolinea
data[which(data$genus == "Potamolinea"), "order"] <- "Coleofasciculales"
data[which(data$genus == "Potamolinea"), "family"] <- "Wilmottiaceae"

# Oxyphoto order issue and unknown family for a Pseudanabaena read
data[which(data$genus == "Pseudanabaena"), "order"] <- "Pseudanabaenales"
data[which(data$genus == "Pseudanabaena"), "family"] <- "Pseudanabaenaceae"

# Incorrect "Cyanobacteriales" order for Richelia and Family
data[which(data$genus == "Richelia"), "order"] <- "Nostocales"
data[which(data$genus == "Richelia"), "family"] <- "Rivulariaceae"

# Incorrect "Cyanobacteriales" order and family for Scytonema
data[which(data$genus == "Scytonema"), "order"] <- "Nostocales"
data[which(data$genus == "Scytonema"), "family"] <- "Scytonemataceae"

# Incorrect "Cyanobacteriales" order and family for Snowella
data[which(data$genus == "Snowella"), "order"] <- "Chroococcales"
data[which(data$genus == "Snowella"), "family"] <- "Microcystaceae"

# Incorrect "Cyanobacteriales" order and family for Sphaerospermopsis
data[which(data$genus == "Sphaerospermopsis"), "order"] <- "Nostocales"
data[which(data$genus == "Sphaerospermopsis"), "family"] <- "Aphanizomenonaceae"

# Incorrect order and family for Synechococcus
data[which(data$genus == "Synechococcus"), "order"] <- "Synechococcales"
data[which(data$genus == "Synechococcus"), "family"] <- "Synechococcaceae"

# Incorrect "Cyanobacteriales" order and family for Synechocystis
data[which(data$genus == "Synechocystis"), "order"] <- "Synechococcales"
data[which(data$genus == "Synechocystis"), "family"] <- "Merismopediaceae"

# Incorrect "Cyanobacteriales" order and family for Tolypothrix
data[which(data$genus == "Tolypothrix"), "order"] <- "Nostocales"
data[which(data$genus == "Tolypothrix"), "family"] <- "Tolypothrichaceae"
# NOTE: highly unresolved (Strunecky et al. 2023)

# Incorrect "Cyanobacteriales" order for Trichodesmium
data[which(data$genus == "Trichodesmium"), "order"] <- "Oscillatoriales"

# Incorrect "Cyanobacteriales" order for Trichormus
data[which(data$genus == "Trichormus"), "order"] <- "Nostocales"

# change the remainder of Oxyphotobacteria and Unknown Family to NA (These have no genus)
data[which(data$order == "Oxyphotobacteria_Incertae_Sedis"), c("order", "family", "genus")] <- NA
data[which(data$family == "Unknown_Family"), c("order", "family", "genus")] <- NA

# lastly, there are some known families without genus under "Cyanobacteriales" order
data[which(data$family == "Nostocaceae"), "order"] <- "Nostocales"
data[which(data$family == "Microcystaceae"), "order"] <- "Chroococcales"
data[which(data$family == "Coleofasciculaceae"), "order"] <- "Coleofasciculales"
data[which(data$family == "Osillatoriaceae"), "order"] <- "Oscillatoriales"
data[which(data$family == "Xenococcaceae"), "order"] <- "Pleurocapsales"
data[which(data$family == "Chroococcidiopsaceae"), "order"] <- "Chroococcidiopsidales"

# now can just remove the remainder of Cyanobacteriales order that aren't given a family or genus
data[which(data$order == "Cyanobacteriales"), "order"] <- NA

# double-check that everything looks correct
names_final <- unique(data[which(data$phylum == "Cyanobacteria"), c("order", "family", "genus")])
view(names_final)

# note: may reassign Trichormus to Anabaena because it's also an equal match in BLAST

#### (10) Fix field date labels ####

# remove TM sample taken on 9/6 at SFE-M-1S (taken inconsistently from other samples; turkey-baster, lots of water)
# & replace with sample taken on 9/8 (taken correctly; mostly mat )
data <- data[-which(data$field_date == "9/6/2022" & data$site_reach == "SFE-M-1S" & data$sample_type == "TM"),]
data[which(data$field_date == "9/8/2022" & data$site_reach == "SFE-M-1S"& data$sample_type == "TM"),]$field_date <- "9/6/2022"

#### (11) Saving ####

# save csv
write.csv(data, "./data/molecular/16s_nochimera_rarefied_95_FINAL.csv", row.names = FALSE)

## save versions of TM and TAC with Microcoleus and Anabaena/Cylindrospermum/Trichormus removed respectively

# for microcoleus
TM_data <- data %>% 
  filter(genus != "Microcoleus") %>% 
  filter(sample_type == "TM")
write.csv(TM_data, "./data/molecular/16s_nochimera_rarefied_95_TM_nomicro.csv", row.names = FALSE)

# for anabaena/cylindrospermum
TAC_data <- data %>% 
  # note: cylindrospermopsis is also reading a highly close match to anabaena in BLAST so will
  # remove that as well
  filter(! genus %in% c("Anabaena","Cylindrospermum","Trichormus",  "Cylindrospermopsis")) %>% 
  filter(sample_type == "TAC")
write.csv(TAC_data, "./data/molecular/16s_nochimera_rarefied_95_TAC_noanacyl.csv", row.names = FALSE)
