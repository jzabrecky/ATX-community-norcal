
library("tidyverse")
library("phyloseq")

data <- read.csv("./data/molecular/16s_nochimera_rarefied_95_filtered_rawreads.csv")

# taxonomy table
taxonomy <- data %>% 
  dplyr::select(feature_ID, domain, phylum, class, order, family, genus, species) %>% 
  unique()
length(unique(taxonomy$feature_ID)) # same length as taxonomy!

# turn dataframe back into a phyloseq object?
# test first across rivers
asv_freq_table_rivers <- data %>% 
  dplyr::select(site, feature_ID, abundance) %>% 
  group_by(site, feature_ID) %>% 
  dplyr::summarize(total_reads = sum(abundance)) %>% 
  pivot_wider(names_from = site, values_from = total_reads, values_fill = 0)

# left join and change names?
together <- left_join(taxonomy, asv_freq_table_rivers, by = c("feature_ID"))
together$feature_ID <- rep(paste0("OTU", 1:nrow(together)))

tax <- as.matrix(together[,2:8])
rownames(tax) <- together$feature_ID
taxphylo <- tax_table(tax)

asv <- as.matrix(together[,9:11])
rownames(asv) <- together$feature_ID
asvphylo <- otu_table(asv, taxa_are_rows = TRUE)

phyloseqobj <- phyloseq(asvphylo, taxphylo)


library(ggVennDiagram)

ggVennDiagram::ggVennDiagram(asv_freq_table_rivers)
