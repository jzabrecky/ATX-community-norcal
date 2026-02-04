# differential abundance analyses walk through

# following: https://sweltr.github.io/high-temp/da.html
# works with both QIIME2 and PICRUST outputs

# install packages
BiocManager::install("microbiomeMarker")
library("microbiomeMarker")


# developer mode
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("yiluheihei/microbiomeMarker")

seq_tab <- readRDS(
  system.file(
    "extdata", "dada2_seqtab.rds",
    package = "microbiomeMarker"
  )
)
tax_tab <- readRDS(
  system.file(
    "extdata", "dada2_taxtab.rds",
    package = "microbiomeMarker"
  )
)
sam_tab <- read.table(
  system.file(
    "extdata", "dada2_samdata.txt",
    package = "microbiomeMarker"
  ),
  sep = "\t",
  header = TRUE,
  row.names = 1
)
ps <- import_dada2(seq_tab = seq_tab, tax_tab = tax_tab, sam_tab = sam_tab)
ps

# normalize data (relative abundances)
normalize(ps, method = "TSS")

# aldex2
data(kostic_crc)
kostic_crc_small <- phyloseq::subset_taxa(
  kostic_crc,
  Phylum %in% c("Firmicutes")
)
mm_lefse <- run_lefse(
  kostic_crc_small,
  group = "DIAGNOSIS"
)
mm_lefse <- run_lefse(
  physeq,
  group = "site"
)
plot_abundance(mm_lefse, group = "DIAGNOSIS")
plot_ef_dot(mm_lefse)
plot_ef_bar(mm_lefse)
plot_cladogram(mm_lefse)

# try to turn my data into a phyloseq object
tm_data <- read.csv("./data/molecular/16s_nochimera_rarefied_95_TM_nomicro.csv") %>% 
  mutate(sample_name = paste(site_reach, ", ", field_date, sep = "")) %>% 
  relocate(sample_name, .before = "site_reach")
seq_df <- tm_data %>% 
  select(sample_name, feature_ID, relative_abundance) %>% 
  pivot_wider(names_from = "feature_ID", values_from = "relative_abundance")
seq_df[is.na(seq_df)] <- 0
seq_matrix <- as.matrix(seq_df[,2:ncol(seq_df)])
rownames(seq_matrix) <- seq_df$sample_name
tax_df <- tm_data %>% 
  select(feature_ID, taxon_full)
sam_df <- tm_data %>% 
  select(sample_name, site_reach, site, field_date) %>% 
  unique()
rownames(sam_df) <- sam_df$sample_name


# to do: make phyloseq object, check transformations
# current transformation should be fine?
# just see what other people do


# converting OTU table to phyloseq object

otumat = matrix(sample(1:100, 100, replace = TRUE), nrow = 10, ncol = 10)
otumat
rownames(otumat) <- paste0("OTU", 1:nrow(otumat))
colnames(otumat) <- paste0("Sample", 1:ncol(otumat))
otumat
myasvs <- tm_data %>% 
  select(feature_ID, sample_name, relative_abundance) %>% 
  pivot_wider(names_from = sample_name, values_from = relative_abundance)
myasvs[is.na(myasvs)] <- 0
myasvmatrix <- as.matrix(myasvs[,2:ncol(myasvs)])
rownames(myasvmatrix) <- myasvs$feature_ID
taxmat = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(otumat), ncol = 7)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxmat
mytax <- tm_data %>% 
  select(feature_ID, domain, phylum, class, order, family, genus, species) %>% 
  unique()
mytaxmatrix <- as.matrix(mytax[,2:ncol(mytax)])
rownames(mytaxmatrix) <- mytax$feature_ID
colnames(mytaxmatrix) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
OTU = otu_table(otumat, taxa_are_rows = TRUE)
OTU = otu_table(myasvmatrix, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
TAX = tax_table(mytaxmatrix)
physeq = phyloseq(OTU, TAX, sam_df)
unique(rownames(TAX))

# checking that rownames for OTUs are unique
any(eval(rownames(myasvmatrix) == unique(rownames(myasvmatrix))) == FALSE)
# they are

sampledata = sample_data(data.frame(
  site = sam_df$site,
  row.names=sample_names(physeq),
  stringsAsFactors=FALSE
))
sampledata

# what about rownames for taxonomy
any(eval(rownames(mytaxmatrix) == unique(rownames(mytaxmatrix))) == FALSE)
rownames(mytaxmatrix)

# finally got a phyloseq object!
plot_bar(physeq, fill = "class")

physeq1 = merge_phyloseq(physeq, sampledata)

library(ANCOMBC)
ancom_da <- ancombc(phyloseq = physeq, formula = "location", 
                    p_adj_method = "fdr", zero_cut = 0.90, lib_cut = 1000, 
                    group = "location", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = FALSE)

ancom_res_df <- data.frame(
  Species = row.names(ancom_da$res$beta),
  beta = unlist(ancom_da$res$beta),
  se = unlist(ancom_da$res$se),
  W = unlist(ancom_da$res$W),
  p_val = unlist(ancom_da$res$p_val),
  q_val = unlist(ancom_da$res$q_val),
  diff_abn = unlist(ancom_da$res$diff_abn))

fdr_ancom <- ancom_res_df %>%
  dplyr::filter(q_val < 0.05)

dim(fdr_ancom)


#subset only the last 400 features for efficiency
data(selex)
selex.sub <- selex[1200:1600,]
conds <- c(rep("NS", 7), rep("S", 7))
x.aldex <- aldex(selex.sub, conds, mc.samples=16, test="t", effect=TRUE, 
                 include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE, gamma=NULL)

mm_lefse <- run_lefse(
  physeq1,
  group = "site"
)
plot_abundance(mm_lefse, group = "site")
plot_ef_dot(mm_lefse)

# how about only filter for cyanobacteria?
cyano <- subset_taxa(physeq1, Phylum=="Cyanobacteria")

mm_lefse <- run_lefse(
  cyano,
  group = "site"
)
plot_ef_dot(mm_lefse)

# How does this information differ from what we currently get form our methods?

ALDEx2::aldex(myASVmatrix, conditions = samp_df)

library(ALDEx2)

data(selex)
selex <- selex[1201:1600,] # subset for efficiency
conds <- c(rep("NS", 7), rep("S", 7))
x <- aldex(selex, conds, mc.samples=2, denom="all",
           test="t", effect=TRUE, paired.test=FALSE)
x
summary(x)
