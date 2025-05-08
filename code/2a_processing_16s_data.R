#### Processing .biom files
### Jordan Zabrecky
## last edited 03.15.2025

## TBD!?!?

#### (1) Loading in libraries and data ####

# loading libraries
lapply(c("tidyverse", "plyr", "biomformat"), require, character.only = T)

## read in .biom files
# only processed version (i.e. chimeras removed) 
biom_norm <- ldply(list.files(path = "./data/molecular/fasta_biom_qza_files", pattern = "nochim"),
                   function(filename) {
  temp <- read_biom(paste("data/molecular/fasta_biom_qza_files/", filename, sep = ""))
  biom <- biom_data(temp)
})

test <- as(biom, 'TsparseMatrix') 
test2 <- data.frame(taxon_id = test @Dimnames[[1]][ test @i + 1],
                   seq_SampleID = test @Dimnames[[2]][ test @j + 1],
                   seq_abund = test @x)


# rarefied version (i.e. where low abundance samples get 
# filtered out and existing features get subsampled.)