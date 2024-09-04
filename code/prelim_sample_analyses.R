#### Assigning sample names and types to sequence reads
### 9.4.2024

# NOTE: our "unknown" IDs from the plates are 104, 115, 208, 1000, and 1001
# our "unknown" IDs from the sample metadata are 108, 801, 552, 554, and 931

# also have sample 27 which we were unsure if it was Microcoleus or not
# (though I don't think it has a corresponding anatoxin concentration measurement bc sample was so small anyways)

## Read in data

# taxonomy data
taxonomy <- read.csv("Bac_counts_taxonomy_cleanedup.csv")

# metadata

### RENAME THIS R CODE LOL