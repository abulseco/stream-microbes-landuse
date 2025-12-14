# Stream landuse project
# Merging sequencing runs
# Last updated: 12/14/2025 by ANB

# PLEASE READ: This script assumes you have run the following:
    # july_2025_16S_dada2.R
    # august_2025_16S_dada2.R

# Setup Environ----
## Load libraries----
library(dada2); library(here)

setwd <- here("MERGE_SEQRUNS_16S/")
list.files(setwd)

# Merging sequencing tables----
# Function to add run-specific suffix to sample names
rename_samples <- function(seqtab, suffix) {
  samples <- paste0(rownames(seqtab), "_", suffix)
  rownames(seqtab) <- samples
  return(seqtab)
}

seqtab_run1 <- readRDS("MERGE_SEQRUNS_16S/seqtab1.rds")
seqtab_run2 <- readRDS("MERGE_SEQRUNS_16S/seqtab2.rds")

# Merge sequence tables
seqtab_all <- mergeSequenceTables(seqtab_run1, seqtab_run2)
dim(seqtab_all)
table(nchar(getSequences(seqtab_all)))
sample_names_check <- rownames(seqtab_all)
print(sample_names_check)
write.table(sample_names_check, "MERGE_SEQRUNS_16S/sample_names_all.txt", sep="\t", quote=F, col.names=NA)
saveRDS(seqtab_all, "MERGE_SEQRUNS_16S/seqtab_all.rds")

# Remove chimeras----
seqtab_nochim_all <- removeBimeraDenovo(seqtab_all, method = "consensus", multithread = TRUE)
saveRDS(seqtab_nochim_all, "MERGE_SEQRUNS_16S/seqtab-nochim-all.rds") 
# seqtab_nochim_all <- readRDS("MERGE_SEQRUNS_16S/seqtab-nochim") 

# Assigning taxonomy----
taxa_all <- assignTaxonomy(seqtab_nochim_all, "DATABASE/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
saveRDS(taxa_all, "MERGE_SEQRUNS_16S/taxa-table-all.rds")
# taxa_allDNA <- readRDS("taxa-table-DNA.rds") # reading back in

# Write sample names
write.table(sample.names1, "sample_names1.txt", sep="\t", quote=F, col.names=NA)
write.table(sample.names2, "sample_names2.txt", sep="\t", quote=F, col.names=NA)
