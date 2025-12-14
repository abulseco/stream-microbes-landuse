# Stream Project
# Summer 2024 collection
# DADA2 Analysis - 16S
# Aug 2025 seq run
# Last updated: 12/14/2025 by ANB

# PLEASE READ: We sequenced the samples for this project on several runs
# November 2024 (1st run, bad reverse reads see below; [no #])
# July 2025 (16S [1] and ITS [2])
# August 2025 (redo for some sediment samples plus water samples; [3])

# The code below is for the August 2025 sequencing run only pre-merge
# August 2025 Run----
path="~/Documents/R-local/2024_pielter-stream-landuse/stream-microbes-landuse/RAW-DATA/2025-AUG/"
list.files(path)

# Inspect reads
# Need to make sure this matches what your sequencs look like

# Add the .gz if your files are still gunzipped
fnFs2 <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs2 <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names
sample.names2 <- sapply(strsplit(basename(fnFs2), "_"), `[`, 1)
list(sample.names2)

# To make sure the forward and reverse reads are the same
count_reads <- function(file) {
  fq <- readFastq(file)
  length(fq)
} # count reads function

# Requires the package "ShortRead"
counts_before <- data.frame(
  sample = basename(fnFs2),
  forward_reads = sapply(fnFs2, count_reads),
  reverse_reads = sapply(fnRs2, count_reads)
)
View(counts_before)
write.csv(counts_before, "read_counts_for_16S-August2025.csv")
# We will need to redo some of these

# Inspect read quality profiles
P1 <- plotQualityProfile(fnFs2[1:15]) # change to your number of samples, max
P2 <- plotQualityProfile(fnRs2[1:15])
P1 # Drops off right at the end
P2 # looks good

## Filter & Trim
# Moving forward, will only be analyzing the forward reads
# Place filtered files in filtered/ subdirectory
filtFs2 <- file.path(path, "filtered", paste0(sample.names2, "_F_filt.fastq.gz"))
filtRs2 <- file.path(path, "filtered", paste0(sample.names2, "_R_filt.fastq.gz"))
names(filtFs2) <- sample.names2
names(filtRs2) <- sample.names2

# Trimming primers
FWD_PRIMER_LEN <-19
REV_PRIMER_LEN <-30

# Here, you are trimming your reads for quality
# Removing primers
# Truncating the ends of the reads
out2 <- filterAndTrim(fnFs2, filtFs2, fnRs2, filtRs2, 
                      truncLen=c(150,150),
                      trimLeft = c(FWD_PRIMER_LEN, REV_PRIMER_LEN),
                      maxN = 0,  # Discard reads with any Ns
                      maxEE = 2,  # Allow max 2 expected errors
                      truncQ = 2,  # Truncate reads at the first quality score <= 2
                      compress = TRUE,  # Compress the output files
                      verbose = TRUE, 
                      multithread = TRUE) 
head(out2)
saveRDS(out2,"DADA2_OUTPUT/out2.rds")

## Learn the error rates
errF2 <- learnErrors(filtFs2, multithread=TRUE)
errR2 <- learnErrors(filtRs2, multithread=TRUE)
saveRDS(errF2,"DADA2_OUTPUT/errF2.rds")
saveRDS(errR2,"DADA2_OUTPUT/errR2.rds")
# errF2 <-readRDS("errF2.rds")
# errR2 <-readRDS("errR2.rds")

# Visualize the errors
plotErrors(errF2, nominalQ=TRUE)
plotErrors(errR2, nominalQ=TRUE)

## Dereplicate----
# This makes downstream analysis a bit quicker
derepFs2 <- derepFastq(filtFs2, verbose=TRUE)
derepRs2 <- derepFastq(filtRs2, verbose=TRUE)
saveRDS(derepFs2, "DADA2_OUTPUT/derepFs2.rds")
saveRDS(derepRs2, "DADA2_OUTPUT/derepRs2.rds")
# To read back in
# derepFs2 <- readRDS("DADA2_OUTPUT/derepFs2.rds")
# derepRs2 <- readRDS("DADA2_OUTPUT/derepRs2.rds")
# Name the derep-class objects by the sample names
names(derepFs2) <- sample.names2
names(derepRs2) <- sample.names2

# Run the DADA2 algorithm
# If you don't replicate, use "filtFs" 
dadaFs2 <- dada(derepFs2, err=errF2, multithread=TRUE)
saveRDS(dadaFs2, "DADA2_OUTPUT/dadaFs2.rds")
dadaRs2 <- dada(derepRs2, err=errR2, multithread=TRUE)
saveRDS(dadaRs2, "DADA2_OUTPUT/dadaRs2.rds")
dadaFs2[[1]]
# To read back in
# dadaFs2 <- readRDS("DADA2_OUTPUT/dadaFs2.rds")
# dadaRs2 <- readRDS("DADA2_OUTPUT/dadaRs2.rds")

# Merge your forward and reverse reads
mergers2 <- mergePairs(dadaFs2, derepFs2, dadaRs2, derepRs2, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers2[[1]])
saveRDS(mergers2, "DADA2_OUTPUT/mergers2.rds")
# To read back in
# mergers2 <- readRDS("DADA2_OUTPUT/mergers2.rds")

# Construct a sequence table
seqtab2 <- makeSequenceTable(mergers2)
dim(seqtab2)
table(nchar(getSequences(seqtab2)))
saveRDS(seqtab2, "DADA2_OUTPUT/seqtab2.rds")

# STOP HERE and proceed to "seq_run_merge.R" to merge seq tables
