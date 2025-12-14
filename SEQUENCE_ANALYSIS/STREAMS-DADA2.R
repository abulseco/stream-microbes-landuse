# Stream Project
# Summer 2024 collection
# DADA2 Analysis
# Last updated: 12/14/2025 by ANB

# PLEASE READ: We sequenced the samples for this project on several runs
    # November 2024 (1st run, bad reverse reads see below; [no #])
    # July 2025 (16S [1] and ITS [2])
    # August 2025 (redo for some sediment samples plus water samples; [3])

# Setup----
## Load libraries----
library(dada2); library(ShortRead); library(patchwork)

# Nov 2024 Run----
path="~/Documents/R-local/2024_pielter-stream-landuse/stream-microbes-landuse/RAW-DATA/2024-NOV/"
list.files(path)

## Inspect reads
# Need to make sure this matches what your sequencs look like

# Add the .gz if your files are still gunzipped
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
list(sample.names)

# To make sure the forward and reverse reads are the same
count_reads <- function(file) {
  fq <- readFastq(file)
  length(fq)
} # count reads function

# Requires the package "ShortRead"
counts_before <- data.frame(
  sample = basename(fnFs),
  forward_reads = sapply(fnFs, count_reads),
  reverse_reads = sapply(fnRs, count_reads)
)
View(counts_before)
write.csv(counts_before, "read_counts_for_16S-July2025.csv")


# Inspect read quality profiles
P1 <- plotQualityProfile(fnFs[1:15]) # change to your number of samples, max
P2 <- plotQualityProfile(fnRs[1:15])
P1
P2 # This doesn't look good. Will do a quick analysis but then may need to return here later

## Filter & Trim
# Moving forward, will only be analyzing the forward reads
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Trimming primers
FWD_PRIMER_LEN <-19
REV_PRIMER_LEN <-30

# Here, you are trimming your reads for quality
# Removing primers
# Truncating the ends of the reads
out <- filterAndTrim(fnFs, filtFs,
                     trimLeft = FWD_PRIMER_LEN,
                     truncLen = 150,  # Example values
                     maxN = 0,  # Discard reads with any Ns
                     maxEE = 2,  # Allow max 2 expected errors
                     truncQ = 2,  # Truncate reads at the first quality score <= 2
                     compress = TRUE,  # Compress the output files
                     verbose = TRUE, 
                     multithread = TRUE) 
head(out)
saveRDS(out,"out.rds")

## Learn the error rates----
# DADA 2's algorithm is based on modeling the error rates in the sequencing process. 
# It leverages a parametric error model (err), which is different for each individual 
# amplicon dataset. The learnErrors method learns this error model from the data, by
# alternating estimation of the error rates and inference of sample composition until 
# they converge on a jointly consistent solution.
errF <- learnErrors(filtFs, multithread=TRUE)
saveRDS(errF,"errF.rds")

# Visualize the errors
plotErrors(errF, nominalQ=TRUE)

## Dereplicate
derepFs <- derepFastq(filtFs, verbose=TRUE)
saveRDS(derepFs, "derepFS.rds")
names(derepFs) <- sample.names

# To read back in
# derepFs <- readRDS("derepFs.rds")
# derepRs <- readRDS("derepRs.rds")
# Name the derep-class objects by the sample names

## Run the DADA2 algorithm----
# If you don't replicate, use "filtFs" 
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
saveRDS(dadaFs, "dadaFs.rds")
dadaFs[[1]]

# To read back in
# dadaFs <- readRDS("dadaFs.rds")

## Construct a sequence table----
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
table(nchar(getSequences(seqtab)))

## Remove chimeras
# Chimeras are common but shouldn't be more than 3-4% 
# If high, then remove primers and redo analysis
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) # This means that ~2.5% were chimeric
saveRDS(seqtab.nochim, "seqtab-nochim.rds")

## Do a quick check in
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)

## Assigning taxonomy----
# Be sure to put the database in your working directory or another known location
# You can download proper database here: https://benjjneb.github.io/dada2/training.html
taxa <- assignTaxonomy(seqtab.nochim, "~/tax/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
saveRDS(taxa, "taxa-table.rds")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# July 2025 Run----
path="~/Documents/R-local/2024_pielter-stream-landuse/stream-microbes-landuse/RAW-DATA/2025-JUL/2025-JUL-16S/"
list.files(path)

# Inspect reads
# Need to make sure this matches what your sequencs look like

# Add the .gz if your files are still gunzipped
fnFs1 <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs1 <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names
sample.names1 <- sapply(strsplit(basename(fnFs1), "_"), `[`, 1)
list(sample.names1)

# To make sure the forward and reverse reads are the same
count_reads <- function(file) {
  fq <- readFastq(file)
  length(fq)
} # count reads function

# Requires the package "ShortRead"
counts_before <- data.frame(
  sample = basename(fnFs1),
  forward_reads = sapply(fnFs1, count_reads),
  reverse_reads = sapply(fnRs1, count_reads)
)
View(counts_before)
write.csv(counts_before, "read_counts_for_16S-July2025.csv")
# We will need to redo some of these

# Inspect read quality profiles
P1 <- plotQualityProfile(fnFs1[1:15]) # change to your number of samples, max
P2 <- plotQualityProfile(fnRs1[1:15])
P1 # Drops off right at the end
P2 # looks good

## Filter & Trim
# Moving forward, will only be analyzing the forward reads
# Place filtered files in filtered/ subdirectory
filtFs1 <- file.path(path, "filtered", paste0(sample.names1, "_F_filt.fastq.gz"))
filtRs1 <- file.path(path, "filtered", paste0(sample.names1, "_R_filt.fastq.gz"))
names(filtFs1) <- sample.names1
names(filtRs1) <- sample.names1

# Trimming primers
FWD_PRIMER_LEN <-19
REV_PRIMER_LEN <-30

# Here, you are trimming your reads for quality
# Removing primers
# Truncating the ends of the reads
out1 <- filterAndTrim(fnFs1, filtFs1, fnRs1, filtRs1, 
                      truncLen=c(145,150),
                      trimLeft = c(FWD_PRIMER_LEN, REV_PRIMER_LEN),
                      maxN = 0,  # Discard reads with any Ns
                      maxEE = 2,  # Allow max 2 expected errors
                      truncQ = 2,  # Truncate reads at the first quality score <= 2
                      compress = TRUE,  # Compress the output files
                      verbose = TRUE, 
                      multithread = TRUE) 
head(out1)
saveRDS(out,"DADA2_OUTPUT/out1.rds")

## Learn the error rates
errF1 <- learnErrors(filtFs1, multithread=TRUE)
errR1 <- learnErrors(filtRs1, multithread=TRUE)
saveRDS(errF1,"DADA2_OUTPUT/errF1.rds")
saveRDS(errR1,"DADA2_OUTPUT/errR1.rds")
# errF1 <-readRDS("errF1.rds")
# errR1 <-readRDS("errR1.rds")

# Visualize the errors
plotErrors(errF1, nominalQ=TRUE)
plotErrors(errR1, nominalQ=TRUE)

## Dereplicate----
# This makes downstream analysis a bit quicker
derepFs1 <- derepFastq(filtFs1, verbose=TRUE)
derepRs1 <- derepFastq(filtRs1, verbose=TRUE)
saveRDS(derepFs1, "DADA2_OUTPUT/derepFS1.rds")
saveRDS(derepRs1, "DADA2_OUTPUT/derepRs1.rds")
# To read back in
# derepFs <- readRDS("derepFs.rds")
# derepRs <- readRDS("derepRs.rds")
# Name the derep-class objects by the sample names
names(derepFs1) <- sample.names1
names(derepRs1) <- sample.names1

# To read back in
# derepFs1 <- readRDS("DADA2_OUTPUT/derepFs1.rds")
# derepRs1 <- readRDS("DADA2_OUTPUT/derepRs1.rds")
# Name the derep-class objects by the sample names

## Run the DADA2 algorithm
# If you don't replicate, use "filtFs" 
# Run the DADA2 algorithm----
# If you don't replicate, use "filtFs" 
dadaFs1 <- dada(derepFs1, err=errF1, multithread=TRUE)
saveRDS(dadaFs1, "DADA2_OUTPUT/dadaFs.rds")
dadaRs1 <- dada(derepRs1, err=errR1, multithread=TRUE)
saveRDS(dadaRs1, "DADA2_OUTPUT/dadaRs.rds")
dadaFs1[[1]]
# To read back in
# dadaFs <- readRDS("dadaFs.rds")
# dadaRs <- readRDS("dadaRs.rds")

# Merge your forward and reverse reads
mergers1 <- mergePairs(dadaFs1, derepFs1, dadaRs1, derepRs1, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers1[[1]])
saveRDS(mergers1, "DADA2_OUTPUT/mergers1.rds")
# To read back in
# mergers1 <- readRDS("DADA2_OUTPUT/mergers.rds")

# Construct a sequence table
seqtab1 <- makeSequenceTable(mergers1)
dim(seqtab1)
table(nchar(getSequences(seqtab1)))
saveRDS(seqtab1, "DADA2_OUTPUT/seqtab1.rds")

# STOP HERE to analyze August sequences, see below to merge

## Remove chimeras
# Chimeras are common but shouldn't be more than 3-4% 
# If high, then remove primers and redo analysis
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) # This means that ~2.5% were chimeric
saveRDS(seqtab.nochim, "seqtab-nochim.rds")

## Do a quick check in
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)

## Assigning taxonomy
# Be sure to put the database in your working directory or another known location
# You can download proper database here: https://benjjneb.github.io/dada2/training.html
taxa <- assignTaxonomy(seqtab.nochim, "~/tax/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
saveRDS(taxa, "taxa-table.rds")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


