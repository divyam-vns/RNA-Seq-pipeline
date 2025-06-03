# Load required libraries
# DESeq2 Analysis in R
library(DESeq2)
library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)

# Ensure that you have a sample information file (sample_info.csv) with columns sample and condition.
# Set paths
gtf_file <- "/path/to/annotation.gtf"
bam_dir <- "/path/to/aligned_bams"
sample_info_file <- "/path/to/sample_info.csv"

# Create TxDb object
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

# Extract exons by gene
exons_by_gene <- exonsBy(txdb, by = "gene")

# List BAM files
bam_files <- list.files(bam_dir, pattern = "*.bam$", full.names = TRUE)

# Read sample information
sample_info <- read.csv(sample_info_file)

# Create BamFileList
bam_list <- BamFileList(bam_files, yieldSize = 2000000)

# Summarize overlaps
se <- summarizeOverlaps(features = exons_by_gene,
                        reads = bam_list,
                        mode = "Union",
                        singleEnd = FALSE,
                        ignore.strand = TRUE,
                        fragments = TRUE)

# Add sample information
colData(se) <- DataFrame(sample_info)

# Create DESeqDataSet
dds <- DESeqDataSet(se, design = ~ condition)

# Run DESeq2
dds <- DESeq(dds)

# Get results
res <- results(dds)

# Save results
write.csv(as.data.frame(res), file = "deseq2_results.csv")
