library(DESeq2)
library(tidyverse)
library(pheatmap)
library(GenomicRanges)
library(rtracklayer)

#args <- commandArgs(trailingOnly = TRUE)
raw_counts_filepath <- "results/feature_counts_table.tsv" #args[1]
meta_data_filepath <- "data/raw_reads/SRR_metadata.csv" #args[2]
gtf_filepath <- "config/refs/Mus_musculus.GRCm39.110.gtf"
#DESeq2_output <- args[3]

# load featureCounts output, metadata adn gtf file
raw_counts <- read.csv(raw_counts_filepath, skip = 1, sep = "\t", row.names="Geneid")
meta_data <- read.csv(meta_data_filepath, row.names="SRR")


# Drop the columns 'Chr', 'Start', 'End', 'Strand', and 'Length'
raw_counts <- raw_counts[, !(colnames(raw_counts) %in% c('Chr', 'Start', 'End', 'Strand', 'Length'))]

# Rename the remaining columns to keep only the part of the column name starting with 'SRR' (included) and ending with '_fastp_' (excluded)
colnames(raw_counts) <- sapply(colnames(raw_counts), function(x) {
  start_index <- regexpr('SRR', x)
  end_index <- regexpr('_fastp_', x)
  if (start_index != -1 & end_index != -1) {
    substr(x, start_index, end_index - 1)
  } else {
    x
  }
})

# keep only relevant columns from metadata and set them as factors
meta_data <- meta_data[,c("cell.type", "treatment")]
meta_data$cell.type <- as.factor(meta_data$cell.type)
meta_data$treatment <- as.factor(meta_data$treatment)

# creating genomic ranges object

gr <- import(gtf_filepath)

#create DESeq2 object

dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = meta_data,
                              design = ~ treatment)


gr_2 <- makeGRangesFromGTF(
  gtf_filepath,
  level = c("genes", "transcripts"),
  ignoreVersion = TRUE,
  extraMcols = TRUE
)

mcols(dds)

rowRanges(dds) <- gr

# keep only row that have more than 10 reads across all samples
dds <- dds[rowSums(counts(dds)) >= 1, ]

# specify factor level (what is treated and what is control)
dds$treatment <- relevel(dds$treatment, ref = "siRNA Control (Luciferase)")

# run DESeq2
dds <- DESeq(dds)

# Perform multiple testing adjustment using Benjamini and Hochberg's approach
res <- results(dds, alpha = 0.05, pAdjustMethod = "BH")

# Extract differentially expressed genes (padj < 0.05)
degs <- subset(res, padj < 0.05)

# Hierarchically cluster the differentially expressed genes based on log10(FPKM + 1)
vsd <- vst(dds)
mat <- assay(vsd)[rownames(degs),]
mat <- mat - rowMeans(mat)
pheatmap(log10(mat + 1))

# Illustrate distance between silenced and control samples of each cell line with PCA
plotPCA(vsd, intgroup = c("treatment", "cell.type"))








fpkm_values <- fpkm(dds)



# Test for significance in gene expression levels between cell lines using Wilcoxon signed-rank test
wilcox.test(log10(mat + 1))

# save the output
saveRDS(dds, file = DESeq2_output)
