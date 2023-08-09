library(DESeq2)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
raw_counts_filepath <- args[1]
meta_data_filepath <- args[1]


raw_counts <- read.csv(raw_counts_filepath, skip = 1, sep = "\t")

# this part might go into snakefile
meta_data <- read.csv(meta_data_filepath)
meta_data <- meta_data[,c("SRR", "Cell type", "Treatement")]


#check that all the columns of featurecounts match the rows of metadata (name and order)
all(colnames(raw_counts) == rownames(meta_data))

dds <- DESeqDataSetFromMatrix(countData = raw_count, colData = meta_data, design = ~ Treatment)

# keep only row that have more than 10 reads across all samples
dds <- dds[rawSums(counts(dds)) >= 10, ]

# specify factor level (what is treated and what is control)
dds$Treatment <- relevel(dds$Treatment, ref = "siRNA Control (Luciferase)")

# run DESeq2
dds <- DESeq(dds)

# Calculate FPKM values
fpkm_values <- fpkm(dds, robust = TRUE)