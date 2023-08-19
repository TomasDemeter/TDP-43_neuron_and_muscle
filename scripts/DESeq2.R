library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggplot2)


#################################
# Nature-like theme for ggplot2 #
#################################

nature_theme <- function(
  base_size = 12,
  base_family = "",
  color_palette = c("#8a3838", "#f07f7e", "#51848a", "#99c0cc")) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      # Set the background color to white
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      
      # Set the size and color of the axis lines and ticks
      axis.line = element_line(size = 1, color = "black"),
      axis.ticks = element_line(size = 1, color = "black"),
      
      # Set the font size and color of the axis text and title
      axis.text = element_text(size = rel(1.5), color = "#5c5a5a"),
      axis.title = element_text(size = rel(1.5), color = "black"),
      
      # Set the font size and color of the legend title and text
      legend.text = element_text(size = rel(1.2), color = "black"),

      axis.text.x = element_text(angle = 45, hjust = 1),

      # Remove the legend title
      legend.title = element_blank()
    )
}


####################
# Loading the data #
####################

#args <- commandArgs(trailingOnly = TRUE)
raw_counts_filepath <- "results/feature_counts_table.tsv" #args[1]
meta_data_filepath <- "data/raw_reads/SRR_metadata.csv" #args[2]
#DESeq2_output <- args[3]

# load featureCounts output, metadata adn gtf file
raw_counts <- read.csv(raw_counts_filepath, skip = 1, sep = "\t", row.names="Geneid")
meta_data <- read.csv(meta_data_filepath, row.names="SRR")


######################################
# Preparing data for FPKM and DESeq2 #
######################################

# save gene length info for later fpmk calculcations
gene_length <- raw_counts["Length"]

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


#################
# DESeq2 analysis
#################

#create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = meta_data,
                              design = ~ cell.type + treatment)


# add gene_length info to dds object
mcols(dds)$basepairs <- gene_length$Length

dds <- dds[which(rowSums(counts(dds)) >= 1),]


# calculate fpkm
fpkm_values <- fpkm(dds)

# Perform PCA on the FPKM values
pca_res <- prcomp(t(fpkm_values))

# Create a data frame with the PCA results and sample information
pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], treatment = colData(dds)$treatment, cell.type = colData(dds)$cell.type)

# Create a new variable that combines treatment and cell type
pca_df$group <- interaction(pca_df$treatment, pca_df$cell.type)

# Generate the PCA plot
pca_plot_fpkm <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 5) +
  xlab(paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100), "%)")) +
  ylab(paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100), "%)")) + 
  scale_color_manual(values = c("#8a3838", "#f07f7e", "#51848a", "#99c0cc"),
  labels = c("C2C12 siLUC", "C2C12 siTDP", "NSC34 siLUC", "NSC34 siTDP")) +
  nature_theme()

print(pca_plot_fpkm)



#____________________________________________________
# Extract the FPKM values for the gene of interest
gene_fpkm <- fpkm_values["ENSMUSG00000041459",]

# Create a data frame with the FPKM values, PCA2 scores, and sample information
fpkm_pca_df <- data.frame(fpkm = gene_fpkm, pca2 = pca_res$x[,2], treatment = colData(dds)$treatment, cell.type = colData(dds)$cell.type)

# Create a new variable that combines treatment and cell type
fpkm_pca_df$group <- interaction(fpkm_pca_df$treatment, fpkm_pca_df$cell.type)

log_gene_fpkm <- log10(gene_fpkm + 1)

# Create a data frame with the log10(FPKM + 1) values and sample information
log_fpkm_df <- data.frame(sample = colnames(fpkm_values), log_fpkm = log_gene_fpkm, treatment = colData(dds)$treatment, cell.type = colData(dds)$cell.type)

# Create a new variable that combines treatment and cell type
log_fpkm_df$group <- interaction(log_fpkm_df$treatment, log_fpkm_df$cell.type)

# Generate the scatter plot
log_fpkm_plot <- ggplot(log_fpkm_df, aes(x = group, y = log_fpkm, color = group)) +
  geom_point(size = 5) +
  ggtitle("") +
  ylab("log10(FPKM + 1)") +
  xlab("") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_x_discrete(labels = c("siLUC", "siTDP", "siLUC", "siTDP")) +
  scale_color_manual(values = c("#8a3838", "#f07f7e", "#51848a", "#99c0cc"),
  labels = c("C2C12 siLUC", "C2C12 siTDP", "NSC34 siLUC", "NSC34 siTDP")) +
  nature_theme()

print(log_fpkm_plot)



#____________________

# specify factor level (what is treated and what is control)
dds$treatment <- relevel(dds$treatment, ref = "siRNA Control (Luciferase)")

# run DESeq2
dds <- DESeq(dds)

# Perform multiple testing adjustment using Benjamini and Hochberg's approach
res <- results(dds, alpha = 0.05, pAdjustMethod = "BH")

# Extract differentially expressed genes (padj < 0.05)
degs <- subset(res, padj < 0.05)


# plot gene expression dispersion
plotDispEsts(dds)


# Hierarchically cluster the differentially expressed genes based on log10(FPKM + 1)

vsd <- vst(dds, blind = FALSE)
#mat <- assay(vsd)[rownames(degs),]
#mat <- mat - rowMeans(mat)
#pheatmap(log10(mat + 1))


# Illustrate distance between silenced and control samples of each cell line with PCA
# generate the PCA plot


# the NSC34 groups might be switched here!!!!
pca_plot <- plotPCA(vsd, intgroup = c("treatment", "cell.type")) +
  geom_point(size = 5) +
  xlab(paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100), "%)")) +
  ylab(paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100), "%)")) +
  coord_cartesian(ylim = c(-4, 4))# +
  #scale_color_manual(values = c("#8a3838", "#f07f7e", "#51848a", "#99c0cc"), labels = c("C2C12 siLUC", "C2C12 siTDP", "NSC34 siLUC", "NSC34 siTDP")) +
  # nature_theme()

print(pca_plot)


#____________

# adjust the y-axis limits
pca_plot <- pca_plot + coord_cartesian(ylim = c(-4, 4))
# display the plot
print(pca_plot)


# save the output
saveRDS(dds, file = DESeq2_output)
