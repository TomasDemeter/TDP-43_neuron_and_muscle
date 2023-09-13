library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ggVennDiagram)
library(biomaRt)
library(ggrepel)
library(dplyr)


#raw_counts_filepath <- "results//featureCounts/feature_counts_table.tsv"
#meta_data_filepath <- "data/raw_reads/SRR_metadata.csv"
#RBPs_path <- "./data/RNAbinding_proteins_GeneIDs.csv"
#output_folder <- "results/DESeq2_output"

####################
# Loading the data #
####################

args <- commandArgs(trailingOnly = TRUE)
raw_counts_filepath <- args[1] 
meta_data_filepath <- args[2] 
RBPs_path <- args[3]
output_folder <- args[4]

# load featureCounts output, metadata adn gtf file
raw_counts <- read.csv(raw_counts_filepath, skip = 1, sep = "\t", row.names="Geneid")
meta_data <- read.csv(meta_data_filepath, row.names="SRR")


################################################
################################################
# DEFINING FUNCTIONS FOR ANALYSIS AND PLOTTING # 
################################################
################################################

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
      axis.line = element_line(linewidth = 1, color = "black"),
      axis.ticks = element_line(linewidth = 1, color = "black"),
      
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


##########################################
# function for preparing data for DESeq2 #
##########################################

# preparing count matrix
prepare_counts <- function(raw_counts) {
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

  return(raw_counts)
}

# preparing sample data
prepare_sample_data <- function(meta_data){
  # rename values in treatment and keep only relevant columns from metadata and set them as factors
  meta_data$treatment <- ifelse(grepl("Control", meta_data$treatment, ignore.case = TRUE), "siLUC", 
  ifelse(grepl("TDP", meta_data$treatment, ignore.case = TRUE), "siTDP", meta_data$treatment))

  sample_data <- meta_data[,c("cell.type", "treatment")]
  sample_data$cell.type <- as.factor(sample_data$cell.type)
  sample_data$treatment <- as.factor(sample_data$treatment)

  return(sample_data)
}


############################
# DESeq2 analysis function #
############################
generate_dds <- function(counts, sample_data, experimental_design, gene_lengths) {
  # create DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = sample_data,
                                design = experimental_design)

  # specify factor level (what is treated and what is control)
  dds$treatment <- relevel(dds$treatment, ref = "siLUC")

  # add gene_length info to dds object
  mcols(dds)$basepairs <- gene_lengths$Length

  # add gene_ids
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  gene_ids <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
                    filters = "ensembl_gene_id",
                    values = rownames(rowData(dds)),
                    mart = ensembl)

  rowData(dds)$gene_id <- gene_ids$mgi_symbol[match(rownames(rowData(dds)), gene_ids$ensembl_gene_id)]

  # filter out genes with less than 1 read in all samples
  dds <- dds[which(rowSums(counts(dds)) >= 1),]

  return(dds)
}


#####################################################
# Function to determine expression levels of Tardbp #
#####################################################
gene_expressions <- function(dds, fpkm, gene) {
  # Extract the FPKM values for the gene of interest (mouse TDP-43)
  gene_fpkm <- fpkm[gene,]

  # Create a data frame with the log10(FPKM + 1) values and sample information
  log_gene_fpkm <- log10(gene_fpkm + 1)
  log_fpkm_df <- data.frame(sample = colnames(fpkm), log_fpkm = log_gene_fpkm, treatment = colData(dds)$treatment, cell.type = colData(dds)$cell.type)

  # Create a new variable that combines treatment and cell type
  log_fpkm_df$group <- paste(log_fpkm_df$cell.type, log_fpkm_df$treatment, sep = " ")
  # Generate the scatter plot
  log_fpkm_plot <- ggplot(log_fpkm_df, aes(x = group, y = log_fpkm, color = group)) +
    geom_point(size = 5) +
    ggtitle("") +
    ylab("log10(FPKM + 1)") +
    xlab("") +
    scale_y_continuous(limits = c(0, NA)) +
    scale_color_manual(values = c("#8a3838", "#f07f7e", "#51848a", "#99c0cc")) +
    nature_theme()

  return(log_fpkm_plot)
}


#############################
# Function for PCA plotting #
#############################
pca_plotting <- function(pca_results) {
  # Create a data frame with the PCA results and sample information
  pca_df <- data.frame(PC1 = pca_results$x[,1], PC2 = pca_results$x[,2], treatment = colData(dds)$treatment, cell.type = colData(dds)$cell.type)
  pca_df$PC2 <- -pca_df$PC2

  # Create a new variable that combines treatment and cell type
  pca_df$group <- paste(pca_df$cell.type, pca_df$treatment, sep = " ")

  # Generate the PCA plot
  pca_plot_fpkm <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 5) +
    xlab(paste0("PC1 (", round(summary(pca_results)$importance[2,1]*100), "%)")) +
    ylab(paste0("PC2 (", round(summary(pca_results)$importance[2,2]*100), "%)")) + 
    scale_color_manual(values = c("#8a3838", "#f07f7e", "#51848a", "#99c0cc")) +
    nature_theme()

  return(pca_plot_fpkm)
}


################################################################
# function for visualising PCA2 and a gene of interest in FPKM #
################################################################
# Create a data frame with the FPKM values, PCA2 scores, and sample information
pc2_vs_gene <- function(dds, fpkm, pca_results, gene) {

  gene_fpkm <- fpkm[gene,]
  fpkm_pca_df <- data.frame(fpkm = gene_fpkm, pca2 = pca_results$x[,2], treatment = colData(dds)$treatment, cell.type = colData(dds)$cell.type)
  fpkm_pca_df$pca2 <- -fpkm_pca_df$pca2 

  # Create a new variable that combines treatment and cell type
  fpkm_pca_df$group <- paste(fpkm_pca_df$cell.type, fpkm_pca_df$treatment, sep = " ") 

  # create plot
  log_fpkm_tdp_plot <- ggplot(fpkm_pca_df, aes(x = fpkm, y = pca2, color = group)) +
    geom_point(size = 5) +
    ggtitle("") +
    ylab(paste0("PC2 (", round(summary(pca_results)$importance[2,2]*100), "%)")) +
    xlab("Tardbp FPKM") +
    scale_color_manual(values = c("#8a3838", "#f07f7e", "#51848a", "#99c0cc")) +
    nature_theme()

  return(log_fpkm_tdp_plot)
}


##########################################################
# function for differential gene expression calculations #
##########################################################
DGE_calculation <- function(dds) {
  # subset dds object
  dds_neuron <- dds[ , dds$cell.type == "NSC34"]
  dds_muscle <- dds[ , dds$cell.type == "C2C12"]

  # remove cell.type from experimental design
  dds_neuron@design <- ~ treatment
  dds_muscle@design <- ~ treatment

  # run DESeq2
  dds_neuron <- DESeq(dds_neuron)
  dds_muscle <- DESeq(dds_muscle)

  # saving the results
  neuronal_res <- results(dds_neuron, alpha = 0.05, contrast = c("treatment", "siTDP", "siLUC"))
  muscle_res <- results(dds_muscle, alpha = 0.05, contrast = c("treatment", "siTDP", "siLUC"))

  # filtering out genes with insignificat expression change 
  neuronal_res_sig <- subset(neuronal_res, padj < 0.05)
  muscle_res_sig <- subset(muscle_res, padj < 0.05)
  
  # Adding gene_id column
  neuronal_res_sig$gene_id <- rowData(dds)[rownames(neuronal_res_sig), "gene_id"]
  muscle_res_sig$gene_id <- rowData(dds)[rownames(muscle_res_sig), "gene_id"]
  
  return(list(neuronal_res_sig = neuronal_res_sig, muscle_res_sig = muscle_res_sig))
}


######################################
# function for plotting venn diagram #
######################################
DGE_venn <- function(data_list, title) {
  # extract the row names of the data frames
  transcript_sets <- lapply(data_list, rownames)
  
  # create a Venn diagram object using ggvenn
  venn <- Venn(transcript_sets)
  
  # process the data for plotting using ggvenn
  data <- process_data(venn)
  
  # define new set names and plot title
  new_set_names <- c("C2C12", "NSC34")
  plot_title <- title
  
  # create a new data frame for set labels with the new set names
  setlabel_data <- venn_setlabel(data)
  setlabel_data$name <- new_set_names
  
  # define new fill colors
  fill_colors <- c("#8a3838", "#6E5E61", "#51848a")
  
  # calculate percentages for region counts
  region_data <- venn_region(data)
  total_count <- sum(region_data$count)
  region_data$percentage <- paste0(round(region_data$count / total_count * 100, 1), "%")
  
  # plot the Venn diagram with the new settings
  DE_venn <- ggplot() +
    # 1. region count layer
    geom_sf(aes(fill = id), data = region_data, alpha = 0.5) +
    scale_fill_manual(values = fill_colors, guide = "none") +
    # 2. set edge layer
    geom_sf(data = venn_setedge(data), show.legend = FALSE, color = NA) +
    # 3. set label layer
    geom_sf_text(aes(label = name), data = setlabel_data, size = 6) +
    # 4. region label layer
    geom_sf_label(aes(label = paste0(count, "\n", percentage)), data = region_data, size = 5) +
    nature_theme() +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    ggtitle(plot_title) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 24))

  return(DE_venn)
}


##############################################################
# Function for plotting expression changes of common targets #
##############################################################
expression_changes_scatter <- function(dds) {
  # Convert the DESeqResults objects to data frames
  neuronal_res_sig_df <- as.data.frame(dds[1])
  muscle_res_sig_df <- as.data.frame(dds[2])

  # Merge the two data frames by row names
  merged_data <- merge(muscle_res_sig_df, neuronal_res_sig_df, by = "row.names")

  # Set the row names of the merged data frame to the gene names
  rownames(merged_data) <- merged_data$Row.names

  # Remove the Row.names column from the merged data frame
  merged_data$Row.names <- NULL
  # calculate Spearman's rank correlation coefficient and p-value
  cor.test <- cor.test(merged_data$muscle_res_sig.log2FoldChange, merged_data$ neuronal_res_sig.log2FoldChange, method = "spearman")
  rho <- cor.test$estimate
  p.value <- cor.test$p.value

  # Format the p-value
  if (p.value < 0.001) {
    sig_label <- "<0.001"
  } else if (p.value < 0.01) {
    sig_label <- "<0.01"
  } else if (p.value < 0.05) {
    sig_label <- "<0.05"
  } else {
    sig_label <- "ns"
  }

  # create the plot
  expression_changes <- ggplot(merged_data) +
    geom_point(aes(x = muscle_res_sig.log2FoldChange, y = neuronal_res_sig.log2FoldChange)) +
    geom_abline(intercept = 0, slope = 1, color = "grey") +
    geom_hline(yintercept = 0, color = "grey") +
    geom_vline(xintercept = 0, color = "grey") +
    geom_smooth(aes(x = muscle_res_sig.log2FoldChange, y = neuronal_res_sig.log2FoldChange), method = "lm", color = "#51848a") +
    xlab("log2 (fold change) in C2C12") +
    ylab("log2 (fold change) in NSC34") + 
    annotate("text", x = 1.9, y = -1.9, label = paste("Spearman's rho =", round(rho, 2), "\np-value =", sig_label), hjust = 1, size = 6) +
    nature_theme() +
    xlim(-2, 2) +
    ylim(-2, 2) +
    theme(axis.text.x = element_text(angle = 0))

  return(expression_changes)
}


##############################
# function for volcano plots #
##############################
volcano_plot <- function(dds, data_frame_1, data_frame_2, cell_names, colour) {
  # Extracting data frames from dds object
  df_1 <- as.data.frame(dds[[data_frame_1]])
  df_2 <- as.data.frame(dds[[data_frame_2]])

  # Add a column to indicate specificity
  df_1$specific <- !rownames(df_1) %in% rownames(df_2)

  # Create the volcano plot
  volcano <- ggplot(df_1,
                    aes(x = log2FoldChange, y = -log10(padj), color = specific)) +
    geom_point() +
    scale_color_manual(values = c("#515050",colour),
                      labels = c("common", paste0(cell_names, "-specific"))) +
    geom_vline(xintercept = c(log2(0.7), log2(1.3)),
              color = "grey") +
    geom_label_repel(data = subset(df_1,
                                    abs(log2FoldChange) > 1.5 & padj < 0.05),
                    aes(label = gene_id), size = 5, force = 10, box.padding = 0.5,
                    direction = "both", show.legend = FALSE) +
    xlab("log2(fold change)") +
    ylab("-log10 (padj)") +
    ggtitle(cell_names) +
    xlim(-3, 3) +
    ylim(0, 100) +
    nature_theme() +
    theme(axis.text.x = element_text(angle = 0),
          plot.title = element_text(hjust = 0.5, size = 24),
          legend.text = element_text(size = 20)) +
    guides(color = guide_legend(override.aes = list(shape = 15, size = 5)))
  
  # Display the plot
  return(volcano)
}

############################
# get RBPs for DE analysis #
############################
get_RBPs <- function(path, fpmks){
  RBPs_IDs <- read.csv(path)
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  RBPs_IDs$ensembl_ids <- getBM(attributes='ensembl_gene_id', filters = 'mgi_symbol', values = RBPs_IDs$GeneID, mart = ensembl)$ensembl_gene_id
  fpkm_values <- as.data.frame(fpkm_values)
  fpkm_values <- fpkm_values %>%
                  rownames_to_column(var = "ensembl_ids")
  fpmk_RBPs <- left_join(RBPs_IDs, fpkm_values, by = "ensembl_ids") 
  fpmk_RBPs[, 3:ncol(fpmk_RBPs)] <- fpmk_RBPs[, 3:ncol(fpmk_RBPs)] + 1

  fpmk_RBPs <- fpmk_RBPs %>%
    pivot_longer(cols = -c(GeneID, ensembl_ids), names_to = "sample", values_to = "expression") %>%
    mutate(cell_type = ifelse(sample %in% colnames(fpmk_RBPs)[1:6], "C2C12", "NSC34")) %>%
    group_by(GeneID, ensembl_ids, cell_type) %>%
    summarize(avg_expression = mean(expression))

  return(fpmk_RBPs)
}

########################
# plot boxplot of RBPs #
########################
RBPs_boxplot <- function(fpmk_RBPs) {
  wilcox_result <- wilcox.test(avg_expression ~ cell_type, data = fpmk_RBPs)
  p_value <- wilcox_result$p.value

  # Create a significance label using asterisks
  if (p_value < 0.001) {
    sig_label <- "***"
  } else if (p_value < 0.01) {
    sig_label <- "**"
  } else if (p_value < 0.05) {
    sig_label <- "*"
  } else {
    sig_label <- "ns"
  }

  RBPs_box <- ggplot(fpmk_RBPs, aes(x = cell_type, y = log10(avg_expression), fill = cell_type)) +
    geom_boxplot(notch = TRUE) +
    scale_fill_manual(values = c("C2C12" = "#8a3838", "NSC34" = "#51848a")) +
    labs(x = "", y = "Log10(FPKM + 1)") +
    ylim(0, 3) +
    nature_theme() +
    theme(legend.position = "none") + 
    annotate("text", x = 1.5, y = 3, label = sig_label, size = 6) +
    geom_segment(aes(x = 1.2, y = 2.9, xend = 1.8, yend = 2.9))
  
  return(RBPs_box)
}


###################################
# Create the scatter plot of RBPs #
###################################
RBP_scatter_plotting <- function(RBPs) {
  
  RBPs <- tidyr::spread(RBPs, cell_type, avg_expression)
  RBPs$cell_type <- ifelse(RBPs$C2C12 > RBPs$NSC34, "C2C12", "NSC34")

  cor.test.result <- cor.test(log10(RBPs$C2C12), log10(RBPs$NSC34), method = "spearman")
  rho <- cor.test.result$estimate
  p.value <- cor.test.result$p.value

    # Format the p-value
  if (p.value < 0.001) {
    sig_label <- "<0.001"
  } else if (p.value < 0.01) {
    sig_label <- "<0.01"
  } else if (p.value < 0.05) {
    sig_label <- "<0.05"
  } else {
    sig_label <- "ns"
  }

  RBP_scatter <- ggplot(RBPs, aes(x = log10(C2C12), y = log10(NSC34))) +
    geom_point(aes(color = cell_type, size = 5), show.legend = FALSE) +
    geom_abline(intercept = 0, slope = 1, color = "grey") +
    scale_color_manual(values = c("C2C12" = "#8a3838", "NSC34" = "#51848a")) +
    labs(x = "Log10(FPKM + 1) in C2C12", y = "Log10(FPKM + 1) in NSC34") +
    geom_label_repel(aes(label = GeneID, color = cell_type), size = 5, force = 10, box.padding = 0.5, direction = "both", show.legend = FALSE) +
    annotate("text", x = 2, y = 0, label = paste("Spearman's rho =", round(rho, 2), "\np-value =", sig_label), hjust = 0, vjust = 0, size = 6) +
    theme() +
    nature_theme()

  return(RBP_scatter)
}



########################
########################
# Running the analysis #
########################
########################

# save gene length info for later fpmk calculcations
gene_lengths <- raw_counts["Length"]

counts <- prepare_counts(raw_counts)
sample_data <- prepare_sample_data(meta_data)

# create experimental design variable
experimental_design_formula <- as.formula("~ cell.type + treatment")

# DESeqDataSet object generation
dds <- generate_dds(counts, sample_data, experimental_design_formula, gene_lengths)

# calculate fpkm
fpkm_values <- fpkm(dds)

# Perform PCA on the FPKM values
pca_res <- prcomp(t(fpkm_values), scale = TRUE)

# tdp-43 expression based on fpkm in all conditions
tdp_expressions <- gene_expressions(dds, fpkm_values, "ENSMUSG00000041459")

# PCA analysis
pca_plot <- pca_plotting(pca_res)

# pca2 vs tdp-43 expression
pc2_tdp <- pc2_vs_gene(dds, fpkm_values, pca_res, "ENSMUSG00000041459")

# DGE of all the genes
DE_full <- DGE_calculation(dds)
DE_full_venn <- DGE_venn(DE_full, "DGE")

# DGE of genes with fpkm > 0.5 in each sample
keep <- apply(fpkm_values, MARGIN = 1, FUN = function(x) all(x > 0.5))
dds_filtered <- dds[keep,]

DE_filtered <- DGE_calculation(dds_filtered)
DE_filtered_venn <- DGE_venn(DE_filtered, "FPKM > 0.5")


# expression changes between cell lines
expression_changes_plot <- expression_changes_scatter(DE_full)


# volcano plots of DGE in each cell line
neuronal_volcano <- volcano_plot(DE_full, "neuronal_res_sig", "muscle_res_sig", "NSC34", "#51848a")
muscle_volcano <- volcano_plot(DE_full, "muscle_res_sig", "neuronal_res_sig", "C2C12", "#8a3838")

# box plot
RBPs <- get_RBPs(RBPs_path, fpkm_values)
RBPs_box <- RBPs_boxplot(RBPs)

# RBPs scatterplot
RBPs_scatter_plot <- RBP_scatter_plotting(RBPs)

# DE RBPs venn 
neurons_DE <- as.data.frame(DE_full[[1]])
muscle_DE <- as.data.frame(DE_full[[2]])
RBPs_neurons_DE <- neurons_DE[rownames(neurons_DE) %in% RBPs$gene, ]
RBPs_muscle_DE <- muscle_DE[rownames(muscle_DE) %in% RBPs$gene, ]
RBPs_list <- list(RBPs_muscle_DE, RBPs_neurons_DE)

RBP_venn <- DGE_venn(RBPs_list, "DE RBPs") # add list of the genes



######################
######################
# Saving the results #
######################
######################

# Save all plots to disk in the PNG format
ggsave(filename = paste0(output_folder, "/tdp_expressions.png"), plot = tdp_expressions, width = 20, height = 15)
ggsave(filename = paste0(output_folder, "/pca_plot.png"), plot = pca_plot, width = 20, height = 15)
ggsave(filename = paste0(output_folder, "/pc2_tdp.png"), plot = pc2_tdp, width = 20, height = 15)
ggsave(filename = paste0(output_folder, "/DE_full_venn.png"), plot = DE_full_venn, width = 20, height = 15)
ggsave(filename = paste0(output_folder, "/DE_filtered_venn.png"), plot = DE_filtered_venn, width = 20, height = 15)
ggsave(filename = paste0(output_folder, "/expression_changes_plot.png"), plot = expression_changes_plot, width = 20, height = 15)
ggsave(filename = paste0(output_folder, "/neuronal_volcano.png"), plot = neuronal_volcano, width = 20, height = 15)
ggsave(filename = paste0(output_folder, "/muscle_volcano.png"), plot = muscle_volcano, width = 20, height = 15)
ggsave(filename = paste0(output_folder, "/RBPs_box.png"), plot = RBPs_box, width = 20, height = 15)
ggsave(filename = paste0(output_folder, "/RBPs_scatter_plot.png"), plot = RBPs_scatter_plot, width = 20, height = 15)
ggsave(filename = paste0(output_folder, "/RBP_venn.png"), plot = RBP_venn, width = 20, height = 15)

# Save DGE results
neurons_DE$ensembl_id <- rownames(neurons_DE)
write.csv(neurons_DE, file = paste0(output_folder, "/neuron_DE.csv"), row.names = FALSE)

muscle_DE$ensembl_id <- rownames(muscle_DE)
write.csv(muscle_DE, file = paste0(output_folder, "/muscle_DE.csv"), row.names = FALSE)


# save DESeqDataSet and fpkm values
saveRDS(dds, file = paste0(output_folder, "/DESeq2_output.rds"))
write.csv(fpkm_values, file = paste0(output_folder, "/fpkm_values.csv"), row.names = TRUE)

# save differentially expressed RBPs
write.csv(RBPs, file = paste0(output_folder, "/RBPs_DE.csv"), row.names = FALSE)
