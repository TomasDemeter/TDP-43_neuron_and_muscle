library(clusterProfiler)
library(ggplot2)
library(ggVennDiagram)
library(EnsDb.Mmusculus.v79)
library(org.Mm.eg.db)
library(enrichplot)
library(dplyr)
library(tidyr)
library(stringr)


DE_path <- "./results/DESeq2_output/"
AS_path <- "./results/DESeq2_output/"
output_folder <- "./results/GO_term_analysis/"
####################
# Loading the data #
####################
args <- commandArgs(trailingOnly = TRUE)
DE_path <-  paste0(args[1], "/")
AS_path <- paste0(args[2], "/")
output_folder <- paste0(args[3], "/") 


neuron_genes <- read.csv(paste0(DE_path, "neuron_DE.csv"))
muscle_genes <- read.csv(paste0(DE_path, "muscle_DE.csv"))

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
            axis.line = element_line(linewidth  = 1, color = "black"),
            axis.ticks = element_line(linewidth  = 1, color = "black"),
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


########################################
# running enrichment GO-terms analysis #
########################################

go_term_analysis <- function(condition_df) {
    ensembl_ids <- condition_df$ensembl_id
    ego <- enrichGO(gene = ensembl_ids,
                    keyType = "ENSEMBL",
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05)
    return(ego)
}


###################################################################
# generating venn diagram of enriched GO-terms in both cell lines #
###################################################################

EGO_venn <- function(ego1, ego2, ego1_name, ego2_name, title) {

    ego_merged <- list(ego1$ID, ego2$ID)

    # create a Venn diagram object using ggvenn
    venn <- Venn(ego_merged)

    # process the data for plotting using ggvenn
    data <- process_data(venn)

    # define new set names and plot title
    new_set_names <- c(ego1_name, ego2_name)
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


#################################################
# plotting erichment GO-terms in each cell line #
#################################################

EGO_terms <- function (ego1, ego2, ego1_name, ego2_name) {
    # Combine the description and counts into a data frame
    ego1_df <- data.frame(description = ego1@result$Description, counts = ego1@result$Count)
    ego2_df <- data.frame(description = ego2@result$Description, counts = ego2@result$Count)

    combined_df <-  inner_join(ego1_df, ego2_df, by = "description", suffix = c(paste0("_", ego1_name), paste0("_", ego2_name)))

    # Use the [[ operator to access the columns by name
    combined_df$total_counts <-  combined_df[[paste0("counts_", ego1_name)]] + combined_df[[paste0("counts_", ego2_name)]]

    top15_combined_df <- head(combined_df[order(-combined_df$total_counts),], 15)

    # Use the paste0() function to create the column names
    df_long <- top15_combined_df %>%
        pivot_longer(cols = c(paste0("counts_", ego1_name), paste0("counts_", ego2_name)), names_to = "Count_Type", values_to = "Counts")

    # Create the horizontal bar chart
    EGO_terms_plot <- ggplot(df_long, aes(x = Counts, y = description, fill = Count_Type)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(x = "Gene no.", y = "") +
        scale_fill_manual(values = setNames(c("#51848a", "#8a3838"), 
                        c(paste0("counts_", ego1_name), paste0("counts_", ego2_name)))) +
        nature_theme() +
        theme(legend.position = "none")
        
    return(EGO_terms_plot)
}
#########################################################
# preparing data with ensembl ids for GO terms download #
#########################################################
ensembl_IDs_AS <- function(path, pattern) {
    files <- list.files(path = path, pattern = paste0(pattern, "$"))
    full_file_paths <- file.path(path, files)
    geneid_data_list <- lapply(full_file_paths, function(x) read.csv(x)["GeneID"])
    data_combined <- do.call(rbind, geneid_data_list)
    colnames(data_combined)[1] <- "ensembl_id"
    return(data_combined)
}
##########################################
# extracting GO terms for df preparation #
##########################################
get_go_terms <- function(cell_type, keywords, group) {
    keywords_pattern <- paste0(keywords, collapse = "|")
    go_terms <- cell_type %>%
                dplyr::filter(str_detect(Description, keywords_pattern)) %>%
                dplyr::filter(-log(p.adjust, base = 10) > 1) %>%
                dplyr::select(Description, p.adjust) %>%
                mutate(GO_group = group)
    return(go_terms)
}
#############################################################
# preparing data frame for GO terms horizontal bar plotting #
#############################################################
get_go_df <- function(cell_type, keywords_list) {
    combined_df <- data.frame()
    for (i in seq_along(keywords_list)) {
        tmp_keywords <- get_go_terms(cell_type, keywords_list[[i]], names(keywords_list)[[i]])
        combined_df <- rbind(combined_df, tmp_keywords)
    }
    
    description_order <- as.vector(combined_df$Description)
    GO_order <- unique(as.vector(combined_df$GO_group))
    
    combined_df$Description <- factor(combined_df$Description, levels = rev(description_order))
    combined_df$GO_group <- factor(combined_df$GO_group, levels = GO_order)
    
    return(combined_df)
}
########################################################################
# generating horizontal bar of GO terms of alternatively spliced genes #
########################################################################
unique_GO_terms_bar <- function(go_terms_dataframe) {
    # Create a named vector of colors
    group_colors <- c("neuron" = "#cb950e", "RNA metabolism" = "#496d88", "DNA" = "#c8919d", "muscle" = "#698969")
    
    # Create the plot
    plot <- ggplot(go_terms_dataframe, aes(x = Description, y = -log(p.adjust, base = 10), fill = GO_group)) +
        geom_bar(stat = "identity") +
        xlab("") +
        ylab("-log10(p.adj)") +
        scale_fill_manual(values = group_colors) + 
        coord_flip() +
        nature_theme()
    
    return(plot)
}


########################
########################
# Running the analysis #
########################
########################

# generating results
DE_go_neuron <- go_term_analysis(neuron_genes)
DE_go_muscle <- go_term_analysis(muscle_genes)


# plotting the resutls
ego_venn_plot <- EGO_venn(DE_go_muscle, DE_go_neuron, "C2C12", "NSC34", "GO terms")
ego_terms_plot <- EGO_terms(DE_go_muscle, DE_go_neuron, "C2C12", "NSC34") 


# getting files from alternative splicing analysis
neuron_data_combined <- ensembl_IDs_AS(AS_path, "_neuron.csv")
muscle_data_combined <- ensembl_IDs_AS(AS_path, "_muscle.csv")


# getting GO terms of alternatively spliced genes
AS_go_neuron <- go_term_analysis(neuron_data_combined)
AS_go_muscle <- go_term_analysis(muscle_data_combined)


# plotting GO terms of alternatively spliced genes
AS_go_venn_plot <- EGO_venn(AS_go_muscle, AS_go_neuron, "C2C12", "NSC34", "GO terms of AS genes")


# extracting cell type specific GO terms
AS_go_neuron_unique <- anti_join(as.data.frame(AS_go_neuron), as.data.frame(AS_go_muscle), by = "Description")
AS_go_muscle_unique <- anti_join(as.data.frame(AS_go_muscle), as.data.frame(AS_go_neuron), by = "Description")


# Filter rows in both data frames where p.adjust > -log10(p.adjust)
AS_go_neuron_filtered <- AS_go_neuron %>% filter(-log(p.adjust, base = 10) > 1)
AS_go_muscle_filtered <- AS_go_muscle %>% filter(-log(p.adjust, base = 10) > 1)


# Merge the filtered data frames
common_go <- inner_join(as.data.frame(AS_go_neuron_filtered), as.data.frame(AS_go_muscle_filtered), by = "Description")
common_go <- common_go %>% mutate(p.adjust = (p.adjust.x + p.adjust.y) / 2)


# defining keywords for GO term clustering and saving them all in a list
neuron_keywords <- list("neuron", "synapse", "synaptic", "dendritic", "axonogenesis", "neurogenesis", "postsynapse", "nervous", "neurotransmitter")
RNA_keywords <- list("RNA", "tRNA", "ncRNA", "mRNA", "ribonucleoprotein")
DNA_keywords <- list("DNA", "telomerase", "chromosome", "histone", "chromatin", "telomere")
muscle_keywords <- list("muscle")

keywords_list <- list(neuron_keywords,
                    RNA_keywords,
                    DNA_keywords,
                    muscle_keywords)

names(keywords_list) <- c("neuron", "RNA metabolism", "DNA", "muscle")


# preparing data frames for plotting of a horizontal bar chart of GO terms
neuronal_go_df <- get_go_df(AS_go_neuron_unique, keywords_list)
muscle_go_df <- get_go_df(AS_go_muscle_unique, keywords_list)
common_go_df <-  get_go_df(common_go, keywords_list)


# for neuronal cells I am plotting only top 15 of each category otherwise wont fit on the plot
top15_neuronal_go_df <- neuronal_go_df %>%
    group_by(GO_group) %>%
    top_n(n = 15, wt = -log(p.adjust, base = 10))


# plotting horizontal bars
neuron_AS_GO_plot <- unique_GO_terms_bar(top15_neuronal_go_df)
muscle_AS_GO_plot <- unique_GO_terms_bar(muscle_go_df)    
common_AS_GO_plot <- unique_GO_terms_bar(common_go_df)



######################
######################
# Saving the reuslts #
######################
######################

# saving the figures
ggsave(filename = paste0(output_folder, "/ego_venn_plot.png"), plot = ego_venn_plot, width = 20, height = 15)
ggsave(filename = paste0(output_folder, "/ego_terms_plot.png"), plot = ego_terms_plot, width = 20, height = 15)
ggsave(filename = paste0(output_folder, "/AS_go_venn_plot.png"), plot = AS_go_venn_plot, width = 20, height = 15)
ggsave(filename = paste0(output_folder, "/neuron_AS_GO_plot.png"), plot = neuron_AS_GO_plot, width = 20, height = 15)
ggsave(filename = paste0(output_folder, "/muscle_AS_GO_plot.png"), plot = muscle_AS_GO_plot, width = 20, height = 15)
ggsave(filename = paste0(output_folder, "/common_AS_GO_plot.png"), plot = common_AS_GO_plot, width = 20, height = 15)


# saving enriched GO-terms analysis
write.csv(DE_go_neuron, file = paste0(output_folder, "/DE_go_neuron.csv"), row.names = FALSE)
write.csv(DE_go_muscle, file = paste0(output_folder, "/DE_go_muscle.csv"), row.names = FALSE)
write.csv(AS_go_neuron, file = paste0(output_folder, "/AS_go_neuron.csv"), row.names = FALSE)
write.csv(AS_go_muscle, file = paste0(output_folder, "/AS_go_muscle.csv"), row.names = FALSE)
