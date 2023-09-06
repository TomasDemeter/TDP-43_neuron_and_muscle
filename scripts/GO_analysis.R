library(clusterProfiler)
library(ggplot2)
library(ggVennDiagram)
library(EnsDb.Mmusculus.v79)
library(org.Mm.eg.db)
library(enrichplot)
library(dplyr)
library(tidyr)
library(stringr)

#input_path <- "results/DESeq2_output/"
#output_folder <- "results/GO_term_analysis"
neuron_genes_filepath <- args[1]
muscle_genes_filepath <- args[2]
output_folder <- args[3] #"results/GO_term_analysis"

neuron_genes <- read.csv(neuron_genes_filepath)
muscle_genes <- read.csv(muscle_genes_filepath)


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
                    pvalueCutoff = 0.05,
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
    dge_venn <- ggplot() +
        # 1. region count layer
        geom_sf(aes(fill = id), data = region_data, alpha = 0.5) +
        scale_fill_manual(values = fill_colors, guide = FALSE) +
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

    return(dge_venn)
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


########################
########################
# Running the analysis #
########################
########################

# generating results
ego_neuron <- go_term_analysis(neuron_genes)
ego_muscle <- go_term_analysis(muscle_genes)

# plotting the resutls
ego_venn_plot <- EGO_venn(ego_muscle, ego_neuron, "C2C12", "NSC34", "GO terms")
ego_terms_plot <- EGO_terms(ego_muscle, ego_neuron, "C2C12", "NSC34") # different pathways compared to the paper but still related to neurodegeneration


######################
######################
# Saving the reuslts #
######################
######################

# saving the figures
ggsave(filename = paste0(output_folder, "/ego_venn_plot.png"), plot = ego_venn_plot)
ggsave(filename = paste0(output_folder, "/ego_terms_plot.png"), plot = ego_terms_plot)

# saving enriched GO-terms analysis
write.csv(ego_neuron, file = paste0(output_folder, "/ego_neuron.csv"), row.names = FALSE)
write.csv(ego_muscle, file = paste0(output_folder, "/ego_muscle.csv"), row.names = FALSE)




AS_path <- "./results/AS_analysis_output"
neuron_files <- list.files(path = AS_path, pattern = "_neuron.csv$")
full_neuron_file_paths <- file.path(AS_path, neuron_files)
neuron_geneid_data_list <- lapply(full_neuron_file_paths, function(x) read.csv(x)["GeneID"])
neuron_data_combined <- do.call(rbind, neuron_geneid_data_list)
colnames(neuron_data_combined)[1] <- "ensembl_id"

muscle_files <- list.files(path = AS_path, pattern = "_muscle.csv$")
full_muscle_file_paths <- file.path(AS_path, muscle_files)
muscle_geneid_data_list <- lapply(full_muscle_file_paths, function(x) read.csv(x)["GeneID"])
muscle_data_combined <- do.call(rbind, muscle_geneid_data_list)
colnames(muscle_data_combined)[1] <- "ensembl_id"




AS_ego_neuron <- go_term_analysis(neuron_data_combined)
AS_ego_muscle <- go_term_analysis(muscle_data_combined)

ego_venn_plot <- EGO_venn(AS_ego_muscle, AS_ego_neuron, "C2C12", "NSC34", "GO terms of AS geens")

AS_ego_neuron_unique <- anti_join(as.data.frame(AS_ego_neuron), as.data.frame(AS_ego_muscle), by = "Description")
AS_ego_muscle_unique <- anti_join(as.data.frame(AS_ego_muscle), as.data.frame(AS_ego_neuron), by = "Description")


EGO_terms_plot <- ggplot(AS_ego_neuron_unique, aes(x = Counts, y = Description, fill = Count_Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Gene no.", y = "") +
    nature_theme() +
    theme(legend.position = "none")




neuron_keywords <- list("neuron", "synapse", "synaptic", "dendritic", "axonogenesis")
neuron_keywords_pattern <- paste0(neuron_keywords, collapse = "|")
neuron_go_terms <- AS_ego_neuron_unique %>%
                    dplyr::filter(str_detect(Description, neuron_keywords_pattern)) %>%
                    dplyr::filter(-log(p.adjust, base = 10) > 2) %>%
                    dplyr::select(Description, p.adjust) %>%
                    mutate(GO_group = "neuron")



RNA_keywords <- list("RNA", "tRNA", "ncRNA", "mRNA", "ribonucleoprotein")
RNA_keywords_pattern <- paste0(RNA_keywords, collapse = "|")
RNA_go_terms <- AS_ego_neuron_unique %>%
                    dplyr::filter(str_detect(Description, RNA_keywords_pattern)) %>%
                    dplyr::filter(-log(p.adjust, base = 10) > 1) %>%
                    dplyr::select(Description, p.adjust) %>%
                    mutate(GO_group = "RNA metabolism")



DNA_keywords <- list("DNA")
DNA_keywords_pattern <- paste0(DNA_keywords, collapse = "|")
DNA_go_terms <- AS_ego_neuron_unique %>%
                    dplyr::filter(str_detect(Description, DNA_keywords_pattern)) %>%
                    dplyr::filter(-log(p.adjust, base = 10) > 1) %>%
                    dplyr::select(Description, p.adjust) %>%
                    mutate(GO_group = "DNA")


neuron_go_df <- rbind(neuron_go_terms, RNA_go_terms, DNA_go_terms)

AS_ego_neuron_unique$Description
