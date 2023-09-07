library(ggplot2)
library(ggVennDiagram)
library(dplyr)
library(tidyr)

########################
# Inputs/ output paths # 
########################
args <- commandArgs(trailingOnly = TRUE)
neuron_as_path <- "./results/rMATS_output/rMATS_neuron/" # args[1]
muscle_as_path <- "./results/rMATS_output/rMATS_muscle/" # args[2]
output_dir <- "./results/AS_analysis_output/" # args[3]

AS_events_types = c("SE", "MXE", "RI", "A3SS", "A5SS")


##############################
# Functions used in analysis # 
##############################
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

splicing_events <- function(cell_path, AS) {

    path <- paste0(cell_path, AS, ".MATS.JCEC.txt")
    file <- read.csv(path, sep = "\t")
    file <- filter(file, FDR < 0.01)
    return(file)
}


get_transcripts <- function(df, col_names) {
    df %>%
        unite(col = "transcript", col_names, sep = "_") %>%
        pull(transcript)
}


horizontal_bar_chart_data <- function(cell_line_1, cell_line_2, i, cell_line_1_name, cell_line_2_name) {

    cell_line_1_length <- length(cell_line_1$GeneID)
    cell_line_2_length <- length(cell_line_2$GeneID)
    common_count <- length(intersect(cell_line_1$GeneID, cell_line_2$GeneID))

    # Create a data frame for the bar chart
    data <- data.frame(
        Category = c(cell_line_1_name, "common", cell_line_2_name),
        Count = c(cell_line_1_length, common_count, cell_line_2_length),
        Splice_type = i)
    
    data <- data %>%
        group_by(Splice_type) %>%
        mutate(Percentage_common = round(sum(Count[Category == 'common']) / sum(Count[Category != 'common']) * 100, 1))
    return(data)
}


horizontal_bar_plotting <- function(plotting_data) {
    ggplot(plotting_data, aes(x = Splice_type, y = Count, fill = Category)) +
        geom_bar(stat = "identity", width = 0.9, color = "black", linewidth = 0.5) +
        coord_flip() +
        xlab("") +
        ylab("No. of events") +
        ylim(0,2000) +
        nature_theme() +
        scale_fill_manual(values = c("#688b68", "#89a164", "#cb940a")) +
        geom_text(data = plotting_data[!duplicated(plotting_data$Splice_type),], 
                    aes(label = paste0(Percentage_common,"%")), size = 7, hjust = -2, nudge_y = 600)
}

AS_venn <- function(transcript_sets, title) {
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
    fill_colors <- c("#B4C5B4", "#A6AA77", "#E6CA84")
    
    # calculate percentages for region counts
    region_data <- venn_region(data)
    total_count <- sum(region_data$count)
    region_data$percentage <- paste0(round(region_data$count / total_count * 100, 1), "%")
    
    # plot the Venn diagram with the new settings
    dge_venn <- ggplot() +
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
    return(dge_venn)
}


############
# Analysis # 
############
# filtering rMATS results by FDR <= 0.01
neuron_df_list <- list()
muscle_df_list <- list()

for (i in AS_events_types){
    tmp <- splicing_events(neuron_as_path, i)
    neuron_df_list <- c(neuron_df_list, list(tmp))
}
names(neuron_df_list) <- c("SE_neuron", "MXE_neuron", "RI_neuron", "A3SS_neuron", "A5SS_neuron")

for (i in AS_events_types){
    tmp <- splicing_events(muscle_as_path, i)
    muscle_df_list <- c(muscle_df_list, list(tmp))
}
names(muscle_df_list) <- c("SE_muscle", "MXE_muscle", "RI_muscle", "A3SS_muscle", "A5SS_muscle")



## VENN DIAGRAM OF TRANSCRIPTS ##
SE_columns <- c("GeneID", "exonStart_0base",  "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")
SE_neuron_transcript <- get_transcripts(neuron_df_list$SE_neuron, SE_columns)
SE_muscle_transcript <- get_transcripts(muscle_df_list$SE_muscle, SE_columns)

MXE_columns <- c("GeneID", "X1stExonStart_0base", "X1stExonEnd", "X2ndExonStart_0base", "X2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")
MXE_neuron_transcript <- get_transcripts(neuron_df_list$MXE_neuron, MXE_columns)
MXE_muscle_transcript <- get_transcripts(muscle_df_list$MXE_muscle, MXE_columns)

RI_columns <- c("GeneID", "riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")
RI_neuron_transcript <- get_transcripts(neuron_df_list$RI_neuron, RI_columns)
RI_muscle_transcript <- get_transcripts(muscle_df_list$RI_muscle, RI_columns)

A3SS_columns <- c("GeneID", "longExonStart_0base", "longExonEnd",   "shortES",   "shortEE", "flankingES", "flankingEE")
A3SS_neuron_transcript <- get_transcripts(neuron_df_list$A3SS_neuron, A3SS_columns)
A3SS_muscle_transcript <- get_transcripts(muscle_df_list$A3SS_muscle, A3SS_columns)

A5SS_columns <- c("GeneID", "longExonStart_0base", "longExonEnd",   "shortES",   "shortEE", "flankingES", "flankingEE")
A5SS_neuron_transcript <- get_transcripts(neuron_df_list$A5SS_neuron, A5SS_columns)
A5SS_muscle_transcript <- get_transcripts(muscle_df_list$A5SS_muscle, A5SS_columns)

# preparing list of alternatively spliced transcripts 
AS_transcripts_neuron <- as.vector(c(SE_neuron_transcript, MXE_neuron_transcript, RI_neuron_transcript, A3SS_neuron_transcript, A5SS_neuron_transcript))
AS_transcript_muscle <- as.vector(c(SE_muscle_transcript, MXE_muscle_transcript, RI_muscle_transcript, A3SS_muscle_transcript, A5SS_muscle_transcript))
AS_transcripts <- list(AS_transcript_muscle, AS_transcripts_neuron)

#  plotting venn diagram of alternatively spliced genes
AS_transcripts_venn <- AS_venn(AS_transcripts, "AS events rMATS")



## VENN DIAGRAM OF GENES ##
# preparing list of alternatively spliced genes
AS_genes_neurons <- c()
for (i in neuron_df_list){
    tmp <- i$GeneID
    AS_genes_neurons <- c(AS_genes_neurons, tmp)
}

AS_genes_muscle <- c()
for (i in muscle_df_list){
    tmp <- i$GeneID
    AS_genes_muscle <- c(AS_genes_muscle, tmp)
}

AS_genes <- list(AS_genes_muscle, AS_genes_neurons)

#  plotting venn diagram of alternatively spliced genes
AS_genes_venn <- AS_venn(AS_genes, "AS Genes")



## BAR CHART OF TRANSCRIPTS ##
# generating data frame for horizontal bar of different alternative splicing events types
plotting_df <- data.frame()
for (k in seq_along(AS_events_types)) {
    i <- AS_events_types[k]
    tmp_df <- horizontal_bar_chart_data(neuron_df_list[[k]], muscle_df_list[[k]], i, "C2C12-specific",  "NSC34-specific")
    plotting_df <- rbind(plotting_df, tmp_df)
}

plotting_df$Splice_type <- factor(plotting_df$Splice_type, levels = rev(AS_events_types))
plotting_df$Category <- factor(plotting_df$Category, levels = c("C2C12-specific", "common", "NSC34-specific"))

# plots splicing event counts 
bar_plot <- horizontal_bar_plotting(plotting_df)


##################
# Saving results #
##################
# create directory for saving results
dir.create(output_dir)

# saving neurons data frames

for (name in names(neuron_df_list)) {
    write.csv(neuron_df_list[[name]], file = paste0(output_dir, name, ".csv"), row.names = FALSE)
}

# saving muscle cell data frames
for (name in names(muscle_df_list)) {
    write.csv(muscle_df_list[[name]], file = paste0(output_dir, name, ".csv"), row.names = FALSE)
}

# saving figures 
ggsave(filename = paste0(output_dir, "splicing_horizontal_bar_plot.png"), plot = bar_plot)
ggsave(filename = paste0(output_dir, "AS_transcripts_venn.png"), plot = AS_transcripts_venn)
ggsave(filename = paste0(output_dir, "AS_genes_venn.png"), plot = AS_genes_venn)
