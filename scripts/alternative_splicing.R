library(ggplot2)
library(ggVennDiagram)
library(dplyr)

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
    file <- filter(file, FDR <= 0.01)

    return(file)
}

bar_chart_data <- function(cell_line_1, cell_line_2, i, cell_line_1_name, cell_line_2_name) {

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

horizontal_plot <- function(plotting_data) {
    ggplot(plotting_data, aes(x = Splice_type, y = Count, fill = Category)) +
        geom_bar(stat = "identity", width = 0.9, color = "black", linewidth = 0.5) +
        coord_flip() +
        xlab("") +
        ylab("No. of events") +
        ylim(0, 1200) +
        nature_theme() +
        scale_fill_manual(values = c("#688b68", "#89a164", "#cb940a")) +
        geom_text(data = plotting_data[!duplicated(plotting_data$Splice_type),], 
                    aes(label = paste0(Percentage_common,"%")), size = 7, hjust = -4)

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
filtered_df_list <- list()
for (i in AS_events_types){
    neuron <- splicing_events(neuron_as_path, i)
    muscle <- splicing_events(muscle_as_path, i)
    tmp_list <- c(list(muscle), list(neuron))
    filtered_df_list <- c(filtered_df_list, tmp_list)
}


# create df for plotting horizontal bar chart
plotting_df <- data.frame()
for (k in seq_along(AS_events_types)) {
    i <- AS_events_types[k]
    j <- seq(1, length(filtered_df_list), by = 2)[k]
    tmp_df <- bar_chart_data(filtered_df_list[[j]], filtered_df_list[[j + 1]], i, "C2C12-specific",  "NSC34-specific")
    plotting_df <- rbind(plotting_df, tmp_df)
}
plotting_df$Splice_type <- factor(plotting_df$Splice_type, levels = rev(AS_events_types))
plotting_df$Category <- factor(plotting_df$Category, levels = c("C2C12-specific", "common", "NSC34-specific"))


# plots splicing event counts and save it into output folder
bar_plot <- horizontal_plot(plotting_df)



#this goes into a function to prepare data for venn diagram
muscle_gene_ids <- c()
for (i in seq(1, length(filtered_df_list), by = 2)) {
    muscle_gene_ids <- c(muscle_gene_ids, filtered_df_list[[i]]$GeneID)
}

neuron_gene_ids <- c()
for (i in seq(2, length(filtered_df_list), by = 2)) {
    neuron_gene_ids <- c(neuron_gene_ids, filtered_df_list[[i]]$GeneID)
}


transcript_sets <- list(muscle_gene_ids, neuron_gene_ids)
venn_diagram_rMATS <- AS_venn(transcript_sets, "AS events rMATS")
venn_diagram_rMATS


# create directory for saving results
dir.create(output_dir)

# saving figures 
ggsave(filename = paste0(output_dir, "splicing_horizontal_bar_plot.png"), plot = bar_plot)
ggsave(filename = paste0(output_dir, "venn_diagram_rMATS.png"), plot = venn_diagram_rMATS)

