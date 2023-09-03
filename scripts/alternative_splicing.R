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


splicing_events <- function(muscle_path, neuron_path, AS) {

    muscle_JCEC_path <- paste0(muscle_path, AS, ".MATS.JCEC.txt")
    neuron_JCEC_path <- paste0(neuron_path, AS, ".MATS.JCEC.txt")
    muscle_JCEC <- read.csv(muscle_JCEC_path, sep = "\t")
    neuron_JCEC <- read.csv(neuron_JCEC_path, sep = "\t")
    muscle_JCEC <- filter(muscle_JCEC, FDR <= 0.01)
    neuron_JCEC <- filter(neuron_JCEC, FDR <= 0.01)

    return(list(muscle_JCEC, neuron_JCEC))
}

bar_chart_data <- function(results, i) {

    muscle_count <- length(results[[1]]$GeneID)
    neuron_count <- length(results[[2]]$GeneID)
    common_count <- length(intersect(results[[1]]$GeneID, results[[2]]$GeneID))

    # Create a data frame for the bar chart
    data <- data.frame(
        Category = c("C2C12-specific", "common", "NSC34-specific"),
        Count = c(muscle_count, common_count, neuron_count),
        Splice_type = i)

    return(data)
}

horizontal_plot <- function(final_plotting_data) {
    ggplot(final_plotting_data, aes(x = Splice_type, y = Count, fill = Category)) +
        geom_bar(stat = "identity", width = 0.9, color = "black", size = 0.5) +
        coord_flip() +
        xlab("") +
        ylab("No. of events") +
        ylim(0, 1200) +
        nature_theme() +
        scale_fill_manual(values = c("#688b68", "#89a164", "#cb940a")) +
        geom_text(data = final_plotting_data[!duplicated(final_plotting_data$Splice_type),], 
                    aes(label = paste0(Percentage_common,"%")), size = 7, hjust = -4)

}


############
# Analysis # 
############
# create directory for saving results
dir.create(output_dir)

# create df to save the results for plotting
final_plotting_data <- data.frame()

# generate results and saves them 
for (i in AS_events_types) {

    # reading input files and generates results per each type of splicing event
    filtered_AS <- splicing_events(muscle_as_path, neuron_as_path, i)

    # saves each data frame as a separate csv file into output dir
    mapply(write.csv, x=filtered_AS, file=c(paste0(output_dir, i, "_muscle.csv"), paste0(output_dir, i, "_neuron.csv")), row.names=FALSE)

    # produces data for horizontal bar chart and saves it into final_plotting_data data frame
    plotting_AS_data <- bar_chart_data(filtered_AS, i)
    final_plotting_data <- rbind(final_plotting_data, plotting_AS_data)

    final_plotting_data <- final_plotting_data %>%
        group_by(Splice_type) %>%
        mutate(Percentage_common = round(sum(Count[Category == 'common']) / sum(Count[Category != 'common']) * 100, 1))

    final_plotting_data$Splice_type <- factor(final_plotting_data$Splice_type, levels = rev(AS_events_types))
    final_plotting_data$Category <- factor(final_plotting_data$Category, levels = c("C2C12-specific", "common", "NSC34-specific"))
}

# plots splicing event counts and save it into output folder
bar_plot <- horizontal_plot(final_plotting_data)
ggsave(filename = paste0(output_dir, "splicing_horizontal_bar_plot.png"), plot = bar_plot)



venn_data <- final_plotting_data %>%
    group_by(Category) %>%
    summarise(total_AS_events = sum(Count)) %>%
    mutate(Percentage_common = round(total_AS_events[Category == 'common'] / sum(total_AS_events[Category != 'common']) * 100, 1))





Venn(list(venn_data$total_AS_events))



t(final_plotting_data)
