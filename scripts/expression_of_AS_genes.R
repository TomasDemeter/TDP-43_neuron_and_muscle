library(tidyverse)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(readr)



fpkm_values_path <- "results/DESeq2_output/fpkm_values.csv"
AS_output_path <- "results/AS_analysis_output/"

###############################
# function for preparing data #
###############################
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

read_column <- function(file, column_name) {
    data <- read.csv(file, stringsAsFactors = FALSE)
    return(data[[column_name]])
}

get_fpkm <- function(path, pattern, fpmk_values) {
    files <- list.files(path = path, pattern = pattern)
    list <- lapply(paste0(path, files), read_column, column_name = "GeneID")
    combined <- data.frame(Col = unlist(list))
    fpkm <- fpkm_values[fpkm_values$X %in% combined$Col, ]

    return(fpkm)
}

get_box_df <- function(common, specific) {
    common <- common %>%
            mutate(average = rowMeans(select(., -X), na.rm = TRUE)) %>%
            mutate(Type = "Common")

    specific <- specific %>%
            mutate(average = rowMeans(select(., -X), na.rm = TRUE)) %>%
            mutate(Type = "Specific")
    df <- data.frame()
    df <- rbind(select(common, average, Type), select(specific, average, Type))
    return(df)
}



AS_expression_boxplot <- function(df) {
    wilcox_result <- wilcox.test(average ~ Type, data = df)
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

    boxplot <- ggplot(df, aes(x = Type, y = log10(average), fill = Type)) +
    geom_boxplot(notch = TRUE) +
    scale_fill_manual(values = c("Common" = "#8aa261", "Specific" = "#cccccc")) +
    labs(x = "", y = "Log10(FPKM + 1)") +
    ylim(0, 3) +
    nature_theme() +
    theme(legend.position = "none") + 
    annotate("text", x = 1.5, y = 3, label = sig_label, size = 6) +
    geom_segment(aes(x = 1.2, y = 2.9, xend = 1.8, yend = 2.9))

    return(boxplot)
}


########################
########################
# Running the analysis #
########################
########################

# read fpkm values and add 1 too all 
fpkm_values <- read.csv(fpkm_values_path)
fpkm_values[,2:ncol(fpkm_values)] <- sapply(fpkm_values[,2:ncol(fpkm_values)], as.numeric)
fpkm_values[,2:ncol(fpkm_values)] <- fpkm_values[,2:ncol(fpkm_values)] + 1


neuron_fpkm <- get_fpkm(AS_output_path, "*_neuron.csv$", fpkm_values)
muscle_fpkm <-  get_fpkm(AS_output_path, "*_muscle.csv$", fpkm_values)


common <- intersect(neuron_fpkm, muscle_fpkm)
specific <- setdiff(union(neuron_fpkm, muscle_fpkm), common)


avg_expressions <- get_box_df(common, specific)
avg_expressions_boxplot <- AS_expression_boxplot(avg_expressions)


ggsave(filename = paste0(AS_output_path, "/avg_expressions_boxplot.png"), plot = avg_expressions_boxplot, width = 20, height = 15)

