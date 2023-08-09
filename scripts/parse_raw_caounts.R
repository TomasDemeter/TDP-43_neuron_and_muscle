args <- commandArgs(trailingOnly = TRUE)

raw_counts_filepath <- args[1]
hisat2_results_directory_name <- args[2]
out_parsed_counts_filepath <- args[3]

# Parse raw counts to have clean column names
raw_counts <- read.csv(raw_counts_filepath, skip=1, sep="\t")
raw_counts_parsed <- raw_counts[, !(names(raw_counts) %in% c("Chr","Start","End","Strand","Length"))]

old_col_names <- colnames(raw_counts_parsed)

new_col_names <- gsub(hisat2_results_directory_name, '', old_col_names)
new_col_names <- gsub('_fastp_Aligned.sortedByCoord.out.bam', '', new_col_names)

colnames(raw_counts_parsed) <- new_col_names

# Write parsed raw counts to .csv file
write.csv(raw_counts_parsed, out_parsed_counts_filepath, row.names=FALSE)
