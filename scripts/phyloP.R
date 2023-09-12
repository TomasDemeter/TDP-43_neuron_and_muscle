library(rtracklayer)
path_to_bigwig <- ".results/bamCoverage/SRR14183630.bigwig"
data <- BigWigFile(path_to_bigwig)

SE_neuron_path <- "results/AS_analysis_output/SE_neuron.csv"
SE_neuron <- read.csv(SE_neuron_path)

head(SE_neuron)
# Assuming 'SE_neuron' is your data frame
bed_df <- SE_neuron[, c("chr", "exonStart_0base", "exonEnd")]
# Convert 0-based start positions to 1-based start positions for BED format
bed_df$exonStart_0base <- bed_df$exonStart_0base + 1
# Write to BED file
write.table(bed_df, file = "exons.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
