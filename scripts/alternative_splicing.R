library(Rmats)

neuron_as_path <- "./results/rMATS_output/rMATS_neuron/summary.txt"
muscle_as_path <- "./results/rMATS_output/rMATS_muscle/summary.txt"

neuron_as <- read.csv(neuron_as_path, sep = "\t")
muscle_as <- read.csv(muscle_as_path, sep = "\t")

neuron_sig_JCEC <- sum(neuron_as$"SignificantEventsJCEC")
muscle_sig_JCEC <- sum(neuron_as$"SignificantEventsJCEC")

'''
function input two file paths 
create dict where keys are different splicing events
each key will associate with a list of length muscle, neuron common
function should also output the actuall filtered dfs?
another function for plotting venn and horizontal bar chart

'''


## SE
muscle_SE_JCEC_path <- "./results/rMATS_output/rMATS_muscle/SE.MATS.JCEC.txt"
muscle_SE_JCEC <- read.csv(muscle_SE_JCEC_path, sep = "\t")
muscle_SE_JCEC <- subset(muscle_SE_JCEC, FDR <= 0.01)
muscle_SE_JCEC <- muscle_SE_JCEC$GeneID

neuron_SE_JCEC_path <- "./results/rMATS_output/rMATS_neuron/SE.MATS.JCEC.txt"
neuron_SE_JCEC <- read.csv(neuron_SE_JCEC_path, sep = "\t")
neuron_SE_JCEC <- subset(neuron_SE_JCEC, FDR <= 0.01)
neuron_SE_JCEC <- neuron_SE_JCEC$GeneID

#SE_NSC34 <-setdiff(neuron_SE_JCEC, muscle_SE_JCEC)
#SE_C2C12 <- setdiff(muscle_SE_JCEC, neuron_SE_JCEC) 
SE_common <- intersect(neuron_SE_JCEC, muscle_SE_JCEC)

length(muscle_SE_JCEC)
length(neuron_SE_JCEC)
length(SE_common)


length(SE_common) / (length(muscle_SE_JCEC) + length(neuron_SE_JCEC)) * 100
