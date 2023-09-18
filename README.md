# Summary 

This repository contains the code and results of my reanalysis of the RNA-seq data from the paper “Cell environment shapes TDP-43 function with implications in neuronal and muscle disease” by Šušnjar et al. (2022). The paper investigated the role of TDP-43, a RNA-binding protein involved in various aspects of mRNA metabolism, in mouse muscle (C2C12) and neuronal (NSC34) cells. The authors identified different sets of transcripts and splicing events that are regulated by TDP-43 in each cell type, and showed that some of them are also altered in human tissues of patients with neurodegenerative and myodegenerative diseases.

I reprocessed the raw sequencing data from scratch using Snakemake, a workflow management system that allows for reproducible and scalable data analysis. Snakemake enables the definition of analysis steps as rules that can be executed in parallel, and automatically handles the dependencies between them. Snakemake also tracks the changes in the input and output files, and only reruns the necessary rules when something is modified.

The main steps of my analysis pipeline are:

- Quality control and trimming of the reads using fastp and multiqc.
- Alignment of the reads to the mouse reference genome (GRCm38) using hisat2, a fast and memory-efficient aligner that can handle spliced alignments.
- Quantification of gene expression and detection of alternative splicing events using rMATS and MAJIQ.
- Differential expression and splicing analysis using DESeq2.
- Functional enrichment analysis using GOseq.

I used hisat2 instead of star, another popular aligner for RNA-seq data, because hisat2 requires less memory and disk space, and has comparable or better accuracy and speed. This makes hisat2 more suitable for large-scale analyses or for running on machines with limited resources.

This is still work in progress.

!rulegraph