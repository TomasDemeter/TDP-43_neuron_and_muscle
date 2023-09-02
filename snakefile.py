###########################################
# Snakemake pipeline for RNA-Seq analysis #
###########################################


#############
# Libraries #
#############
import os
import glob
import pandas as pd


#################
# Configuration #
#################
configfile: "config/config.yaml" # where to find parameters
WORKING_DIR = config["working_dir"]
RESULT_DIR  = config["result_dir"]
samplefile  = config["samples"]


##########################
# Samples and conditions #
##########################
# create lists containing the sample names and conditions
samples = pd.read_csv(config["samples"], dtype = str, index_col = 0)
SAMPLES = samples.index.tolist()


###########################
# Input functions for rules
###########################
def sample_is_single_end(sample_name):
    """This function checks samples metadata for Library layout column 
    and evaulate whether samples are single or paired end"""
    library_layout_col = None
    for col in samples.columns:
        if col.lower().replace(" ", "") == "librarylayout":
            library_layout_col = col
            break
    if library_layout_col is None:
        raise ValueError("Column 'librarylayout' not found in DataFrame")
    sample = samples.loc[sample_name]
    return sample[library_layout_col].lower() != "paired"

def get_fastq(wildcards):
    """This function checks if the sample has paired end or single end reads 
    and returns either the value of the sample attribute of the wildcards object if single end 
    or with _1 or _2 attached to the string if paired"""
    if sample_is_single_end(wildcards.SRR):
        return WORKING_DIR + "raw_reads/" + f"{wildcards.SRR}" + ".fastq.gz"
    else:
        return (
            WORKING_DIR + "raw_reads/" + f"{wildcards.SRR}_1" + ".fastq.gz",
            WORKING_DIR + "raw_reads/" + f"{wildcards.SRR}_2" + ".fastq.gz",
        )

def get_trim_names(wildcards):
    """
    This function:
        1. Checks if the sample is paired end or single end
        2. Returns the correct input and output trimmed file names. 
    """
    inFile = get_fastq(wildcards)

    if sample_is_single_end(wildcards.SRR):
        return "--in1 " + inFile[0] + " --out1 " + WORKING_DIR + "trimmed/" + wildcards.SRR + "_trimmed.fq.gz" 
    else:
        return (
            "--in1 " + inFile[0] + " --in2 " + inFile[1] + " --out1 " +
            WORKING_DIR + "trimmed/" + wildcards.SRR + "_R1_trimmed.fq.gz --out2 " +
            WORKING_DIR + "trimmed/" + wildcards.SRR + "_R2_trimmed.fq.gz"
        )

def get_hisat2_names(wildcards):
    """
    This function:
        1. Checks if the sample is paired end or single end.
        2. Returns the correct input file names for Hisat2 mapping step.
    """
    if sample_is_single_end(wildcards.SRR):
        return "-1 " + WORKING_DIR + "trimmed/" + wildcards.SRR + "_R1_trimmed.fq.gz"     
    else:
        return "-1 " + WORKING_DIR + "trimmed/" + wildcards.SRR + "_R1_trimmed.fq.gz" + " -2 " + WORKING_DIR + "trimmed/" + wildcards.SRR + "_R2_trimmed.fq.gz"


def rmart_inputs(samplefile):
    snakemake_file = pd.read_csv(samplefile)
    bam_files = glob.glob(os.path.abspath(RESULT_DIR + "hisat2_aligned/*.bam"))
    
    controls = [snakemake_file.iloc[i, 0] for i, j in enumerate(snakemake_file.treatment) if "control" in j.lower()]
    experiments = [snakemake_file.iloc[i, 0] for i, j in enumerate(snakemake_file.treatment) if "control" not in j.lower()]
    
    control_bams = [i for i in bam_files if any(c in i for c in controls)]
    experiments_bams = [i for i in bam_files if any(c in i for c in experiments)]
    
    controls_txt_path = os.path.abspath(RESULT_DIR + 'hisat2_aligned/control_bams.txt')
    experiments_txt_path = os.path.abspath(RESULT_DIR + 'hisat2_aligned/experiments_bams.txt')

    with open(controls_txt_path, 'w') as f:
        f.write(','.join(control_bams))
    
    with open(experiments_txt_path, 'w') as f:
        f.write(','.join(experiments_bams))
        
    return "--b1 " + controls_txt_path + " --b2 " + experiments_txt_path

###################
# Desired outputs #
###################
BAM_FILES   = expand(RESULT_DIR + "hisat2_aligned/{SRR}_fastp_Aligned.sortedByCoord.out.bam", SRR = SAMPLES)
MULTIQC     = RESULT_DIR + "MultiQC/multiqc_report.html"
COUNTS      = RESULT_DIR + "featureCounts/feature_counts_table.tsv"
RMATS       = RESULT_DIR + "rMATS_output/"
DESEQ       = RESULT_DIR + "DESeq2_output/DESeq2_output.rds"
GO_ANALYSIS = RESULT_DIR + "GO_term_analysis/"

############
# Pipeline #
############
rule all:
    input:
        BAM_FILES,
        COUNTS,
        MULTIQC,
        RMATS,
        DESEQ,
        GO_ANALYSIS
    message:
        "RNA-seq pipeline run complete!"


#############################
# Genome reference indexing #
#############################
rule hisat2_index:
    input:
        fasta   = config["refs"]["genome"],
        gtf     = config["refs"]["gtf"]
    output:
        genome_index = expand("{index}/GRCm39_index.{n}.ht2", index = config["refs"]["index"], n = range(1, 9)),
        splice_sites = config["refs"]["splice_sites"]
    message:
        "generating Hisat2 genome index"
    threads:
        14
    params:
        dirs = config["refs"]["directories"]
    shell:
        "mkdir -p {params.dirs}; "
        "extract_splice_sites.py {input.gtf} > {output.splice_sites}; "
        "hisat2-build -f {input.fasta} "
        "{config[refs][index]}/GRCm39_index "
        "-p {threads}"


#########################
# RNA-seq read trimming #
#########################
rule fastp:
    input:
        get_fastq
    output:
        html = WORKING_DIR + "fastp/{SRR}_fastp.html",
        json = WORKING_DIR + "fastp/{SRR}_fastp.json"
    message:
        "trimming {wildcards.SRR} reads"
    threads: 7
    log:
        log_file = WORKING_DIR + "fastp/{SRR}.log.txt"
    params:
        in_and_out_files    = get_trim_names,
        phread_quality      = config["fastp"]["phread_quality"],
        base_limit          = config["fastp"]["base_limit"],
        percent_limit       = config["fastp"]["percent_limit"]
    shell:
        "mkdir -p {WORKING_DIR}fastp; "
        "mkdir -p {WORKING_DIR}trimmed; "
        "fastp --thread {threads} "
        "--qualified_quality_phred {params.phread_quality} "
        "--unqualified_percent_limit {params.percent_limit} "
        "--n_base_limit {params.base_limit} "
        "--html {output.html} "
        "--json {output.json} "
        "{params.in_and_out_files} "
        "2>{log.log_file}"
        

###########################
# RNA-Seq read alignement #
###########################
rule hisat2_samtools:
    input:
        genome_index    = rules.hisat2_index.output.genome_index,
        splice_sites    = rules.hisat2_index.output.splice_sites,
        fastp_json      = WORKING_DIR + "fastp/{SRR}_fastp.json"
    output:
        bam     = RESULT_DIR + "hisat2_aligned/{SRR}_fastp_Aligned.sortedByCoord.out.bam",
        log     = RESULT_DIR + "hisat2_aligned/{SRR}_fastp_Log.final.out"
    params:
        hisat2_input_file_names =  get_hisat2_names,
    message:
        "Mapping {wildcards.SRR} reads to genome"
    threads: 7
    shell:
        "hisat2 -x {config[refs][index]}/GRCm39_index "
        "{params.hisat2_input_file_names} "
        "2> {output.log} "
        "-p {threads} "
        "--known-splicesite-infile {input.splice_sites} "
        "--new-summary "
        "| samtools sort -o {output.bam}"

####################################
# Produce table of raw gene counts #
####################################
rule featureCounts:
    input:
        bams    = expand(RESULT_DIR + "hisat2_aligned/{SRR}_fastp_Aligned.sortedByCoord.out.bam", SRR = SAMPLES),
        gtf     = config["refs"]["gtf"]
    output:
        raw_feature_couts = RESULT_DIR + "featureCounts/feature_counts_table.tsv",
    message: "Producing the table of raw counts (counting read multimappers)"
    threads: 14
    shell:
        "featureCounts -T {threads} "
        "-s 0 "
        "-t exon "
        "-g gene_id "
        "-F 'GTF' "
        "-a {input.gtf} "
        "-o {output} "
        "-p {input.bams}"


##################
# MultiQC report #
##################
rule multiqc:
    input:
        fastp_input     = expand(WORKING_DIR + "fastp/{SRR}_fastp.json", SRR = SAMPLES),
        hisat2_input    = expand(RESULT_DIR + "hisat2_aligned/{SRR}_fastp_Log.final.out", SRR = SAMPLES),
        feature_counts  = rules.featureCounts.output.raw_feature_couts
    output:
        RESULT_DIR + "MultiQC/multiqc_report.html"
    params:
        fastp_directory         = WORKING_DIR + "fastp/",
        hisat2_directory        = RESULT_DIR + "hisat2_aligned/",
        feature_counts_input    = RESULT_DIR + "featureCounts/feature_counts_table.tsv.summary",
        outdir                  = RESULT_DIR + "MultiQC/"
    message: "Summarising fastp, hisat2 and featureCounts reports with multiqc"
    shell:
        "multiqc --force "
        "--outdir {params.outdir} "
        "{params.fastp_directory} "
        "{params.hisat2_directory} "
        "{params.feature_counts_input} "
        "--module fastp "
        "--module hisat2 "
        "--module featureCounts"   


#############################################
# Alternative splicing analysis using rMATS #
#############################################       
rule rMATS:
    output:
        rmats_output = directory(RESULT_DIR + "rMATS_output/") #try summary.txt later
    params:
        read_length = config["rmats"]["read_length"],
        FDR_cutoff  = config["rmats"]["FDR_cutoff"],
        read_type   = config["rmats"]["read_type"],
        gtf_file    = config["refs"]["gtf"],
        inputs      = rmart_inputs(samplefile),
        rmats_exe   = config["rmats"]['rmats_executable'],
        temp_output = RESULT_DIR + "rMATS_output/tmp/"
    threads: 14
    message: "Running rMATS"
    shell:
        "mkdir -p {output.rmats_output}; "
        "{params.rmats_exe} "
        "{params.inputs} "
        "--gtf {params.gtf_file} "
        "-t {params.read_type} "
        "--nthread {threads} "
        "--readLength {params.read_length} "
        "--cstat {params.FDR_cutoff} "
        "--od {output.rmats_output} "
        "--tmp {params.temp_output}"


###########################################
# Produce table of normalised gene counts #
###########################################
rule DESeq2:
    input:
        raw       = RESULT_DIR + "feature_counts_table.tsv",
        metadata  = config["samples"]
    output:
        neuron_dge  = RESULT_DIR + "DESeq2_output/neuron_dge.csv",
        muscle_dge  = RESULT_DIR + "DESeq2_output/muscle_dge.csv",
        fpkm_values = RESULT_DIR + "DESeq2_output/fpkm_values.csv",
        DESeq2_dds  = RESULT_DIR + "DESeq2_output/DESeq2_output.rds",
        output_dir  = directory(RESULT_DIR + "DESeq2_output/")
    message:
        "Running DESeq2"
    shell:
        "mkdir -p {RESULT_DIR}DESeq2_output; "
        "Rscript --vanilla scripts/DESeq2.R {input.raw} {input.metadata} {output.output_dir}"


###############
# GO analysis # 
###############
rule GOanalysis:
    input:
        neuronal_dge = rules.DESeq2.output.neuron_dge,
        muscular_dge = rules.DESeq2.output.muscle_dge
    output:
        ego_neuron  = RESULT_DIR + "GO_term_analysis/ego_neuron.csv",
        ego_muscle  = RESULT_DIR + "GO_term_analysis/ego_muscle.csv",
        output_dir  = RESULT_DIR + "GO_term_analysis"
    message:
        "Running GO_analysis"
    shell:
        "mkdir -p {output.output_dir}; "
        "Rscript --vanilla scripts/GO_analysis {input.neuronal_dge} {input.muscular_dge} {output.output_dir}"