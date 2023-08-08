#########################################
# Snakemake pipeline for RNA-Seq analysis
#########################################


###########
# Libraries
###########
import os
import glob
import pandas as pd


###############
# Configuration
###############
configfile: "config/config.yaml" # where to find parameters
WORKING_DIR = config["working_dir"]
RESULT_DIR = config["result_dir"]
samplefile = config["samples"]


########################
# Samples and conditions
########################
# create lists containing the sample names and conditions
samples = pd.read_csv(config["samples"], dtype=str, index_col=0)
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
        return (WORKING_DIR + "raw_reads/" + f"{wildcards.SRR}_1" + ".fastq.gz", WORKING_DIR + "raw_reads/" + f"{wildcards.SRR}_2" + ".fastq.gz")

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
        return "--in1 " + inFile[0] + " --in2 " + inFile[1] + " --out1 " + WORKING_DIR + "trimmed/" + wildcards.SRR + "_R1_trimmed.fq.gz --out2 "  + WORKING_DIR + "trimmed/" + wildcards.SRR + "_R2_trimmed.fq.gz"

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


#################
# Desired outputs
#################
MULTIQC = RESULT_DIR + "multiqc_report.html"
BAM_FILES = expand(RESULT_DIR + "hisat2_aligned/{SRR}_fastp_Aligned.sortedByCoord.out.bam", SRR = SAMPLES)
#MAPPING_REPORT = RESULT_DIR + "mapping_summary.csv"

rule all:
    input:
        MULTIQC,
        BAM_FILES
    message:
        "RNA-seq pipeline run complete!"


###########################
# Genome reference indexing
###########################
rule hisat2_index:
    input:
        fasta = config["refs"]["genome"],
        gtf =   config["refs"]["gtf"]
    output:
        genome_index = expand("{index}/GRCm39_index.{n}.ht2", index=config["refs"]["index"], n=range(1,9)),
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

#######################
# RNA-seq read trimming
#######################
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
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"],
        in_and_out_files =  get_trim_names
    resources: cpus=5
    shell:
        "mkdir -p {WORKING_DIR}fastp; "
        "mkdir -p {WORKING_DIR}trimmed; "
        "fastp --thread {threads}  "
        "--html {output.html} "
        "--json {output.json} "
        "--qualified_quality_phred {params.qualified_quality_phred} "
        "{params.in_and_out_files} "
        "2>{log.log_file}"

rule multiqc:
    input:
        fastp_input = expand(WORKING_DIR + "fastp/{SRR}_fastp.json", SRR = SAMPLES),
        hisat2_input = expand(RESULT_DIR + "hisat2_aligned/{SRR}_fastp_Log.final.out", SRR = SAMPLES)
    output:
        RESULT_DIR + "multiqc_report.html"
    params:
        fastp_directory = WORKING_DIR + "fastp/",
        hisat2_directory = RESULT_DIR + "hisat2_aligned/",
        outdir = RESULT_DIR
    message: "Summarising fastp and hisat2 reports with multiqc"
    shell:
        "multiqc --force "
        "--outdir {params.outdir} "
        "{params.fastp_directory} "
        "{params.hisat2_directory} "
        "--module fastp "
        "--module hisat2"
        

#########################
# RNA-Seq read alignement
#########################
rule hisat2_samtools:
    input:
        genome_index = expand("{index}/GRCm39_index.{n}.ht2", index=config["refs"]["index"], n=range(1,9))
    output:
        bam = RESULT_DIR + "hisat2_aligned/{SRR}_fastp_Aligned.sortedByCoord.out.bam",
        log = RESULT_DIR + "hisat2_aligned/{SRR}_fastp_Log.final.out"
    params:
        hisat2_input_file_names =  get_hisat2_names,
        splice_sites = config["refs"]["splice_sites"]
    message:
        "Mapping {wildcards.SRR} reads to genome"
    threads: 7
    shell:
        "hisat2 -x {config[refs][index]}/GRCm39_index "
        "{params.hisat2_input_file_names} "
        "2> {output.log} "
        "-p {threads} "
        "--known-splicesite-infile {params.splice_sites} "
        "--new-summary "
        "| samtools sort -o {output.bam}"

'''
rule generate_mapping_summary:
    input:
        expand(RESULT_DIR + "hisat2_aligned/{sample}_fastp_Log.final.out", sample = SAMPLES)
    output:
        RESULT_DIR + "mapping_summary.csv"
    message:
        "Concatenating Hisat2 mapping report and generating .csv mapping summary."
    params:
        directory_with_mapping_reports = RESULT_DIR + "hisat2_aligned/",
        star_directory_name = RESULT_DIR + "hisat2_aligned/"
    shell:
        "python scripts/generate_mapping_summary.py {params.directory_with_mapping_reports} {params.star_directory_name} {output}"
'''