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


########################
# Samples and conditions
########################
# read the table containing the sample, condition and fastq file information
units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)

# create lists containing the sample names and conditions
SAMPLES = units.index.get_level_values('sample').unique().tolist()
samples = pd.read_csv(config["units"], dtype=str,index_col=0,sep="\t")
samplefile = config["units"]


###########################
# Input functions for rules
###########################
def sample_is_single_end(sample):
    """This function detect missing value in the column 2 of the units.tsv"""
    if "fq2" not in samples.columns:
        return True
    else:
        return pd.isnull(samples.loc[(sample), "fq2"])

def get_fastq(wildcards):
    """This function checks if the sample has paired end or single end reads and returns 1 or 2 names of the fastq files"""
    if sample_is_single_end(wildcards.sample):
        return samples.loc[(wildcards.sample), ["fq1"]].dropna()
    else:
        return samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def get_trim_names(wildcards):
    """
    This function:
      1. Checks if the sample is paired end or single end
      2. Returns the correct input and output trimmed file names. 
    """
    if sample_is_single_end(wildcards.sample):
        inFile = samples.loc[(wildcards.sample), ["fq1"]].dropna()
        return "--in1 " + inFile[0] + " --out1 " + WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz" 
    else:
        inFile = samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
        return "--in1 " + inFile[0] + " --in2 " + inFile[1] + " --out1 " + WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz --out2 "  + WORKING_DIR + "trimmed/" + wildcards.sample + "_R2_trimmed.fq.gz"

def get_hisat2_names(wildcards):
    """
    This function:
      1. Checks if the sample is paired end or single end.
      2. Returns the correct input file names for Hisat2 mapping step.
    """
    if sample_is_single_end(wildcards.sample):
        return WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz"     
    else:
        return WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz " + WORKING_DIR + "trimmed/" + wildcards.sample + "_R2_trimmed.fq.gz"


#################
# Desired outputs
#################
MULTIQC = RESULT_DIR + "multiqc_report.html"
BAM_FILES = expand(RESULT_DIR + "hisat2/{sample}_Aligned.sortedByCoord.out.bam", sample = SAMPLES)
MAPPING_REPORT = RESULT_DIR + "mapping_summary.csv"

if config["keep_working_dir"] == False:
    rule all:
        input:
            MULTIQC,
            BAM_FILES, 
            MAPPING_REPORT,
            RESULT_DIR + "raw_counts.parsed.tsv",
            RESULT_DIR + "scaled_counts.tsv"
        message:
            "RNA-seq pipeline run complete!"
        shell:
            "cp config/config.yaml {RESULT_DIR};"
            "cp config/samples.tsv {RESULT_DIR};"
            "rm -r {WORKING_DIR}"
elif config["keep_working_dir"] == True:
    rule all:
        input:
            MULTIQC,
            BAM_FILES, 
            MAPPING_REPORT,
            RESULT_DIR + "raw_counts.parsed.tsv",
            RESULT_DIR + "scaled_counts.tsv"
        message:
            "RNA-seq pipeline run complete!"
        shell:
            "cp config/config.yaml {RESULT_DIR};"
            "cp config/samples.tsv {RESULT_DIR};"
else:
    raise ValueError('Please specify only "True" or "False" for the "keep_working_dir" parameter in the config file.')

'''
###########################
# Genome reference indexing
###########################
rule hisat2_index:
    input:
        fasta = config["refs"]["genome"],
        gtf =   config["refs"]["gtf"]
    output:
        genome_index = [WORKING_DIR + "GRCm39_index/"]
    message:
        "generating Hisat2 genome index"
    params:
        genome_dir = WORKING_DIR + "genome/"
    threads:
        16
    resources:
        mem_mb=50000
    shell:
        "mkdir -p {params.genome_dir}; " # create dir if not existing already
        "extract_splice_sites.py {input.gtf} > splice_sites.txt; " # 
        "hisat2-build -p {threads} "
        "-ss splice_sites.txt "
        "-f {input.fasta}"
'''        


#######################
# RNA-seq read trimming
#######################
rule fastp:
    input:
        get_fastq
    output:
        fq1  = WORKING_DIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz",
        fq2  = WORKING_DIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz",
        html = WORKING_DIR + "fastp/{sample}_fastp.html",
        json = WORKING_DIR + "fastp/{sample}_fastp.json"
    message:"trimming {wildcards.sample} reads"
    threads: 5
    log:
        RESULT_DIR + "fastp/{sample}.log.txt"
    params:
        sampleName = "{sample}",
        in_and_out_files =  get_trim_names,
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"]
    resources: cpus=5
    shell:
        "touch {output.fq2}; "
        "fastp --thread {threads}  "
        "--html {output.html} "
        "--json {output.json} "
        "--qualified_quality_phred {params.qualified_quality_phred} "
        "{params.in_and_out_files} "
        "2>{log}"

rule multiqc:
    input:
        expand(WORKING_DIR + "fastp/{sample}_fastp.json", sample = SAMPLES)
    output:
        RESULT_DIR + "multiqc_report.html"
    params:
        fastp_directory = WORKING_DIR + "fastp/",
        outdir = RESULT_DIR
    message: "Summarising fastp reports with multiqc"
    shell:
        "multiqc --force "
        "--outdir {params.outdir} "
        "{params.fastp_directory} "


#########################
# RNA-Seq read alignement
#########################
rule hisat2_samtools:
    input:
        read1= WORKING_DIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz",
        read2= WORKING_DIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz"
    output:
        bam = RESULT_DIR + "hisat2/{sample}_Aligned.sortedByCoord.out.bam",
        log = RESULT_DIR + "hisat2/{sample}_Log.final.out"
    message:
        "Mapping {wildcards.sample} reads to genome"
    params:
        sample_name             =  "{sample}",
        prefix                  =  RESULT_DIR + "hisat2/{sample}_"
        hisat2_input_file_names =  get_hisat2_names,
        prefix                  =  RESULT_DIR + "hisat2/{sample}_",
        genome_index            =  WORKING_DIR + "GRCm39_index/"
    threads: 10
    resources: cpus=10
    shell:
        "hisat2 -x {params.genome_index} "
        "-p {threads} "
        "-1 {input.read1} "
        "-2 {input.read2} "
        "| samtools sort -o {output.bam}"










# add esearch 

# Read the accession numbers from the file
with open('data/accession_numbers.txt') as f:
    accessions = [line.strip() for line in f]

SRA, FRR = glob_wildcards("data/raw/raw_reads/{sra}_{frr}.fastq.gz")

rule all:
    input:
        expand("data/raw/raw_reads/{accession}.fastq", accession=accessions),
        expand("data/raw/raw_QC/{sra}_{frr}_fastqc.{extension}", sra=SRA, frr=FRR, extension=["zip", "html"]),
        expand("data/hisat2_aligned/{sra}_aligned.bam", sra=SRA)

rule prefetch_sra:
    output:
        temp("{accession}.sra")
    shell:
        "prefetch {wildcards.accession}"

rule fasterq_dump:
    input:
        rules.prefetch_sra.output
    params:
        outdir = 'data/raw/raw_reads'
    output:
        "{params.outdir}/{accession}.fastq"
    shell:
        "fasterq-dump {input} --outdir {params.outdir} --split-files --skip-technical --progress"

rule rawFastqc:
    input:
        raw_reads="data/raw/raw_reads/{sra}_{frr}.fastq.gz"
    output:
        zip="data/raw/raw_QC/{sra}_{frr}_fastqc.zip",
        html="data/raw/raw_QC/{sra}_{frr}_fastqc.html"
    threads:
        1
    params:
        path="data/raw/raw_QC"
    shell:
        "fastqc {input.raw_reads} --threads {threads} -o {params.path}"

rule trimmomatic:
    input:
        read1="data/trimmed_reads/{sra}_1.fastq",
        read2="data/trimmed_reads/{sra}_2.fastq"
    output:
        forwardPaired="data/trimmed_reads/{sra}_1.fastq",
        reversePaired="data/trimmed_reads/{sra}_2.fastq"
    threads:
        1
    params:
        basename="data/trimmed_reads/{sra}",
        log="data/trimmed_reads/{sra}.log"
    shell:
        "trimmomatic PE -threads {threads} {input.read1} {input.read2} {output.forwardPaired} {output.reversePaired} ILLUMINACLIP:data/NEBNext_Ultra_II_RNA_Library_Prep_Kit.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36 2>{params.log}"
