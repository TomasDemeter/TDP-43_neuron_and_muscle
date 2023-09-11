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
samples = pd.read_csv(samplefile, dtype = str, index_col = 0)
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


def rMATS_inputs():
    samples = pd.read_csv(samplefile)
    bam_files = glob.glob(os.path.abspath("results/hisat2_aligned/*.bam"))

    # replace values in treatment column
    for i in samples.treatment:
        if "control" in i.lower():
            samples.treatment = samples.treatment.replace(i, "control")
        else:
            samples.treatment = samples.treatment.replace(i, "experimental")
            
    # create a dictionary with treatment and cell type as keys and SRR numbers as values
    d = {f"{treatment}_{cell_type}": [] for treatment, cell_type in samples[['treatment', 'cell type']].drop_duplicates().itertuples(index=False)}
    
    # map SRR numbers to bam file paths
    for treatment, cell_type, srr in samples[['treatment', 'cell type', 'SRR']].itertuples(index=False):
        d[f"{treatment}_{cell_type}"].append(srr)
    srr_to_bam = {bam.split('/')[-1].split('_')[0]: bam for bam in bam_files}

    # Replace the values in group_dict with the corresponding bam file paths
    for key, value in d.items():
        d[key] = [srr_to_bam[srr] for srr in value]
        
    srr_to_bam = {bam.split('/')[-1].split('_')[0]: bam for bam in bam_files}
    
    # create txt files
    os.makedirs("./results/rMATS_output", exist_ok=True)
    for key, values in d.items():
        with open(f"./results/rMATS_output/{key}.txt", "w") as f:
            for value in values:
                f.write(value + ',')


###################
# Desired outputs #
###################
BAM_FILES   = expand(RESULT_DIR + "hisat2_aligned/{SRR}_fastp_Aligned.sortedByCoord.out.bam", SRR = SAMPLES)
MULTIQC     = RESULT_DIR + "MultiQC/multiqc_report.html"
COUNTS      = RESULT_DIR + "featureCounts/feature_counts_table.tsv"
RMATS_N     = RESULT_DIR + "rMATS_output/rMATS_neuron/summary.txt"
RMATS_M     = RESULT_DIR + "rMATS_output/rMATS_muscle/summary.txt"
RMATS_BAMS  = RESULT_DIR + "rMATS_output/control_C2C12.txt"
DESEQ       = RESULT_DIR + "DESeq2_output/DESeq2_output.rds"
AS_SPLICING = RESULT_DIR + "AS_analysis_output/"
GO_ANALYSIS = RESULT_DIR + "GO_term_analysis/"
BIGWIG      = expand(RESULT_DIR + "bamCoverage/{SRR}.bigwig", SRR = SAMPLES)

############
# Pipeline #
############
rule all:
    input:
        BAM_FILES,
        COUNTS,
        MULTIQC,
        RMATS_N,
        RMATS_M,
        DESEQ,
        AS_SPLICING,
        GO_ANALYSIS,
        BIGWIG
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


#################################### 
# Generating input files for rMATS #
####################################    
rule rMATS_inputs:
    input:
        bam = expand(rules.hisat2_samtools.output.bam, SRR=SAMPLES)
    output:
        rmats_output_dir            = directory(config["rmats"]["output_dir"]),
        experimental_NSC34_paths    = RESULT_DIR + "rMATS_output/experimental_NSC34.txt",
        control_NSC34_paths         = RESULT_DIR + "rMATS_output/control_NSC34.txt",
        experimental_C2C12_paths    = RESULT_DIR + "rMATS_output/experimental_C2C12.txt",
        control_C2C12_paths         = RESULT_DIR + "rMATS_output/control_C2C12.txt"
    message: "Generating bam txt files with bam paths"
    run:
        rMATS_inputs()


############################################# 
# Alternative splicing analysis using rMATS #
#############################################
# Neurons       
rule rMATS_neuron:
    input:
        NSC34_experimental_paths    = rules.rMATS_inputs.output.experimental_NSC34_paths,
        NSC34_control_paths         = rules.rMATS_inputs.output.control_NSC34_paths
    output:
        rmats_summary       = RESULT_DIR + "rMATS_output/rMATS_neuron/summary.txt",
        rmats_output_dir    = directory(RESULT_DIR + "rMATS_output/rMATS_neuron/"),
        temp_output_dir     = directory(RESULT_DIR + "rMATS_output/rMATS_neuron/tmp/")
    params:
        read_length         = config["rmats"]["read_length"],
        read_type           = config["rmats"]["read_type"],
        FDR_cutoff          = config["rmats"]["FDR_cutoff"],
        gtf_file            = config["refs"]["gtf"],
        rmats_exe           = config["rmats"]['rmats_executable']
    threads: 14
    message: "Running rMATS on neuronal cells"
    shell:
        "{params.rmats_exe} "
        "--b1 {input.NSC34_experimental_paths} "
        "--b2 {input.NSC34_control_paths} "
        "--gtf {params.gtf_file} "
        "--cstat {params.FDR_cutoff} "
        "-t {params.read_type} "
        "--nthread {threads} "
        "--readLength {params.read_length} "
        "--od {output.rmats_output_dir} "
        "--tmp {output.temp_output_dir}"

# Muscle cells
rule rMATS_muscle:
    input:
        C2C12_experimental_paths    = rules.rMATS_inputs.output.experimental_C2C12_paths,
        C2C12_control_paths         = rules.rMATS_inputs.output.control_C2C12_paths
    output:
        rmats_summary       = RESULT_DIR + "rMATS_output/rMATS_muscle/summary.txt",
        rmats_output_dir    = directory(RESULT_DIR + "rMATS_output/rMATS_muscle/"),
        temp_output_dir     = directory(RESULT_DIR + "rMATS_output/rMATS_muscle/tmp/")
    params:
        read_length         = config["rmats"]["read_length"],
        read_type           = config["rmats"]["read_type"],
        FDR_cutoff          = config["rmats"]["FDR_cutoff"],
        gtf_file            = config["refs"]["gtf"],
        rmats_exe           = config["rmats"]['rmats_executable']
    threads: 14
    message: "Running rMATS on muscle cells"
    shell:
        "{params.rmats_exe} "
        "--b1 {input.C2C12_experimental_paths} "
        "--b2 {input.C2C12_control_paths} "
        "--gtf {params.gtf_file} "
        "--cstat {params.FDR_cutoff} "
        "-t {params.read_type} "
        "--nthread {threads} "
        "--readLength {params.read_length} "
        "--od {output.rmats_output_dir} "
        "--tmp {output.temp_output_dir}"

###########################################
# Produce table of normalised gene counts #
###########################################
rule DESeq2:
    input:
        raw         = RESULT_DIR + "featureCounts/feature_counts_table.tsv",
        metadata    = config["samples"],
        RBPs        = config["DESeq2"]["RNAbinding_proteins"]
    output:
        neuron_DE   = RESULT_DIR + "DESeq2_output/neuron_DE.csv",
        muscle_DE   = RESULT_DIR + "DESeq2_output/muscle_DE.csv",
        fpkm_values = RESULT_DIR + "DESeq2_output/fpkm_values.csv",
        DESeq2_dds  = RESULT_DIR + "DESeq2_output/DESeq2_output.rds",
        output_dir  = directory(RESULT_DIR + "DESeq2_output/")
    message:
        "Running DESeq2"
    shell:
        "mkdir -p {RESULT_DIR}DESeq2_output; "
        "Rscript --vanilla scripts/DESeq2.R {input.raw} {input.metadata} {input.RBPs} {output.output_dir}"

#################################
# Alternative splicing analysis # 
#################################
rule alternative_splicing:
    input:
        neuron_rMATS = rules.rMATS_neuron.output.rmats_output_dir,
        muscle_rMATS = rules.rMATS_muscle.output.rmats_output_dir
    output:
        output_dir  = directory(RESULT_DIR + "AS_analysis_output/")
    message:
        "Running alternative_splicing R script"
    shell:
        "Rscript --vanilla scripts/alternative_splicing.R {input.neuron_rMATS} {input.muscle_rMATS} {output.output_dir}"

###############
# GO analysis # 
###############
rule GO_analysis:
    input:
        differential_expression_results  = rules.DESeq2.output.output_dir,
        alternative_splicing_results = rules.alternative_splicing.output.output_dir
    output:
        output_dir  = directory(RESULT_DIR + "GO_term_analysis"),
        DE_go_neuron   = RESULT_DIR + "GO_term_analysis/DE_go_neuron.csv",
        DE_go_muscle   = RESULT_DIR + "GO_term_analysis/DE_go_muscle.csv",
        AS_go_neuron   = RESULT_DIR + "GO_term_analysis/AS_go_neuron.csv",
        AS_go_muscle   = RESULT_DIR + "GO_term_analysis/AS_go_muscle.csv"
    message:
        "Running GO_analysis R script"
    shell:
        "mkdir -p {output.output_dir}; "
        "Rscript --vanilla scripts/GO_analysis.R {input.differential_expression_results} {input.alternative_splicing_results} {output.output_dir}"

###############
# bamCoverage # 
###############
rule bamCoverage:
    input:
        bam = rules.hisat2_samtools.output.bam
    output:
        bigwig      = RESULT_DIR + "bamCoverage/{SRR}.bigwig"
    params:
        output_dir  = directory(RESULT_DIR + "bamCoverage"),
    threads: 7
    message:
        "Running bamCoverage"
    shell:
        "mkdir -p {params.output_dir}; "
        "samtools index {input.bam}; "
        "bamCoverage --bam {input.bam} --outFileName {output.bigwig} --numberOfProcessors {threads}"