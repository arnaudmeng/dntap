#!/usr/bin/python

################################################################################
# 
# DNTAP has been designed to perform both de novo transcriptome assembly and 
# transcriptome analysis. It includes:
#    - RNA-seq data filtering and trimming (Trimmomatic)
#    - De novo assembly (Trinity)
#    - Transcriptome assembly evaluation (Transrate) 
#    - Protein coding domains prediction (Transdescoder) 
#    - Protein functional annotation (InterProScan 5)
#
# This project has been implemented in Python - Snakemake
# Snakemake link : https://snakemake.readthedocs.io/en/stable/index.html
#
# Example command line : 
# > snakemake \
#   --snakefile dntap.py \
#   --configfile <config.yaml> \
#   --cores <max_n_cores>
# 
################################################################################

# Imports
import os 

# Get current working directory  
dir_path = os.getcwd()

# Configuration file
configfile: "/home/meng/PIPELINE/pipeline_v4/config.yaml"
    
# User defined samples
SAMPLES = config["samples"] if config["samples"] is not None else []

# User defined ouput directory
OUT_DIR = config["directories"]["outdir"]

# Relative output directories
FASTQC_RAW_DIR = OUT_DIR + "fastqc_raw"
FASTQC_TRIMMED_DIR = OUT_DIR + "fastqc_trimmed"
TRIMMOMATIC_DIR = OUT_DIR + "trimmomatic_out/"
TRINITY_DIR = OUT_DIR + "trinity_out"
TRANSRATE_DIR = OUT_DIR + "transrate_out"
TRANSDECODER_DIR = OUT_DIR + "transdecoder_out"
INTERPROSCAN_DIR = OUT_DIR + "interproscan_out"

# Software launch command
fastqc = config["software"]["fastqc"]
trimmomatic = config["software"]["trimmomatic"]
trinity = config["software"]["trinity"]
transrate = config["software"]["transrate"]
transdecoder_longorfs = config["software"]["transdecoder_longorfs"]
transdecoder_predict = config["software"]["transdecoder_predict"]
interproscan = config["software"]["interproscan"]

# ALL
rule all:
    input:
        fastqc_raw_out = FASTQC_RAW_DIR,
        fastqc_trimmed_out = FASTQC_TRIMMED_DIR,
        transrate_out = TRANSRATE_DIR,
        interproscan_out = INTERPROSCAN_DIR + "/interproscan.out"    
        ,
# FASTQC : evaluate quality of raw FASTQ files
rule raw_fastqc:
    input:
        r1 = SAMPLES["forward"],
        r2 = SAMPLES["reverse"]
    output:
        fastqc_raw_out = FASTQC_RAW_DIR
    log:
        OUT_DIR + "logs/fastqc/raw_fastqc.log"
    threads: 
        config["threads"]["fastqc"]
    shell:
        """
        mkdir {output.fastqc_raw_out}
        
        {fastqc} \
        {input.r1} \
        {input.r2} \
        --outdir {output.fastqc_raw_out} \
        --threads {threads} &> {log}
        """

# TRIMMOMATIC : designed to trim and clean FASTQ files
rule trimmomatic:
    input:
        r1 = SAMPLES["forward"],
        r2 = SAMPLES["reverse"]
    output:
        r1_paired = TRIMMOMATIC_DIR + "forward.trimmomatic.paired.fastq",
        r1_unpaired = TRIMMOMATIC_DIR + "forward.trimmomatic.unpaired.fastq",
        r2_paired = TRIMMOMATIC_DIR + "reverse.trimmomatic.paired.fastq",
        r2_unpaired = TRIMMOMATIC_DIR + "reverse.trimmomatic.unpaired.fastq"
    log:
        OUT_DIR + "logs/trimmomatic/trimmomatic.log"
    threads: 
        config["threads"]["fastqc"]
    params:
        trimmomatic_params = config["trimmomatic_params"]
    shell:
        """
        {trimmomatic} PE \
            -threads {threads} \
            {input.r1} \
            {input.r2} \
            {output.r1_paired} \
            {output.r1_unpaired} \
            {output.r2_paired} \
            {output.r2_unpaired} \
            {params.trimmomatic_params} 2> {log}
        """

# FASTQC : evaluate quality of trimmed FASTQ files
rule trim_fastqc:
    input:
        r1 = TRIMMOMATIC_DIR + "forward.trimmomatic.paired.fastq",
        r2 = TRIMMOMATIC_DIR + "reverse.trimmomatic.paired.fastq"
    output:
        fastqc_trimmed_out = FASTQC_TRIMMED_DIR
    log:
        OUT_DIR + "logs/fastqc/trimmed_fastqc.log"
    threads: 
        config["threads"]["fastqc"]
    shell:
        """
        mkdir {output.fastqc_trimmed_out}
                
        {fastqc} \
        {input.r1} \
        {input.r2} \
        --outdir {output.fastqc_trimmed_out} \
        --threads {threads} &> {log}
        """

# TRINITY : designed to assemble de novo RNA-seq data
rule trinity:
    input:
        left = TRIMMOMATIC_DIR + "forward.trimmomatic.paired.fastq",
        right = TRIMMOMATIC_DIR + "reverse.trimmomatic.paired.fastq"
    output:
        trinity_out = TRINITY_DIR + "/Trinity.fasta"
    log:
        OUT_DIR + "logs/trinity/trinity.log"
    params:
        max_memory = config["trinity_params"]["max_memory"],
        trinity_dir = TRINITY_DIR
    threads: 
        config["threads"]["trinity"]
    shell:
        """
        {trinity} \
        --seqType fq \
        --left {input.left} \
        --right {input.right} \
        --output {params.trinity_dir} \
        --CPU {threads} \
        --max_memory {params.max_memory} > {log}
        """

# TRANSRATE : evaluate de novo assembled contigs
rule transrate:
    input:
        assembly = TRINITY_DIR + "/Trinity.fasta",
        left = TRIMMOMATIC_DIR + "forward.trimmomatic.paired.fastq",
        right = TRIMMOMATIC_DIR + "reverse.trimmomatic.paired.fastq"
    output:
        transrate_out = TRANSRATE_DIR
    log:
        OUT_DIR + "logs/transrate/transrate.log"
    threads: 
        config["threads"]["transrate"]
    shell:
        """
        {transrate} \
        --assembly {input.assembly} \
        --left {input.left} \
        --right {input.right} \
        --output {output.transrate_out} \
        --threads {threads} > {log}
        """

# TRANSDECODER : identifies candidate coding regions within transcript sequences
rule transdecoder:
    input:
        assembly = TRINITY_DIR + "/Trinity.fasta",
    output:
        transdecoder_out = TRANSDECODER_DIR
    log:
        OUT_DIR + "logs/transdecoder/transdecoder.log"
    params:
        min_protein_len = config["transdecoder_params"]["min_protein_len"],
    threads: 
        config["threads"]["transdecoder"]
    shell:
        """
        mkdir {output.transdecoder_out}
        
        cd {output.transdecoder_out}
        
        {transdecoder_longorfs} \
        -t {input.assembly} \
        -m {params.min_protein_len} &> {log}
        
        {transdecoder_predict} \
        -t {input.assembly} \
        --cpu {threads} &>> {log}
        
        cd {dir_path}
        """
