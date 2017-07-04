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

# Software executable
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
        fastqc_raw_out = FASTQC_RAW_DIR,            # FASTQC on raw FASTQ
        fastqc_trimmed_out = FASTQC_TRIMMED_DIR,    # FASTQC on filtered FASTQ
        transrate_out = TRANSRATE_DIR,              # assembly evaluation
        interproscan_out = INTERPROSCAN_DIR         # final results
        

# FASTQC: This rule is use to generate an evaluation report raw FASTQ files 
# provided by the user.
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


# TRIMMOMATIC: This rule is use to filter raw FASTQ files. 
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


# FILTERED FASTQC: This rule is use to generate an evaluation report on 
# filtered FASTQ files previously processed by Trimmomatic.
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


# TRINITY: This rule is use to de novo assemble filtered FASTQ files into 
# contigs. 
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


# TRANSRATE: This rule is use to generate an evaluation report on previously
# de novo assembled contigs.
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


# TRANSDECODER: This rule is use to predict protein coding domains from 
# previoudly de novo assemble contigs.
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


# INTERPROSCAN: This rule is use to search for functional annotation of 
# previously predicted protein coding domains.
rule interproscan:
    input:
        prediction = TRANSDECODER_DIR
    output:
        interproscan_out = INTERPROSCAN_DIR
    log:
        OUT_DIR + "logs/interproscan/interproscan.log"
    params:
        transdecoder_out = TRANSDECODER_DIR + "/Trinity.fasta.transdecoder.pep",
        out_format = config["interproscan_params"]["out_format"],
        db = config["interproscan_params"]["db"]
    threads: 
        config["threads"]["interproscan"]
    shell:
        """
        mkdir {output.interproscan_out}
        
        sed -i 's/*//g' {params.transdecoder_out} 
        
        {interproscan} \
        -i {params.transdecoder_out} \
        -d {output.interproscan_out} \
        -f {params.out_format} \
        -appl {params.db} \
        -cpu {threads} \
        -dp \
        --goterms &> {log}
        """
