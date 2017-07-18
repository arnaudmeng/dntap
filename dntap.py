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
# To generate pipeline diagram: 
# > snakemake \
#   --snakefile dntap.py \
#   --configfile dntap_config.yaml \
#   --dag | dot -Tpng > diag1.png
#
################################################################################

# Imports
import os 

# Get current working directory  
dir_path = os.getcwd()
    
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

# function to create fake inputs
def make_fake_inputs(index):
        
    if (index == "pe"): 
        
        sample_dir = os.path.dirname(SAMPLES["forward"])
        fake_file = sample_dir + "/none"
        os.system("touch " + fake_file)
        SAMPLES["single"] = fake_file
        
    elif (index == "se"): 
        
        sample_dir = os.path.dirname(SAMPLES["single"])
        fake_file = sample_dir + "/none"
        os.system("touch " + fake_file)
        SAMPLES["forward"] = fake_file
        SAMPLES["reverse"] = fake_file

# function to define inputs to RAW_FASTQC rule
def define_raw_fastqc_inputs(wildcards):
    data_type = config["data_type"]["type"]
    
    if (data_type == "pe"):
        input = [SAMPLES["forward"],SAMPLES["reverse"]]
        
    elif (data_type == "se"):
        input = SAMPLES["single"]
        
    return input

# creating fake files
make_fake_inputs(config["data_type"]["type"])

# ALL
rule all:
    input:
        fastqc_raw_out = FASTQC_RAW_DIR,            # FASTQC on raw FASTQ
        interproscan_out = INTERPROSCAN_DIR,        # final results
        transrate_out = TRANSRATE_DIR,              # assembly evaluation
        fastqc_trimmed_out = FASTQC_TRIMMED_DIR,    # FASTQC on filtered FASTQ


# FASTQC: This rule is use to generate an evaluation report raw FASTQ files 
# provided by the user.
rule raw_fastqc:
    input:
        fastq = define_raw_fastqc_inputs
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
        {input.fastq} \
        --outdir {output.fastqc_raw_out} \
        --threads {threads} &> {log}
        """


# TRIMMOMATIC: This rule is use to filter raw FASTQ files. 
rule trimmomatic:
    input:
        forward = SAMPLES["forward"],
        reverse = SAMPLES["reverse"],
        single = SAMPLES["single"]
    output:
        out = TRIMMOMATIC_DIR
    log:
        OUT_DIR + "logs/trimmomatic/trimmomatic.log"
    threads: 
        config["threads"]["fastqc"]
    params:
        trimmomatic_params = config["trimmomatic_params"]
    run:
        if (config["data_type"]["type"] == "pe"):
            
            shell("""
            {trimmomatic} PE \
            -threads {threads} \
            {input.forward} \
            {input.reverse} \
            {output.out}forward.trimmomatic.paired.fastq  \
            {output.out}forward.trimmomatic.unpaired.fastq \
            {output.out}reverse.trimmomatic.paired.fastq \
            {output.out}reverse.trimmomatic.unpaired.fastq \
            {params.trimmomatic_params} 2> {log}
            """)
            
        elif (config["data_type"]["type"] == "se"):
            
            shell("""
            {trimmomatic} SE \
            -threads {threads} \
            {input.single} \
            {output.out}single.trimmomatic.fastq \
            {params.trimmomatic_params} 2> {log}
            """)


# FILTERED FASTQC: This rule is use to generate an evaluation report on 
# filtered FASTQ files previously processed by Trimmomatic.
rule trim_fastqc:
    input:
        TRIMMOMATIC_DIR,
    output:
        fastqc_trimmed_out = FASTQC_TRIMMED_DIR
    log:
        OUT_DIR + "logs/fastqc/trimmed_fastqc.log"
    threads: 
        config["threads"]["fastqc"]
    run:
        if (config["data_type"]["type"] == "pe"):
            
            shell("""
            mkdir {output.fastqc_trimmed_out}
            
            {fastqc} \
            {input}/forward.trimmomatic.paired.fastq \
            {input}/reverse.trimmomatic.paired.fastq \
            --outdir {output.fastqc_trimmed_out} \
            --threads {threads} &> {log}
            """)
            
        elif (config["data_type"]["type"] == "se"):
            
            shell("""
            mkdir {output.fastqc_trimmed_out}
            
            {fastqc} \
            {input}/single.trimmomatic.fastq \
            --outdir {output.fastqc_trimmed_out} \
            --threads {threads} &> {log}
            """)


# TRINITY: This rule is use to de novo assemble filtered FASTQ files into 
# contigs. 
rule trinity:
    input:
        TRIMMOMATIC_DIR
    output:
        trinity_out = TRINITY_DIR + "/Trinity.fasta"
    log:
        OUT_DIR + "logs/trinity/trinity.log"
    params:
        max_memory = config["trinity_params"]["max_memory"],
        trinity_dir = TRINITY_DIR
    threads: 
        config["threads"]["trinity"]
    run:
        if (config["data_type"]["type"] == "pe"):
            
            shell("""
            {trinity} \
            --seqType fq \
            --left {input}forward.trimmomatic.paired.fastq \
            --right {input}reverse.trimmomatic.paired.fastq \
            --output {params.trinity_dir} \
            --CPU {threads} \
            --max_memory {params.max_memory} > {log}
            """)
            
        elif (config["data_type"]["type"] == "se"):
            
            shell("""
            {trinity} \
            --seqType fq \
            --single {input}single.trimmomatic.fastq \
            --output {params.trinity_dir} \
            --CPU {threads} \
            --max_memory {params.max_memory} > {log}
            """)


# TRANSRATE: This rule is use to generate an evaluation report on previously
# de novo assembled contigs.
rule transrate:
    input:
        TRINITY_DIR + "/Trinity.fasta",
    output:
        transrate_out = TRANSRATE_DIR
    log:
        OUT_DIR + "logs/transrate/transrate.log"
    params:
        trimmomatic_dir = TRIMMOMATIC_DIR
    threads: 
        config["threads"]["transrate"]
    run:
        if (config["data_type"]["type"] == "pe"):
            
            shell("""
            {transrate} \
            --assembly {input} \
            --left {params.trimmomatic_dir}/forward.trimmomatic.paired.fastq \
            --right {params.trimmomatic_dir}/reverse.trimmomatic.paired.fastq \
            --output {output.transrate_out} \
            --threads {threads} > {log}
            """)
            
        elif (config["data_type"]["type"] == "se"):
            
            shell("""
            {transrate} \
            --assembly {input} \
            --output {output.transrate_out} \
            --threads {threads} > {log}
            """)


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
