# The dntap (*de novo* transcriptome analysis pipeline)

## Description
This pipeline has been designed to process **RNA-seq data**. It allows to obtain results and reports from your data very easily by simply fill a configuration file and launch a command line. 
It takes advantages of **[Snakemake](https://snakemake.readthedocs.io/en/stable/)** workflow engine. You can refer to the Snakemake **[publication](https://academic.oup.com/bioinformatics/article/28/19/2520/290322/Snakemake-a-scalable-bioinformatics-workflow)** for more details.

The pipeline is composed of 6 steps:
1. `.fastq` evaluation (**[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**)
2. Quality filtering `.fastq` (**[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)**)
3. Newly filtered `.fastq` evaluation (**[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**)
4. *De novo* transcriptome assembly (**[Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)**)
5. Assembly metrics evaluation (**[Transrate](http://hibberdlab.com/transrate/)**)
6. Protein coding domains prediction (**[Transdecoder](https://transdecoder.github.io/)**)
7. Functional annotation of prediction protein coding domains (**[InterProScan 5](https://github.com/ebi-pf-team/interproscan/wiki/HowToRun)**)

Note that this version of the pipeline only supported paired-end `.fastq`

![pipeline diagram](https://github.com/arnaudmeng/dntap/extra/diag1.png)

## Quick example

**dnatp** takes 2 arguments: 
1. The snakefile `dntap.py` 
2. The configuration file `dntap_config.yaml`

```bash
snakemake --snakefile dntap.py --configfile dntap_config.yaml --cores 20
```
+ It is possible to provided a naximum of cores to be used by `dntap.py` steps, here: `--cores 20`.

## Pre-requirement and installation

The pipeline is provided with ready-to-use binairies and sources for:

+ **FastQC** (v0.11.5)
+ **Trimmomatic** (v0.36)
+ **Transrate** (v1.0.3)
+ **TransDecoder** (v3.0.1)

However it requires installation for:

+ **[Snakemake](http://snakemake.readthedocs.io/en/stable/getting_started/installation.html)**
+ **[Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Installing-Trinity)** (2.4.0)
+ **[InterProScan](https://github.com/ebi-pf-team/interproscan/wiki/HowToDownload)** (v5.24-63)

Note that the current version of the pipeline do not take advantages of the `ipr_lookup service` of **InterProScan**. 
However you need to download and install the `Panther` database if you want InterProScan to perform research on it.       

## Configuration file

### Software location

To use the pipeline proprely you must start with filling software location in the `dntap_config.yaml`
```h
software:
    # You must provide absolute path to the executable file
    # e.g. fastqc: /path/to/fastqc/install/directory/fastqc
    
    fastqc: path/to/src/FastQC/fastqc
    trimmomatic: path/to/src/Trimmomatic-0.36/trimmomatic-0.36.jar
    trinity: /path/to/src/trinityrnaseq-Trinity-v2.4.0/Trinity
    transrate: /path/to/src/transrate-1.0.3-linux-x86_64/transrate
    transdecoder_longorfs: /path/to/src/TransDecoder-3.0.1/TransDecoder.LongOrfs
    transdecoder_predict: /path/to/src/TransDecoder-3.0.1/TransDecoder.Predict
    interproscan: /path/to/src/interproscan-5.24-63.0/interproscan.sh
```

If **Trinity** install location is `/usr/local/trinityrnaseq-Trinity-v2.4.0/Trinity`
Simply remplace the absolute path location in the `dntap_config.yaml` as follow:
```h
trinity: /usr/local/trinityrnaseq-Trinity-v2.4.0/Trinity
```

### Inputs and parameters

Finally you can provide `.fastq` files location to be process by the pipeline.
```h
# You must provide absolute path to paired-end RNA-seq file (.fastq / .fq).

forward: /path/to/sample/reads.left.fq
reverse: /path/to/sample/reads.right.fq
```

Change parameters if needed:
```h
trimmomatic_params:
    MINLEN:32 SLIDINGWINDOW:10:20 LEADING:5 TRAILING:5
    
trinity_params:
    max_memory: 20G
    
transdecoder_params:
    min_protein_len: 100
    
interproscan_params:
    out_format: tsv
    db: TIGRFAM, SFLD, ProDom, Hamap, SMART, CDD, ProSiteProfiles, ProSitePatterns, SUPERFAMILY, PRINTS, PANTHER, Gene3D, PIRSF, Pfam, Coils
```

You also have the possibility to set maximum threads to be use by each step of the pipeline:
```h
threads:
    # You can set maximum threads to be use be each step.

    fastqc: 20
    trimmomatic: 6
    trinity: 20
    transrate: 20
    transdecoder: 20
    interproscan: 20
```

Note that it cannot exceed the maximum number provided in the command line `--cores 20`.
