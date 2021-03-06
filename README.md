[![Build Status](https://www.travis-ci.com/maxibor/anonymap.svg?token=pwT9AgYi4qJY4LTp9WUy&branch=master)](https://www.travis-ci.com/maxibor/anonymap) [![DOI](https://zenodo.org/badge/174518757.svg)](https://zenodo.org/badge/latestdoi/174518757)


# Anonymap
- Trimming ([AdapterRemoval](https://github.com/MikkelSchubert/adapterremoval)) 
- Mapping ([Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)) 
- Anonymizing (Python + [Pysam](https://pysam.readthedocs.io/en/latest/))


## Introduction

Anonymap is a simple nextflow pipeline for trimming, aligning, and anonymizing the sam file to lift off data access restrictions.

## Dependencies
- [Conda](https://conda.io/en/latest/miniconda.html)
    - For Mac 
    ```
    $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    $ chmod +x https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    $ ./https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    ```
    - For Linux
    ```
    $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $ chmod +x https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $ ./https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    ```

- [Nextflow](https://www.nextflow.io/)  
    ```
    $ conda install -c bioconda nextflow
    ```



## Help 

```
$ nextflow run maxibor/anonymap --h
N E X T F L O W  ~  version 19.04.0
Launching `main.nf` [romantic_woese] - revision: dca0e16c86
=========================================
 Anonymap
 Homepage: https://github.com/maxibor/anonymap
 Author: Maxime Borry <borry@shh.mpg.de>
=========================================
Usage:
The typical command for running the pipeline is as follows:
nextflow run maxibor/anonymap --btindex "/path/to/bowtieIndex/*.bt2" --reads "/path/to/reads/*_R{1,2}.fastq.gz"
Mandatory arguments:
  --reads                       Path to input sequencing read file(s) (must be surrounded with double quotes). Example: "/path/to/reads/*_R{1,2}.fastq.gz"
  --btindex                     Path to Bowtie2 Index (must be surrounded with dpouble quotes). Example: "/path/to/bowtieIndex/genome_basename*"

Options:
  --pairedEnd                   To specify wheter reads are single or paired end. Default = true
  --phred                       Specifies the fastq quality encoding (33 | 64). Defaults to 33
  --mode                        Read anonymization mode. All reads are anonymized ('all') or only the mapped reads ('mapped'). Defaults to mapped

Other options:
  --results                     Path of result directory. Defaults to ./results
  --help  --h                   Shows this help page
```

## Output
- `result_directory/<sample_name>.flagstat.txt`: Samtools flagstat alignment statistics
- `result_directory/<sample_name>.anonym.bam`: Anonymized alignment bam file.
- `result_directory/pipeline_info/anonymap_log.txt`: Anonymap log file

## Anonymisation

Anonymap replaces the query sequence in the aligment `SAM` files by `NNNNNN...`  

By default only the mapped reads are anonymized.  
You can anonymize all reads by using `--mode 'all'`

## Getting the Bowtie2 Index of reference genome

Pre-indexed genomes for Bowtie2 are provided by Illumina on their [iGenome server](https://support.illumina.com/sequencing/sequencing_software/igenome.html)