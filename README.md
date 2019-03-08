[![Build Status](https://www.travis-ci.com/maxibor/anonymap.svg?token=pwT9AgYi4qJY4LTp9WUy&branch=master)](https://www.travis-ci.com/maxibor/anonymap)

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
N E X T F L O W  ~  version 0.32.0
Launching `main.nf` [prickly_blackwell] - revision: 524a5b2fa5

=========================================
 Anonymap
 Homepage: https://github.com/maxibor/anonymap
 Author: Maxime Borry <borry[at]shh.mpg.de>
=========================================
Usage:
The typical command for running the pipeline is as follows:
nextflow run maxibor/anonymap --btindex "/path/to/bowtieIndex/*.bt2" --reads "/path/to/reads/*_R{1,2}.fastq.gz"
Mandatory arguments:
  --reads                       Path to input sequencing read file(s) (must be surrounded with double quotes). Example: "/path/to/reads/*_R{1,2}.fastq.gz"
  --btindex                     Path to Bowtie2 Index (must be surrounded with dpouble quotes). Example: "/path/to/bowtieIndex/*.bt2"

Options:
  --pairedEnd                   To specify wheter reads are single or paired end. Default = true
  --phred                       Specifies the fastq quality encoding (33 | 64). Defaults to 33
  --mode                        Read anonymization mode. All reads are anonymized ('all') or only the mapped reads ('mapped'). Defaults to mapped

Other options:
  --results                     Path of result directory. Defaults to ./results
  --help  --h                   Shows this help page
```

## Anonymisation

- Anonymap replaces the query sequence in the aligment `SAM` files by `NNNNNN...`
- Anonymap replaces the reference start position of every alignment by 1

These two steps combined impede anyone to have acess to the original sequence.  
By default only the mapped reads are anonymized.  
You can anonymize all reads by using `--mode 'all'`