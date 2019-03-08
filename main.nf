#!/usr/bin/env nextflow

params.reads = ""
params.pairedEnd = true
params.phred = 33
params.btindex = ''
params.bowtie = 'very-sensitive'
params.results = "./results"
params.mode = 'mapped'


bowtie_setting = ''


def helpMessage() {
    log.info"""
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
      --btindex                     Path to Bowtie2 Index (must be surrounded with dpouble quotes). Example: "/path/to/bowtieIndex/*.bt2"

    Options:
      --pairedEnd                   To specify wheter reads are single or paired end. Default = ${params.pairedEnd}
      --phred                       Specifies the fastq quality encoding (33 | 64). Defaults to ${params.phred}
      --mode                        Read anonymization mode. All reads are anonymized ('all') or only the mapped reads ('mapped'). Defaults to ${params.mode}

    Other options:
      --results                     Path of result directory. Defaults to ${params.results}
      --help  --h                   Shows this help page

    """.stripIndent()
}

// Show help message
params.help = false
params.h = false
if (params.help || params.h){
    helpMessage()
    exit 0
}

// Bowtie setting check
if (params.bowtie == 'very-fast'){
    bowtie_setting = '--very-fast'
} else if (params.bowtie == 'very-sensitive'){
    bowtie_setting = '--very-sensitive -N 1'
} else {
    println "Problem with --bowtie. Make sure to choose between 'very-fast' or 'very-sensitive'"
    exit(1)
}

Channel
    .fromPath(params.btindex)
    .ifEmpty {exit 1, 'Cannot find any bowtie2 index matching : ${params.btindex}\nPlease make sure to specify the index with "/path/to/bowtieIndexDir/*.bt2"'}
    .set {bt_index_genome}


Channel
    .fromFilePairs( params.reads, followLinks: true, size: params.pairedEnd ? 2 : 1,  )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\n" }
	.set {reads_to_trim}

process AdapterRemoval {
    tag "$name"

    conda 'bioconda::adapterremoval'

    errorStrategy 'ignore'

    label 'expresso'

    input:
        set val(name), file(reads) from reads_to_trim

    output:
        set val(name), file('*.trimmed.fastq') into trimmed_reads
        file("*.settings") into adapter_removal_results

    script:
        settings = name+".settings"
        if (params.pairedEnd){
            out1 = name+".pair1.trimmed.fastq"
            out2 = name+".pair2.trimmed.fastq"
            """
            AdapterRemoval --basename $name --file1 ${reads[0]} --file2 ${reads[1]} --trimns --trimqualities --minquality 20 --minlength 30 --output1 $out1 --output2 $out2 --threads ${task.cpus} --qualitybase ${params.phred} --settings $settings
            """
        } else {
            se_out = name+".trimmed.fastq"
            """
            AdapterRemoval --basename $name --file1 ${reads[0]} --trimns --trimqualities --minquality 20 --minlength 30 --output1 $se_out --threads ${task.cpus} --qualitybase ${params.phred} --settings $settings
            """
        }      
}

process Bowtie2Align {
    tag "$name"

    conda 'bioconda::bowtie2 bioconda::samtools'

    errorStrategy 'ignore'

    label 'intenso'

    publishDir "${params.results}/alignments", pattern: "*.flagstat.txt", mode: 'copy'

    input:
        set val(name), file(reads) from trimmed_reads
        file(index) from bt_index_genome.collect()
    output:
        set val(name), file("*.aligned.sorted.bam") into alignment_genome
        set val(name), file("*.flagstat.txt") into align1_multiqc
    script:
        index_name = index.toString().tokenize(' ')[0].tokenize('.')[0]
        samfile = name+".aligned.sam"
        fstat = name+".flagstat.txt"
        outfile = name+".aligned.sorted.bam"
        if (params.pairedEnd) {
            """
            bowtie2 -x $index_name -1 ${reads[0]} -2 ${reads[1]} $bowtie_setting --threads ${task.cpus} > $samfile
            samtools view -S -b -@ ${task.cpus} $samfile | samtools sort -@ ${task.cpus} -o $outfile
            samtools flagstat $samfile > $fstat
            """
        } else {
            """
            bowtie2 -x $index_name -U ${reads[0]} $bowtie_setting --threads ${task.cpus} > $samfile
            samtools view -S -b -@ ${task.cpus} $samfile | samtools sort -@ ${task.cpus} -o $outfile
            samtools flagstat $samfile > $fstat
            """
        }
}

process anonymize {
    tag "$name"

    conda 'python=3.6 bioconda::pysam'

    label 'ristretto'

    echo true

    publishDir "${params.results}/alignments", mode: 'copy'

    input:
        set val(name), file(sam) from alignment_genome
    output:
        set val(name), file("*.anonym.sam") into anonysam
    script:
        outfile = name+".anonym.sam" 
        """
        anonymize -m ${params.mode} -o $outfile $sam
        """
}