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
      --btindex                     Path to Bowtie2 Index (must be surrounded with dpouble quotes). Example: "/path/to/bowtieIndex/genome_basename*"

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

def summary = [:]
summary["reads"] = params.reads
summary["bowtie index"] = params.btindex
summary["Read mode"] = params.pairedEnd
summary["PHRED"] = params.phred
summary["Anonymization mode"] = params.mode
summary["result directory"] = params.results
log.info summary.collect { k,v -> "${k.padRight(25)}: $v" }.join("\n")
log.info "\033[2m----------------------------------------------------\033[0m"

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

    publishDir "${params.results}", pattern: "*.flagstat.txt", mode: 'copy'

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
            samtools flagstat $outfile > $fstat
            """
        } else {
            """
            bowtie2 -x $index_name -U ${reads[0]} $bowtie_setting --threads ${task.cpus} > $samfile
            samtools view -S -b -@ ${task.cpus} $samfile | samtools sort -@ ${task.cpus} -o $outfile
            samtools flagstat $outfile > $fstat
            """
        }
}

process anonymize {
    tag "$name"

    conda 'python=3.6 bioconda::pysam'

    label 'ristretto'

    echo true

    publishDir "${params.results}", mode: 'copy'

    input:
        set val(name), file(sam) from alignment_genome
    output:
        set val(name), file("*.anonym.bam") into anonysam
    script:
        outfile = name+".anonym.bam" 
        """
        anonymize -m ${params.mode} -o $outfile $sam
        """
}

workflow.onComplete {

    def report_fields = [:]
    report_fields['runName'] = custom_runName ?: workflow.runName
    report_fields['success'] = workflow.success
    report_fields['dateComplete'] = workflow.complete
    report_fields['duration'] = workflow.duration
    report_fields['exitStatus'] = workflow.exitStatus
    report_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    report_fields['errorReport'] = (workflow.errorReport ?: 'None')
    report_fields['commandLine'] = workflow.commandLine
    report_fields['projectDir'] = workflow.projectDir
    report_fields['summary'] = summary
    report_fields['summary']['Date Started'] = workflow.start
    report_fields['summary']['Date Completed'] = workflow.complete
    report_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    report_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) report_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) report_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) report_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) report_fields['summary']['Docker image'] = workflow.container
    report_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    report_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    report_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/report_template.txt")
    def txt_template = engine.createTemplate(tf).make(report_fields)
    def report_txt = txt_template.toString()

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.results}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_tf = new File( output_d, "anonymap_log.txt" )
    output_tf.withWriter { w -> w << report_txt }
}