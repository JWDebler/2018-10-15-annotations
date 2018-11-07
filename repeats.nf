#!/usr/bin/env nextflow

//This pipeline uses a reference genome (of high quality ideally)
//to create a repeat library which it then uses to mask
//repeats in the genomes we feed it through the input directory

//+++++++++++++++++ SETUP++++++++++++++++++++++++
params.workdir = "/home/johannes/rdrive/PPG_SEQ_DATA-LICHTJ-SE00182/johannes/notebook/2018-10-15-annotations"

//Reference genome fasta file//params.reference = "${params.workdir}/input/206-15.contigs.nucl.fasta"

//Folder containing fasta files of genomes to be masked
//single input:
//params.input = "${params.workdir}/input/206-15.contigs.nucl.fasta"
//multiple inputs:
params.input = "${params.workdir}/input/*contigs.nucl.dos2unix.fasta"
//Output
params.outdir = "${params.workdir}/output/"

//Minimum repeat length
params.minRepLength = 250
//+++++++++++++++++++++++++++++++++++++++++++++++

//reference = file(params.reference)

//Creating a channel containing all the sequences we want masked
sequences = Channel
.fromPath(params.input)
.map{ file -> tuple(file.simpleName, file)}
.set{genomes}


log.info "====================================================================="
log.info "Repeat annotation using a reference genome"
//log.info "Sample : ${reference}"
log.info "Output    : ${params.outdir}"
log.info "====================================================================="


//Creating a repeat library from the reference genome
process RepeatModeler {
  container 'robsyme/nf-repeatmasking'
  publishDir "${params.outdir}/${sampleID}/repeatlibrary", mode: 'copy', pattern: '*.classified'
  tag {sampleID}
  cpus 2

  input:
  //file 'genome.fasta' from reference
  set sampleID, 'genome.fasta' from genomes

  output:
  set sampleID, 'genome.fasta', 'RM*/consensi.fa.classified' into repeatlibrary

  """
BuildDatabase \
 -name reference \
 genome.fasta

RepeatModeler \
 -database reference \
 -pa ${task.cpus}
  """
}


//Masking repeats using the library created above
process RepeatMasker {
  container 'robsyme/nf-repeatmasking'
  tag {sampleID}
  publishDir "${params.outdir}/${sampleID}", mode: 'copy', pattern: '*.tbl'
  publishDir "${params.outdir}/${sampleID}", mode: 'copy', pattern: '*.gff'
  
  cpus 8

  input:
  set sampleID, "genome.fasta", "library.fasta" from repeatlibrary

  output:
  set sampleID, '*.fasta.out', 'genome.fasta'  into repeatmaskerout
  set "${sampleID}.repeats.repeatmodeler.tbl", "${sampleID}.repeats.repeatmodeler.gff"

  """
RepeatMasker \
 -no_is \
 -nolow \
 -gff \
 -pa ${task.cpus} \
 -lib library.fasta \
 genome.fasta
 mv *.tbl ${sampleID}.repeats.repeatmodeler.tbl
 mv *.gff ${sampleID}.repeats.repeatmodeler.gff
  """
}

process removeShortMatches {
  tag {sampleID}
  container 'robsyme/nf-repeatmasking'
  publishDir "${params.outdir}/${sampleID}", mode: 'copy'
  
  input:
  set sampleID, "rm.out", 'genome.fasta' from repeatmaskerout
  
  output:
  set sampleID, "${sampleID}.repeats.repeatmodeler.out", "${sampleID}*fasta" into repeatMaskerKnowns


  """
head -n 3 rm.out > ${sampleID}.repeats.repeatmodeler.out
tail -n +4 rm.out | awk '\$7 - \$6 > ${params.minRepLength}' >> ${sampleID}.repeats.repeatmodeler.out
tail -n +4 rm.out | awk 'BEGIN{OFS="\\t"} \$7 - \$6 > ${params.minRepLength} {print \$5, \$6, \$7}' >> tmp.bed
maskFastaFromBed -fi genome.fasta -bed tmp.bed -fo ${sampleID}.repeats.repeatmodeler.masked.soft.fasta -soft
maskFastaFromBed -fi genome.fasta -bed tmp.bed -fo ${sampleID}.repeats.repeatmodeler.masked.hard.fasta
  """
}

workflow.onComplete {
    log.info "========================================================"
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'Failed' }"
    log.info "========================================================"
}