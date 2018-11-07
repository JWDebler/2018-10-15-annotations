//+++++++++++++++++ SETUP++++++++++++++++++++++++
params.workdir = "/home/johannes/rdrive/PPG_SEQ_DATA-LICHTJ-SE00182/johannes/notebook/2018-10-15-annotations"

//Folder containing fasta files of genomes to be annotated
params.repeatmasker = "${params.workdir}/output/*/*.repeats.masked.soft.fasta"
params.repeatmodeer = "${params.workdir}/output/*/*.repeatmodeler.masked.soft.fasta"

params.outputdir = "${params.workdir}/output"

//point this to the scripts directory of this repository
params.scripts = "/home/johannes/rdrive/PPG_SEQ_DATA-LICHTJ-SE00182/johannes/notebook/2018-10-15-annotations/scripts"

//Minimum repeat length
params.minRepLength = 250
//+++++++++++++++++++++++++++++++++++++++++++++++

//Creating a channel containing all the sequences we want annotated
sequencesRepeatmasker = Channel
.fromPath(params.input)
.map{ file -> tuple(file.simpleName, file)}
.tap{repeatmaskerMasked}


sequencesRepeatmodeler = Channel
.fromPath(params.input)
.map{ file -> tuple(file.simpleName, file)}
.tap{repeatmodelerMasked}


//GenemarkES annotation
process annotation_genemark_repeatmasker {
    tag {sampleID}
    cpus 10
    publishDir "${params.outputdir}/${sampleID}/", mode: 'copy', pattern: '*.gtf'

    input:
	set sampleID, "${sampleID}.*.fasta" from repeatmaskerMasked

    output:
	set sampleID, "${sampleID}.genemark.repeatmasker.gtf", "${sampleID}..fasta" into genemarkOutputRepeatmasker

    """
    /opt/genemark-ES/gmes_petap.pl --ES --fungus --soft_mask 1 --cores ${task.cpus} --sequence ${sampleID}.*.fasta
    mv genemark.gtf ${sampleID}.genemark.repeatmasker.gtf
    """
}

//GenemarkES annotation
process annotation_genemark_repeatmodeler {
    tag {sampleID}
    cpus 10
    publishDir "${params.outputdir}/${sampleID}/", mode: 'copy', pattern: '*.gtf'

    input:
	set sampleID, "${sampleID}.*.fasta" from repeatmodelerMasked

    output:
	set sampleID, "${sampleID}.genemark.repeatmodeler.gtf", "${sampleID}..fasta" into genemarkOutputRepeatmodeler

    """
    /opt/genemark-ES/gmes_petap.pl --ES --fungus --soft_mask 1 --cores ${task.cpus} --sequence ${sampleID}.*.fasta
    mv genemark.gtf ${sampleID}.genemark.repeatmodeler.gtf
    """
}

//pull proteins out of genemark annotation
process extractProteinsFromGenemarkRepeatmasker {
  publishDir "${params.outputdir}/${sampleID}", mode: 'copy'
  tag {sampleID}

  input:
  set sampleID, "${sampleID}.genemark.gtf", "input.fasta" from genemarkOutputRepeatmasker

  output:
  set sampleID, "${sampleID}.genemark.proteins.repeatmasker.fasta" into proteinsFromGenemarkRepeatmasker

  """
  /opt/genemark-ES/get_sequence_from_GTF.pl ${sampleID}.genemark.gtf input.fasta
  mv prot_seq.faa ${sampleID}.genemark.proteins.repeatmasker.fasta
  """
}

//pull proteins out of genemark annotation
process extractProteinsFromGenemarkRepeatmodeler {
  publishDir "${params.outputdir}/${sampleID}", mode: 'copy'
  tag {sampleID}

  input:
  set sampleID, "${sampleID}.genemark.gtf", "input.fasta" from genemarkOutputRepeatmodeler

  output:
  set sampleID, "${sampleID}.genemark.proteins.repeatmodeler.fasta" into proteinsFromGenemarkRepeatmodeler

  """
  /opt/genemark-ES/get_sequence_from_GTF.pl ${sampleID}.genemark.gtf input.fasta
  mv prot_seq.faa ${sampleID}.genemark.proteins.repeatmodeler.fasta
  """
}


workflow.onComplete {
    log.info "========================================================"
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'Failed' }"
    log.info "========================================================"
}