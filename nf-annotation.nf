#!/usr/bin/env nextflow

//+++++++++++++++++ SETUP++++++++++++++++++++++++
params.workdir = "/home/johannes/rdrive/PPG_SEQ_DATA-LICHTJ-SE00182/johannes/notebook/2018-10-15-annotations"

//Folder containing fasta files of genomes to be annotated
params.input = "${params.workdir}/input/*contigs.nucl.fasta"

params.outputdir = "${params.workdir}/output"

//point this to the scripts directory of this repository
params.scripts = "/home/johannes/rdrive/PPG_SEQ_DATA-LICHTJ-SE00182/johannes/notebook/2018-10-15-annotations/scripts"

//Minimum repeat length
params.minRepLength = 250
//+++++++++++++++++++++++++++++++++++++++++++++++


log.info "==================================================="
log.info "|| Annotation of a folder full of FASTA files    ||"      
log.info "|| using GeneMarkES, infernal and tRNAscanSE     ||"
log.info "|| RepeatMasker ...";
log.info "|| Output  : ${params.outputdir}/"
log.info "=============================================================================="


//Creating a channel containing all the sequences we want annotated
sequences = Channel
.fromPath(params.input)
.map{ file -> tuple(file.simpleName, file)}
.tap{genomesFordos2unix}

//Making sure there are no strange invisible characters left that would mess with programs down the road
//and also cleaning up fasta headers with 'sed', because they contain spaces wich mess
//with other tools down the road
process dos2unix {
    tag {sampleID}
    publishDir "${params.workdir}/input", mode: 'copy'

    input:
    set sampleID, 'genome.fasta' from genomesFordos2unix

    output:
    set sampleID, "${sampleID}.contigs.nucl.dos2unix.fasta" into dos2unixed

    """
    dos2unix -n genome.fasta dos2unix.fasta
    sed 's, .*\$,,g' -i dos2unix.fasta
    mv dos2unix.fasta ${sampleID}.contigs.nucl.dos2unix.fasta
    """
}


dos2unixed
.tap{genomesForGenemark}
.tap{genomesForTRNAscanSE}
.tap{genomesForInfernal}
.tap{genomesForRepeats}
.tap{genomesForAugustus}

//Repeatmasking
process RepeatMasker {
  container 'robsyme/nf-repeatmasking'
  tag {sampleID}
  publishDir "${params.outputdir}/${sampleID}/", mode: 'copy', pattern: '*.tbl'
  publishDir "${params.outputdir}/${sampleID}/", mode: 'copy', pattern: '*.gff'
  cpus 8

  input:
  set sampleID, "${sampleID}.fasta" from genomesForRepeats

  output:
  set "${sampleID}.repeats.tbl", "${sampleID}.repeats.gff"
  set sampleID, '*.fasta.out', "${sampleID}.fasta"  into repeatmaskerout

  """
RepeatMasker \
 -no_is \
 -nolow \
 -gff \
 -pa ${task.cpus} \
 -species fungi \
 ${sampleID}.fasta
 mv *.tbl ${sampleID}.repeats.tbl
 mv *.gff ${sampleID}.repeats.gff
  """
}

process cleanUpRempeatMaskerOutput {
  tag {sampleID}
  container 'robsyme/nf-repeatmasking'
  publishDir "${params.outputdir}/${sampleID}", mode: 'copy'
  
  input:
  set sampleID, "rm.out", 'genome.fasta' from repeatmaskerout
  
  output:
  set sampleID, "${sampleID}.repeats.out", "${sampleID}*fasta" into repeatMaskerKnowns


  """
head -n 3 rm.out > ${sampleID}.repeats.out
tail -n +4 rm.out | awk '\$7 - \$6 > ${params.minRepLength}' >> ${sampleID}.repeats.out
tail -n +4 rm.out | awk 'BEGIN{OFS="\\t"} \$7 - \$6 > ${params.minRepLength} {print \$5, \$6, \$7}' >> tmp.bed
maskFastaFromBed -fi genome.fasta -bed tmp.bed -fo ${sampleID}.repeats.masked.soft.fasta -soft
maskFastaFromBed -fi genome.fasta -bed tmp.bed -fo ${sampleID}.repeats.masked.hard.fasta
  """
}

//GenemarkES annotation
process annotation_genemark {
    tag {sampleID}
    cpus 10
    publishDir "${params.outputdir}/${sampleID}/", mode: 'copy', pattern: '*.gtf'

    input:
	set sampleID, "${sampleID}.*.fasta" from genomesForGenemark

    output:
	set sampleID, "${sampleID}.genemark.gtf", "${sampleID}..fasta" into genemarkOutput

    """
    /opt/genemark-ES/gmes_petap.pl --ES --fungus --cores ${task.cpus} --sequence ${sampleID}.*.fasta
    mv genemark.gtf ${sampleID}.genemark.gtf
    """
}

//pull proteins out of genemark annotation
process extractProteinsFromGenemark {
  publishDir "${params.outputdir}/${sampleID}", mode: 'copy'
  tag {sampleID}

  input:
  set sampleID, "${sampleID}.genemark.gtf", "input.fasta" from genemarkOutput

  output:
  set sampleID, "${sampleID}.genemark.proteins.fasta" into proteinsFromGenemark

  """
  /opt/genemark-ES/get_sequence_from_GTF.pl ${sampleID}.genemark.gtf input.fasta
  mv prot_seq.faa ${sampleID}.genemark.proteins.fasta
  """
}

proteinsFromGenemark
.tap {proteinsForEffectorP}
.tap {proteinsForDeepsig}
.tap {proteinsForInterproscan}
.tap {proteinsForX}

//tRNA annotation with tRNAscanSE
process annotation_trnascan {
    tag {sampleID}
    publishDir "${params.outputdir}/${sampleID}/", mode: 'copy'

    input:
	set sampleID, "${sampleID}.*.fasta" from genomesForTRNAscanSE

    output:
	set sampleID, "${sampleID}.trnascanSE.gff3"
    """
    /home/johannes/local/bin/tRNAscan-SE -o trnascanoutput.out ${sampleID}.*.fasta 
    ${params.scripts}/convert_tRNAScanSE_to_gff3.pl --input=trnascanoutput.out > ${sampleID}.trnascanSE.gff3
    """
}

//RNA annotation with infernal
process annotaton_infernal {
    tag {sampleID}
    publishDir "${params.outputdir}/${sampleID}/", mode: 'copy'
    cpus 10
    
    input:
	set sampleID, "${sampleID}.*.fasta" from genomesForInfernal

    output:
	set sampleID, "${sampleID}.infernal.gff3"
    """
    wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/13.0/Rfam.clanin
    wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/13.0/Rfam.cm.gz
    gzip -d Rfam.cm.gz
    cmpress Rfam.cm
    cmscan --rfam --cpu ${task.cpus} --cut_ga --nohmmonly --tblout contigs.cmscan.tbl --fmt 2 --clanin Rfam.clanin Rfam.cm ${sampleID}.*.fasta > infernal.cmscan
    grep -v ^# contigs.cmscan.tbl > contigs.cmscan.clean.tbl && awk '{printf "%s\tinfernal\t%s\t%d\t%d\t%s\t%s\t.\tNote=RfamID-%s\\n" ,\$4,\$2,\$10,\$11,\$17,\$12,\$3}'  contigs.cmscan.clean.tbl > ${sampleID}.infernal.gff3
    """
}


//effector prediction with effectorP
process effectorP {
  tag { sampleID }
  publishDir "${params.outputdir}/${sampleID}", mode: 'copy'

  input:
  set sampleID, "input.fasta" from proteinsForEffectorP

  output:
  file "${sampleID}.effectorP.tsv"

  """
  /opt/effectorP/EffectorP_2.0/Scripts/EffectorP.py \
  -s \
  -o ${sampleID}.effectorP.tsv \
  -i input.fasta
  """
}

//signal peptide prediction with deepsig
process deepsig {
  tag { sampleID }
  publishDir "${params.outputdir}/${sampleID}", mode: 'copy'

  input:
  set sampleID, "input.fasta" from proteinsForDeepsig

  output:
  file "${sampleID}.deepsig.out"

  """
  runDeepSig.sh input.fasta euk ${sampleID}.deepsig.out 8
  """

}

//protein property annotation with interproscan including SignalP 
process interproscan {
  tag { sampleID }
  publishDir "${params.outputdir}/${sampleID}", mode: 'copy'
  cpus 12

  input:
  set sampleID, "proteins.fasta" from proteinsForInterproscan

  output:
  file "${sampleID}.interproscan.tsv"

  """
  interproscan.sh \
  --applications SignalP_EUK,Pfam,TMHMM,PANTHER,PRINTS,ProDom,ProSitePatterns,ProSiteProfiles,MobiDBLite\
  --cpu ${task.cpus} \
  --seqtype p \
  --disable-precalc \
  --goterms \
  --pathways \
  --iprlookup\
  --input proteins.fasta \
  --output-file-base ${sampleID}.interproscan \
  --format tsv
  """

}

workflow.onComplete {
    log.info "========================================================"
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'Failed' }"
    log.info "========================================================"
}
