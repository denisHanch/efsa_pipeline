#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { fastqc }    from './modules/fastqc/fastqc.nf'
include { trimming }  from './modules/trimming/trimming.nf'
include { multiqc }   from './modules/multiqc/multiqc.nf'
include { validate }  from './modules/validate/validate.nf'

workflow {

  // parse user’s chosen steps
  def stepsList = params.steps.tokenize(',')

  // load all FASTQ files
  Channel
    .fromPath("${params.in_dir}/*.fastq")
    .set { samples_ch }

  // chain steps conditionally
  def chA = stepsList.contains('fastqc')    ? fastqc(samples_ch)   : samples_ch
  def chB = stepsList.contains('trimming')  ? trimming(chA)        : chA
  def chC = stepsList.contains('multiqc')   ? multiqc(chB)         : chB

  // publish whatever’s left into data/outputs
  chC.publishDir params.out_dir, mode: 'copy'
}
