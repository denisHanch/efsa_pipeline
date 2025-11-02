#!/usr/bin/env nextflow

include { nanoplot; multiqc } from '../modules/qc.nf'
include { sv_long; mapping_long; mapping_long as mapping_long_plasmid; sv_long as sv_long_plasmid }  from '../modules/subworkflow.nf'
include { logUnmapped; logWorkflowCompletion; loadFastqFiles } from '../modules/logs.nf'
include { calc_unmapped; calc_unmapped as calc_unmapped_plasmid; get_unmapped_reads;get_unmapped_reads as get_unmapped_reads_plasmid } from '../modules/mapping.nf'


workflow long_ref {

    take:
        fastqs
        fasta
        mapping_tag
        plasmid_fasta
        out_folder_name

    main:
        // qc
        nanoplot(fastqs, out_folder_name)
        
        // mapping to the reference
        mapping_long(fastqs, fasta, mapping_tag, out_folder_name) | set { indexed_bam }

         // printout % unmapped reads
        calc_unmapped(indexed_bam) | set { pct }
        logUnmapped(pct, params.long_threshold,  "${out_folder_name}-${mapping_tag}")
        get_unmapped_reads(indexed_bam, out_folder_name) | set { unmapped_fastq }
        
        // mapping reads to plasmid & variant calling
        if (plasmid_fasta) {
            Channel.from(plasmid_fasta) | set { plasmid_fasta }

            mapping_long_plasmid(unmapped_fastq, plasmid_fasta, mapping_tag, "${out_folder_name}-plasmid") | set { unmapped_bam }
            get_unmapped_reads_plasmid(unmapped_bam, "${out_folder_name}-plasmid") | set { unmapped_fastq }
        }

        // SV calling against the reference
        if (out_folder_name == "long-ref") { 
            sv_long(fasta, indexed_bam, mapping_tag, out_folder_name) | set { sv_vcf }
        } else {
            sv_vcf = Channel.empty()
        }

    emit:
        sv_vcf
        unmapped_fastq
}

out_folder_name = "long-ref"

workflow {
    // Processing inputs
    log.info  "Processing files in directory: ${params.in_dir}"

    pacbio_fastqs = loadFastqFiles("${params.in_dir}/pacbio/*.fastq.gz")
    ont_fastqs = loadFastqFiles("${params.in_dir}/ont/*.fastq.gz")

    def ref_plasmid = file("$params.in_dir").listFiles()?.findAll { it.name =~ /ref_plasmid\.(fa|fna|fasta)$/ } ?: []
    def mod_plasmid = file("$params.in_dir").listFiles()?.findAll { it.name =~ /mod_plasmid\.(fa|fna|fasta)$/ } ?: []

    Channel.fromPath("$params.in_dir/*ref.{fa,fna,fasta}", checkIfExists: true) | set { ref_fasta }
    
    if (pacbio_fastqs) {
        long_ref(pacbio_fastqs, ref_fasta, "map-pb", ref_plasmid, out_folder_name)
    }
}

logWorkflowCompletion(out_folder_name)