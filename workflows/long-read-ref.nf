#!/usr/bin/env nextflow

include { nanoplot; multiqc } from '../modules/qc.nf'
include { sv_long; mapping_long; mapping_long as mapping_long_plasmid; sv_long as sv_long_plasmid }  from '../modules/subworkflow.nf'
include { logUnmapped; logUnmapped as logUnmapped_plasmid; logWorkflowCompletion; loadFastqFiles } from '../modules/logs.nf'
include { calc_unmapped; calc_unmapped as calc_unmapped_plasmid; calc_total_reads; get_unmapped_reads;get_unmapped_reads as get_unmapped_reads_plasmid } from '../modules/mapping.nf'


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

        get_unmapped_reads(indexed_bam, out_folder_name) | set { unmapped_fastq }
        
        // printout % unmapped reads
        calc_total_reads(indexed_bam) | set { total_reads }
        calc_unmapped(unmapped_fastq) | set { nreads }
        logUnmapped(nreads, total_reads, out_folder_name, "")

        
        // mapping reads to plasmid & variant calling
        if (plasmid_fasta) {
            Channel.from(plasmid_fasta) | set { plasmid_fasta }

            mapping_long_plasmid(unmapped_fastq, plasmid_fasta, mapping_tag, "${out_folder_name}-plasmid") | set { unmapped_bam }
            get_unmapped_reads_plasmid(unmapped_bam, "${out_folder_name}-plasmid") | set { unmapped_fastq }

            calc_unmapped_plasmid(unmapped_fastq) | set { nreads }
            logUnmapped_plasmid(nreads, total_reads, "${out_folder_name}-plasmid", "against plasmid")
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

    Channel.fromPath("$params.in_dir/*{ref,reference_genome}.{fa,fna,fasta}", checkIfExists: true) | set { ref_fasta }
    
    if (pacbio_fastqs) {
        long_ref(pacbio_fastqs, ref_fasta, "map-pb", ref_plasmid, out_folder_name)
    }
}

logWorkflowCompletion(out_folder_name)