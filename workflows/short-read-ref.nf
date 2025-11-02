#!/usr/bin/env nextflow

include { multiqc } from '../modules/qc.nf'
include { calc_unmapped; get_unmapped_reads; bwa_index; bwa_index as bwa_index_plasmid; get_unmapped_reads as get_unmapped_reads_plasmid } from '../modules/mapping.nf'
include { freebayes; bcftools_stats; bcftools_stats as bcftools_stats_plasmid } from '../modules/variant_calling.nf'
include { qc; mapping; sv; annotate_vcf; mapping as mapping_plasmid } from '../modules/subworkflow.nf'
include { logUnmapped; logWorkflowCompletion; describePipeline } from '../modules/logs.nf'


workflow short_ref {
    take:
        trimmed
        fasta
        out_folder_name
        plasmid_fasta

    main:
        // mapping to the reference
        bwa_index(fasta, out_folder_name) | set { fasta_index } 
        mapping(fasta, fasta_index, trimmed, out_folder_name) | set { indexed_bam }

        // printout % unmapped reads
        calc_unmapped(indexed_bam) | set { pct }
        logUnmapped(pct, params.short_threshold, out_folder_name)

        get_unmapped_reads(indexed_bam, out_folder_name) | set { unmapped_fastq }

        // mapping reads to plasmid & variant calling
        if (plasmid_fasta) {
            Channel.from(plasmid_fasta) | set { mod_plasmid_fasta }

            bwa_index_plasmid(mod_plasmid_fasta, "${out_folder_name}-plasmid") | set { unmapped_fasta_index } 
            mapping_plasmid(mod_plasmid_fasta, unmapped_fasta_index, unmapped_fastq, "${out_folder_name}-plasmid") | set { unmapped_bam }
            get_unmapped_reads_plasmid(unmapped_bam, "${out_folder_name}-plasmid") | set { unmapped_fastq }

        }

         if (out_folder_name == "short-ref") { 
            // SNP & variant calling
            freebayes(fasta, fasta_index, indexed_bam, out_folder_name) | set { vcf }
            bcftools_stats(vcf, out_folder_name) | set { bcftools_out }
            
            // Optional annotation of vcf files
            def gtf_files = file("$params.in_dir").listFiles()?.findAll { it.name =~ /ref\.gtf$/ } ?: []
            def gff_files = file("$params.in_dir").listFiles()?.findAll { it.name =~ /ref\.gff$/ } ?: []

            if (gtf_files) {
                Channel.fromPath(gtf_files) | set { gtf }
                annotate_vcf(fasta, gtf, vcf, "gtf", "gtf22") | set {qc_vcf}
            
                qc_vcf.mix(bcftools_out).collect() | set { qc_out }
                multiqc(qc_out, out_folder_name, 'varint_calling')
            } else if (gff_files) {
                Channel.fromPath(gff_files) | set { gff }
                annotate_vcf(fasta, gff, vcf, "gff", "gff3") | set {qc_vcf}
            
                qc_vcf.mix(bcftools_out).collect() | set { qc_out }
                multiqc(qc_out, out_folder_name, 'varint_calling')
            } else {
                multiqc(bcftools_out, out_folder_name, 'varint_calling')
            }
        
            // SVs variant calling
            sv(fasta, indexed_bam, out_folder_name) | set { sv_vcf }
        } else {
            sv_vcf = Channel.empty()
        }
    emit:
        sv_vcf
        unmapped_fastq

}

out_folder_name = "short-ref"

workflow { 

    log.info  "Processing files in directory: ${params.in_dir}"
    
    def plasmid_files = file("$params.in_dir").listFiles()?.findAll { it.name =~ /ref_plasmid\.(fa|fna|fasta)$/ } ?: []

    Channel.fromPath("$params.in_dir/illumina/*.fastq.gz") | set { fastqs }
    Channel.fromPath("$params.in_dir/*ref*.{fa,fna,fasta}") | set { fasta }
    
    
    Channel.from(short_read_files).map { file ->
    def matcher = file.name =~ /^(.+?)(?:[_\.](S[0-9]+_L[0-9]+_)?(R[12]|[12]))?\.f(ast)?q\.gz$/
    if( matcher.matches() ) {
        [ matcher[0][1], file ]
        }
    }
    | filter { it }
    | groupTuple(sort: true)
    | set { fastqs }


    // QC and trimming module
    qc(fastqs, out_folder_name) | set { trimmed }

    short_ref(trimmed, fasta, "short-ref", plasmid_files)
}


logWorkflowCompletion(out_folder_name, !params.map_to_mod_fa)