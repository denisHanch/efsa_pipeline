#!/usr/bin/env nextflow

include { multiqc } from '../modules/qc.nf'
include { calc_unmapped; bwa_index; bwa_index as bwa_index_plasmid; get_unmapped_reads } from '../modules/mapping.nf'
include { freebayes; bcftools_stats; freebayes as freebayes_plasmid; bcftools_stats as bcftools_stats_plasmid } from '../modules/variant_calling.nf'
include { qc; mapping; sv; annotate_vcf; mapping as mapping_plasmid; sv as sv_plasmid } from '../modules/subworkflow.nf'
include { logUnmapped; logWorkflowCompletion } from '../modules/logs.nf'


out_folder_name = "short-ref"

workflow short_ref {
    take:
        trimmed
        fasta

    main:        
        // mapping to the reference
        bwa_index(fasta, out_folder_name) | set { fasta_index } 
        mapping(fasta, fasta_index, trimmed, out_folder_name) | set { indexed_bam }

        // printout % unmapped reads
        calc_unmapped(indexed_bam) | set { pct }
        logUnmapped(pct, params.short_threshold, out_folder_name)

        // mapping reads to plasmid & variant calling
        def plasmid_files = file("$params.in_dir").listFiles()?.findAll { it.name =~ /plasmid\.(fa|fna|fasta)$/ } ?: []

        if (plasmid_files) {
            Channel.from(plasmid_files) | set { mod_plasmid_fasta }

            get_unmapped_reads(indexed_bam, "short-ref-plasmid") | set { unmapped_fastq }
            bwa_index_plasmid(mod_plasmid_fasta, out_folder_name) | set { unmapped_fasta_index } 
            mapping_plasmid(mod_plasmid_fasta, unmapped_fasta_index, unmapped_fastq, "short-ref-plasmid") | set { unmapped_bam }
            sv_plasmid(mod_plasmid_fasta, unmapped_bam, "short-ref-plasmid")
            freebayes_plasmid(mod_plasmid_fasta, unmapped_fasta_index, unmapped_bam, "short-ref-plasmid") | set { vcf_plasmid }
            bcftools_stats_plasmid(vcf_plasmid, "short-ref-plasmid")
        }

        // SNPs variant calling against the reference
        freebayes(fasta, fasta_index, indexed_bam, out_folder_name) | set { vcf }
        bcftools_stats(vcf, out_folder_name) | set { bcftools_out }
        
        // Annotate SNPs & QC

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
            multiqc(qc_vcf, out_folder_name, 'varint_calling')
        }
    
        // SVs variant calling against the reference
        sv(fasta, indexed_bam, out_folder_name) | set { sv_vcf }
    
    emit:
        sv_vcf

}


workflow { 
    log.info  "Processing files in directory: ${params.in_dir}"
    
    Channel.fromPath("$params.in_dir/illumina/*.fastq.gz") | set { fastqs }
    Channel.fromPath("$params.in_dir/*ref*.{fa,fna,fasta}") | set { fasta }
    
    fastqs
    | map {[(it.name =~ /^([^_]+)(_((S[0-9]+_L[0-9]+_)?R[12]_001|[12]))?\.fastq.gz/)[0][1], it]}
    | groupTuple(sort: true)
    | set { fastqs }

    // QC and trimming module
    qc(fastqs, out_folder_name) | set { trimmed }

    short_ref(trimmed, fasta)
}


logWorkflowCompletion(out_folder_name, !params.map_to_mod_fa)