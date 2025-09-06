#!/usr/bin/env nextflow

include { nanoplot; multiqc } from '../modules/qc.nf'
include { samtool_stats; minimap2; samtools_sort; samtool_index_bam } from '../modules/mapping.nf'
include { sniffles; debreak; cute_sv; survivor } from '../modules/sv_calling.nf'
include { bcftools_stats } from '../modules/variant_calling.nf'

out_folder_name = "long-ref"


workflow long_ref {
    // Processing inputs
    Channel.fromPath("$params.in_dir/*ref.{fa,fna,fasta}") | set { fasta }

    println("Processing files in directory: ${params.in_dir}")

    Channel.fromPath("${params.in_dir}/tmp2/*_subreads.fastq.gz")
        .map { file -> 
            def name = file.baseName.replaceFirst('.fastq', '')
            return [name, file]
        }
        .set { fastqs }


    // fastqs.view()
    
    // QC
     nanoplot(fastqs)

    // mapping
    minimap2(fastqs, fasta) | set { sam }
    samtools_sort(sam, out_folder_name) | samtool_index_bam | set { indexed_bam }

    // variant calling
    cute_sv(fasta, indexed_bam) | set { cute_vcf }
    debreak(fasta, indexed_bam) | set { debreak_vcf }
    sniffles(indexed_bam) | set { sniffles_vcf }

    survivor(cute_vcf, debreak_vcf, sniffles_vcf) | set { merged_vcf }
    bcftools_stats(merged_vcf, out_folder_name) | set { bcftools_out }
    multiqc(out_folder_name, bcftools_out)
}

workflow {
 long_ref()
}