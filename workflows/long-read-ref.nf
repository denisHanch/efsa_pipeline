#!/usr/bin/env nextflow

include { nanoplot; multiqc } from '../modules/qc.nf'
include { sv_long; mapping_long; mapping_long as mapping_long_plasmid; sv_long as sv_long_plasmid }  from '../modules/subworkflow.nf'
include { logUnmapped; logWorkflowCompletion } from '../modules/logs.nf'
include { calc_unmapped; get_unmapped_reads } from '../modules/mapping.nf'


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
        logUnmapped(pct, params.long_threshold,  "long-ref-${mapping_tag}")
        get_unmapped_reads(indexed_bam, out_folder_name) | set { unmapped_fastq }
        
        // mapping reads to plasmid & variant calling
        if (plasmid_fasta) {
            Channel.from(plasmid_files) | set { mod_plasmid_fasta }

            mapping_long_plasmid(unmapped_fastq, mod_plasmid_fasta, mapping_tag, "long-ref-plasmid") | set { unmapped_bam }
            get_unmapped_reads(unmapped_bam, "long-ref-plasmid") | set { unmapped_fastq }
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

    Channel.fromPath("${params.in_dir}/pacbio/*.fastq.gz")
        .map { file -> 
            def name = file.baseName.replaceFirst('.fastq', '')
            return [name, file]
        }
        .set { pacbio_fastqs }

    Channel.fromPath("${params.in_dir}/ont/*.fastq.gz")
        .map { file -> 
            def name = file.baseName.replaceFirst('.fastq', '')
            return [name, file]
        }
        .set { ont_fastqs }

    def ref_plasmid = file("$params.in_dir").listFiles()?.findAll { it.name =~ /ref_plasmid\.(fa|fna|fasta)$/ } ?: []
    def mod_plasmid = file("$params.in_dir").listFiles()?.findAll { it.name =~ /mod_plasmid\.(fa|fna|fasta)$/ } ?: []

    
    Channel.fromPath("$params.in_dir/*ref.{fa,fna,fasta}", checkIfExists: true) | set { ref_fasta }
    
    if (pacbio_fastqs) {
        long_ref(pacbio_fastqs, ref_fasta, "map-pb", ref_plasmid, out_folder_name)
    }
}

logWorkflowCompletion(out_folder_name, params.map_to_mod_fa)