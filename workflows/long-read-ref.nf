#!/usr/bin/env nextflow

include { nanoplot; multiqc } from '../modules/qc.nf'
include { sv_long; mapping_long; mapping_long as mapping_long_plasmid; sv_long as sv_long_plasmid }  from '../modules/subworkflow.nf'
include { logUnmapped } from '../modules/logs.nf'
include { calc_unmapped; get_unmapped_reads } from '../modules/mapping.nf'

out_folder_name = "long-ref"

workflow long_ref {

    take:
        fastqs
        fasta
        mapping_tag

    main:
        // qc
        nanoplot(fastqs, out_folder_name)
        
        // mapping to the reference
        mapping_long(fastqs, fasta, mapping_tag, out_folder_name) | set { indexed_bam }

         // printout % unmapped reads
        calc_unmapped(indexed_bam) | set { pct }
        logUnmapped(pct, params.long_threshold,  "long-ref-${mapping_tag}")
        
        // mapping reads to plasmid & variant calling
        Channel.fromPath("$params.in_dir/*plasmid.{fa,fna,fasta}") | set { mod_plasmid_fasta }

        if (mod_plasmid_fasta) {
            get_unmapped_reads(indexed_bam, "long-ref-plasmid") | set { unmapped_fastq }
            mapping_long_plasmid(unmapped_fastq, mod_plasmid_fasta, mapping_tag, "long-ref-plasmid") | set { unmapped_bam }
            sv_long_plasmid(mod_plasmid_fasta, unmapped_bam,  "long-ref-plasmid")
        }

        // SV calling against the reference
        sv_long(fasta, indexed_bam, out_folder_name)

    emit:
        log.info "▶ The long read processing pipeline completed successfully."
    }


workflow {
    // Processing inputs
    log.info  "Processing files in directory: ${params.in_dir}"

    Channel.fromPath("${params.in_dir}/pacbio/*_subreads.fastq.gz")
        .map { file -> 
            def name = file.baseName.replaceFirst('.fastq', '')
            return [name, file]
        }
        .set { pacbio_fastqs }

    Channel.fromPath("${params.in_dir}/ont/*_subreads.fastq.gz")
        .map { file -> 
            def name = file.baseName.replaceFirst('.fastq', '')
            return [name, file]
        }
        .set { ont_fastqs }


    Channel.fromPath("$params.in_dir/*ref.{fa,fna,fasta}", checkIfExists: true) | set { ref_fasta }
    
     if (pacbio_fastqs) {
        long_ref(pacbio_fastqs, ref_fasta, "map-pb")
    }

    if (ont_fastqs) {        
        long_ref(ont_fastqs, ref_fasta,  "map-ont")
    }
}
