#!/usr/bin/env nextflow

include { nucmer; deltaFilter; showCoords; syri } from '../modules/ref_x_mod.nf'
include { logWorkflowCompletion } from '../modules/logs.nf'


workflow ref_mod {
    take:
        ref_fasta
        mod_fasta
    main:
        def prefix_name = "ref_x_mod"

        ref_mod_fasta = ref_fasta
            .combine(mod_fasta)
            .map { ref, mod -> tuple(prefix_name, ref, mod) }

        ref_mod_fasta | nucmer | set { delta }
        deltaFilter(prefix_name, delta) | set { filtered_delta }
        showCoords(prefix_name, filtered_delta) | set { coords }
        syri(ref_mod_fasta, coords, filtered_delta) | set { sv_vcf }
    emit: 
        sv_vcf
}

workflow {
    Channel.fromPath("$params.in_dir/*ref*.{fa,fna,fasta}", checkIfExists: true) | set { ref_fasta }
    Channel.fromPath("$params.in_dir/*{assembled_genome,mod}.{fa,fna,fasta}", checkIfExists: true) | set { mod_fasta }

    ref_mod(ref_fasta, mod_fasta)
}


logWorkflowCompletion("referece to modified fasta comparision", true)