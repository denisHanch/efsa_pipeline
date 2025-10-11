#!/usr/bin/env nextflow

include { nucmer; deltaFilter; showCoords; syri } from '../modules/ref_x_mod.nf'
include { logWorkflowCompletion } from '../modules/logs.nf'


workflow ref_mod {
    main:
        def prefix_name = "ref_x_mod"

        ref_mod_fasta = Channel.of(tuple(prefix_name, file("$params.in_dir/tmp/*ref.{fa,fna,fasta}"), file("$params.in_dir/tmp/*mod.{fa,fna,fasta}")))
        ref_mod_fasta | nucmer | set { delta }
        deltaFilter(prefix_name, delta) | set { filtered_delta }
        showCoords(prefix_name, filtered_delta) | set { coords }
        syri(ref_mod_fasta, coords, filtered_delta)
}

workflow {
    ref_mod()
}

logWorkflowCompletion("referece to modified fasta comparision", true)