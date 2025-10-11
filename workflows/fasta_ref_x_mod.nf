#!/usr/bin/env nextflow

include { nucmer; deltaFilter; showCoords; syri } from '../modules/ref_x_mod.nf'

workflow ref_mod {
    main:
        def prefix_name = "ref_x_mod"
        ref_mod_fasta = Channel.of(tuple(prefix_name, file("$params.in_dir/*ref*.{fa,fna,fasta}", deep: true), file("$params.in_dir/*{assembled_genome,mod}.{fa,fna,fasta}", deep: true)))
        ref_mod_fasta.view()
        ref_mod_fasta | nucmer | set { delta }
        deltaFilter(prefix_name, delta) | set { filtered_delta }
        showCoords(prefix_name, filtered_delta) | set { coords }
        syri(ref_mod_fasta, coords, filtered_delta)

    emit:
        log.info "▶ The reference to modified fasta comparison pipeline completed successfully."
}

workflow {
    ref_mod()
}