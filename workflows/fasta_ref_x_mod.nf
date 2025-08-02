#!/usr/bin/env nextflow

include { nucmer; deltaFilter; showCoords; syri } from '../modules/ref_x_mod.nf'


workflow {
    def prefix_name = "ref_x_mod"

    Channel
        .of(tuple(prefix_name, file(params.fastqDir + "/ref.fa"), file(params.fastqDir + "/mod.fa"))) set { ref_mod_fasta }
        ref_mod_fasta | nucmer | set { delta }
        deltaFilter(prefix_name, delta) | set { filtered_delta }
        showCoords(prefix_name, filtered_delta) | set { coords }
        syri(ref_mod_fasta, coords, filtered_delta)
}