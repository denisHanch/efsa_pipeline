include { sortVcf; indexVcf; truvari } from './variant_calling.nf'

workflow truvari_comparison {

    // Inputs
    take:
        ref_fasta
        vcfs_ch

    main:
        // Sorting and indexing VCFs
        sortVcf(vcfs_ch) | indexVcf | set { indexed_vcfs }

        // Preprocessing channel for Truvari input
        split_ch = indexed_vcfs.branch {
            ref_mod: it[0] == "ref_x_modsyri"
            others:  it[0] != "ref_x_modsyri"
        }

        ref_mod_ch = split_ch.ref_mod.collect()
        others_ch  = split_ch.others

        vcf_pairs_ch = ref_mod_ch.combine(others_ch)

        // Run Truvari comparison
        truvari(ref_fasta, vcf_pairs_ch)
}

workflow.onComplete {
    if (workflow.success) {
        log.info "✅ Truvari stage finished successfully after comparing pipelines.\n"
    } else {
        log.err " ❌ Truvari failed to finish the comparision."
    }
}