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
            assembly: it[0] == "assembly"
            others:  it[0] != "assembly"
        }

        assembly_vcfs = split_ch.assembly.collect()
        vcf_pairs_ch  = assembly_vcfs.combine(split_ch.others)

        // Run Truvari comparison
        truvari(ref_fasta, vcf_pairs_ch)
}

workflow.onComplete {
    if (workflow.success) {
        log.info "✅ Truvari: the comparison of vcf files finished successfully.\n"
    } else {
        log.err "❌ Truvari: the comparison of vcf files failed.\n"
    }
}