
include { fastqc; multiqc; trimgalore } from './qc.nf'
include { bwa_mapping; samtool_index_bam; samtools_sort; samtool_stats; picard; calc_unmapped } from './mapping.nf'
include { convert_bcf_to_vcf; delly; samtools_index; picard_dict } from '../modules/sv_calling.nf'
include { snpeff; build_config } from '../modules/variant_calling.nf'


workflow qc {
    take:
        fastqs
        out_folder_name

    main:
        trimgalore(fastqs, out_folder_name) | set { trimmed }
        fastqc(trimmed, out_folder_name) | set {fastqc_out}
        multiqc(fastqc_out, out_folder_name, 'qc')

    emit:
        trimmed
}


workflow mapping {
    take:
        fasta
        fasta_index
        trimmed
        out_folder_name

    main:
        bwa_mapping(fasta, fasta_index, trimmed, out_folder_name) | set { sam } 
        samtools_sort(sam, out_folder_name) | set { bam }
        samtool_stats(bam, out_folder_name) | set { stats_out }
        samtool_index_bam(bam, out_folder_name) | set { indexed_bam }
        picard(fasta, indexed_bam, out_folder_name) | set { picard_out }
        
        picard_out.collect() | set { qc_out }
        multiqc(qc_out, out_folder_name, 'mapping')


    emit:
        indexed_bam
}

workflow sv {
    take:
        fasta
        indexed_bam
        out_folder_name


    main:
        samtools_index(fasta, out_folder_name) | set { fai }
        picard_dict(fasta, out_folder_name) | set { dict }
        delly(indexed_bam, fasta, fai, dict, out_folder_name) | set { bcf }
        convert_bcf_to_vcf(bcf, out_folder_name) | set { sv_vcf }

    emit:
        sv_vcf
}


workflow annotate_vcf {
    take:
        fasta
        gtf
        vcf


    main:
        build_config(fasta, gtf) | set { snpeff_out }
        genome_id = snpeff_out.map { genome_id, snpeff_config -> genome_id }
        snpeff_config = snpeff_out.map { genome_id, snpeff_config -> snpeff_config }
        snpeff(vcf, genome_id, snpeff_config) | set { snpeff_output }
        annotated_vcfs = snpeff_output.map { id, vcf, html -> tuple(id, vcf) }
        qc_vcf = snpeff_output.map { id, vcf, html -> html }

    emit:
        qc_vcf
}