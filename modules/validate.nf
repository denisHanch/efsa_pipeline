#!/usr/bin/env nextflow

process validate {

    tag "validate"
    publishDir "${params.valid_dir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        filename == 'validated_params.json' ? filename : null
    }

    input:
    path config_json

    output:
    path 'validated_params.json', emit: params_json
    path '{*.fasta,*.fasta.gz,*.gff,*.gff3,*.tsv,*/*.fastq.gz,*/*.bam}', optional: true, emit: validated_files

    script:
    def val_level_arg  = params.validation_level     ? "--validation-level ${params.validation_level}"  : ""
    def log_level_arg  = params.logging_level        ? "--logging-level ${params.logging_level}"         : ""
    def org_type_arg   = params.organism_type        ? "--type ${params.organism_type}"                  : ""
    def defrag_arg     = params.force_defragment_ref ? "--force-defragment-ref"                          : ""
    """
    validation.sh --config ${config_json} --threads ${params.max_cpu} ${val_level_arg} ${log_level_arg} ${org_type_arg} ${defrag_arg}
    """
}
