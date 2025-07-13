// Import modules
include { multiqc }   from '../modules/multiqc/multiqc.nf'
include { validation }  from '../modules/validation/validation.nf'

workflow ADVANCED_ANALYSIS {
    take:
    input_ch
    
    main:
    // Validation on trimmed reads
    validation(input_ch, file("${projectDir}/modules/validation/validation.py"))
    def val_ch = validation.out.report
    def val_log = validation.out.log
    val_log
        .map { line -> "[VALIDATION]: $line" }
        .view()

    // MultiQC aggregates validation & QC reports
    multiqc(input_ch)
    def mq_ch = multiqc.out.report
    def mq_log = multiqc.out.log
    mq_log
        .map { line -> "[MULTIQC]: $line" }
        .view()

    
    emit:
    val_out = val_ch
    mq_out = mq_ch
}
