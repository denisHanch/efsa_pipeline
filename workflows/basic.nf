// Import modules
include { fastqc }    from '../modules/fastqc/fastqc.nf'
include { trimming }  from '../modules/trimming/trimming.nf'

workflow BASIC_ANALYSIS {
    take:
    samples_ch
    
    main:
    // 2) FASTQC: access named outputs using .out.emit_name
    fastqc(samples_ch)
    def qc_ch = fastqc.out.report
    def qc_log = fastqc.out.log
    qc_log
        .map { line -> "[FASTQC]: $line" }
        .view()

    // 3) Trimming on QCed reads
    trimming(qc_ch)
    def trim_ch = trimming.out.report
    def trim_log = trimming.out.log
    trim_log
        .map { line -> "[TRIMMING]: $line" }
        .view()
    
    emit:
    fastqc_out = qc_ch
    trimm_report = trim_ch

}