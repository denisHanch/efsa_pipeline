
def logUnmapped(pct, out_folder_name) {
    pct.collect{
    def percentage = it.toDouble()
    def logLine = "The percentage of unmapped reads in ${out_folder_name} pipeline: ${percentage} %"

    if (percentage == 0.0)
        log.info logLine
    else if (percentage <= params.threshold)
        log.warn logLine
    else
        log.error logLine + "\nChange parameter map_to_mod_fa in nextflow.config to true."
    }
}
