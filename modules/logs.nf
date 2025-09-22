
def logUnmapped(pct, threshold, out_folder_name) {
    pct.collect{
    def percentage = it.toDouble()
    def log_line = "The percentage of unmapped reads in ${out_folder_name} pipeline: ${percentage} %"

    if (percentage == 0.0)
        log.info log_line
    else if (percentage <= threshold)
        log.warn log_line
    else
        log.error log_line + "\nChange parameter map_to_mod_fa in nextflow.config to true."
    }
}