#!/usr/bin/env nextflow

process validation {
    tag { sample.baseName }

    // Use the local "pipeline-validation:dev"
    container 'pipeline-validation:dev'

    input:
      path sample
      path validation_script

    output:
      path "*.${task.process}.txt", emit: report
      stdout emit: log

    script:
    """
    echo "Starting validation"
    python $validation_script $sample ${sample.baseName}.${task.process}.txt
    """
}