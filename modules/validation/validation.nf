#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Only use a fallback if the user didn't provide --input_dir
params.input_dir = params.input_dir ?: "$baseDir/data"  // or remove this line entirely

// Safety: expand relative/tilde and normalize
def input_dir = file(params.input_dir).toAbsolutePath()

process validation {
  container 'pipeline-validation:dev'
  errorStrategy 'ignore'

  input:
  path config_json       
  path validation_script 
  path all_files         

  output:
  path "validation.log", emit: log_file
  path "exit_code.txt", emit: exit_code
  stdout emit: log

  """
  echo "Starting validation"
  python "$validation_script" "$config_json" > validation.log 2>&1
  """
}


workflow {
def all_files_ch = Channel
    .fromPath("${input_dir}/*", type: 'file')
    .filter { it.name != 'config.json' }
    .ifEmpty { error "No files (except config.json) found in: ${input_dir}" }
    .collect()

  def config_ch = Channel
      .fromPath("${input_dir}/config.json", type: 'file', followLinks: true)
      .ifEmpty { error "Missing config.json in: ${input_dir}" }

  def script_ch = Channel
      .fromPath('/EFSA_workspace/modules/validation/main.py', type: 'file')
      .ifEmpty { error "Script not found: /EFSA_workspace/modules/validation/main.py" }

    
  validation(config_ch, script_ch, all_files_ch)
}