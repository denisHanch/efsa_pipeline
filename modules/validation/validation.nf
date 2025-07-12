process validation {
  tag { sample.baseName }

  // Use the local "pipeline-validation:dev" image you built
  container 'pipeline-validation:dev'

  input:
    path sample

  output:
    path "${sample.baseName}.validationd.txt"

  script:
  """
  python modules/validation/validation.py $sample ${sample.baseName}.validation.txt
  """
}
