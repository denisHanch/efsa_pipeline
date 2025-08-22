# EFSA Pipeline

From within dev container:  

Run pipeline for short-read processing:
`nextflow run workflows/short-read-ref.nf -resume`

Run pipeline for long-read processing:
`nextflow run workflows/long-read-ref.nf -resume`

Run pipeline for ref & mod fasta comparisons:
`nextflow run workflows/fasta_ref_x_mod.nf -resume`

