# Directory Structures

## Validation Outputs

Validation outputs are not published to a persistent top-level directory. They are generated within the `validate` process task workspace and consumed directly by downstream channels.

The validation stage emits:

- `validated_params.json` (runtime pipeline switches and file-path metadata)
- validated FASTA/FASTQ/BAM/GFF/TSV files used by subsequent workflows

## `data/outputs` Directory Structure

After successful pipeline execution, the outputs are organized as follows:

```
data/outputs
├── fasta_ref_mod       → Results from reference vs modified FASTA comparison (if run_ref_x_mod is true)
├── illumina            → Short-read (Illumina) mapping results
├── logs/               → Pipeline logs, Nextflow reports, trace data, and process manifest
├── ont                 → Long-read (Oxford Nanopore) mapping results
├── pacbio              → Long-read (PacBio) mapping results
├── tables              → Per-SV csv tables
└── unmapped_stats      → Summary statistics of unmapped reads for each workflow
```

## Output Documentation

A detailed description of each output subfolder is available in the **[Output Documentation](../outputs/index.md)**:

- [Reference vs Modified FASTA Pipeline](../outputs/fasta-ref-mod.md)
- [Short-Read Processing Pipeline (Illumina)](../outputs/illumina.md)
- [Long-Read Processing Pipeline (PacBio & Oxford Nanopore)](../outputs/long-reads.md)
- [Unmapped Reads Statistics](../outputs/unmapped-stats.md)
- [Logs](../outputs/logs.md)

## See Also

- [Running the Pipeline](running-pipeline.md) - How to execute the pipeline
- [Runtime Messages](runtime-messages.md) - Understanding pipeline progress
