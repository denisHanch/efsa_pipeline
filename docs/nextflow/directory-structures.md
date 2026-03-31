# Directory Structures

## `data/valid` Directory Structure

Each validation run creates a timestamped subdirectory. Previous runs are preserved.

```
data/valid/
├── validated_params.json          # Fixed path — always here for Nextflow -params-file
│
└── run_YYYYMMDD_HHMMSS/           # Per-run output directory
    ├── SampleName_ref.fasta       # Validated reference genome FASTA
    ├── SampleName_mod.fasta       # Validated modified genome FASTA (if provided)
    ├── SampleName_ref_plasmid.fasta  # Reference plasmid (if identified)
    ├── SampleName_mod_plasmid.fasta  # Modified plasmid (if identified)
    ├── SampleName_contig_0.fasta  # Contig files from fragmented assembly (if applicable)
    ├── SampleName_ref.gff         # Validated GFF annotation (if provided)
    ├── validation.log             # Structured JSON log for this run
    ├── report.txt                 # Human-readable validation statistics
    │
    ├── illumina/
    │   ├── SampleName_R1.fastq.gz
    │   └── SampleName_R2.fastq.gz
    │
    ├── ont/
    │   ├── SampleName.fastq.gz    # FASTQ input (converted/copied)
    │   └── SampleName.bam         # BAM input (copied as-is, if provided)
    │
    └── pacbio/
        ├── SampleName.fastq.gz    # FASTQ input (converted/copied)
        └── SampleName.bam         # BAM input (copied as-is, if provided)
```

### File Descriptions

| File / Folder              | Description                                                                   |
| -------------------------- | ----------------------------------------------------------------------------- |
| `validated_params.json`    | Produced by the validation step; loaded by Nextflow via `-params-file`. Contains all validated file paths and pipeline flags. Always written at `data/valid/` (top level) so the Nextflow command never needs updating between runs. |
| `run_YYYYMMDD_HHMMSS/`     | Timestamped directory for one validation run. All validated files, the log, and the report go here. Previous runs are never overwritten. |
| `*_ref.fasta`              | Validated reference genome FASTA.                                             |
| `*_mod.fasta`              | Validated modified genome FASTA.                                              |
| `*_ref_plasmid.fasta`      | Reference plasmid sequences (if identified).                                  |
| `*_mod_plasmid.fasta`      | Modified plasmid sequences (if identified).                                   |
| `*_contig_N.fasta`         | Individual contig files from a fragmented modified assembly.                  |
| `*_ref.gff`                | Validated GFF/GTF feature annotation (if provided).                          |
| `validation.log`           | Structured JSON log for the run (auto-incremented if re-run in same second).  |
| `report.txt`               | Human-readable statistics for all validated files.                            |
| `illumina/`                | Paired-end Illumina FASTQ reads.                                              |
| `ont/`                     | Nanopore reads — FASTQ and/or BAM depending on input.                         |
| `pacbio/`                  | PacBio reads — FASTQ and/or BAM depending on input.                           |

## `data/outputs` Directory Structure

After successful pipeline execution, the outputs are organized as follows:

```
data/outputs
├── fasta_ref_mod       → Results from reference vs modified FASTA comparison (if run_ref_x_mod is set to true in `data/valid/validated_params.json`)
├── illumina            → Short-read (Illumina) mapping results
├── logs/               → Pipeline logs, Nextflow reports, trace data, and process manifest
├── ont                 → Long-read (Oxford Nanopore) mapping results
├── pacbio              → Long-read (PacBio) mapping results
├── tables              → Per-SV csv tables
├── truvari             → Variant comparison results from Truvari (if run_truvari is set to true in `data/valid/validated_params.json`)
└── unmapped_stats      → Summary statistics of unmapped reads for each workflow
```

## Output Documentation

A detailed description of each output subfolder is available in the **[Output Documentation](../outputs/index.md)**:

- [Reference vs Modified FASTA Pipeline](../outputs/fasta-ref-mod.md)
- [Short-Read Processing Pipeline (Illumina)](../outputs/illumina.md)
- [Long-Read Processing Pipeline (PacBio & Oxford Nanopore)](../outputs/long-reads.md)
- [VCF Comparison Pipeline with Truvari Output](../outputs/truvari.md)
- [Unmapped Reads Statistics](../outputs/unmapped-stats.md)
- [Logs](../outputs/logs.md)

## See Also

- [Running the Pipeline](running-pipeline.md) - How to execute the pipeline
- [Runtime Messages](runtime-messages.md) - Understanding pipeline progress