# Directory Structures

## `data/valid` Directory Structure

This directory contains all input data used by the Nextflow pipeline.

```
data/valid/
├── assembled_genome.fasta
├── reference_genome.fasta
├── ref_plasmid.fa             # Reference plasmid sequences (if used)
├── mod_plasmid.fa             # Modified/assembled plasmid sequences (if used)
├── ref_feature.gff            # Genome annotation file GTF/GFF (if used)
│
├── illumina/
│   ├── SampleName_1.fastq.gz
│   ├── SampleName_2.fastq.gz
│
├── ont/
│   └── SampleName.fastq.gz
│
└── pacbio/
    └── SampleName.fastq.gz
```

### File Descriptions

| File / Folder            | Description                                           |
| ------------------------ | ----------------------------------------------------- |
| `reference_genome.fasta` | The primary reference genome sequence.                |
| `assembled_genome.fasta` | Assembled or modified genome for comparison/analysis. |
| `ref_plasmid.fa`         | Reference plasmid sequences.                          |
| `mod_plasmid.fa`         | Modified or assembled plasmid sequences.              |
| `ref_feature.gff`        | GFF feature file for annotations.                     |
| `illumina/`              | Paired-end Illumina short reads.                      |
| `ont/`                   | Oxford Nanopore long reads.                           |
| `pacbio/`                | PacBio long reads.                                    |


## `data/outputs` Directory Structure

After successful pipeline execution, the outputs are organized as follows:

```
data/outputs
├── fasta_ref_mod       → Results from reference vs modified FASTA comparison
├── illumina            → Short-read (Illumina) mapping results
├── logs                → Pipeline logs and Nextflow reports
├── ont                 → Long-read (Oxford Nanopore) mapping results
├── pacbio              → Long-read (PacBio) mapping results
├── tables              → Per-SV csv tables
├── truvari             → Variant comparison results from Truvari (if --run_truvari set to true)
└── unmapped_stats      → Summary statistics of unmapped reads for each workflow
```

## Output Documentation

A detailed description of each output subfolder is available in the **[Output Documentation](../outputs/index.md)**:

- [Reference vs Modified FASTA Pipeline Output](../outputs/fasta-ref-mod.md)
- [Short-Read Processing Pipeline Output (Illumina)](../outputs/illumina.md)
- [Long-Read Processing Pipeline Output (PacBio & Oxford Nanopore)](../outputs/long-reads.md)
- [VCF Comparison Pipeline with Truvari Output](../outputs/truvari.md)
- [Unmapped Reads Statistics](../outputs/unmapped-stats.md)
- [Logs](../outputs/logs.md)

## See Also

- [Running the Pipeline](running-pipeline.md) - How to execute the pipeline
- [Runtime Messages](runtime-messages.md) - Understanding pipeline progress
