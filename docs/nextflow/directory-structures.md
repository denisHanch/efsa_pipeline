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

### Supported Extensions

| Data Type       | Supported Extensions                   |
| --------------- | -------------------------------------- |
| FASTA sequences | `.fa`, `.fna`, `.fasta`                |
| GFF annotations | `.gff`, `.gtf`                         |
| FASTQ reads     | `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz` |

## `data/outputs` Directory Structure

After successful pipeline execution, the outputs are organized as follows:

```
data/outputs
├── fasta_ref_mod       → Results from reference vs modified FASTA comparison
├── illumina            → Short-read (Illumina) mapping results
├── logs                → Pipeline logs and Nextflow reports
├── ont                 → Long-read (Oxford Nanopore) mapping results
├── pacbio              → Long-read (PacBio) mapping results
├── truvari             → Variant comparison results from Truvari
└── unmapped_stats      → Summary statistics of unmapped reads for each workflow
```

A detailed description of each output subfolder is available in the **[Outputs Documentation](../outputs/)** section.

## Quick Links to Output Documentation

- [Reference vs Modified FASTA Output](../outputs/fasta-ref-mod.md)
- [Illumina Output](../outputs/illumina.md)
- [Long-Read Output (PacBio & ONT)](../outputs/long-reads.md)
- [Truvari Comparison Output](../outputs/truvari.md)
- [Unmapped Reads Statistics](../outputs/unmapped-stats.md)
- [Logs](../outputs/logs.md)

## See Also

- [Running the Pipeline](running-pipeline.md) - How to execute the pipeline
- [Runtime Messages](runtime-messages.md) - Understanding pipeline progress
