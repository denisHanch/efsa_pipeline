# EFSA Pipeline Documentation

Welcome to the EFSA Pipeline documentation. This pipeline provides comprehensive tools for genomic analysis, including quality control, validation, and variant detection.

## Quick Navigation

### Getting Started

- [Quick Start Guide](getting-started/quick-start.md)
- [Docker Setup](getting-started/docker-setup.md)

### Input Validation

- [Validation Overview](validation/OVERVIEW.md)
- [Supported File Formats](validation/FILE_FORMATS.md)
- [Configuration Guide](validation/CONFIG_GUIDE.md)
- [Validation Settings](validation/SETTINGS.md)
- [Performance Tips](validation/PERFORMANCE.md)

### Nextflow Pipeline

- [Running the Pipeline](nextflow/running-pipeline.md)
- [Runtime Messages](nextflow/runtime-messages.md)
- [Directory Structures](nextflow/directory-structures.md)
- [Pipeline Visualization](nextflow/pipeline-visualization.md)

### Output Documentation

- [Reference vs Modified FASTA Pipeline Output](outputs/fasta-ref-mod.md)
- [Short-Read Pipeline Output (Illumina)](outputs/illumina.md)
- [Long-Read Pipeline Output (PacBio & Oxford Nanopore)](outputs/long-reads.md)
- [VCF Comparison Pipeline Output](outputs/truvari.md)
- [Unmapped Reads Statistics](outputs/unmapped-stats.md)
- [Logs](outputs/logs.md)
- [Generation of per structral variation type CSV tables](outputs/sv-tables.md)

## About

The EFSA Pipeline is designed to process genomic data from multiple sequencing platforms (Illumina, PacBio, and Oxford Nanopore) and perform comprehensive variant analysis including structural variant detection and comparison.
