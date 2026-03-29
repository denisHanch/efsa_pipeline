# validation_pkg Technical Documentation

Technical documentation for the `validation_pkg` bioinformatics validation package.

## Overview

The `validation_pkg` is a comprehensive Python package for validating and processing genomic data files. It provides:

- **Multi-format support**: FASTA, GenBank, FASTQ, BAM, GFF, GTF, BED
- **Automatic format detection and conversion**
- **Three validation levels**: strict (thorough), trust (fast), minimal (copy only)
- **Parallel processing**: Multi-threaded compression and validation
- **Comprehensive reporting**: Detailed validation reports and logging
- **Inter-file validation**: Check consistency across multiple files

**Package Location:** `modules/validation/validation-pkg/`

**Main Entry Point:** Used by `modules/validation/main.py`

---

## Documentation Index

### Getting Started

- **[API_REFERENCE.md](API_REFERENCE.md)** - Core API documentation
  - Functional API (`validate_genome`, `validate_reads`, `validate_feature`)
  - Configuration management (`ConfigManager`, `Config`)
  - Output metadata objects
  - Quick start examples

### Validators

- **[VALIDATORS.md](VALIDATORS.md)** - Validator class documentation
  - GenomeValidator - Genome and plasmid files
  - ReadValidator - Sequencing reads (FASTQ, BAM)
  - FeatureValidator - Annotations (GFF, GTF, BED)
  - Validation levels comparison
  - Common patterns

### Configuration

- **[SETTINGS.md](SETTINGS.md)** - Complete settings reference
  - Global options (threads, validation_level, logging_level)
  - GenomeValidator settings (15+ options)
  - ReadValidator settings
  - FeatureValidator settings
  - Settings inheritance and priorities
  - Performance impact analysis

### Advanced Features

- **[INTERFILE_VALIDATION.md](INTERFILE_VALIDATION.md)** - Cross-file validation
  - Genome ↔ Genome validation (sequence structure comparison)
  - Read ↔ Read validation (paired-end completeness)
  - Settings and configuration
  - Complete workflow examples
  - Troubleshooting guide

- **[ADVANCED_FEATURES.md](ADVANCED_FEATURES.md)** - Advanced functionality
  - Plasmid handling (detection, splitting, merging)
  - Format conversion (GenBank→FASTA, BED→GFF, BAM→FASTQ)
  - Sequence ID management
  - Paired-end pattern detection
  - Coordinate system conversion (0-based ↔ 1-based)
  - Parallel compression

### Error Handling & Logging

- **[ERROR_HANDLING.md](ERROR_HANDLING.md)** - Exception handling
  - Exception hierarchy
  - Common exceptions and solutions
  - Error handling patterns
  - Troubleshooting guide

- **[LOGGING_REPORTING.md](LOGGING_REPORTING.md)** - Logging and reports
  - Logging system (`setup_logging`, `get_logger`)
  - Validation reports (`ValidationReport`)
  - Log file structure
  - Report formats (text, JSON)

### Examples

- **[EXAMPLES.md](EXAMPLES.md)** - Practical code examples
  - Quick start examples
  - Common use cases (bacterial genomes, paired-end reads, etc.)
  - Complete workflows
  - Integration examples (Nextflow, Snakemake)

---

## Quick Reference

### Installation

```bash
# Package is part of EFSA pipeline
# For standalone use (production):
pip install -e /path/to/validation-pkg

# For development (includes testing dependencies):
pip install -r /path/to/validation-pkg/requirements-dev.txt
# Or using extras:
pip install -e "/path/to/validation-pkg[dev]"
```

### Basic Usage

```python
from validation_pkg import ConfigManager, validate_genome, validate_reads

# Load configuration
config = ConfigManager.load("config.json")

# Validate files
ref_result = validate_genome(config.ref_genome)
reads_results = validate_reads(config.reads)

print(f"Validated {ref_result.num_sequences} sequences")
print(f"Validated {len(reads_results)} read files")
```

### Common Tasks

| Task | Documentation |
|------|---------------|
| Validate genome files | [API_REFERENCE.md](API_REFERENCE.md#genome-validation) |
| Validate read files | [API_REFERENCE.md](API_REFERENCE.md#read-validation) |
| Validate feature files | [API_REFERENCE.md](API_REFERENCE.md#feature-validation) |
| Compare genomes | [INTERFILE_VALIDATION.md](INTERFILE_VALIDATION.md#genome-genome-validation) |
| Check paired-end reads | [INTERFILE_VALIDATION.md](INTERFILE_VALIDATION.md#read-read-validation) |
| Split plasmids | [ADVANCED_FEATURES.md](ADVANCED_FEATURES.md#plasmid-handling) |
| Convert formats | [ADVANCED_FEATURES.md](ADVANCED_FEATURES.md#format-conversion) |
| Handle errors | [ERROR_HANDLING.md](ERROR_HANDLING.md) |
| Generate reports | [LOGGING_REPORTING.md](LOGGING_REPORTING.md#validation-reports) |
| Custom settings | [SETTINGS.md](SETTINGS.md) |

---

## Package Architecture

### Main Components

```
validation_pkg/
├── __init__.py           # Public API exports
├── config_manager.py     # Configuration loading and parsing
├── exceptions.py         # Exception hierarchy
├── logger.py            # Structured logging system
├── report.py            # Validation report generation
├── validators/
│   ├── genome_validator.py      # Genome file validation
│   ├── read_validator.py        # Read file validation
│   ├── feature_validator.py     # Feature file validation
│   ├── interfile_genome.py      # Genome comparison
│   └── interfile_read.py        # Paired-end checking
└── utils/
    ├── base_settings.py  # BaseSettings, BaseOutputMetadata, BaseValidatorSettings
    ├── base_validator.py # BaseValidator abstract class
    ├── file_handler.py   # File I/O and compression utilities
    ├── formats.py        # CodingType, GenomeFormat, ReadFormat, FeatureFormat enums
    ├── path_utils.py     # Path resolution and directory-traversal security
    └── sequence_stats.py # N50 calculation
```

### Data Flow

```
1. Config File (JSON)
   ↓
2. ConfigManager.load()
   ↓
3. Config Object (GenomeConfig, ReadConfig, FeatureConfig)
   ↓
4. Validator (GenomeValidator, ReadValidator, FeatureValidator)
   ↓
5. OutputMetadata (validation results and statistics)
   ↓
6. Inter-File Validation (optional)
   ↓
7. ValidationReport (final report)
```

---

## Key Concepts

### Validation Levels

All validators support three validation levels:

| Level | Speed | Thoroughness | Use Case |
|-------|-------|--------------|----------|
| **strict** | 1x (baseline) | Full validation + statistics | First run, quality checks |
| **trust** | 10-15x faster | Basic validation only | Trusted data, format conversion |
| **minimal** | 100x+ faster | No validation (copy only) | Pre-validated files |

See: [VALIDATORS.md](VALIDATORS.md#validation-levels)

### Settings Hierarchy

Settings are merged from multiple sources:

1. Default values (lowest priority)
2. Config file global options
3. Config file per-file options
4. Validator Settings object (highest priority)

See: [SETTINGS.md](SETTINGS.md#settings-inheritance)

### Output Organization

Validated files are organized by type:

```
data/valid/
├── reference_genome.fasta      # Genomes
├── modified_genome.fasta
├── ref_feature.gff3           # Features
├── illumina/                  # Reads by NGS type
│   ├── sample_R1.fastq.gz
│   └── sample_R2.fastq.gz
├── ont/
│   └── nanopore.fastq.gz
└── pacbio/
    └── pacbio.fastq.gz
```

---

## Common Workflows

### 1. Basic Validation

```python
from validation_pkg import ConfigManager, validate_genome, validate_reads

config = ConfigManager.load("config.json")
ref = validate_genome(config.ref_genome)
reads = validate_reads(config.reads)
```

See: [EXAMPLES.md](EXAMPLES.md#quick-start)

### 2. Validation with Custom Settings

```python
from validation_pkg import GenomeValidator

settings = GenomeValidator.Settings(
    plasmid_split=True,
    validation_level='trust'
)
result = validate_genome(config.ref_genome, settings)
```

See: [VALIDATORS.md](VALIDATORS.md#basic-usage)

### 3. Complete Pipeline with Reports

```python
from validation_pkg import (
    setup_logging,
    ValidationReport,
    validate_genome,
    validate_reads,
    genomexgenome_validation
)

setup_logging(log_file='./logs/validation.log')
report = ValidationReport('./logs/report.txt')

# Validate
ref = validate_genome(config.ref_genome)
mod = validate_genome(config.mod_genome)
reads = validate_reads(config.reads)

# Report
report.write(ref, 'genome')
report.write(mod, 'genome')
report.write(reads, 'read')

# Inter-file check
check = genomexgenome_validation(ref, mod)
report.write(check, 'genomexgenome')

report.flush(format='text')
```

See: [EXAMPLES.md](EXAMPLES.md#complete-workflows)

---

## Supported Formats

### Input Formats

| File Type | Formats | Extensions |
|-----------|---------|------------|
| **Genome** | FASTA, GenBank | `.fasta`, `.fa`, `.fna`, `.gb`, `.gbk`, `.genbank` |
| **Reads** | FASTQ, BAM | `.fastq`, `.fq`, `.bam` |
| **Features** | GFF, GTF, BED | `.gff`, `.gff3`, `.gtf`, `.bed` |
| **Compression** | gzip, bzip2 | `.gz`, `.bz2` |

### Output Formats

| File Type | Output Format |
|-----------|---------------|
| **Genome** | FASTA (uncompressed) |
| **Reads** | FASTQ (gzip compressed) |
| **Features** | GFF3 (uncompressed) |

---

## Performance Optimization

### Speed Optimization

1. **Use trust mode** for faster validation:
   ```python
   settings = GenomeValidator.Settings(validation_level='trust')
   ```

2. **Install parallel compression tools**:
   ```bash
   sudo apt-get install pigz pbzip2  # 2-6x faster
   ```

3. **Increase threads**:
   ```json
   {"options": {"threads": 16}}
   ```

See: [SETTINGS.md](SETTINGS.md#performance-impact)

### Memory Optimization

- Use streaming parsers (automatic)
- Minimal mode for large files
- Process files individually, not all at once

---

## Troubleshooting

### Common Issues

| Problem | Solution | Documentation |
|---------|----------|---------------|
| "Missing required field" | Fix config.json structure | [ERROR_HANDLING.md](ERROR_HANDLING.md#problem-missing-required-field-ref_genome_filename) |
| "File not found" | Check file paths | [ERROR_HANDLING.md](ERROR_HANDLING.md#problem-file-not-found-pathtofile) |
| "Duplicate sequence IDs" | Use `replace_id_with` setting | [ERROR_HANDLING.md](ERROR_HANDLING.md#problem-duplicate-sequence-ids-found) |
| "Invalid coordinates" | Fix GFF/BED file or use relaxed settings | [ERROR_HANDLING.md](ERROR_HANDLING.md#problem-invalid-coordinates-start-end) |
| "Genome sequence count mismatch" | Adjust inter-file validation settings | [INTERFILE_VALIDATION.md](INTERFILE_VALIDATION.md#problem-genome-sequence-count-mismatch) |
| "Missing R1 for R2" | Add missing file or allow incomplete pairs | [INTERFILE_VALIDATION.md](INTERFILE_VALIDATION.md#problem-found-r2-files-without-matching-r1) |

See: [ERROR_HANDLING.md](ERROR_HANDLING.md#troubleshooting-guide)

---

## Related Documentation

### User Documentation

Located in `docs/validation/`:

- [OVERVIEW.md](../OVERVIEW.md) - User-facing validation overview
- [FILE_FORMATS.md](../FILE_FORMATS.md) - Supported file formats
- [CONFIG_GUIDE.md](../CONFIG_GUIDE.md) - Configuration file guide
- [SETTINGS.md](../SETTINGS.md) - Validation levels (user guide)
- [PERFORMANCE.md](../PERFORMANCE.md) - Performance tips

### Package Documentation

This technical documentation (located in `docs/validation/validation_pkg/`):

- API reference and programmatic usage
- Detailed settings and options
- Advanced features
- Error handling and troubleshooting

---

## Version Information

**Package Version:** 0.1.0
**Author:** Dominika Bohuslavova
**License:** EUPL-1.2

**Python Requirements:**
- Python 3.10+
- BioPython
- structlog
- pysam (optional, for BAM support)

**External Tools (optional):**
- pigz - Parallel gzip (2-4x faster compression)
- pbzip2 - Parallel bzip2 (3-6x faster compression)
- samtools - BAM file handling

---

## Getting Help

1. **Check documentation**: Start with [API_REFERENCE.md](API_REFERENCE.md) and [EXAMPLES.md](EXAMPLES.md)
2. **Review error messages**: See [ERROR_HANDLING.md](ERROR_HANDLING.md)
3. **Enable debug logging**:
   ```python
   setup_logging(console_level='DEBUG', log_file='debug.log')
   ```
4. **Check log files**: `./logs/validation_*.log`

---

## Contributing

For issues or contributions to the validation package, contact the EFSA pipeline development team.

---

**Last Updated:** 2024-12-06
**Documentation Version:** 1.0
