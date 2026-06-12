# Validator Classes

Detailed documentation for `GenomeValidator` and `ReadValidator` classes.

## Table of Contents

- [Overview](#overview)
- [GenomeValidator](#genomevalidator)
- [ReadValidator](#readvalidator)
- [Validation Levels](#validation-levels)
- [Common Patterns](#common-patterns)

---

## Overview

The validation package provides three main validator classes, each specialized for a specific file type:

| Validator | File Types | Input Formats | Output Format |
|-----------|------------|---------------|---------------|
| **GenomeValidator** | Genome, Plasmid | FASTA, GenBank | FASTA |
| **ReadValidator** | Sequencing Reads | FASTQ, BAM | FASTQ (compressed) |

All validators support:
- Three validation levels (strict, trust, minimal)
- Automatic compression/decompression
- Parallel processing (via threading)
- Detailed logging and error reporting

---

## GenomeValidator

Validates and processes genome/plasmid files in FASTA and GenBank formats.

### Basic Usage

```python
from validation_pkg import GenomeValidator, ConfigManager

config = ConfigManager.load("config.json")

# Use default settings
validator = GenomeValidator(config.ref_genome)
result = validator.run()

print(f"Sequences: {result.num_sequences}")
print(f"Total size: {result.total_genome_size} bp")
```

### Constructor

```python
GenomeValidator(genome_config, settings=None)
```

**Parameters:**
- `genome_config` (GenomeConfig): Configuration from ConfigManager
- `settings` (GenomeValidator.Settings, optional): Custom validation settings

**Attributes:**
- `genome_config`: Input genome configuration
- `output_dir`: Output directory path
- `validation_level`: Validation level ('strict', 'trust', or 'minimal')
- `threads`: Number of threads for compression
- `input_path`: Input file path
- `settings`: Validator settings object
- `sequences`: List of parsed SeqRecord objects (after validation)

### Methods

#### `run()`

Execute the validation workflow.

**Returns:** `GenomeOutputMetadata`

**Workflow:**
1. Parse genome file using BioPython
2. Validate sequences based on validation level
3. Apply filters (length, ID modification)
4. Handle plasmids (split or merge)
5. Save to FASTA format with optional compression

**Example:**
```python
validator = GenomeValidator(config.ref_genome)
result = validator.run()

# Access results
print(f"Output file: {result.output_file}")
print(f"Processing time: {result.elapsed_time:.2f}s")
```

### Validation Checks

#### Strict Mode
- ✓ Parse all sequences
- ✓ Validate all sequence IDs (no empty IDs)
- ✓ Check for empty sequences
- ✓ Calculate total genome size
- ✓ Apply all edits (filtering, ID replacement, plasmid handling)

#### Trust Mode
- ✓ Parse all sequences
- ✓ Validate only first sequence
- ✓ Apply all edits
- ○ Skip statistics calculation
- ○ Skip extensive validation

#### Minimal Mode
- ○ No parsing/validation
- ○ Direct file copy (rename and move)
- ○ Requires FASTA format input

### Special Features

#### Plasmid Handling

GenomeValidator can automatically detect and handle plasmid sequences:

```python
# Split plasmids into separate files
settings = GenomeValidator.Settings(plasmid_split=True)
validator = GenomeValidator(config.ref_genome, settings)
result = validator.run()

print(f"Main genome: {result.output_file}")
print(f"Plasmids: {result.plasmid_filenames}")
```

**Plasmid Detection:**
- Longest sequence = main chromosome
- All other sequences = plasmids
- Or use `main_first=True` to select first sequence as main

**Plasmid Options:**
- `plasmid_split=True`: Save each plasmid to separate file
- `plasmids_to_one=True`: Merge all plasmids into one file
- `is_plasmid=True`: Treat all sequences as plasmids (no main chromosome)

#### Sequence Filtering

Filter sequences by minimum length:

```python
# Remove sequences shorter than 1000 bp
settings = GenomeValidator.Settings(min_sequence_length=1000)
validator = GenomeValidator(config.ref_genome, settings)
result = validator.run()

print(f"Sequences filtered: {result.num_sequences_filtered}")
```

#### Sequence ID Replacement

Replace sequence IDs with custom values:

```python
# Replace IDs with 'chr' (auto-increments for multiple sequences)
settings = GenomeValidator.Settings(replace_id_with='chr')
validator = GenomeValidator(config.ref_genome, settings)
result = validator.run()

# Output sequence IDs: 'chr', 'chr1', 'chr2', ...
# Original IDs stored in description field
```

### Format Conversion

GenBank → FASTA conversion happens automatically:

```python
# Input: genome.gbk (GenBank format)
# Output: genome.fasta (FASTA format)
validator = GenomeValidator(config.ref_genome)
result = validator.run()

# All GenBank annotations converted to FASTA descriptions
# Original GenBank file preserved
```

---

## ReadValidator

Validates and processes sequencing read files in FASTQ and BAM formats.

### Basic Usage

```python
from validation_pkg import ReadValidator, ConfigManager

config = ConfigManager.load("config.json")

# Validate single read file
validator = ReadValidator(config.reads[0])
result = validator.run()

print(f"NGS type: {result.ngs_type}")
print(f"Pairing detected: {result.illumina_pairing_detected}")
```

### Constructor

```python
ReadValidator(read_config, settings=None)
```

**Parameters:**
- `read_config` (ReadConfig): Configuration from ConfigManager
- `settings` (ReadValidator.Settings, optional): Custom validation settings

**Attributes:**
- `read_config`: Input read configuration
- `output_dir`: Output directory path
- `validation_level`: Validation level
- `threads`: Number of threads for compression
- `input_path`: Input file path
- `ngs_type`: Sequencing platform (illumina, ont, pacbio)
- `settings`: Validator settings object

### Methods

#### `run()`

Execute the validation workflow.

**Returns:** `ReadOutputMetadata`

**Workflow:**
1. Parse read file (FASTQ or BAM)
2. Validate reads based on validation level
3. Detect paired-end patterns (Illumina)
4. Save to compressed FASTQ format

**Example:**
```python
validator = ReadValidator(config.reads[0])
result = validator.run()

# Check if paired-end
if result.read_number:
    print(f"Paired-end detected: R{result.read_number}")
    print(f"Base name: {result.base_name}")
```

### Validation Checks

#### Strict Mode
- ✓ Parse all reads
- ✓ Validate all read IDs and sequences
- ✓ Compress output with gzip

#### Trust Mode
- ✓ Parse first 10 reads for validation
- ✓ Validate first 10 reads only
- ✓ Copy original file to output with compression conversion
- ○ Skip extensive validation

#### Minimal Mode
- ○ No parsing/validation
- ✓ Copy and compress file
- ○ Requires FASTQ format input

### Special Features

#### Paired-End Detection

Automatically detects Illumina paired-end patterns:

```python
validator = ReadValidator(config.reads[0])
result = validator.run()

if result.read_number:
    print(f"Pattern detected!")
    print(f"  Base name: {result.base_name}")
    print(f"  Read number: R{result.read_number}")

# Supported patterns:
# - sample_R1_001.fastq, sample_R2_001.fastq
# - sample_1.fastq, sample_2.fastq
# - sample_R1.fastq, sample_R2.fastq
# - sample.R1.fastq, sample.R2.fastq
# - sampleR1.fastq, sampleR2.fastq
```

Use inter-file validation to check R1↔R2 completeness (see [INTERFILE_VALIDATION.md](INTERFILE_VALIDATION.md)).

#### BAM Handling

BAM files are **copied to the output directory by default** (`ignore_bam=True`, `keep_bam=True`). Downstream tools must handle BAM directly. To enable BAM→FASTQ conversion, set `ignore_bam=False`:

```python
# Default: BAM is copied as-is, no FASTQ conversion
validator = ReadValidator(config.reads[0])
result = validator.run()  # Output: reads.bam (copy)

# Enable FASTQ conversion
settings = ReadValidator.Settings(ignore_bam=False, keep_bam=False)
validator = ReadValidator(config.reads[0], settings)
result = validator.run()  # Output: reads.fastq.gz

# Convert AND keep original BAM
settings = ReadValidator.Settings(ignore_bam=False, keep_bam=True)
```

**Requirements:** `samtools` in PATH. BAM must contain sequence and quality scores.

**Limitations:** Secondary/supplementary alignments are skipped; paired-end info may be lost.

### Output Organization

By default all reads are written to the base output directory. To organize reads into subdirectories by NGS type, enable the `outdir_by_ngs_type` setting:

```python
settings = ReadValidator.Settings(outdir_by_ngs_type=True)

# Output layout:
# output_dir/
#   illumina/
#     sample_R1.fastq.gz
#   ont/
#     nanopore_reads.fastq.gz
#   pacbio/
#     pacbio_reads.fastq.gz
```

---

## Validation Levels

All validators support three validation levels that control thoroughness vs performance.

### Comparison Table

| Aspect | Strict | Trust | Minimal |
|--------|--------|-------|---------|
| **Parse file** | ✓ All | ✓ All | ✗ None |
| **Validate** | ✓ All records | ✓ All records | ✗ None |
| **Statistics** | ✓ Full | ○ Limited | ✗ None |
| **Edits** | ✓ Apply | ✓ Apply | ✗ None |
| **Speed** | 1x (baseline) | 10-15x faster | 100x+ faster |
| **Use case** | First run, full validation | Trusted data, need edits | Pre-validated files |

**Notes:**
- **GenomeValidator**: Trust mode validates first sequence only (faster)
- **ReadValidator**: Trust mode parses only first 10 reads, then copies original file (faster)

### Setting Validation Level

#### Via Config File (Global)

```json
{
  "options": {
    "validation_level": "trust"
  }
}
```

#### Via Config File (Per-File)

```json
{
  "ref_genome_filename": {
    "filename": "genome.fasta",
    "validation_level": "strict"
  }
}
```

#### Via Settings Object

```python
settings = GenomeValidator.Settings(validation_level='trust')
validator = GenomeValidator(config.ref_genome, settings)
```

**Note:** Config file settings override Settings object settings.

---

## Common Patterns

### Pattern 1: Batch Processing with Same Settings

```python
from validation_pkg import ConfigManager, GenomeValidator

config = ConfigManager.load("config.json")

# Create settings once
settings = GenomeValidator.Settings(
    plasmid_split=True,
    min_sequence_length=500
)

# Apply to multiple genomes
for genome_config in [config.ref_genome, config.mod_genome]:
    validator = GenomeValidator(genome_config, settings)
    result = validator.run()
    print(f"{result.output_filename}: {result.num_sequences} sequences")
```

### Pattern 2: Progressive Validation

```python
# First run: strict validation
strict_settings = GenomeValidator.Settings(validation_level='strict')
validator = GenomeValidator(config.ref_genome, strict_settings)
result = validator.run()

# Subsequent runs: trust mode for speed
trust_settings = GenomeValidator.Settings(validation_level='trust')
validator2 = GenomeValidator(config.mod_genome, trust_settings)
result2 = validator2.run()
```

### Pattern 3: Error Recovery

```python
from validation_pkg.exceptions import ValidationError

try:
    validator = GenomeValidator(config.ref_genome)
    result = validator.run()
except ValidationError as e:
    print(f"Validation failed: {e}")

    # Retry with more lenient settings
    relaxed_settings = GenomeValidator.Settings(
        allow_empty_id=True,
        min_sequence_length=0
    )
    validator = GenomeValidator(config.ref_genome, relaxed_settings)
    result = validator.run()
```

### Pattern 4: Performance Optimization

```python
# Optimize for large files
settings = ReadValidator.Settings(
    validation_level='trust'  # Skip extensive validation
)

# Use maximum threads
config.options['threads'] = 16

validator = ReadValidator(config.reads[0], settings)
result = validator.run()
```

---

## See Also

- [SETTINGS.md](SETTINGS.md) - Complete settings reference
- [API_REFERENCE.md](API_REFERENCE.md) - API documentation
- [EXAMPLES.md](EXAMPLES.md) - Code examples
- [ERROR_HANDLING.md](ERROR_HANDLING.md) - Exception handling
