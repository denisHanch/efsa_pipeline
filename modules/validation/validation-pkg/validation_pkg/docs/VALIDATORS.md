# Validator Classes

Detailed documentation for `GenomeValidator`, `ReadValidator`, and `FeatureValidator` classes.

## Table of Contents

- [Overview](#overview)
- [GenomeValidator](#genomevalidator)
- [ReadValidator](#readvalidator)
- [FeatureValidator](#featurevalidator)
- [Validation Levels](#validation-levels)
- [Common Patterns](#common-patterns)

---

## Overview

The validation package provides three main validator classes, each specialized for a specific file type:

| Validator | File Types | Input Formats | Output Format |
|-----------|------------|---------------|---------------|
| **GenomeValidator** | Genome, Plasmid | FASTA, GenBank | FASTA |
| **ReadValidator** | Sequencing Reads | FASTQ, BAM | FASTQ (compressed) |
| **FeatureValidator** | Annotations | GFF, GTF, BED | GFF3 |

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
print(f"GC content: {result.gc_content:.2f}%")
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
- ✓ Calculate statistics (GC content, N50, total size)
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

print(f"Reads: {result.num_reads:,}")
print(f"N50: {result.n50:,} bp")
print(f"Mean length: {result.mean_read_length:.1f} bp")
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
4. Calculate statistics (N50, mean length, etc.)
5. Save to compressed FASTQ format

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
- ✓ Check for invalid characters (if enabled)
- ✓ Calculate statistics (N50, total bases, length distribution)
- ✓ Compress output with gzip

#### Trust Mode
- ✓ Parse first 10 reads for validation
- ✓ Validate first 10 reads only
- ✓ Copy original file to output with compression conversion
- ○ Skip full read count
- ○ Skip statistics calculation
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

**Requirements:** `pysam` installed or `samtools` in PATH. BAM must contain sequence and quality scores.

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

## FeatureValidator

Validates and processes feature annotation files in GFF, GTF, and BED formats using **gffread** as the primary processing engine.

### Basic Usage

```python
from validation_pkg import FeatureValidator, ConfigManager

config = ConfigManager.load("config.json")

# Validate feature file
validator = FeatureValidator(config.ref_feature)
result = validator.run()

print(f"Features: {result.num_features}")
print(f"Types: {result.feature_types}")
print(f"Sequences: {result.sequence_ids}")
```

### Constructor

```python
FeatureValidator(feature_config, settings=None)
```

**Parameters:**
- `feature_config` (FeatureConfig): Configuration from ConfigManager
- `settings` (FeatureValidator.Settings, optional): Custom validation settings

**Attributes:**
- `feature_config`: Input feature configuration
- `output_dir`: Output directory path
- `validation_level`: Validation level
- `threads`: Number of threads for compression
- `input_path`: Input file path
- `settings`: Validator settings object
- `features`: List of parsed Feature objects (after validation)

### Architecture

FeatureValidator uses a **simplified 3-step workflow**:

```
1. _parse_input()   → Parse using gffread (strict/trust) or skip (minimal)
2. _edit_features() → Apply edits in Python (strict mode only)
3. _save_output()   → Write GFF3 output
```

**Key Design Principles:**
- **Trust gffread**: Validation handled by gffread's `-v -E` flags
- **Standards compliance**: gffread follows GFF3 specification exactly
- **Separation of concerns**: gffread for parsing, Python for editing

### Methods

#### `run()`

Execute the validation workflow.

**Returns:** `OutputMetadata`

**Workflow by Mode:**

**Strict Mode:**
```
Input → gffread (parse + validate) → Edit (sort, ID replace) → Save GFF3
```

**Trust Mode:**
```
Input → gffread (parse + validate) → Save GFF3 (no edits)
```

**Minimal Mode:**
```
Input → Copy file (no gffread, no parsing)
```

**Example:**
```python
# Strict mode with sorting and ID replacement
settings = FeatureValidator.Settings(
    sort_by_position=True,
    replace_id_with='chr1'
)
validator = FeatureValidator(config.ref_feature, settings)
result = validator.run()

print(f"Output: {result.output_file}")
print(f"Features: {result.num_features}")
```

### Validation Modes

#### Strict Mode
- ✓ Parse all features via **gffread**
- ✓ Validate with gffread's `-v -E` flags (structure, syntax, coordinates)
- ✓ Convert formats (GTF/BED → GFF3)
- ✓ Apply edits (sorting, ID replacement)
- ✓ Trust gffread for standards compliance

**Use when:** Processing new/untrusted files, need editing, require strict GFF3 compliance

#### Trust Mode
- ✓ Parse all features via **gffread**
- ✓ Validate with gffread's `-v -E` flags
- ✓ Convert formats (GTF/BED → GFF3)
- ○ Skip edits (sorting and ID replacement not applied)
- ✓ Trust gffread for standards compliance

**Use when:** Processing trusted files, want fast validation and format conversion without edits

#### Minimal Mode
- ✗ No parsing (gffread not used)
- ✗ No validation
- ✗ No format conversion
- ✓ Direct file copy
- ⚠️ **Requires:** GFF format + matching compression

**Use when:** Pre-validated files, need fastest processing, no changes required

### Special Features

#### Format Conversion (Automatic)

All format conversions handled by **gffread**:

**GTF → GFF3:**
```python
# Input: genes.gtf
# Output: genes.gff3
validator = FeatureValidator(config.ref_feature)
result = validator.run()

# gffread handles:
# - Attribute conversion (gene_id → ID, transcript_id → Parent)
# - Feature hierarchy (gene → transcript → exon)
# - Coordinate validation
```

**BED → GFF3:**
```python
# Input: features.bed (0-based, half-open)
# Output: features.gff3 (1-based, closed)
validator = FeatureValidator(config.ref_feature)
result = validator.run()

# gffread handles:
# - Coordinate system transformation
# - BED columns → GFF3 attributes
# - Feature type inference
```

#### Feature Sorting (Strict and Trust Modes)

Sort features by genomic position in Python:

```python
settings = FeatureValidator.Settings(sort_by_position=True)
validator = FeatureValidator(config.ref_feature, settings)
result = validator.run()

# Features sorted by: seqname → start → end
# Applied AFTER gffread parsing
```

**Note:** Only applied in strict mode. Skipped in trust and minimal modes.

#### Sequence ID Replacement (Strict Mode Only)

Replace sequence IDs in Python:

```python
# Replace all seqnames with 'chr1'
settings = FeatureValidator.Settings(replace_id_with='chr1')
validator = FeatureValidator(config.ref_feature, settings)
result = validator.run()

# All features in column 1 get 'chr1' as seqname
# Applied AFTER gffread parsing
```

**Note:** Only applied in strict mode. Skipped in trust and minimal modes.

### gffread Integration

**What gffread does:**
- Parse and validate GFF/GTF/BED files
- Convert formats to GFF3
- Check syntax and structure (`-v` verbose, `-E` expose errors)
- Adjust parent features to match children coordinates
- Infer missing features (e.g., mRNA from CDS)

**What Python does:**
- Sort features by position (**strict mode only**)
- Replace sequence IDs (**strict mode only**)
- Write compressed output

**Installation:**
```bash
conda install -c bioconda gffread
```

### Settings Reference

```python
class Settings:
    sort_by_position: bool = True        # Sort features — strict mode only
    check_coordinates: bool = True       # Python-level coordinate validation (start≥1, start≤end)
    replace_id_with: str = None          # Replace seqnames — strict mode only
    coding_type: str = None              # Output compression ('gz', 'bz2')
    output_filename_suffix: str = None   # Add suffix to output filename
    output_subdir_name: str = None       # Save to subdirectory
```

**Notes:**
- `check_coordinates`: validates `start >= 1` and `start <= end` after parsing; strict mode checks all features (parallel when threads > 1 and ≥ 1 000 features), trust mode checks first 10 only; issues are logged as warnings and do not stop processing
- Editing settings (`sort_by_position`, `replace_id_with`) apply in **strict mode only**
- Trust and minimal modes skip all editing

---

## Validation Levels

All validators support three validation levels that control thoroughness vs performance.

### Comparison Table

| Aspect | Strict | Trust | Minimal |
|--------|--------|-------|---------|
| **Parse file** | ✓ All | ✓ All | ✗ None |
| **Validate** | ✓ All records | ✓ All records* | ✗ None |
| **Statistics** | ✓ Full | ○ Limited | ✗ None |
| **Edits** | ✓ Apply | ✓ Apply | ✗ None |
| **Speed** | 1x (baseline) | ~1x (FeatureValidator)† | 100x+ faster |
| **Use case** | First run, full validation | Trusted data, need edits | Pre-validated files |

**Notes:**
- \* **FeatureValidator**: Both strict and trust modes use gffread for validation (same thoroughness), but only strict mode applies Python edits (sort, ID replace)
- † **FeatureValidator**: Trust mode is similar speed to strict (both use gffread); strict mode is slightly slower due to Python edits
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
    validation_level='trust',  # Skip extensive validation
    check_invalid_chars=False  # Skip character checking
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
