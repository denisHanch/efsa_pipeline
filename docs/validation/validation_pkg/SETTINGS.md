# Settings Reference

Complete reference for all validator settings in the validation_pkg package.

## Table of Contents

- [Overview](#overview)
- [Global Options](#global-options)
- [GenomeValidator Settings](#genomevalidator-settings)
- [ReadValidator Settings](#readvalidator-settings)
- [FeatureValidator Settings](#featurevalidator-settings)
- [Settings Inheritance](#settings-inheritance)
- [Performance Impact](#performance-impact)

---

## Overview

The validation package uses a **two-level settings system**:

1. **Global Options**: Apply to all files (set in config.json `options` section)
2. **Validator Settings**: Specific to each validator type (set programmatically or per-file)

### Settings Hierarchy

```
Config File Global Options (lowest priority)
    ↓
Config File Per-File Options
    ↓
Validator Settings Object (highest priority)
```

**Example:**
```python
# Config file has validation_level='trust'
# But Settings object overrides to 'strict'
settings = GenomeValidator.Settings(validation_level='strict')
validator = GenomeValidator(config.ref_genome, settings)
# Result: Uses 'strict' mode
```

---

## Global Options

Global options can be set in the config.json file and apply to ALL validators.

### Allowed Global Options

Only these fields are allowed in the `options` section:

```json
{
  "options": {
    "threads": 8,
    "validation_level": "strict",
    "logging_level": "INFO"
  }
}
```

### `threads`

Number of threads for parallel compression.

**Type:** `int` or `null`
**Default:** Auto-detect (uses CPU count)
**Range:** 1-64 (warning if > system cores)
**Recommended:** 8-16 for strict mode, 4 for trust/minimal

**Example:**
```json
{"threads": 16}
```

**Validation:**
- Must be positive integer or null
- Warning if exceeds system CPU cores
- Warning if exceeds 16 (diminishing returns)

**Performance Impact:**
- **Strict mode**: 3-7x speedup with 8+ threads
- **Trust/Minimal mode**: Minimal benefit (mainly compression)

### `validation_level`

Validation thoroughness level.

**Type:** `str`
**Default:** `"strict"` (if not specified)
**Options:** `"strict"`, `"trust"`, `"minimal"`

| Level | Speed | Use Case |
|-------|-------|----------|
| `strict` | 1x (baseline) | First run, quality checks, statistics |
| `trust` | 10-15x faster | Trusted data, format conversion |
| `minimal` | 100x+ faster | Pre-validated files, rename/move only |

**Example:**
```json
{"validation_level": "trust"}
```

**Behavior by Validator:**
- **Genome**: See [GenomeValidator Settings](#genomevalidator-settings)
- **Read**: See [ReadValidator Settings](#readvalidator-settings)
- **Feature**: See [FeatureValidator Settings](#featurevalidator-settings)

### `logging_level`

Console logging verbosity.

**Type:** `str`
**Default:** `"INFO"`
**Options:** `"DEBUG"`, `"INFO"`, `"WARNING"`, `"ERROR"`

| Level | Output |
|-------|--------|
| `DEBUG` | All messages (very verbose) |
| `INFO` | Progress, statistics, warnings, errors |
| `WARNING` | Warnings and errors only |
| `ERROR` | Errors only |

**Example:**
```json
{"logging_level": "WARNING"}
```

**Note:** File logs always use DEBUG level regardless of this setting.

---

## GenomeValidator Settings

Settings for `GenomeValidator` class. Create using `GenomeValidator.Settings()`.

### Creating Settings

```python
from validation_pkg import GenomeValidator

# Default settings
settings = GenomeValidator.Settings()

# Custom settings
settings = GenomeValidator.Settings(
    plasmid_split=True,
    min_sequence_length=1000,
    replace_id_with='chr'
)

# Update existing settings
settings = settings.update(validation_level='trust')
```

### Validation Options

#### `validation_level`

Validation thoroughness (overrides global option if set).

**Type:** `str`
**Default:** Inherited from global options or `"strict"`
**Options:** `"strict"`, `"trust"`, `"minimal"`

**Strict Mode:**
- Parse all sequences with BioPython
- Validate all sequence IDs (no duplicates, no empty IDs)
- Check for empty sequences
- Calculate statistics (GC content, N50, total size)
- Apply all edits

**Trust Mode:**
- Parse all sequences
- Validate only first sequence
- Apply all edits
- Skip statistics calculation
- ~10-15x faster than strict

**Minimal Mode:**
- No parsing/validation
- Direct file copy
- Only renames and moves file
- Requires FASTA input format
- ~100x+ faster than strict

#### `allow_empty_sequences`

Allow sequences with empty sequence data.

**Type:** `bool`
**Default:** `False`
**Raises:** `GenomeValidationError` if False and empty sequences found

**Example:**
```python
settings = GenomeValidator.Settings(allow_empty_sequences=True)
```

**Use Case:** Some GenBank files may have sequences without sequence data (e.g., reference records).

#### `allow_empty_id`

Allow sequences with empty IDs.

**Type:** `bool`
**Default:** `False`
**Raises:** `GenomeValidationError` if False and empty IDs found

**Example:**
```python
settings = GenomeValidator.Settings(allow_empty_id=True)
```

#### `warn_n_sequences`

Warn if number of sequences exceeds this threshold.

**Type:** `int`
**Default:** `2`
**Raises:** Warning (not error)

**Example:**
```python
# Warn if more than 5 sequences
settings = GenomeValidator.Settings(warn_n_sequences=5)
```

**Use Case:** Bacterial genomes typically have 1-2 sequences (chromosome + plasmid). More sequences may indicate assembly issues.

### Plasmid Handling

#### `is_plasmid`

Treat all sequences as plasmids (no main chromosome).

**Type:** `bool`
**Default:** `False`

**Example:**
```python
# For plasmid-only files
settings = GenomeValidator.Settings(is_plasmid=True)
```

**Behavior:**
- If `plasmid_split=True`: Each sequence saved to separate file
- If `plasmids_to_one=True`: All sequences merged into one file
- If both False: All sequences saved to main output file

#### `plasmid_split`

Separate plasmid sequences into different files.

**Type:** `bool`
**Default:** `False`
**Mutually Exclusive:** Cannot be True if `plasmids_to_one=True`

**Example:**
```python
settings = GenomeValidator.Settings(plasmid_split=True)
validator = GenomeValidator(config.ref_genome, settings)
result = validator.run()

print(f"Main: {result.output_file}")
print(f"Plasmids: {result.plasmid_filenames}")
# Output:
# Main: reference_genome.fasta
# Plasmids: ['ref_plasmid_1.fasta', 'ref_plasmid_2.fasta']
```

**File Naming:**
- Main chromosome: `{basename}_genome.fasta`
- Plasmids: `{basename}_plasmid_1.fasta`, `{basename}_plasmid_2.fasta`, ...

#### `plasmids_to_one`

Merge all plasmid sequences into a single file.

**Type:** `bool`
**Default:** `False`
**Mutually Exclusive:** Cannot be True if `plasmid_split=True`

**Example:**
```python
settings = GenomeValidator.Settings(plasmids_to_one=True)
validator = GenomeValidator(config.ref_genome, settings)
result = validator.run()

print(f"Main: {result.output_file}")
print(f"Merged plasmids: {result.plasmid_filenames[0]}")
# Output:
# Main: reference_genome.fasta (chromosome only)
# Merged plasmids: ref_plasmid.fasta (all plasmids)
```

#### `main_longest`

Select longest sequence as main chromosome.

**Type:** `bool`
**Default:** `True`
**Mutually Exclusive:** Cannot be True if `main_first=True`

**Example:**
```python
settings = GenomeValidator.Settings(main_longest=True)
```

**Use Case:** Default behavior. Assumes chromosome is longer than plasmids.

#### `main_first`

Select first sequence as main chromosome.

**Type:** `bool`
**Default:** `False`
**Mutually Exclusive:** Cannot be True if `main_longest=True`

**Example:**
```python
settings = GenomeValidator.Settings(main_first=True, main_longest=False)
```

**Use Case:** When file structure is known and chromosome is always first.

### Editing Options

#### `replace_id_with`

Replace all sequence IDs with this value.

**Type:** `str` or `None`
**Default:** `None` (no replacement)

**Behavior:**
- Single sequence: Uses exact value
- Multiple sequences: Auto-appends increments
- Original IDs stored in description field

**Example:**
```python
settings = GenomeValidator.Settings(replace_id_with='chr')
validator = GenomeValidator(config.ref_genome, settings)
result = validator.run()

# Input IDs: NC_000913.3, plasmid_pA
# Output IDs: chr, chr1
# Descriptions: chr description="Original_id:NC_000913.3"
```

#### `min_sequence_length`

Minimum sequence length to keep (in base pairs).

**Type:** `int`
**Default:** `100`
**Range:** 0-∞

**Example:**
```python
# Remove sequences shorter than 1000 bp
settings = GenomeValidator.Settings(min_sequence_length=1000)
validator = GenomeValidator(config.ref_genome, settings)
result = validator.run()

print(f"Filtered: {result.num_sequences_filtered} sequences")
```

**Use Case:** Remove small contigs, filter assembly artifacts.

### Output Format

#### `coding_type`

Output compression format.

**Type:** `str` or `None`
**Default:** `None` (uncompressed)
**Options:** `"gz"`, `"bz2"`, `None`

**Example:**
```python
# Output compressed with gzip
settings = GenomeValidator.Settings(coding_type='gz')
```

**Performance:**
- Uses `pigz` if available (2-4x faster than gzip)
- Uses `pbzip2` if available (3-6x faster than bzip2)
- Falls back to standard gzip/bzip2

#### `output_filename_suffix`

Add suffix to output filename.

**Type:** `str` or `None`
**Default:** `None`

**Example:**
```python
settings = GenomeValidator.Settings(output_filename_suffix='_filtered')
# Input: genome.fasta
# Output: genome_filtered.fasta
```

#### `output_subdir_name`

Create subdirectory for output files.

**Type:** `str` or `None`
**Default:** `None`

**Example:**
```python
settings = GenomeValidator.Settings(output_subdir_name='genomes')
# Output: {output_dir}/genomes/genome.fasta
```

---

## ReadValidator Settings

Settings for `ReadValidator` class.

### Creating Settings

```python
from validation_pkg import ReadValidator

# Default settings
settings = ReadValidator.Settings()

# Custom settings
settings = ReadValidator.Settings(
    validation_level='trust',
    check_invalid_chars=False
)
```

### Validation Options

#### `validation_level`

Validation thoroughness.

**Type:** `str`
**Default:** Inherited from global options or `"strict"`
**Options:** `"strict"`, `"trust"`, `"minimal"`

**Strict Mode:**
- Parse all reads with BioPython
- Validate all read IDs and sequences
- Check for invalid characters (if enabled)
- Calculate statistics (N50, total bases, length distribution)
- Output: gzip compressed FASTQ

**Trust Mode:**
- Parse all reads (fast line counting)
- Validate only first 10 reads
- Count total reads
- Output: gzip compressed FASTQ
- ~10-15x faster than strict

**Minimal Mode:**
- No parsing/validation
- Copy and compress file
- Requires FASTQ input
- ~100x+ faster than strict

#### `check_invalid_chars`

Check for invalid characters in read sequences.

**Type:** `bool`
**Default:** `False`

**Valid Characters:** A, T, C, G, N (case-insensitive)

**Example:**
```python
# Enable character checking
settings = ReadValidator.Settings(check_invalid_chars=True)
```

**Performance Impact:** Moderate (5-10% slower when enabled in strict mode)

#### `allow_empty_id`

Allow reads with empty IDs.

**Type:** `bool`
**Default:** `False`
**Raises:** `ReadValidationError` if False and empty IDs found

**Example:**
```python
settings = ReadValidator.Settings(allow_empty_id=True)
```

### Output Format

#### `coding_type`

Output compression format.

**Type:** `str` or `None`
**Default:** `None` (uncompressed)
**Options:** `"gz"`, `"bz2"`, `None`

**Example:**
```python
settings = ReadValidator.Settings(coding_type='bz2')
```

#### `output_filename_suffix`

Add suffix to output filename.

**Type:** `str` or `None`
**Default:** `None`

**Example:**
```python
settings = ReadValidator.Settings(output_filename_suffix='_qc')
# Input: reads.fastq
# Output: reads_qc.fastq.gz
```

#### `output_subdir_name`

Save output to a named subdirectory.

**Type:** `str` or `None`
**Default:** `None` (writes directly to output_dir)

**Example:**
```python
settings = ReadValidator.Settings(output_subdir_name='illumina_reads')
# Output: {output_dir}/illumina_reads/reads.fastq.gz
```

**Note:** To automatically organize reads by NGS type into subdirectories, use `outdir_by_ngs_type=True` instead.

#### `outdir_by_ngs_type`

Automatically create a subdirectory named after the NGS type.

**Type:** `bool`
**Default:** `False`

**Example:**
```python
settings = ReadValidator.Settings(outdir_by_ngs_type=True)
# Output: {output_dir}/illumina/reads.fastq.gz
#         {output_dir}/ont/reads.fastq.gz
```

---

## FeatureValidator Settings

Settings for `FeatureValidator` class.

### Creating Settings

```python
from validation_pkg import FeatureValidator

# Default settings
settings = FeatureValidator.Settings()

# Custom settings
settings = FeatureValidator.Settings(
    sort_by_position=True,
    replace_id_with='chr1'
)
```

### Validation Options

#### `validation_level`

Validation thoroughness.

**Type:** `str`
**Default:** Inherited from global options or `"strict"`
**Options:** `"strict"`, `"trust"`, `"minimal"`

**Strict Mode:**
- Parse all features via gffread (validates syntax, coordinates, structure)
- Apply Python edits (sorting, ID replacement)
- Convert to GFF3 format

**Trust Mode:**
- Parse all features via gffread (same validation as strict)
- Skip Python edits (no sorting or ID replacement)
- Similar speed to strict (gffread runs in both modes)

**Minimal Mode:**
- No parsing/validation (gffread not used)
- Direct file copy
- Requires GFF format input

#### `sort_by_position`

Sort features by genomic position.

**Type:** `bool`
**Default:** `True`

**Sort Order:** seqname → start → end

**Example:**
```python
settings = FeatureValidator.Settings(sort_by_position=True)
```

**Use Case:** Required by many downstream tools (e.g., tabix, bedtools).

#### `check_coordinates`

Legacy setting preserved for backward compatibility. Coordinate validation is now handled entirely by gffread (`-v -E` flags) and this field is not used internally.

**Type:** `bool`
**Default:** `True`

### Editing Options

#### `replace_id_with`

Replace sequence IDs (column 1) for all features.

**Type:** `str` or `None`
**Default:** `None` (no replacement)

**Behavior:**
- All features get the same sequence ID in column 1
- Only applied in strict mode
- Original sequence IDs are not preserved

**Example:**
```python
settings = FeatureValidator.Settings(replace_id_with='chr1')

# Input GFF:
# NC_000913.3  source  gene  100  200  .  +  .  ID=gene1
#
# Output GFF:
# chr1  source  gene  100  200  .  +  .  ID=gene1
```

### Output Format

#### `coding_type`

Output compression format.

**Type:** `str` or `None`
**Default:** `None` (uncompressed)
**Options:** `"gz"`, `"bz2"`, `None`

**Example:**
```python
settings = FeatureValidator.Settings(coding_type='gz')
```

#### `output_filename_suffix`

Add suffix to output filename.

**Type:** `str` or `None`
**Default:** `None`

**Example:**
```python
settings = FeatureValidator.Settings(output_filename_suffix='_sorted')
# Input: features.gff
# Output: features_sorted.gff3
```

#### `output_subdir_name`

Create subdirectory for output files.

**Type:** `str` or `None`
**Default:** `None`

**Example:**
```python
settings = FeatureValidator.Settings(output_subdir_name='annotations')
# Output: {output_dir}/annotations/features.gff3
```

---

## Settings Inheritance

Understanding how settings are merged from different sources.

### Priority Order (Highest to Lowest)

1. **Validator Settings Object** (programmatic)
2. **Config File Per-File Options**
3. **Config File Global Options**
4. **Default Values**

### Example: Complete Override Chain

```json
{
  "options": {
    "validation_level": "minimal",
    "threads": 8
  },
  "ref_genome_filename": {
    "filename": "genome.fasta",
    "validation_level": "trust"
  }
}
```

```python
# Config: global=minimal, file=trust
settings = GenomeValidator.Settings(validation_level='strict')
validator = GenomeValidator(config.ref_genome, settings)

# Final validation_level: 'strict' (Settings object wins)
# Final threads: 8 (from config, not overridden)
```

### Warning Messages

When file-level options override global options:

```
WARNING: File-level option 'validation_level=strict' overrides
         global option 'validation_level=trust'
```

---

## Performance Impact

### Validation Level Performance

| Validator | Strict | Trust | Minimal |
|-----------|--------|-------|---------|
| **Genome** | 1x | ~10x faster (first seq only) | 100x+ faster |
| **Read** | 1x | ~10x faster (first 10 reads) | 100x+ faster |
| **Feature** | 1x | Similar (gffread, no Python edits) | 100x+ faster |

### Threading Performance

#### Genome Validation (Strict Mode)

| Threads | Speed | Use Case |
|---------|-------|----------|
| 1 | 1x | Single-core systems |
| 4 | 3x | Small genomes |
| 8 | 5x | **Recommended** |
| 16 | 7x | Large genomes |
| 32+ | 7-8x | Diminishing returns |

#### Read Validation (All Modes)

Threading mainly benefits compression:
- Strict mode: 2-3x speedup with 8 threads
- Trust/Minimal mode: 1-2x speedup (compression only)

### Option-Specific Performance

| Option | Performance Impact | Notes |
|--------|-------------------|-------|
| `check_invalid_chars=False` | 5-10% faster | Reads only |
| `sort_by_position=False` | 10-20% faster | Features only |
| `plasmid_split=True` | 5-10% slower | Genomes only (extra I/O) |
| `min_sequence_length>0` | Negligible | Filtering is fast |
| `replace_id_with` | Negligible | String replacement is fast |

---

## See Also

- [VALIDATORS.md](VALIDATORS.md) - Validator class documentation
- [API_REFERENCE.md](API_REFERENCE.md) - API reference
- [ADVANCED_FEATURES.md](ADVANCED_FEATURES.md) - Advanced feature guides
- [PERFORMANCE.md](../PERFORMANCE.md) - Performance optimization guide
