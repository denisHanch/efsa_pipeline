# Inter-File Validation

Documentation for cross-file consistency validation in the validation_pkg package.

## Table of Contents

- [Overview](#overview)
- [Genome ↔ Genome Validation](#genome-genome-validation)
- [Read ↔ Read Validation](#read-read-validation)
- [Complete Workflow](#complete-workflow)
- [Best Practices](#best-practices)
- [Troubleshooting](#troubleshooting)

---

## Overview

Inter-file validation checks consistency **across multiple validated files**. This is essential for ensuring that reference and modified genomes are compatible, or that paired-end reads are complete.

### Validation Types

| Validation | Purpose | Use Case |
|------------|---------|----------|
| **Genome ↔ Genome** | Compare ref vs mod genomes | Ensure consistent sequence structure |
| **Read ↔ Read** | Check paired-end completeness | Verify R1↔R2 matching |

### When to Use

Inter-file validation should be run **after** individual file validation:

```python
# Step 1: Validate individual files
ref_result = validate_genome(config.ref_genome)
mod_result = validate_genome(config.mod_genome)

# Step 2: Validate inter-file consistency
from validation_pkg import genomexgenome_validation
result = genomexgenome_validation(ref_result, mod_result)
```

---

## Genome ↔ Genome Validation

Compare reference and modified genomes for structural consistency.

### Function: `genomexgenome_validation()`

```python
from validation_pkg import genomexgenome_validation, GenomeXGenomeSettings

result = genomexgenome_validation(
    ref_genome_result,
    mod_genome_result,
    settings=GenomeXGenomeSettings()
)
```

**Parameters:**
- `ref_genome_result` (GenomeOutputMetadata): Result from reference genome validation
- `mod_genome_result` (GenomeOutputMetadata): Result from modified genome validation
- `settings` (GenomeXGenomeSettings, optional): Validation settings

**Returns:** `dict` with keys:
- `passed` (bool): True if all checks passed
- `warnings` (List[str]): List of warning messages
- `errors` (List[str]): List of error messages
- `metadata` (dict): Validation metadata

**Raises:**
- `InterFileValidationError`: If critical validation fails (non-empty `errors` list in result)

---

### GenomeXGenomeSettings

Settings for genome-to-genome validation.

```python
from validation_pkg import GenomeXGenomeSettings

settings = GenomeXGenomeSettings(
    same_number_of_sequences=True,
    same_sequence_ids=False,
    same_sequence_lengths=False
)
```

#### Options

##### `same_number_of_sequences`

Require both genomes to have the same number of sequences.

**Type:** `bool`
**Default:** `True`

**Example:**
```python
settings = GenomeXGenomeSettings(same_number_of_sequences=True)
result = genomexgenome_validation(ref_result, mod_result, settings)

# Fails if:
# - ref_genome has 2 sequences (chromosome + plasmid)
# - mod_genome has 1 sequence (chromosome only)
```

**Use Case:** Ensure structural similarity between genomes.

##### `same_sequence_ids`

Require both genomes to have matching sequence IDs (order-independent).

**Type:** `bool`
**Default:** `False`

**Example:**
```python
settings = GenomeXGenomeSettings(same_sequence_ids=True)
result = genomexgenome_validation(ref_result, mod_result, settings)

# Fails if:
# - ref_genome has IDs: ['chr1', 'plasmid_A']
# - mod_genome has IDs: ['chr1', 'plasmid_B']
```

**Use Case:** Ensure genomes reference the same sequences (important for feature comparison).

##### `same_sequence_lengths`

Require matching sequence lengths for common sequence IDs.

**Type:** `bool`
**Default:** `False`
**Requires:** `same_sequence_ids=True`

**Example:**
```python
settings = GenomeXGenomeSettings(
    same_sequence_ids=True,
    same_sequence_lengths=True
)
result = genomexgenome_validation(ref_result, mod_result, settings)

# Fails if:
# - Both have 'chr1' sequence
# - ref chr1: 4,000,000 bp
# - mod chr1: 4,001,000 bp (different length)
```

**Use Case:** Verify genomes have identical sequence lengths (no large insertions/deletions).

---

### Basic Example

```python
from validation_pkg import (
    ConfigManager,
    validate_genome,
    genomexgenome_validation,
    GenomeXGenomeSettings
)

# Load configuration
config = ConfigManager.load("config.json")

# Validate genomes
ref_result = validate_genome(config.ref_genome)
mod_result = validate_genome(config.mod_genome)

# Inter-file validation with default settings
result = genomexgenome_validation(ref_result, mod_result)

if result['passed']:
    print("✓ Genomes are compatible")
else:
    print("✗ Genome validation failed:")
    for error in result['errors']:
        print(f"  - {error}")
```

---

### Advanced Example

```python
from validation_pkg import GenomeXGenomeSettings

# Strict validation: require everything to match
strict_settings = GenomeXGenomeSettings(
    same_number_of_sequences=True,
    same_sequence_ids=True,
    same_sequence_lengths=True
)

result = genomexgenome_validation(
    ref_result,
    mod_result,
    settings=strict_settings
)

# Check results
print(f"Status: {'PASSED' if result['passed'] else 'FAILED'}")

# Metadata details
metadata = result['metadata']
print(f"Ref sequences: {metadata['ref_num_sequences']}")
print(f"Mod sequences: {metadata['mod_num_sequences']}")
print(f"Common IDs: {metadata['common_sequence_ids']}")
print(f"Ref-only IDs: {metadata['ref_only_ids']}")
print(f"Mod-only IDs: {metadata['mod_only_ids']}")
print(f"Length mismatches: {metadata['length_mismatches']}")
```

---

### Result Metadata

The `metadata` dict contains detailed validation information:

```python
{
    'ref_num_sequences': 2,
    'mod_num_sequences': 2,
    'common_sequence_ids': ['chr1', 'plasmid_A'],
    'ref_only_ids': [],
    'mod_only_ids': [],
    'length_mismatches': {
        # sequence_id: {'ref_length': int, 'mod_length': int, 'difference': int}
    },
    # Populated when settings.characterize=True (requires minimap2):
    'contigs_found': True,           # whether mapped sequences were identified
    'plasmids_found': True,          # whether unmapped sequences were identified
    'contig_files': ['/path/to/contig_0.fasta'],  # individual contig FASTA files
    'plasmid_file': '/path/to/plasmid.fasta'      # merged plasmid FASTA (None if none found)
}
```

---

## Read ↔ Read Validation

Check paired-end read completeness (R1↔R2 matching based on filename regex).

### Function: `readxread_validation()`

```python
from validation_pkg import readxread_validation, ReadXReadSettings

result = readxread_validation(
    reads_results,
    settings=ReadXReadSettings()
)
```

**Parameters:**
- `reads_results` (List[ReadOutputMetadata]): Results from read validation
- `settings` (ReadXReadSettings, optional): Validation settings

**Returns:** `dict` with keys:
- `passed` (bool): True if all checks passed
- `warnings` (List[str]): List of warning messages
- `errors` (List[str]): List of error messages
- `metadata` (dict): Validation metadata

**Raises:**
- `ReadValidationError`: If critical validation fails

---

### ReadXReadSettings

Settings for read-to-read validation.

```python
from validation_pkg import ReadXReadSettings

settings = ReadXReadSettings(
    pair_end_basename=True,
    allow_missing_r1=False
)
```

#### Options

##### `pair_end_basename`

Check that R2 files have matching R1 files.

**Type:** `bool`
**Default:** `True`

**Example:**
```python
settings = ReadXReadSettings(pair_end_basename=True)
result = readxread_validation(reads_results, settings)

# Checks:
# - sample_R1.fastq.gz exists → OK
# - sample_R2.fastq.gz exists → OK
# - orphan_R2.fastq.gz exists but no orphan_R1 → ERROR
```

**Use Case:** Ensure paired-end datasets are complete before analysis.

##### `allow_missing_r1`

Allow R2 files without matching R1 files.

**Type:** `bool`
**Default:** `False`

**Example:**
```python
# Allow incomplete pairs (warning instead of error)
settings = ReadXReadSettings(
    pair_end_basename=True,
    allow_missing_r1=True
)
result = readxread_validation(reads_results, settings)

# Missing R1 for R2 → WARNING (not error)
```

**Use Case:** Handling incomplete datasets, QC filtering that removed some R1 files.

---

### Basic Example

```python
from validation_pkg import (
    ConfigManager,
    validate_reads,
    readxread_validation
)

# Load configuration
config = ConfigManager.load("config.json")

# Validate all read files
reads_results = validate_reads(config.reads)

# Inter-file validation
result = readxread_validation(reads_results)

if result['passed']:
    print("✓ All paired-end reads are complete")
else:
    print("✗ Paired-end validation failed:")
    for error in result['errors']:
        print(f"  - {error}")
```

---

### Advanced Example

```python
from validation_pkg import ReadXReadSettings

# Validate reads
reads_results = validate_reads(config.reads)

# Check paired-end completeness
settings = ReadXReadSettings(
    pair_end_basename=True,
    allow_missing_r1=False
)

result = readxread_validation(reads_results, settings)

# Detailed results
metadata = result['metadata']
print(f"Pairs checked: {metadata['pairs_checked']}")
print(f"Complete pairs: {metadata['complete_pairs']}")
print(f"Missing R1: {metadata['missing_r1']}")
print(f"Duplicate R1: {metadata['duplicate_r1']}")
print(f"Duplicate R2: {metadata['duplicate_r2']}")
```

---

### Supported Paired-End Patterns

ReadValidator automatically detects these Illumina paired-end patterns:

| Pattern | Example R1 | Example R2 |
|---------|------------|------------|
| Lane numbers | `sample_R1_001.fastq` | `sample_R2_001.fastq` |
| Simple numbers | `sample_1.fastq` | `sample_2.fastq` |
| Standard | `sample_R1.fastq` | `sample_R2.fastq` |
| Dot separator | `sample.R1.fastq` | `sample.R2.fastq` |
| No separator | `sampleR1.fastq` | `sampleR2.fastq` |
| With suffix | `sample_R1_combined.fastq` | `sample_R2_combined.fastq` |

**Base Name Extraction:**
- `sample_R1_001.fastq` → base: `sample`, read: `1`
- `sample_R2_001.fastq` → base: `sample`, read: `2`
- Validates that both R1 and R2 exist for the same base name

---

### Result Metadata

The `metadata` dict contains detailed validation information:

```python
{
    'pairs_checked': 3,
    'complete_pairs': ['sample1', 'sample2'],
    'missing_r1': ['orphan_R2.fastq.gz'],
    'duplicate_r1': [],
    'duplicate_r2': []
}
```

---

## Complete Workflow

Complete example with all validation steps.

### Full Pipeline

```python
from validation_pkg import (
    ConfigManager,
    validate_genome,
    validate_reads,
    validate_feature,
    genomexgenome_validation,
    readxread_validation,
    GenomeXGenomeSettings,
    ReadXReadSettings,
    ValidationReport
)

# 1. Load configuration
config = ConfigManager.load("config.json")

# 2. Setup report
report = ValidationReport("./logs/validation_report.txt")

# 3. Validate genomes
ref_genome_result = validate_genome(config.ref_genome)
mod_genome_result = validate_genome(config.mod_genome)
report.write(ref_genome_result, "genome")
report.write(mod_genome_result, "genome")

# 4. Validate reads
reads_results = validate_reads(config.reads)
report.write(reads_results, "read")

# 5. Validate features
if config.ref_feature:
    feature_result = validate_feature(config.ref_feature)
    report.write(feature_result, "feature")

# 6. Inter-file validation: Genome ↔ Genome
genome_settings = GenomeXGenomeSettings(
    same_number_of_sequences=True,
    same_sequence_ids=False
)
genome_check = genomexgenome_validation(
    ref_genome_result,
    mod_genome_result,
    settings=genome_settings
)
report.write(genome_check, "genomexgenome")

if not genome_check['passed']:
    print("WARNING: Genome validation failed:")
    for error in genome_check['errors']:
        print(f"  {error}")

# 7. Inter-file validation: Read ↔ Read
read_settings = ReadXReadSettings(
    pair_end_basename=True,
    allow_missing_r1=False
)
read_check = readxread_validation(reads_results, settings=read_settings)
report.write(read_check, "readxread")

if not read_check['passed']:
    print("WARNING: Read validation failed:")
    for error in read_check['errors']:
        print(f"  {error}")

# 8. Generate final report
report.flush(format="text")

# Summary
print("\n=== Validation Summary ===")
print(f"Genomes: {len([ref_genome_result, mod_genome_result])} validated")
print(f"Reads: {len(reads_results)} validated")
print(f"Genome check: {'PASSED' if genome_check['passed'] else 'FAILED'}")
print(f"Read check: {'PASSED' if read_check['passed'] else 'FAILED'}")
```

---

## Best Practices

### 1. Run Inter-File Validation After Individual Validation

```python
# ✓ GOOD: Validate files first, then check consistency
ref_result = validate_genome(config.ref_genome)
mod_result = validate_genome(config.mod_genome)
check = genomexgenome_validation(ref_result, mod_result)

# ✗ BAD: Cannot run inter-file validation without results
check = genomexgenome_validation(None, None)  # Error!
```

### 2. Use Appropriate Settings for Your Use Case

```python
# Comparing genomes from same species/strain
settings = GenomeXGenomeSettings(
    same_number_of_sequences=True,  # Should have same structure
    same_sequence_ids=False,        # IDs may differ
    same_sequence_lengths=False     # Lengths may differ slightly
)

# Comparing nearly identical genomes (ref vs mutant)
strict_settings = GenomeXGenomeSettings(
    same_number_of_sequences=True,
    same_sequence_ids=True,         # Should have matching IDs
    same_sequence_lengths=True      # Lengths should match exactly
)
```

### 3. Handle Errors Gracefully

```python
result = genomexgenome_validation(ref_result, mod_result)

if not result['passed']:
    # Log errors but continue pipeline if acceptable
    for error in result['errors']:
        logger.warning(f"Genome compatibility issue: {error}")

    # Decide whether to continue or abort
    if result['metadata']['ref_num_sequences'] != result['metadata']['mod_num_sequences']:
        raise ValueError("Cannot proceed with different sequence counts")
```

### 4. Check Metadata for Details

```python
result = readxread_validation(reads_results)

# Don't just check passed/failed
if result['metadata']['missing_r1']:
    print(f"Warning: {len(result['metadata']['missing_r1'])} R2 files without R1")
    for missing in result['metadata']['missing_r1']:
        print(f"  - {missing}")
```

### 5. Use ValidationReport for Complete Tracking

```python
report = ValidationReport("./logs/report.txt")

# Add all results
report.write(ref_result, "genome")
report.write(mod_result, "genome")
report.write(genome_check, "genomexgenome")

# Generate comprehensive report
report.flush(format="text")  # or "json"
```

---

## Troubleshooting

### Problem: "Genome sequence count mismatch"

**Error:**
```
Genome sequence count mismatch: reference has 2 sequence(s), modified has 1 sequence(s)
```

**Cause:** Different number of sequences in genomes (e.g., plasmid present in one genome but not the other).

**Solutions:**
1. **Disable check** if acceptable:
   ```python
   settings = GenomeXGenomeSettings(same_number_of_sequences=False)
   ```

2. **Split plasmids** before validation:
   ```python
   genome_settings = GenomeValidator.Settings(plasmid_split=True)
   ref_result = validate_genome(config.ref_genome, genome_settings)
   ```

3. **Filter sequences** by length:
   ```python
   settings = GenomeValidator.Settings(min_sequence_length=10000)
   ```

---

### Problem: "Genome sequence ID mismatch"

**Error:**
```
Genome sequence ID mismatch: reference-only: ['chr1'], modified-only: ['chromosome1']
```

**Cause:** Different sequence IDs between genomes.

**Solutions:**
1. **Standardize IDs** before validation:
   ```python
   settings = GenomeValidator.Settings(replace_id_with='chr')
   ref_result = validate_genome(config.ref_genome, settings)
   mod_result = validate_genome(config.mod_genome, settings)
   ```

2. **Disable ID check** if acceptable:
   ```python
   settings = GenomeXGenomeSettings(same_sequence_ids=False)
   ```

---

### Problem: "Found R2 file(s) without matching R1"

**Error:**
```
Found R2 file(s) without matching R1 for base name 'sample': ['sample_R2.fastq.gz']
```

**Cause:** Paired-end reads are incomplete (R1 file missing).

**Solutions:**
1. **Add missing R1 file** to config.json:
   ```json
   {
     "reads": [
       {"filename": "sample_R1.fastq.gz", "ngs_type": "illumina"},
       {"filename": "sample_R2.fastq.gz", "ngs_type": "illumina"}
     ]
   }
   ```

2. **Allow missing R1** if acceptable:
   ```python
   settings = ReadXReadSettings(allow_missing_r1=True)
   result = readxread_validation(reads_results, settings)
   # This will produce WARNING instead of ERROR
   ```

3. **Remove orphan R2 file** from config if not needed.

---

### Problem: "No paired-end patterns detected"

**Message:**
```
No paired-end patterns detected in results, validation passed (nothing to check)
```

**Cause:** Read files don't follow Illumina paired-end naming patterns.

**Solutions:**
1. **Rename files** to follow supported patterns (see [Supported Patterns](#supported-paired-end-patterns))
2. **This is not an error** - validation passes if no patterns detected (nothing to validate)
3. **Single-end reads** don't need paired-end validation

---

## See Also

- [API_REFERENCE.md](API_REFERENCE.md) - API documentation
- [VALIDATORS.md](VALIDATORS.md) - Validator documentation
- [LOGGING_REPORTING.md](LOGGING_REPORTING.md) - Report generation
- [EXAMPLES.md](EXAMPLES.md) - Complete workflow examples
