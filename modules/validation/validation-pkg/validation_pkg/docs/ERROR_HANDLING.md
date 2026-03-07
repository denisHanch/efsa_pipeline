# Error Handling

Documentation for exception handling in the validation_pkg package.

## Table of Contents

- [Exception Hierarchy](#exception-hierarchy)
- [Common Exceptions](#common-exceptions)
- [Exception Details](#exception-details)
- [Error Handling Patterns](#error-handling-patterns)
- [Troubleshooting Guide](#troubleshooting-guide)

---

## Exception Hierarchy

All validation exceptions inherit from `ValidationError` base class.

```
ValidationError (base)
├── ConfigurationError
├── FileNotFoundError
├── CompressionError
├── FileFormatError
│   ├── FastaFormatError
│   ├── GenBankFormatError
│   ├── BedFormatError
│   ├── GffFormatError
│   ├── FastqFormatError
│   └── BamFormatError
├── GenomeValidationError
├── ReadValidationError
├── FeatureValidationError
└── InterFileValidationError
```

---

## Common Exceptions

### ValidationError

Base exception for all validation errors.

**Import:**
```python
from validation_pkg.exceptions import ValidationError
```

**Catching All Validation Errors:**
```python
try:
    result = validate_genome(config.ref_genome)
except ValidationError as e:
    print(f"Validation failed: {e}")
```

**Use Case:** Catch all package-specific errors in one handler.

---

### ConfigurationError

Raised when configuration file is invalid.

**Common Causes:**
- Missing required fields in config.json
- Invalid JSON syntax
- Invalid option values
- Invalid file paths

**Example:**
```python
from validation_pkg import ConfigManager
from validation_pkg.exceptions import ConfigurationError

try:
    config = ConfigManager.load("config.json")
except ConfigurationError as e:
    print(f"Configuration error: {e}")
    # Fix config.json and retry
```

**Common Error Messages:**
```
Missing required field: ref_genome_filename
Missing required field: reads
Invalid validation_level 'invalid'. Must be one of: strict, trust, minimal
'threads' must be a positive integer, got -1
Invalid global options: {'unknown_option'}. Only 'threads', 'validation_level', 'logging_level' are allowed
```

**Solutions:**
- Check config.json structure against [Configuration Guide](../CONFIG_GUIDE.md)
- Validate JSON syntax with a JSON validator
- Ensure all file paths exist
- Check option values match allowed values

---

### FileNotFoundError

Raised when required files are missing.

**Note:** This is `validation_pkg.exceptions.FileNotFoundError`, not built-in `FileNotFoundError`.

**Import:**
```python
from validation_pkg.exceptions import FileNotFoundError as ValidationFileNotFoundError
```

**Common Causes:**
- Incorrect file path in config.json
- File doesn't exist
- Path is relative but config is in different directory

**Example:**
```python
try:
    config = ConfigManager.load("config.json")
except ValidationFileNotFoundError as e:
    print(f"File not found: {e}")
    # Check file paths in config.json
```

**Solutions:**
- Verify file paths are correct
- Use absolute paths or paths relative to config.json location
- Check file permissions

---

### FileFormatError

Base exception for file format errors.

**Subclasses:**
- `FastaFormatError`: FASTA format issues
- `GenBankFormatError`: GenBank format issues
- `BedFormatError`: BED format issues
- `GffFormatError`: GFF/GTF format issues
- `FastqFormatError`: FASTQ format issues
- `BamFormatError`: BAM format issues

**Example:**
```python
from validation_pkg.exceptions import FastaFormatError

try:
    result = validate_genome(config.ref_genome)
except FastaFormatError as e:
    print(f"FASTA format error: {e}")
    # Check input file format
```

---

### GenomeValidationError

Raised when genome validation fails.

**Common Causes:**
- Empty sequence IDs (when `allow_empty_id=False`)
- Empty sequences (when `allow_empty_sequences=False`)
- No sequences found in file
- File parsing failure

**Example:**
```python
from validation_pkg.exceptions import GenomeValidationError

try:
    result = validate_genome(config.ref_genome)
except GenomeValidationError as e:
    print(f"Genome validation failed: {e}")
    # Check genome file quality
```

**Common Error Messages:**
```
Sequence at index 0 has no ID
Sequence 'chr1' has zero length
No sequences found in FASTA file
```

**Solutions:**
- Use `allow_empty_id=True` if empty IDs are acceptable
- Use `allow_empty_sequences=True` if empty sequences are acceptable
- Verify the file contains valid sequence data

---

### ReadValidationError

Raised when read validation fails.

**Common Causes:**
- Empty read IDs
- Invalid characters in sequences
- Empty read file
- Inconsistent quality scores

**Example:**
```python
from validation_pkg.exceptions import ReadValidationError

try:
    result = validate_read(config.reads[0])
except ReadValidationError as e:
    print(f"Read validation failed: {e}")
```

**Common Error Messages:**
```
Sequence has no ID
Contains 150 invalid character(s): ['X', 'Y']
Empty read file
Paired-end validation failed: Missing R1 for R2
```

**Solutions:**
- Check read file quality
- Use `allow_empty_id=True` if acceptable
- Use `check_invalid_chars=False` to skip character validation
- Verify paired-end completeness

---

### FeatureValidationError

Raised when feature validation fails.

**Common Causes:**
- gffread tool not available
- Invalid GFF/GTF/BED format that gffread cannot parse

**Example:**
```python
from validation_pkg.exceptions import FeatureValidationError

try:
    result = validate_feature(config.ref_feature)
except FeatureValidationError as e:
    print(f"Feature validation failed: {e}")
```

**Common Error Messages:**
```
gffread tool required.
gffread failed: <stderr from gffread>
```

**Solutions:**
- Install gffread: `conda install -c bioconda gffread`
- Verify GFF/GTF/BED format is valid
- Check gffread output for specific coordinate or syntax errors

---

### CompressionError

Raised when file decompression fails.

**Common Causes:**
- Corrupted compressed file
- Unsupported compression format
- Incomplete download
- Wrong file extension

**Example:**
```python
from validation_pkg.exceptions import CompressionError

try:
    result = validate_genome(config.ref_genome)
except CompressionError as e:
    print(f"Compression error: {e}")
```

**Solutions:**
- Verify file integrity (checksum)
- Re-download file
- Check file extension matches actual compression
- Decompress manually and use uncompressed file

---

### InterFileValidationError

Raised when inter-file consistency checks fail.

**Common Causes:**
- Sequence count mismatch between genomes
- Sequence ID mismatch
- Missing R1 for R2 reads

**Example:**
```python
from validation_pkg.exceptions import InterFileValidationError

try:
    result = genomexgenome_validation(ref_result, mod_result)
except InterFileValidationError as e:
    print(f"Inter-file validation failed: {e}")
```

**Solutions:**
- See [INTERFILE_VALIDATION.md](INTERFILE_VALIDATION.md) troubleshooting section

---

## Exception Details

### Accessing Exception Information

All exceptions provide detailed error messages:

```python
try:
    result = validate_genome(config.ref_genome)
except GenomeValidationError as e:
    # Error message
    print(f"Error: {e}")

    # Exception type
    print(f"Type: {type(e).__name__}")

    # Stack trace (for debugging)
    import traceback
    traceback.print_exc()
```

### Logger Integration

All exceptions are automatically logged:

```python
from validation_pkg import setup_logging, get_logger

setup_logging(log_file="./logs/validation.log")

try:
    result = validate_genome(config.ref_genome)
except ValidationError as e:
    # Error is already logged to file
    # Additional context can be added
    logger = get_logger()
    logger.error(f"Additional context: {e}")
```

---

## Error Handling Patterns

### Pattern 1: Catch All Validation Errors

```python
from validation_pkg import ConfigManager, validate_genome
from validation_pkg.exceptions import ValidationError

try:
    config = ConfigManager.load("config.json")
    ref_result = validate_genome(config.ref_genome)
    mod_result = validate_genome(config.mod_genome)

except ValidationError as e:
    print(f"Validation failed: {e}")
    exit(1)
```

**Use Case:** Simple scripts where any validation error should stop execution.

---

### Pattern 2: Specific Exception Handling

```python
from validation_pkg.exceptions import (
    ConfigurationError,
    FileNotFoundError,
    GenomeValidationError
)

try:
    config = ConfigManager.load("config.json")
    result = validate_genome(config.ref_genome)

except ConfigurationError as e:
    print(f"Fix your config.json: {e}")
    exit(1)

except FileNotFoundError as e:
    print(f"File missing: {e}")
    exit(1)

except GenomeValidationError as e:
    print(f"Genome quality issue: {e}")
    # Continue with warnings
    logger.warning(f"Genome validation issue: {e}")

except ValidationError as e:
    print(f"Other validation error: {e}")
    exit(1)
```

**Use Case:** Different error types need different handling strategies.

---

### Pattern 3: Retry with Relaxed Settings

```python
from validation_pkg import GenomeValidator
from validation_pkg.exceptions import GenomeValidationError

# First try: strict validation
try:
    settings = GenomeValidator.Settings(
        allow_empty_id=False,
        allow_empty_sequences=False
    )
    result = validate_genome(config.ref_genome, settings)

except GenomeValidationError as e:
    print(f"Strict validation failed: {e}")
    print("Retrying with relaxed settings...")

    # Retry with relaxed settings
    relaxed_settings = GenomeValidator.Settings(
        allow_empty_id=True,
        allow_empty_sequences=True
    )
    result = validate_genome(config.ref_genome, relaxed_settings)
    print("Validation succeeded with relaxed settings")
```

**Use Case:** Handle borderline quality data gracefully.

---

### Pattern 4: Continue on Error

```python
from validation_pkg import validate_reads
from validation_pkg.exceptions import ReadValidationError

successful_results = []
failed_files = []

for read_config in config.reads:
    try:
        result = validate_read(read_config)
        successful_results.append(result)

    except ReadValidationError as e:
        print(f"Skipping {read_config.filename}: {e}")
        failed_files.append((read_config.filename, str(e)))

print(f"Successfully validated: {len(successful_results)}")
print(f"Failed: {len(failed_files)}")

if failed_files:
    print("\nFailed files:")
    for filename, error in failed_files:
        print(f"  {filename}: {error}")
```

**Use Case:** Batch processing where some files may fail but others can continue.

---

### Pattern 5: Logging and Re-raising

```python
from validation_pkg import get_logger
from validation_pkg.exceptions import ValidationError

logger = get_logger()

try:
    result = validate_genome(config.ref_genome)

except ValidationError as e:
    # Log with additional context
    logger.error(f"Validation failed for {config.ref_genome.filename}")
    logger.error(f"Error details: {e}")

    # Add to validation report
    logger.add_validation_issue(
        level='ERROR',
        category='genome',
        message=f'Validation failed: {e}',
        details={'file': str(config.ref_genome.filepath)}
    )

    # Re-raise for caller to handle
    raise
```

**Use Case:** Add context to errors while preserving exception propagation.

---

### Pattern 6: Graceful Degradation

```python
from validation_pkg import GenomeValidator
from validation_pkg.exceptions import ValidationError

def validate_with_fallback(genome_config):
    """Try strict, fall back to trust, then minimal."""

    # Try strict
    try:
        settings = GenomeValidator.Settings(validation_level='strict')
        return validate_genome(genome_config, settings)
    except ValidationError as e:
        print(f"Strict failed: {e}, trying trust mode...")

    # Try trust
    try:
        settings = GenomeValidator.Settings(validation_level='trust')
        return validate_genome(genome_config, settings)
    except ValidationError as e:
        print(f"Trust failed: {e}, trying minimal mode...")

    # Try minimal
    try:
        settings = GenomeValidator.Settings(validation_level='minimal')
        return validate_genome(genome_config, settings)
    except ValidationError as e:
        print(f"All validation levels failed: {e}")
        raise

result = validate_with_fallback(config.ref_genome)
```

**Use Case:** Maximize validation when possible, but accept lower quality if necessary.

---

## Troubleshooting Guide

### Problem: "Missing required field: ref_genome_filename"

**Error Type:** `ConfigurationError`

**Cause:** config.json doesn't have required field.

**Solution:**
```json
{
  "ref_genome_filename": {"filename": "reference.fasta"},
  "reads": [
    {"filename": "reads.fastq", "ngs_type": "illumina"}
  ]
}
```

---

### Problem: "File not found: /path/to/file"

**Error Type:** `FileNotFoundError`

**Cause:** File doesn't exist at specified path.

**Solutions:**
1. Check file path is correct
2. Use absolute paths or paths relative to config.json
3. Verify file exists: `ls -la /path/to/file`

---

### Problem: "Invalid characters in sequence"

**Error Type:** `ReadValidationError`

**Cause:** Read sequences contain non-ATCGN characters. Note: character validation is only supported for reads (`check_invalid_chars`), not for genome sequences.

**Solutions:**
1. Fix input file (remove invalid characters)
2. Disable check (disabled by default, `check_invalid_chars=False`):

---

### Problem: Feature file parsing fails

**Error Type:** `FeatureValidationError`

**Cause:** gffread cannot parse the input file (invalid format or missing tool).

**Solutions:**
1. Install gffread: `conda install -c bioconda gffread`
2. Verify the file format is valid GFF3, GTF, or BED
3. Check gffread output for specific errors

---

### Problem: "Genome sequence count mismatch"

**Error Type:** `InterFileValidationError`

**Cause:** Different number of sequences in genomes.

**Solutions:**
- See [INTERFILE_VALIDATION.md](INTERFILE_VALIDATION.md#problem-genome-sequence-count-mismatch)

---

### Problem: "Found R2 file(s) without matching R1"

**Error Type:** `ReadValidationError` or `InterFileValidationError`

**Cause:** Paired-end reads incomplete.

**Solutions:**
- See [INTERFILE_VALIDATION.md](INTERFILE_VALIDATION.md#problem-found-r2-files-without-matching-r1)

---

### Problem: "Corrupted gzip file"

**Error Type:** `CompressionError`

**Cause:** File is damaged or incomplete.

**Solutions:**
1. Re-download file
2. Verify checksum (if available)
3. Decompress manually:
   ```bash
   gunzip -t file.fastq.gz  # Test integrity
   gunzip file.fastq.gz     # Decompress
   ```

---

### General Debugging Tips

1. **Enable DEBUG logging:**
   ```python
   from validation_pkg import setup_logging
   setup_logging(console_level='DEBUG', log_file='debug.log')
   ```

2. **Check log files:**
   ```bash
   tail -f ./logs/validation_*.log
   ```

3. **Validate input files separately:**
   ```bash
   # Test FASTA format
   grep "^>" genome.fasta | head

   # Test FASTQ format
   head -n 4 reads.fastq

   # Test GFF format
   head features.gff
   ```

4. **Use minimal mode for format testing:**
   ```python
   # Minimal mode fails fast if format is wrong
   settings = GenomeValidator.Settings(validation_level='minimal')
   result = validate_genome(config.ref_genome, settings)
   ```

5. **Check BioPython compatibility:**
   ```python
   from Bio import SeqIO
   records = list(SeqIO.parse("genome.fasta", "fasta"))
   print(f"BioPython parsed {len(records)} sequences")
   ```

---

## See Also

- [API_REFERENCE.md](API_REFERENCE.md) - API documentation
- [VALIDATORS.md](VALIDATORS.md) - Validator documentation
- [INTERFILE_VALIDATION.md](INTERFILE_VALIDATION.md) - Inter-file validation
- [LOGGING_REPORTING.md](LOGGING_REPORTING.md) - Logging and reporting
