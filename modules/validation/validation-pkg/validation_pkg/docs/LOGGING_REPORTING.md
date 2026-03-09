# Logging and Reporting

Documentation for logging and report generation in the validation_pkg package.

## Table of Contents

- [Logging System](#logging-system)
- [Validation Reports](#validation-reports)
- [Log File Structure](#log-file-structure)
- [Examples](#examples)

---

## Logging System

The validation package uses `structlog` for structured, colored logging.

### Setup Logging

```python
from validation_pkg import setup_logging

# Basic setup (console only)
setup_logging()

# With log file
setup_logging(
    console_level='INFO',
    log_file='./logs/validation.log'
)

# Debug mode
setup_logging(
    console_level='DEBUG',
    log_file='./logs/debug.log'
)
```

**Parameters:**
- `console_level` (str): Console verbosity - 'DEBUG', 'INFO', 'WARNING', 'ERROR'
- `log_file` (Path or str, optional): Path to log file (auto-creates directory)

### Get Logger

```python
from validation_pkg import get_logger

logger = get_logger()

# Standard logging
logger.info("Processing genome...")
logger.warning("Large number of sequences detected")
logger.error("Validation failed")
logger.debug("Detailed debugging information")

# Structured logging with context
logger.info("File validated",
    file_name="genome.fasta",
    num_sequences=2,
    total_size=4000000
)
```

### Logging Levels

| Level | Use Case | Example |
|-------|----------|---------|
| **DEBUG** | Detailed diagnostics | "Parsing sequence 1 of 2..." |
| **INFO** | Progress updates | "✓ Genome validated successfully" |
| **WARNING** | Issues that don't stop execution | "Warning: More than 2 sequences detected" |
| **ERROR** | Validation failures | "ERROR: Duplicate sequence IDs found" |
| **CRITICAL** | System failures | "CRITICAL: Cannot write output file" |

### Colored Output

Console output is automatically colored:
- **DEBUG**: Gray
- **INFO**: Blue
- **WARNING**: Yellow
- **ERROR**: Red
- **CRITICAL**: Bold Red

### Context-Aware Logging

Add context to log messages:

```python
logger = get_logger()

# File context
logger.info("Starting validation", file_context="genome.fasta")

# Worker context (for parallel processing)
logger.debug("Processing chunk", worker_id=1)

# Category context
logger.info("Validating features", category="feature")
```

### Validation Issues

Track specific validation issues:

```python
logger = get_logger()

# Add validation issue
logger.add_validation_issue(
    level='ERROR',
    category='genome',
    message='Duplicate sequence IDs found',
    details={'duplicate_ids': ['chr1', 'chr1']}
)

# Retrieve issues (direct attribute access)
for issue in logger.validation_issues:
    print(f"{issue['level']}: {issue['message']}")
```

### Timing Measurements

Track operation timing:

```python
logger = get_logger()

# Start timer
logger.start_timer('genome_validation')

# ... do work ...

# Stop timer and log
elapsed = logger.stop_timer('genome_validation')
logger.info(f"Validation completed in {elapsed:.2f}s")

# Get all timings
timings = logger.get_timers()
```

---

## Validation Reports

Generate comprehensive validation reports.

### ValidationReport Class

```python
from validation_pkg import ValidationReport

# Create report
report = ValidationReport('./logs/validation_report.txt')

# Add validation results
report.write(genome_result, 'genome')
report.write(reads_results, 'read')  # Can be list
report.write(feature_result, 'feature')

# Add inter-file validation
report.write(genome_check, 'genomexgenome')
report.write(read_check, 'readxread')

# Generate report file
report.flush(format='text')  # or 'json'
```

### Report Formats

#### Text Format (Default)

Human-readable report with sections:

```python
report.flush(format='text')
```

**Output Structure:**
```
==================================================================================================
  VALIDATION PIPELINE REPORT
==================================================================================================
  Generated:       2024-12-05 15:30:45
  Total Duration:  45.23s

SUMMARY
==================================================================================================
  Overall Status: ✓ PASSED
  Files Processed: 5
    ├─ Genomes:  2
    ├─ Reads:    2
    └─ Features: 1
  Inter-file Validations: 2
    ├─ Passed:  2 ✓

==================================================================================================
FILE VALIDATION RESULTS
==================================================================================================

[1] GENOME FILE
  Input:  /data/inputs/reference.fasta
  Output: /data/valid/reference_genome.fasta
  Time:   12.34s

  Statistics:
    Sequences: 1
    Total Length: 4,641,652 bp
    Sequence IDs: NC_000913.3

...
```

#### JSON Format

Machine-readable report:

```python
report.flush(format='json')
```

**Output Structure:**
```json
{
  "report_metadata": {
    "generated": "2024-12-05T15:30:45",
    "duration_seconds": 45.23,
    "version": "1.0.0"
  },
  "summary": {
    "total_files": 5,
    "genome_files": 2,
    "read_files": 2,
    "feature_files": 1,
    "overall_status": "PASSED"
  },
  "file_validations": [
    {
      "validator_type": "genome",
      "output_data": {...},
      "input_settings": {...}
    }
  ],
  "inter_file_validations": [...]
}
```

### Auto-Incrementing Reports

Reports auto-increment if file exists:

```python
# First run: validation_report.txt
# Second run: validation_report_1.txt
# Third run: validation_report_2.txt
report = ValidationReport('./logs/validation_report.txt')
```

---

## Log File Structure

### Log File Location

Default: Auto-generated in `./logs/` directory

```python
setup_logging(log_file='./logs/validation_20241205_153045.log')
```

### Log File Format

Each log entry includes:
- Timestamp
- Log level
- Message
- Context (file, worker, category)
- Structured data (if provided)

**Example Log Entries:**
```
2024-12-05 15:30:45 [INFO    ] Loading configuration from: config.json
2024-12-05 15:30:46 [INFO    ] [genome] Validating genome file: reference.fasta
2024-12-05 15:30:47 [DEBUG   ] [genome] Parsing FASTA format...
2024-12-05 15:30:48 [DEBUG   ] [genome] Found 1 sequences
2024-12-05 15:30:49 [INFO    ] [genome] ✓ Genome validated successfully
2024-12-05 15:30:50 [WARNING ] [read] More than 1000000 reads detected (may be slow)
2024-12-05 15:30:55 [ERROR   ] [feature] Validation failed: Invalid coordinates
```

### Log File Rotation

For long-running processes, implement log rotation:

```python
from pathlib import Path
from datetime import datetime

# Create timestamped log file
timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
log_file = Path(f'./logs/validation_{timestamp}.log')

setup_logging(log_file=log_file)
```

---

## Examples

### Example 1: Basic Logging

```python
from validation_pkg import setup_logging, get_logger, ConfigManager, validate_genome

# Setup logging
setup_logging(
    console_level='INFO',
    log_file='./logs/validation.log'
)

logger = get_logger()

# Load and validate
try:
    logger.info("Starting validation workflow")

    config = ConfigManager.load("config.json")
    logger.info(f"Loaded config with {len(config.reads)} read files")

    result = validate_genome(config.ref_genome)
    logger.info(f"✓ Validated {result.num_sequences} sequences")

except Exception as e:
    logger.error(f"Validation failed: {e}")
    raise
```

### Example 2: Complete Workflow with Report

```python
from validation_pkg import (
    setup_logging,
    get_logger,
    ConfigManager,
    validate_genome,
    validate_reads,
    genomexgenome_validation,
    ValidationReport
)

# Setup
setup_logging(console_level='INFO', log_file='./logs/validation.log')
logger = get_logger()

# Create report
report = ValidationReport('./logs/validation_report.txt')

# Load config
logger.info("Loading configuration...")
config = ConfigManager.load("config.json")

# Validate genomes
logger.info("Validating genomes...")
ref_result = validate_genome(config.ref_genome)
mod_result = validate_genome(config.mod_genome)
report.write(ref_result, 'genome')
report.write(mod_result, 'genome')

# Validate reads
logger.info("Validating reads...")
reads_results = validate_reads(config.reads)
report.write(reads_results, 'read')

# Inter-file validation
logger.info("Running inter-file validation...")
genome_check = genomexgenome_validation(ref_result, mod_result)
report.write(genome_check, 'genomexgenome')

if genome_check['passed']:
    logger.info("✓ All validations passed")
else:
    logger.warning("⚠ Some validations failed")

# Generate report
report.flush(format='text')
logger.info(f"Report written to: {report.report_path}")
```

### Example 3: Debug Mode with Timing

```python
from validation_pkg import setup_logging, get_logger, validate_genome

# Enable debug logging
setup_logging(console_level='DEBUG', log_file='./logs/debug.log')
logger = get_logger()

# Time validation
logger.start_timer('genome_validation')

result = validate_genome(config.ref_genome)

elapsed = logger.stop_timer('genome_validation')
logger.info(f"Validation completed in {elapsed:.2f}s")

# Check all timings
for name, time in logger.get_timers().items():
    logger.debug(f"{name}: {time:.2f}s")
```

### Example 4: Custom Validation Issues

```python
from validation_pkg import get_logger

logger = get_logger()

# Track custom issues
def validate_custom_rules(genome_result):
    if genome_result.num_sequences > 10:
        logger.add_validation_issue(
            level='WARNING',
            category='genome',
            message='Unusually high sequence count',
            details={
                'num_sequences': genome_result.num_sequences,
                'threshold': 10
            }
        )

    if genome_result.gc_content < 30 or genome_result.gc_content > 70:
        logger.add_validation_issue(
            level='WARNING',
            category='genome',
            message='GC content outside typical range',
            details={
                'gc_content': genome_result.gc_content,
                'normal_range': '30-70%'
            }
        )

# Run custom validation
result = validate_genome(config.ref_genome)
validate_custom_rules(result)

# Review issues
issues = logger.validation_issues
if issues:
    logger.warning(f"Found {len(issues)} validation issues")
    for issue in issues:
        logger.warning(f"{issue['category']}: {issue['message']}")
```

### Example 5: Multi-Format Reports

```python
from validation_pkg import ValidationReport

report = ValidationReport('./logs/validation_report.txt')

# Add all results
report.write(ref_result, 'genome')
report.write(mod_result, 'genome')
report.write(reads_results, 'read')

# Generate both text and JSON
report.flush(format='text')  # Creates validation_report.txt
report.flush(format='json')  # Creates validation_report.json

print("Reports generated:")
print(f"  Text: {report.report_path}")
print(f"  JSON: {report.report_path.with_suffix('.json')}")
```

---

## See Also

- [API_REFERENCE.md](API_REFERENCE.md) - API documentation
- [ERROR_HANDLING.md](ERROR_HANDLING.md) - Exception handling
- [EXAMPLES.md](EXAMPLES.md) - Complete examples
