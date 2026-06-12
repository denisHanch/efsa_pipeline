# Logging

Documentation for the logging system in the `validation_pkg` package.

## Table of Contents

- [Logging System](#logging-system)
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

### Log File Auto-Increment

If `validation.log` already exists, a new run creates `validation_001.log`, `validation_002.log`, etc. — previous logs are never overwritten.

---

## Examples

### Example 1: Basic Logging

```python
from validation_pkg import setup_logging, get_logger, ConfigManager, GenomeValidator

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

    result = GenomeValidator(config.ref_genome).run()
    logger.info(f"✓ Validated {result.num_sequences} sequences")

except Exception as e:
    logger.error(f"Validation failed: {e}")
    raise
```

### Example 2: Debug Mode with Timing

```python
from validation_pkg import setup_logging, get_logger, GenomeValidator

# Enable debug logging
setup_logging(console_level='DEBUG', log_file='./logs/debug.log')
logger = get_logger()

# Time validation
logger.start_timer('genome_validation')

result = GenomeValidator(config.ref_genome).run()

elapsed = logger.stop_timer('genome_validation')
logger.info(f"Validation completed in {elapsed:.2f}s")

# Check all timings
for name, time in logger.get_timers().items():
    logger.debug(f"{name}: {time:.2f}s")
```

### Example 3: Custom Validation Issues

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

# Run custom validation
result = GenomeValidator(config.ref_genome).run()
validate_custom_rules(result)

# Review issues
issues = logger.validation_issues
if issues:
    logger.warning(f"Found {len(issues)} validation issues")
    for issue in issues:
        logger.warning(f"{issue['category']}: {issue['message']}")
```

---

## See Also

- [API_REFERENCE.md](API_REFERENCE.md) - API documentation
- [ERROR_HANDLING.md](ERROR_HANDLING.md) - Exception handling
