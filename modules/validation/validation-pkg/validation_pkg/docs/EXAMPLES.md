# Code Examples

Practical code examples for using the validation_pkg package.

## Table of Contents

- [Quick Start](#quick-start)
- [Common Use Cases](#common-use-cases)
- [Complete Workflows](#complete-workflows)
- [Integration Examples](#integration-examples)

---

## Quick Start

### Minimal Example

```python
from validation_pkg import ConfigManager, validate_genome, validate_reads

# Load configuration
config = ConfigManager.load("./data/inputs/config.json")

# Validate genomes
ref_result = validate_genome(config.ref_genome)
mod_result = validate_genome(config.mod_genome)

# Validate reads
reads_results = validate_reads(config.reads)

print(f"✓ Validated {ref_result.num_sequences} reference sequences")
print(f"✓ Validated {len(reads_results)} read files")
```

### With Custom Settings

```python
from validation_pkg import GenomeValidator, ReadValidator

# Custom genome settings
genome_settings = GenomeValidator.Settings(
    plasmid_split=True,
    min_sequence_length=1000
)

# Custom read settings
read_settings = ReadValidator.Settings(
    validation_level='trust',
    check_invalid_chars=False
)

# Validate with settings
ref_result = validate_genome(config.ref_genome, genome_settings)
reads_results = validate_reads(config.reads, read_settings)
```

---

## Common Use Cases

### 1. Bacterial Genome with Plasmids

```python
from validation_pkg import (
    ConfigManager,
    GenomeValidator,
    validate_genome
)

config = ConfigManager.load("config.json")

# Split chromosome and plasmids
settings = GenomeValidator.Settings(
    plasmid_split=True,
    replace_id_with='seq',
    min_sequence_length=1000
)

# Validate reference genome
ref_result = validate_genome(config.ref_genome, settings)

print(f"Main chromosome: {ref_result.output_file}")
print(f"Number of plasmids: {ref_result.plasmid_count}")
if ref_result.plasmid_filenames:
    for plasmid in ref_result.plasmid_filenames:
        print(f"  - {plasmid}")

# Validate modified genome with same settings
mod_result = validate_genome(config.mod_genome, settings)

# Compare structure
print(f"\nStructure comparison:")
print(f"  Ref sequences: {ref_result.num_sequences}")
print(f"  Mod sequences: {mod_result.num_sequences}")
```

### 2. Paired-End Illumina Reads

```python
from validation_pkg import (
    ConfigManager,
    validate_reads,
    readxread_validation,
    ReadXReadSettings
)

config = ConfigManager.load("config.json")

# Validate all read files
reads_results = validate_reads(config.reads)

# Check paired-end completeness
settings = ReadXReadSettings(pair_end_basename=True)
check = readxread_validation(reads_results, settings)

# Display results
for result in reads_results:
    print(f"{result.output_filename}:")
    print(f"  Reads: {result.num_reads:,}")
    if result.read_number:
        print(f"  Paired-end: R{result.read_number} (base: {result.base_name})")
    if result.n50:
        print(f"  N50: {result.n50:,} bp")

# Check completeness
if not check['passed']:
    print("\n⚠ Warning: Paired-end validation issues:")
    for error in check['errors']:
        print(f"  {error}")
```

### 3. Format Conversion (GenBank → FASTA)

```python
from validation_pkg import ConfigManager, GenomeValidator

config = ConfigManager.load("config.json")

# Input: genome.gbk (GenBank)
# Output: genome.fasta (FASTA)
settings = GenomeValidator.Settings(
    replace_id_with='chr',
    coding_type=None  # Uncompressed output
)

result = validate_genome(config.ref_genome, settings)

print(f"Converted GenBank → FASTA:")
print(f"  Input: {result.input_file}")
print(f"  Output: {result.output_file}")
print(f"  Sequences: {result.num_sequences}")
```

### 4. Feature File Processing

```python
from validation_pkg import ConfigManager, FeatureValidator

config = ConfigManager.load("config.json")

# Sort features and standardize chromosome names
settings = FeatureValidator.Settings(
    sort_by_position=True,
    replace_id_with='chr1',
    coding_type='gz'  # Compressed output
)

result = validate_feature(config.ref_feature, settings)

print(f"Feature validation:")
print(f"  Features: {result.num_features:,}")
print(f"  Types: {', '.join(result.feature_types)}")
print(f"  Sequences: {', '.join(result.sequence_ids)}")
```

### 5. Fast Validation (Trust Mode)

```python
from validation_pkg import (
    ConfigManager,
    GenomeValidator,
    ReadValidator,
    validate_genome,
    validate_reads
)

config = ConfigManager.load("config.json")

# Use trust mode for speed
genome_settings = GenomeValidator.Settings(validation_level='trust')
read_settings = ReadValidator.Settings(validation_level='trust')

# Validate quickly
ref_result = validate_genome(config.ref_genome, genome_settings)
reads_results = validate_reads(config.reads, read_settings)

print(f"Fast validation completed:")
print(f"  Genome time: {ref_result.elapsed_time:.2f}s")
print(f"  Reads time: {sum(r.elapsed_time for r in reads_results):.2f}s")
```

---

## Complete Workflows

### Workflow 1: Full Validation Pipeline

```python
from validation_pkg import (
    setup_logging,
    get_logger,
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

# 1. Setup logging
setup_logging(
    console_level='INFO',
    log_file='./logs/validation.log'
)
logger = get_logger()

# 2. Create validation report
report = ValidationReport('./logs/validation_report.txt')

# 3. Load configuration
logger.info("Loading configuration...")
config = ConfigManager.load("./data/inputs/config.json")

# 4. Validate genomes
logger.info("Validating genomes...")
ref_genome = validate_genome(config.ref_genome)
mod_genome = validate_genome(config.mod_genome)
report.write(ref_genome, 'genome')
report.write(mod_genome, 'genome')

logger.info(f"✓ Reference: {ref_genome.num_sequences} sequences, "
           f"{ref_genome.total_genome_size:,} bp")
logger.info(f"✓ Modified: {mod_genome.num_sequences} sequences, "
           f"{mod_genome.total_genome_size:,} bp")

# 5. Validate reads
logger.info("Validating reads...")
reads_results = validate_reads(config.reads)
report.write(reads_results, 'read')

total_reads = sum(r.num_reads for r in reads_results if r.num_reads)
logger.info(f"✓ Validated {len(reads_results)} read files, "
           f"{total_reads:,} total reads")

# 6. Validate features (if present)
if config.ref_feature:
    logger.info("Validating features...")
    feature_result = validate_feature(config.ref_feature)
    report.write(feature_result, 'feature')
    logger.info(f"✓ Validated {feature_result.num_features:,} features")

# 7. Inter-file validation: Genome comparison
logger.info("Checking genome compatibility...")
genome_settings = GenomeXGenomeSettings(
    same_number_of_sequences=True,
    same_sequence_ids=False
)
genome_check = genomexgenome_validation(
    ref_genome,
    mod_genome,
    settings=genome_settings
)
report.write(genome_check, 'genomexgenome')

if genome_check['passed']:
    logger.info("✓ Genomes are compatible")
else:
    logger.warning("⚠ Genome compatibility issues detected")

# 8. Inter-file validation: Paired-end reads
logger.info("Checking paired-end completeness...")
read_settings = ReadXReadSettings(pair_end_basename=True)
read_check = readxread_validation(reads_results, settings=read_settings)
report.write(read_check, 'readxread')

if read_check['passed']:
    logger.info("✓ Paired-end reads are complete")
else:
    logger.warning("⚠ Paired-end validation issues detected")

# 9. Generate final report
report.flush(format='text')
logger.info(f"Report generated: {report.report_path}")

# 10. Summary
print("\n" + "="*60)
print("VALIDATION SUMMARY")
print("="*60)
print(f"Genomes validated: 2")
print(f"Reads validated: {len(reads_results)}")
print(f"Genome check: {'PASSED' if genome_check['passed'] else 'FAILED'}")
print(f"Read check: {'PASSED' if read_check['passed'] else 'FAILED'}")
print(f"\nReport: {report.report_path}")
print("="*60)
```

### Workflow 2: Batch Processing with Error Handling

```python
from validation_pkg import (
    ConfigManager,
    validate_genome,
    GenomeValidator
)
from validation_pkg.exceptions import ValidationError

# Multiple configuration files
config_files = [
    "./data/sample1/config.json",
    "./data/sample2/config.json",
    "./data/sample3/config.json"
]

results = []
errors = []

for config_file in config_files:
    try:
        print(f"\nProcessing {config_file}...")

        # Load configuration
        config = ConfigManager.load(config_file)

        # Custom settings
        settings = GenomeValidator.Settings(
            plasmid_split=True,
            min_sequence_length=500
        )

        # Validate
        ref_result = validate_genome(config.ref_genome, settings)
        mod_result = validate_genome(config.mod_genome, settings)

        results.append({
            'config': config_file,
            'ref': ref_result,
            'mod': mod_result
        })

        print(f"✓ Success: {ref_result.num_sequences} ref sequences, "
              f"{mod_result.num_sequences} mod sequences")

    except ValidationError as e:
        print(f"✗ Failed: {e}")
        errors.append({
            'config': config_file,
            'error': str(e)
        })

# Summary
print(f"\n{'='*60}")
print(f"Processed {len(config_files)} samples:")
print(f"  Successful: {len(results)}")
print(f"  Failed: {len(errors)}")

if errors:
    print(f"\nFailed samples:")
    for err in errors:
        print(f"  {err['config']}: {err['error']}")
```

### Workflow 3: Progressive Validation

```python
from validation_pkg import (
    ConfigManager,
    GenomeValidator,
    validate_genome
)

config = ConfigManager.load("config.json")

# Step 1: Minimal validation (fastest, check format only)
print("Step 1: Quick format check...")
minimal_settings = GenomeValidator.Settings(validation_level='minimal')
try:
    result = validate_genome(config.ref_genome, minimal_settings)
    print(f"✓ Format is valid")
except Exception as e:
    print(f"✗ Format error: {e}")
    exit(1)

# Step 2: Trust validation (medium speed, basic checks)
print("\nStep 2: Basic validation...")
trust_settings = GenomeValidator.Settings(validation_level='trust')
result = validate_genome(config.ref_genome, trust_settings)
print(f"✓ Basic validation passed: {result.num_sequences} sequences")

# Step 3: Strict validation (thorough, with statistics)
print("\nStep 3: Comprehensive validation...")
strict_settings = GenomeValidator.Settings(
    validation_level='strict',
    plasmid_split=True,
    min_sequence_length=1000
)
result = validate_genome(config.ref_genome, strict_settings)

print(f"✓ Comprehensive validation complete:")
print(f"  Sequences: {result.num_sequences}")
print(f"  Total size: {result.total_genome_size:,} bp")
print(f"  GC content: {result.gc_content:.2f}%")
print(f"  N50: {result.n50:,} bp")
```

---

## Integration Examples

### Integration 1: Custom Validation Script

```python
#!/usr/bin/env python3
"""
Custom validation script for EFSA pipeline.
"""

import sys
import argparse
from pathlib import Path
from validation_pkg import (
    setup_logging,
    ConfigManager,
    validate_genome,
    validate_reads,
    GenomeValidator,
    ReadValidator,
    ValidationReport
)

def main():
    parser = argparse.ArgumentParser(description='Validate genomic data')
    parser.add_argument('config', help='Path to config.json')
    parser.add_argument('--trust', action='store_true',
                       help='Use trust mode for speed')
    parser.add_argument('--split-plasmids', action='store_true',
                       help='Split plasmids to separate files')
    args = parser.parse_args()

    # Setup
    setup_logging(console_level='INFO', log_file='./logs/validation.log')

    # Settings
    genome_settings = GenomeValidator.Settings(
        validation_level='trust' if args.trust else 'strict',
        plasmid_split=args.split_plasmids
    )
    read_settings = ReadValidator.Settings(
        validation_level='trust' if args.trust else 'strict'
    )

    # Load and validate
    config = ConfigManager.load(args.config)

    ref_result = validate_genome(config.ref_genome, genome_settings)
    mod_result = validate_genome(config.mod_genome, genome_settings)
    reads_results = validate_reads(config.reads, read_settings)

    # Report
    report = ValidationReport('./logs/validation_report.txt')
    report.write(ref_result, 'genome')
    report.write(mod_result, 'genome')
    report.write(reads_results, 'read')
    report.flush(format='text')

    print(f"✓ Validation complete. Report: {report.report_path}")
    return 0

if __name__ == '__main__':
    sys.exit(main())
```

**Usage:**
```bash
# Strict validation
python validate.py ./data/inputs/config.json

# Fast validation with plasmid splitting
python validate.py ./data/inputs/config.json --trust --split-plasmids
```

### Integration 2: Nextflow Integration

```python
#!/usr/bin/env python3
"""
Validation module for Nextflow pipeline.
Outputs JSON for downstream processing.
"""

import json
import sys
from validation_pkg import (
    ConfigManager,
    validate_genome,
    validate_reads,
    genomexgenome_validation
)

def main(config_path, output_json):
    # Validate
    config = ConfigManager.load(config_path)
    ref_result = validate_genome(config.ref_genome)
    mod_result = validate_genome(config.mod_genome)
    reads_results = validate_reads(config.reads)

    # Inter-file check
    genome_check = genomexgenome_validation(ref_result, mod_result)

    # Prepare output for Nextflow
    output = {
        'ref_genome': {
            'file': ref_result.output_file,
            'sequences': ref_result.num_sequences,
            'size': ref_result.total_genome_size
        },
        'mod_genome': {
            'file': mod_result.output_file,
            'sequences': mod_result.num_sequences,
            'size': mod_result.total_genome_size
        },
        'reads': [
            {
                'file': r.output_file,
                'reads': r.num_reads,
                'ngs_type': r.ngs_type_detected
            }
            for r in reads_results
        ],
        'validation': {
            'genome_check_passed': genome_check['passed'],
            'errors': genome_check.get('errors', [])
        }
    }

    # Write JSON
    with open(output_json, 'w') as f:
        json.dump(output, f, indent=2)

    # Exit with error if validation failed
    if not genome_check['passed']:
        sys.exit(1)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
```

### Integration 3: Snakemake Integration

```python
# Snakemake rule for validation

rule validate_genomes:
    input:
        config="data/inputs/config.json"
    output:
        ref="data/valid/reference_genome.fasta",
        mod="data/valid/modified_genome.fasta",
        report="logs/validation_report.txt"
    log:
        "logs/validation.log"
    run:
        from validation_pkg import (
            setup_logging,
            ConfigManager,
            validate_genome,
            ValidationReport
        )

        setup_logging(log_file=log[0])
        config = ConfigManager.load(input.config)

        ref_result = validate_genome(config.ref_genome)
        mod_result = validate_genome(config.mod_genome)

        report = ValidationReport(output.report)
        report.write(ref_result, 'genome')
        report.write(mod_result, 'genome')
        report.flush(format='text')
```

---

## See Also

- [API_REFERENCE.md](API_REFERENCE.md) - API documentation
- [VALIDATORS.md](VALIDATORS.md) - Validator documentation
- [ADVANCED_FEATURES.md](ADVANCED_FEATURES.md) - Advanced features
- [ERROR_HANDLING.md](ERROR_HANDLING.md) - Error handling
