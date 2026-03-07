# validation_pkg API Reference

Technical API documentation for the `validation_pkg` bioinformatics validation package.

## Table of Contents

- [Quick Start](#quick-start)
- [Functional API](#functional-api)
- [Configuration Management](#configuration-management)
- [Validator Classes](#validator-classes)
- [Output Metadata](#output-metadata)
- [See Also](#see-also)

---

## Quick Start

### Installation

```bash
# The package is installed as part of the EFSA pipeline
# For standalone use (production):
pip install -e /path/to/validation-pkg

# For development (includes testing dependencies):
pip install -r /path/to/validation-pkg/requirements-dev.txt
# Or using extras:
pip install -e "/path/to/validation-pkg[dev]"
```

**Note:** As of the latest optimization, test dependencies (pytest, etc.) have been moved to `requirements-dev.txt` to reduce container initialization time. Production installations only require the core dependencies: biopython, numpy, pysam, structlog, and typing_extensions.

### Basic Usage

```python
from validation_pkg import ConfigManager, validate_genome, validate_reads, validate_feature

# Load configuration
config = ConfigManager.load("config.json")

# Validate genome files
ref_result = validate_genome(config.ref_genome)
mod_result = validate_genome(config.mod_genome)

# Validate read files
reads_results = validate_reads(config.reads)

# Validate feature file (if present)
if config.ref_feature:
    feature_result = validate_feature(config.ref_feature)
```

---

## Functional API

The functional API provides simplified wrapper functions for common validation tasks.

### Genome Validation

#### `validate_genome(genome_config, settings=None)`

Validate a single genome file.

**Parameters:**
- `genome_config` (GenomeConfig): Genome configuration from ConfigManager
- `settings` (GenomeValidator.Settings, optional): Custom validation settings

**Returns:**
- `GenomeOutputMetadata`: Metadata object with validation results

**Example:**
```python
from validation_pkg import ConfigManager, validate_genome, GenomeValidator

config = ConfigManager.load("config.json")

# Use default settings
result = validate_genome(config.ref_genome)
print(f"Validated {result.num_sequences} sequences")

# Use custom settings
settings = GenomeValidator.Settings()
settings = settings.update(
    plasmid_split=True,
    min_sequence_length=1000
)
result = validate_genome(config.ref_genome, settings)
```

#### `validate_genomes(genome_configs, settings=None)`

Validate multiple genome files with the same settings.

**Parameters:**
- `genome_configs` (List[GenomeConfig]): List of genome configurations
- `settings` (GenomeValidator.Settings, optional): Settings applied to all genomes

**Returns:**
- `List[GenomeOutputMetadata]`: List of metadata objects

**Example:**
```python
# Validate both ref and mod genomes with same settings
genome_list = [config.ref_genome, config.mod_genome]
results = validate_genomes(genome_list, settings)

for result in results:
    print(f"{result.output_filename}: {result.num_sequences} sequences")
```

---

### Read Validation

#### `validate_read(read_config, settings=None)`

Validate a single read file.

**Parameters:**
- `read_config` (ReadConfig): Read configuration from ConfigManager
- `settings` (ReadValidator.Settings, optional): Custom validation settings

**Returns:**
- `ReadOutputMetadata`: Metadata object with validation results

**Example:**
```python
from validation_pkg import validate_read, ReadValidator

# Validate single read file
result = validate_read(config.reads[0])
print(f"Reads: {result.num_reads}")
print(f"N50: {result.n50} bp")
```

#### `validate_reads(read_configs, settings=None)`

Validate multiple read files.

**Parameters:**
- `read_configs` (List[ReadConfig]): List of read configurations
- `settings` (ReadValidator.Settings, optional): Settings applied to all reads

**Returns:**
- `List[ReadOutputMetadata]`: List of metadata objects

**Example:**
```python
# Validate all read files
results = validate_reads(config.reads)

for result in results:
    print(f"{result.output_filename}:")
    print(f"  Reads: {result.num_reads:,}")
    print(f"  Mean length: {result.mean_read_length:.1f} bp")
    print(f"  N50: {result.n50:,} bp")
```

---

### Feature Validation

#### `validate_feature(feature_config, settings=None)`

Validate a single feature annotation file.

**Parameters:**
- `feature_config` (FeatureConfig): Feature configuration from ConfigManager
- `settings` (FeatureValidator.Settings, optional): Custom validation settings

**Returns:**
- `FeatureOutputMetadata`: Metadata object with validation results

**Example:**
```python
from validation_pkg import validate_feature, FeatureValidator

# Validate with custom settings
settings = FeatureValidator.Settings()
settings = settings.update(
    replace_id_with='chr1',
    sort_by_position=True
)
result = validate_feature(config.ref_feature, settings)
print(f"Features: {result.num_features}")
print(f"Types: {result.feature_types}")
```

#### `validate_features(feature_configs, settings=None)`

Validate multiple feature files.

**Parameters:**
- `feature_configs` (List[FeatureConfig]): List of feature configurations
- `settings` (FeatureValidator.Settings, optional): Settings applied to all features

**Returns:**
- `List[FeatureOutputMetadata]`: List of metadata objects

---

## Configuration Management

### ConfigManager

Configuration file loader and parser.

#### `ConfigManager.load(config_path)`

Load and validate a JSON configuration file.

**Parameters:**
- `config_path` (str): Path to config.json file

**Returns:**
- `Config`: Validated configuration object

**Raises:**
- `FileNotFoundError`: If config file doesn't exist
- `ConfigurationError`: If configuration is invalid
- `json.JSONDecodeError`: If JSON is malformed

**Example:**
```python
from validation_pkg import ConfigManager

# Load configuration
config = ConfigManager.load("./data/inputs/config.json")

# Access configuration components
print(f"Reference genome: {config.ref_genome.filepath}")
print(f"Number of read files: {len(config.reads)}")
print(f"Output directory: {config.output_dir}")

# Access global options
threads = config.threads          # Returns int or None
logging_level = config.logging_level   # Returns str (default: 'INFO')
validation_level = config.validation_level  # Returns str (default: 'strict')
```

---

### Config Object

Main configuration container returned by `ConfigManager.load()`.

**Attributes:**
- `ref_genome` (GenomeConfig): Reference genome configuration (required)
- `mod_genome` (GenomeConfig, optional): Modified genome configuration
- `reads` (List[ReadConfig]): Read file configurations (required, non-empty)
- `ref_plasmid` (GenomeConfig, optional): Reference plasmid configuration
- `mod_plasmid` (GenomeConfig, optional): Modified plasmid configuration
- `ref_feature` (FeatureConfig, optional): Reference feature configuration
- `mod_feature` (FeatureConfig, optional): Modified feature configuration
- `options` (dict): Global options (threads, validation_level, logging_level)
- `config_dir` (Path): Directory containing config file
- `output_dir` (Path): Base output directory for validated files

**Properties:**
- `threads`: Thread count from options (returns int or None)
- `logging_level`: Logging level from options (returns str, default: 'INFO')
- `validation_level`: Validation level from options (returns str, default: 'strict')

---

### GenomeConfig

Configuration for genome or plasmid files.

**Attributes:**
- `filename` (str): Original filename
- `filepath` (Path): Absolute resolved path
- `basename` (str): Filename without extension
- `coding_type` (CodingType): Compression format (GZIP, BZIP2, or NONE)
- `detected_format` (GenomeFormat): File format (FASTA or GENBANK)
- `output_dir` (Path): Base output directory
- `global_options` (dict): Merged global + file-level options

---

### ReadConfig

Configuration for sequencing read files.

**Attributes:**
- `filename` (str): Original filename
- `filepath` (Path): Absolute resolved path
- `basename` (str): Filename without extension
- `ngs_type` (str): Sequencing platform ("illumina", "ont", or "pacbio")
- `coding_type` (CodingType): Compression format
- `detected_format` (ReadFormat): File format (FASTQ or BAM)
- `output_dir` (Path): Base output directory
- `global_options` (dict): Merged global + file-level options

---

### FeatureConfig

Configuration for feature annotation files.

**Attributes:**
- `filename` (str): Original filename
- `filepath` (Path): Absolute resolved path
- `basename` (str): Filename without extension
- `coding_type` (CodingType): Compression format
- `detected_format` (FeatureFormat): File format (GFF, GTF, or BED)
- `output_dir` (Path): Base output directory
- `global_options` (dict): Merged global + file-level options

---

## Validator Classes

For advanced use cases, you can instantiate validator classes directly.

### GenomeValidator

```python
from validation_pkg.validators import GenomeValidator

# Create validator with custom settings
settings = GenomeValidator.Settings(
    plasmid_split=True,
    min_sequence_length=1000,
    replace_id_with='chr'
)

validator = GenomeValidator(config.ref_genome, settings)
result = validator.run()
```

See [VALIDATORS.md](VALIDATORS.md) for detailed validator documentation.

### ReadValidator

```python
from validation_pkg.validators import ReadValidator

settings = ReadValidator.Settings(
    validation_level='trust',
    check_invalid_chars=False
)

validator = ReadValidator(config.reads[0], settings)
result = validator.run()
```

### FeatureValidator

```python
from validation_pkg.validators import FeatureValidator

settings = FeatureValidator.Settings(
    sort_by_position=True,
    replace_id_with='chr1'
)

validator = FeatureValidator(config.ref_feature, settings)
result = validator.run()
```

---

## Output Metadata

All validators return metadata objects with validation results.

### GenomeOutputMetadata

**Key Attributes:**
- `output_file` (str): Path to validated genome file
- `num_sequences` (int): Number of sequences
- `total_genome_size` (int): Total length in bp (strict mode)
- `longest_sequence_id` (str): ID of longest sequence
- `longest_sequence_length` (int): Length of longest sequence
- `gc_content` (float): GC content percentage (strict mode)
- `n50` (int): N50 metric in bp (strict mode)
- `plasmid_count` (int): Number of plasmids detected
- `plasmid_filenames` (List[str]): Plasmid output files (if split)
- `sequence_ids` (List[str]): Sequence IDs for inter-file validation
- `sequence_lengths` (dict): Mapping of sequence IDs to lengths
- `validation_level` (str): Validation level used
- `elapsed_time` (float): Processing time in seconds

**Methods:**
- `to_dict()`: Convert to dictionary
- `__str__()`: Human-readable string representation

### ReadOutputMetadata

**Key Attributes:**
- `output_file` (str): Path to validated read file
- `num_reads` (int): Total number of reads
- `n50` (int): N50 read length in bp (strict mode)
- `total_bases` (int): Sum of all read lengths (strict mode)
- `mean_read_length` (float): Average read length (strict mode)
- `longest_read_length` (int): Longest read length (strict mode)
- `shortest_read_length` (int): Shortest read length (strict mode)
- `ngs_type` (str): Configured NGS platform from config (illumina, ont, pacbio)
- `illumina_pairing_detected` (str): Set to 'illumina' if paired-end pattern detected
- `base_name` (str): Illumina paired-end base name extracted from filename
- `read_number` (int): Read number (1 or 2) from pattern detection
- `validation_level` (str): Validation level used
- `elapsed_time` (float): Processing time in seconds

### FeatureOutputMetadata

**Key Attributes:**
- `output_file` (str): Path to validated feature file
- `num_features` (int): Total number of features
- `feature_types` (List[str]): Unique feature types (gene, CDS, exon, etc.)
- `sequence_ids` (List[str]): Sequence IDs referenced by features
- `validation_level` (str): Validation level used
- `elapsed_time` (float): Processing time in seconds

---

## See Also

- [VALIDATORS.md](VALIDATORS.md) - Detailed validator class documentation
- [SETTINGS.md](SETTINGS.md) - Complete settings reference for all validators
- [INTERFILE_VALIDATION.md](INTERFILE_VALIDATION.md) - Inter-file validation API
- [EXAMPLES.md](EXAMPLES.md) - Code examples and use cases
- [ERROR_HANDLING.md](ERROR_HANDLING.md) - Exception handling guide
