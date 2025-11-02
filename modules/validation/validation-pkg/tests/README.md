# Testing Guide

This document explains the testing infrastructure for the validation package. It covers test organization, running tests, creating fixtures, and debugging test failures.

## Table of Contents

- [Quick Start](#quick-start)
- [Test Organization](#test-organization)
- [Running Tests](#running-tests)
- [Test Structure](#test-structure)
- [Creating Test Fixtures](#creating-test-fixtures)
- [Test Naming Conventions](#test-naming-conventions)
- [Debugging Tests](#debugging-tests)
- [Coverage Reports](#coverage-reports)

---

## Quick Start

**Install test dependencies:**

```bash
pip install -e .
# or
pip install -r requirements.txt
```

**Run all tests:**

```bash
# From project root
python -m pytest tests/

# With verbose output
python -m pytest tests/ -v

# With coverage
python -m pytest tests/ --cov=validation_pkg
```

**Run specific test file:**

```bash
python -m pytest tests/test_genome_validator.py -v
```

**Run specific test class:**

```bash
python -m pytest tests/test_settings.py::TestUpdateMethod -v
```

**Run specific test:**

```bash
python -m pytest tests/test_formats.py::TestGenomeFormat::test_fasta_extensions -v
```

---

## Test Organization

### Test Suite Overview

The package has **460 tests** organized by module:

```
tests/
├── test_genome_validator.py      # Genome validation (FASTA/GenBank)
├── test_read_validator.py        # Read validation (FASTQ/BAM)
├── test_feature_validator.py     # Feature validation (GFF/GTF/BED)
├── test_config_manager.py        # Configuration parsing
├── test_logger.py                # Logging system
├── test_exceptions.py            # Exception hierarchy
├── test_file_handler.py          # File I/O and compression (54 tests)
├── test_formats.py               # Format enums (51 tests)
├── test_settings.py              # BaseSettings pattern (42 tests)
├── test_parallel.py              # Parallel processing (9 tests)
└── fixtures/                     # Test data files
    ├── genome/
    ├── reads/
    └── features/
```

### Test Categories

**Validator Tests:**
- Genome validation (FASTA/GenBank format validation, plasmid handling)
- Read validation (FASTQ/BAM validation, quality scores, paired-end)
- Feature validation (GFF/GTF/BED coordinate checking, sorting)

**Core System Tests:**
- Configuration management (JSON parsing, path resolution, format detection)
- Logging (structured logging, issue tracking, report generation)
- Exception handling (error hierarchy, error messages)

**Utility Tests:**
- File handling (compression detection, format conversion, I/O operations)
- Format enums (extension matching, BioPython mapping, flexible input)
- Settings (immutable pattern, serialization, field validation)
- Parallel processing (multi-file validation, thread safety, error handling)

---

## Running Tests

### Basic Commands

**All tests:**
```bash
python -m pytest tests/
```

**Specific test file:**
```bash
python -m pytest tests/test_genome_validator.py
```

**Specific test class:**
```bash
python -m pytest tests/test_settings.py::TestUpdateMethod
```

**Specific test:**
```bash
python -m pytest tests/test_formats.py::TestGenomeFormat::test_fasta_extensions
```

### Useful Options

**Verbose output** (show each test name):
```bash
python -m pytest tests/ -v
```

**Very verbose** (show print statements):
```bash
python -m pytest tests/ -vv -s
```

**Stop on first failure:**
```bash
python -m pytest tests/ -x
```

**Run last failed tests only:**
```bash
python -m pytest tests/ --lf
```

**Run tests matching pattern:**
```bash
python -m pytest tests/ -k "genome"       # All genome-related tests
python -m pytest tests/ -k "not slow"     # Skip slow tests
```

**Parallel test execution** (requires pytest-xdist):
```bash
pip install pytest-xdist
python -m pytest tests/ -n auto  # Use all CPU cores
python -m pytest tests/ -n 4     # Use 4 workers
```

### Test Markers

Tests can be marked for selective running:

```python
import pytest

@pytest.mark.slow
def test_large_file_validation():
    """Test with 10GB file (slow)."""
    pass

@pytest.mark.integration
def test_full_pipeline():
    """Integration test of complete workflow."""
    pass
```

**Run only marked tests:**
```bash
python -m pytest tests/ -m slow          # Run slow tests only
python -m pytest tests/ -m "not slow"    # Skip slow tests
python -m pytest tests/ -m integration   # Run integration tests
```

---

## Test Structure

### Test Class Organization

Tests are organized into classes by functionality:

```python
import pytest
from validation_pkg.validators.genome_validator import GenomeValidator

class TestGenomeValidator:
    """Tests for genome validator."""

    class TestSettings:
        """Tests for genome validator settings."""

        def test_default_values(self):
            """Test default settings values."""
            settings = GenomeValidator.Settings()
            assert settings.validation_level == 'strict'
            assert settings.min_sequence_length == 0

        def test_immutable_update(self):
            """Test settings immutable pattern."""
            settings = GenomeValidator.Settings()
            new_settings = settings.update(validation_level='trust')

            assert settings.validation_level == 'strict'  # Original unchanged
            assert new_settings.validation_level == 'trust'  # New instance

    class TestValidation:
        """Tests for validation methods."""

        def test_fasta_validation(self, tmp_path):
            """Test FASTA file validation."""
            # Test implementation
            pass

        def test_genbank_validation(self, tmp_path):
            """Test GenBank file validation."""
            # Test implementation
            pass
```

### Fixtures

**Built-in pytest fixtures:**
- `tmp_path`: Temporary directory (Path object)
- `tmp_path_factory`: Factory for multiple temporary directories
- `capsys`: Capture stdout/stderr
- `monkeypatch`: Modify objects/environment variables

**Custom fixtures** (defined in `conftest.py` or test files):

```python
import pytest
from pathlib import Path

@pytest.fixture
def sample_fasta(tmp_path):
    """Create a sample FASTA file for testing."""
    fasta_file = tmp_path / "sample.fasta"
    fasta_file.write_text(
        ">seq1\n"
        "ATCGATCGATCG\n"
        ">seq2\n"
        "GCTAGCTAGCTA\n"
    )
    return fasta_file

@pytest.fixture
def genome_config(sample_fasta):
    """Create a genome config object."""
    from validation_pkg.config_manager import GenomeConfig
    from validation_pkg.utils.formats import GenomeFormat, CodingType

    return GenomeConfig(
        filename=str(sample_fasta),
        file_format=GenomeFormat.FASTA,
        coding_type=CodingType.NONE,
        settings_dict={},
        extra={}
    )

# Usage
def test_genome_validation(genome_config):
    """Test genome validation with fixture."""
    validator = GenomeValidator(genome_config, "/tmp/output")
    output_path = validator.validate()
    assert Path(output_path).exists()
```

### Test Patterns

**Pattern 1: Arrange-Act-Assert (AAA)**

```python
def test_settings_update():
    """Test settings update method."""
    # Arrange
    settings = GenomeValidator.Settings()

    # Act
    new_settings = settings.update(min_sequence_length=100)

    # Assert
    assert settings.min_sequence_length == 0      # Original unchanged
    assert new_settings.min_sequence_length == 100  # New instance
```

**Pattern 2: Parameterized Tests**

```python
@pytest.mark.parametrize("extension,expected_format", [
    ("genome.fasta", GenomeFormat.FASTA),
    ("genome.fa", GenomeFormat.FASTA),
    ("genome.fna", GenomeFormat.FASTA),
    ("genome.gbk", GenomeFormat.GENBANK),
    ("genome.gb", GenomeFormat.GENBANK),
])
def test_format_detection(extension, expected_format):
    """Test format detection from various extensions."""
    detected = GenomeFormat(Path(extension).suffix)
    assert detected == expected_format
```

**Pattern 3: Exception Testing**

```python
import pytest
from validation_pkg.exceptions import ValidationError

def test_invalid_validation_level():
    """Test that invalid validation level raises error."""
    with pytest.raises(ValueError, match="validation_level must be"):
        GenomeValidator.Settings(validation_level='invalid')

def test_negative_min_length():
    """Test that negative min_sequence_length raises error."""
    settings = GenomeValidator.Settings()
    with pytest.raises(ValueError, match="min_sequence_length must be >= 0"):
        settings.update(min_sequence_length=-1)
```

---

## Creating Test Fixtures

### File Fixtures

Test data files are stored in `tests/fixtures/`:

```
tests/fixtures/
├── genome/
│   ├── valid_fasta.fa
│   ├── valid_genbank.gbk
│   ├── invalid_fasta.fa
│   └── compressed.fasta.gz
├── reads/
│   ├── valid_illumina_R1.fastq
│   ├── valid_illumina_R2.fastq
│   ├── valid_ont.fastq
│   └── invalid_quality.fastq
└── features/
    ├── valid.gff
    ├── valid.gtf
    ├── valid.bed
    └── invalid_coords.gff
```

### Creating FASTA Fixtures

```python
from pathlib import Path

def create_fasta_fixture(output_path: Path, num_sequences: int = 3):
    """Create a valid FASTA file for testing."""
    content = []
    for i in range(num_sequences):
        content.append(f">sequence_{i}")
        content.append("ATCGATCGATCGATCGATCG")

    output_path.write_text("\n".join(content) + "\n")
    return output_path

# Usage in test
def test_fasta_parsing(tmp_path):
    fasta_file = create_fasta_fixture(tmp_path / "test.fasta")
    # Test validation
```

### Creating FASTQ Fixtures

```python
def create_fastq_fixture(output_path: Path, num_records: int = 10):
    """Create a valid FASTQ file for testing."""
    content = []
    for i in range(num_records):
        content.append(f"@read_{i}")
        content.append("ATCGATCGATCG")
        content.append("+")
        content.append("IIIIIIIIIIII")  # Quality scores

    output_path.write_text("\n".join(content) + "\n")
    return output_path
```

### Creating GFF Fixtures

```python
def create_gff_fixture(output_path: Path, num_features: int = 5):
    """Create a valid GFF file for testing."""
    content = ["##gff-version 3"]
    for i in range(num_features):
        start = i * 1000 + 1
        end = start + 500
        content.append(
            f"chr1\ttest\tgene\t{start}\t{end}\t.\t+\t.\tID=gene_{i}"
        )

    output_path.write_text("\n".join(content) + "\n")
    return output_path
```

### Compressed Fixtures

```python
import gzip

def create_compressed_fasta(output_path: Path):
    """Create a gzip-compressed FASTA file."""
    content = ">seq1\nATCGATCG\n>seq2\nGCTAGCTA\n"

    with gzip.open(output_path, 'wt') as f:
        f.write(content)

    return output_path

# Usage
def test_compressed_genome(tmp_path):
    gz_file = create_compressed_fasta(tmp_path / "genome.fasta.gz")
    # Test validation
```

---

## Test Naming Conventions

### File Names

- Prefix with `test_`: `test_genome_validator.py`
- Match module name: `genome_validator.py` → `test_genome_validator.py`
- One test file per module (exceptions for large modules)

### Class Names

- Prefix with `Test`: `TestGenomeValidator`
- Describe what's being tested: `TestSettings`, `TestValidation`
- Nest classes for organization: `TestGenomeValidator.TestSettings`

### Function Names

- Prefix with `test_`: `test_fasta_validation()`
- Be descriptive: `test_settings_immutable_update()`
- Use underscores for readability
- State what's being tested and expected result:
  - `test_invalid_level_raises_error()`
  - `test_compression_detection_from_gz_extension()`

**Good examples:**
```python
def test_settings_update_returns_new_instance()
def test_invalid_validation_level_raises_value_error()
def test_fasta_parsing_with_multiple_sequences()
def test_parallel_processing_with_four_workers()
```

**Bad examples:**
```python
def test_update()  # Too vague
def test_1()       # No description
def testFasta()    # Wrong naming convention
```

---

## Debugging Tests

### Print Debugging

Use `-s` flag to see print statements:

```bash
python -m pytest tests/test_genome_validator.py -v -s
```

```python
def test_validation(genome_config):
    print(f"Config: {genome_config}")  # Will be shown with -s
    validator = GenomeValidator(genome_config, "/tmp")
    print(f"Validator created: {validator}")
    result = validator.validate()
    print(f"Result: {result}")
    assert result is not None
```

### Using pdb Debugger

**Drop into debugger on failure:**
```bash
python -m pytest tests/test_genome_validator.py --pdb
```

**Set breakpoint in test:**
```python
def test_validation(genome_config):
    validator = GenomeValidator(genome_config, "/tmp")

    import pdb; pdb.set_trace()  # Debugger stops here

    result = validator.validate()
    assert result is not None
```

### Verbose Assertion Output

```python
def test_validation():
    expected = "genome.fasta"
    actual = validator.get_filename()

    # Good: Shows both values on failure
    assert actual == expected, f"Expected {expected}, got {actual}"

    # Better: pytest shows diff automatically
    assert actual == expected
```

### Capturing Logs

```python
import logging

def test_with_logs(caplog):
    """Test that captures log messages."""
    with caplog.at_level(logging.INFO):
        validator.validate()

    # Check log messages
    assert "Validating genome" in caplog.text
    assert any("Success" in record.message for record in caplog.records)
```

### Common Issues

**Issue: Test passes locally but fails in CI**
- Check file paths (use `tmp_path`, not hardcoded paths)
- Check platform differences (Linux vs. macOS vs. Windows)
- Check environment variables

**Issue: Test is flaky (passes sometimes, fails sometimes)**
- Check for timing issues (add delays or retries)
- Check for leftover state (clean up in teardown)
- Check for parallel test interference (use unique temp paths)

**Issue: Test fixture not found**
- Check fixture location (`conftest.py` or same file)
- Check fixture scope (`function`, `class`, `module`, `session`)
- Check fixture name matches parameter name

---

## Coverage Reports

### Generate Coverage Report

**Text summary:**
```bash
python -m pytest tests/ --cov=validation_pkg
```

**HTML report** (detailed, browsable):
```bash
python -m pytest tests/ --cov=validation_pkg --cov-report=html
open htmlcov/index.html  # macOS
xdg-open htmlcov/index.html  # Linux
```

**Terminal report with missing lines:**
```bash
python -m pytest tests/ --cov=validation_pkg --cov-report=term-missing
```

**XML report** (for CI):
```bash
python -m pytest tests/ --cov=validation_pkg --cov-report=xml
```

### Interpreting Coverage

**Good coverage targets:**
- **Overall**: >80% coverage
- **Core validators**: >90% coverage
- **Utilities**: >85% coverage
- **Error handling**: 100% coverage for exception paths

**Coverage report example:**
```
Name                                Stmts   Miss  Cover   Missing
-----------------------------------------------------------------
validation_pkg/__init__.py             45      2    96%   78-79
validation_pkg/config_manager.py      312      8    97%   245, 301-307
validation_pkg/validators/genome.py   428     15    97%   89, 234-241, 567
-----------------------------------------------------------------
TOTAL                                3847     98    97%
```

**Focus on:**
- **Stmts**: Total statements
- **Miss**: Uncovered statements
- **Cover**: Coverage percentage
- **Missing**: Line numbers not covered

---

## Best Practices

### Writing Tests

1. **Test one thing per test**
   ```python
   # Good
   def test_settings_update_validation_level()
   def test_settings_update_min_sequence_length()

   # Bad
   def test_settings_update()  # Tests everything at once
   ```

2. **Use descriptive assertions**
   ```python
   # Good
   assert output_path.exists(), f"Output file not created: {output_path}"

   # Bad
   assert output_path.exists()  # No context on failure
   ```

3. **Clean up after tests**
   ```python
   def test_creates_file(tmp_path):
       """Test uses tmp_path - automatic cleanup."""
       output = tmp_path / "output.txt"
       # No manual cleanup needed
   ```

4. **Test error cases**
   ```python
   def test_invalid_input_raises_error():
       """Test error handling."""
       with pytest.raises(ValueError):
           validator.validate_invalid_data()
   ```

5. **Use fixtures for common setup**
   ```python
   @pytest.fixture
   def validator(genome_config):
       return GenomeValidator(genome_config, "/tmp")

   def test_validation(validator):
       result = validator.validate()
       # ...
   ```

### Running Tests in CI

**GitHub Actions example:**

```yaml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - uses: actions/setup-python@v2
        with:
          python-version: '3.9'
      - name: Install dependencies
        run: |
          pip install -e .
          pip install pytest pytest-cov
      - name: Run tests
        run: pytest tests/ --cov=validation_pkg --cov-report=xml
      - name: Upload coverage
        uses: codecov/codecov-action@v2
```

---

## See Also

- [VALIDATOR_DESIGN.md](../docs/VALIDATOR_DESIGN.md) - Validator architecture
- [README.md](../README.md) - Main package documentation
- [examples/](../examples/) - Usage examples
- [pytest documentation](https://docs.pytest.org/) - Official pytest docs
