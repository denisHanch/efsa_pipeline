# Validation Module Manual

## Overview

The **Validation Module** is designed to check the integrity, format, and encoding of genomic data files used in all pipelines. It supports multiple file types (FASTA, FASTQ, BAM, GTF, GFF, GBK) and can handle compressed files (gzip, bzip2). The module is extensible, allowing for custom validators for each file type.

---

## Directory Structure

```
modules/validation/
    main.py                     # Entry point for validation
    requirements.txt            # Python dependencies
    validation/
        validate.py             # Core Validator classes and subclasses
        utils.py                # Utility functions for file management and validation
        format.py               # Functions for encoding/decoding and format conversion
        README.md               # This manual
    tests/
        run_tests.py            # Pytest-based test runner
        validator.py            # Unit tests for Validator classes
        expected_results.json   # Expected outcomes for test cases
        <test_cases>/           # Individual test folders with config.json and test files
```

---

## Main Components
