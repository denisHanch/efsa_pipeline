# Configuration File Guide

This document defines the **JSON** configuration file that `ConfigManager` uses to locate and process your files (genomes, reads, features).

## Table of Contents

- [Quick Overview](#quick-overview)
- [Configuration Structure](#configuration-structure)
- [File Format Support](#file-format-support)
- [Compression Support](#compression-support)
- [Configuration Objects](#configuration-objects)
- [Validator Settings in Config](#validator-settings-in-config)
  - [Global Options](#global-options-options-field)
  - [File-Level Settings](#file-level-settings)
- [Examples](#examples)
- [Validation Rules](#validation-rules)

---

## Quick Overview

- **Filename:** `config.json` (standard)
- **Location:** Same directory as your input files (paths are relative to config location)
- **Format:** JSON (UTF-8)
- **Purpose:** Specify reference/modified genomes, reads, and optional plasmids & features
---

## Configuration Structure

### Top-Level Keys

| Key | Type | Required | Description |
|---|---|:---:|---|
| `ref_genome_filename` | GenomeConfig | ✅ | Reference genome (FASTA or GenBank) |
| `mod_genome_filename` | GenomeConfig | ✅ | Modified genome (FASTA or GenBank) |
| `ref_plasmid_filename` | GenomeConfig | ❌ | Reference plasmid (FASTA or GenBank) |
| `mod_plasmid_filename` | GenomeConfig | ❌ | Modified plasmid (FASTA or GenBank) |
| `reads` | List[ReadConfig] | ✅ | Read files (FASTQ/BAM), minimum one |
| `ref_feature_filename` | FeatureConfig | ❌ | Features for reference genome (BED, GFF, GTF) |
| `mod_feature_filename` | FeatureConfig | ❌ | Features for modified genome (BED, GFF, GTF) |
| `options` | dict | ❌ | Additional options (e.g., `{"threads": 8}`) |

**Important:**
- All file paths must be **relative to the config file directory**
- At minimum, provide `ref_genome_filename`,`mod_genome_filename` and `reads`

---

## File Format Support

### Genome and Plasmid Files

| Format | Extensions | 
|--------|-----------|
| FASTA | `.fa`, `.fasta`, `.fna` | 
| GenBank | `.gb`, `.gbk`, `.genbank` |

### Feature Files

| Format | Extensions | 
|--------|-----------|
| GFF | `.gff`, `.gff3` |
| GTF | `.gtf` | 
| BED | `.bed` | 

### Read Files

| Format | Extensions | 
|--------|-----------|
| FASTQ | `.fastq`, `.fq` |
| BAM | `.bam` | 

---

## Compression Support

All file types support transparent compression. The package automatically detects and handles compressed files.

### Supported Compression Formats

| Type | Extensions | Detection | Performance |
|------|-----------|-----------|-------------|
| **gzip** | `.gz`, `.gzip` | Automatic | Faster with `pigz`|
| **bzip2** | `.bz2`, `.bzip2` | Automatic | Faster with `pbzip2`|

**Performance tip:** Install `pigz` and `pbzip2` for 3-4x faster compression/decompression:
```bash
# Ubuntu/Debian
sudo apt-get install pigz pbzip2

# macOS
brew install pigz pbzip2
```

### Compression Examples

```json
{
  "ref_genome_filename": {"filename": "genome.fasta.gz"},
  "reads": [
    {"filename": "reads.fastq.bz2", "ngs_type": "illumina"}
  ]
}
```

---

## Configuration Objects

### GenomeConfig Object

Specifies a genome or plasmid file.

**Structure:**
```json
{
  "filename": "path/to/genome.fasta"
}
```

**Optional fields (validator settings):**
```json
{
  "filename": "path/to/genome.fasta.gz",
  "validation_level": "trust",
}
```

**Example:**
```json
"ref_genome_filename": {"filename": "data/reference.gbk"}
```

### ReadConfig Object

Specifies sequencing read files. Supports individual files or directories.

#### Option 1: Individual Files

```json
"reads": [
  {
    "filename": "samples/reads_R1.fastq.gz",
    "ngs_type": "illumina"
  },
  {
    "filename": "samples/reads_R2.fastq.gz",
    "ngs_type": "illumina"
  }
]
```

#### Option 2: Directory of Files

All files in the directory inherit the same settings:

```json
"reads": [
  {
    "directory": "ont_reads/",
    "ngs_type": "ont",
    "validation_level": "trust"
  }
]
```

#### NGS Type Values

| Value |  
|-------|
| `"illumina"` |
| `"ont"` | 
| `"pacbio"` | 

<!-- **Note:** BAM files are automatically detected and converted to FASTQ. Default NGS type for BAM: `pacbio`. -->

### FeatureConfig Object

Specifies feature annotation files (BED, GFF, GTF).

**Structure:**
```json
{
  "filename": "path/to/features.gff"
}
```

**Optional fields (validator settings):**
```json
{
  "filename": "path/to/features.gff",
  "validation_level": "strict",
}
```

**Example:**
```json
"ref_feature_filename": {"filename": "annotations/reference.gff"}
```

---

## Validator Settings in Config

You can specify validator-specific settings directly in your config file at **two levels**:

1. **Global options** (applies to ALL files)
2. **File-level settings** (applies to specific files)

These settings customize validation behavior without modifying code.

### Global Options (`options` field)

**Allowed global options:**
- `threads`: Number of threads for parallel processing (positive integer, default: 8)
  - System automatically detects CPU cores and warns if threads exceed available cores
- `validation_level`: `"strict"`, `"trust"`, or `"minimal"` (default: `"strict"`)
- `logging_level`: `"DEBUG"`, `"INFO"`, `"WARNING"`, `"ERROR"`, or `"CRITICAL"` (default: `"INFO"`)
  - Controls console logging verbosity

**Example:**
```json
{
  "ref_genome_filename": {"filename": "genome.fasta"},
  "reads": [
    {"filename": "reads.fastq", "ngs_type": "illumina"}
  ],
  "options": {
    "threads": 8,
    "validation_level": "trust",
    "logging_level": "DEBUG"
  }
}
```

**Result:** ALL files will use `threads=8`, `validation_level='trust'`, and `logging_level='DEBUG'` by default.

**Validation:**
- Invalid option names → `ConfigurationError` (e.g., `"abc"` not allowed)
- Invalid threads → `ConfigurationError` (e.g., negative numbers)
- Invalid validation_level → `ConfigurationError` (e.g., `"invalid_level"`)
- Invalid logging_level → `ConfigurationError` (e.g., `"verbose"`, must be DEBUG/INFO/WARNING/ERROR/CRITICAL)

### File-Level Settings

Override global options for specific files by adding settings to individual file entries. Same options as for global ones:
- `validation_level`: `"strict"`, `"trust"`, or `"minimal"` (overrides global)
- `threads`: Number of threads for compression (int, overrides global)

**Warnings:** When a file-level setting overrides a global option, a WARNING is logged:
```
WARNING: File-level setting 'validation_level=strict' overrides global option 'validation_level=trust'
```

---

## Examples

### Minimal Configuration

Simplest valid configuration with required fields only:

```json
{
  "ref_genome_filename": {"filename": "ref.fasta"},
  "mod_genome_filename": {"filename": "mod.fasta"},
  "reads": [
    {"filename": "reads.fastq", "ngs_type": "illumina"}
  ]
}
```

### Full Configuration

Complete example with all optional fields and global options:

```json
{
  "ref_genome_filename": {
    "filename": "ref.gbk",
    "validation_level": "strict",
    "threads": 8
  },
  "mod_genome_filename": {
    "filename": "mod.fasta.gz",
    "validation_level": "strict",
    "threads": 8
  },
  "ref_plasmid_filename": {
    "filename": "plasmid_ref.gbk",
    "validation_level": "strict",
    "threads": 8
    },
  "mod_plasmid_filename": {
    "filename": "plasmid_mod.fasta",
    "validation_level": "strict",
    "threads": 8
    },
  "reads": [
    {
      "filename": "illumina_R1.fastq.gz",
      "ngs_type": "illumina",
      "validation_level": "strict",
      "threads": 8
    },
    {
      "filename": "illumina_R2.fastq.gz",
      "ngs_type": "illumina",
      "validation_level": "strict",
      "threads": 8
    },
    {
      "directory": "ont_reads/",
      "ngs_type": "ont",
      "validation_level": "strict",
      "threads": 8
    }
  ],
  "ref_feature_filename": {
    "filename": "features_ref.gff3",
    "validation_level": "strict",
    "threads": 8
  },
  "mod_feature_filename": {
    "filename": "features_mod.bed",
    "validation_level": "strict",
    "threads": 8
    },
  "options": {
    "threads": 8,
    "validation_level": "strict",
    "logging_level":"INFO"
  }
}
```
