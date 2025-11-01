# Configuration File Guide

This document defines the **JSON** configuration file that `ConfigManager` uses to locate and process your bioinformatics files (genomes, reads, features).

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
- **Created by:** Users or analysis pipelines
- **Processed by:** `ConfigManager.load("config.json")`

---

## Configuration Structure

### Top-Level Keys

| Key | Type | Required | Description |
|---|---|:---:|---|
| `ref_genome_filename` | GenomeConfig | ✅ | Reference genome (FASTA or GenBank) |
| `mod_genome_filename` | GenomeConfig | ❌ | Modified genome (FASTA or GenBank) |
| `ref_plasmid_filename` | GenomeConfig | ❌ | Reference plasmid (FASTA or GenBank) |
| `mod_plasmid_filename` | GenomeConfig | ❌ | Modified plasmid (FASTA or GenBank) |
| `reads` | List[ReadConfig] | ❌ | Read files (FASTQ/BAM), minimum one if provided |
| `ref_feature_filename` | FeatureConfig | ❌ | Features for reference genome (BED, GFF, GTF) |
| `mod_feature_filename` | FeatureConfig | ❌ | Features for modified genome (BED, GFF, GTF) |
| `options` | dict | ❌ | Additional options (e.g., `{"threads": 8}`) |

**Important:**
- All file paths must be **relative to the config file directory**
- At minimum, provide `ref_genome_filename` or include `reads`
- Genome validation requires at least one genome file

---

## File Format Support

### Genome and Plasmid Files

| Format | Extensions | Notes |
|--------|-----------|-------|
| FASTA | `.fa`, `.fasta`, `.fna` | Most common genome format |
| GenBank | `.gb`, `.gbk`, `.genbank` | Includes annotations |

### Feature Files

| Format | Extensions | Notes |
|--------|-----------|-------|
| GFF | `.gff`, `.gff3` | General Feature Format |
| GTF | `.gtf` | Gene Transfer Format |
| BED | `.bed` | Browser Extensible Data |

### Read Files

| Format | Extensions | Notes |
|--------|-----------|-------|
| FASTQ | `.fastq`, `.fq` | Quality scores included |
| BAM | `.bam` | Binary alignment format (converted to FASTQ) - **IGNORED NOW** |

---

## Compression Support

All file types support transparent compression. The package automatically detects and handles compressed files.

### Supported Compression Formats

| Type | Extensions | Detection | Performance |
|------|-----------|-----------|-------------|
| **gzip** | `.gz`, `.gzip` | Automatic | Fast with `pigz` (parallel) |
| **bzip2** | `.bz2`, `.bzip2` | Automatic | Fast with `pbzip2` (parallel) |

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
  "min_sequence_length": 500,
  "plasmid_split": true
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

| Value | Description | Typical Read Length |
|-------|-------------|-------------------|
| `"illumina"` | Illumina short reads (SE or PE) | 50-300 bp |
| `"ont"` | Oxford Nanopore long reads | 1-100+ kb |
| `"pacbio"` | PacBio long reads | 10-50+ kb |

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

**New in v0.2.0:** Set global defaults that apply to **all files** in your config automatically.

**Allowed global options:**
- `threads`: Number of threads for parallel processing (positive integer)
- `validation_level`: `"strict"`, `"trust"`, or `"minimal"`

**Why only these two?** These are the only settings that make sense globally. Other settings like `plasmid_split` or `sort_by_position` are file-specific and should be set per-file.

**Example:**
```json
{
  "ref_genome_filename": {"filename": "genome.fasta"},
  "reads": [
    {"filename": "reads.fastq", "ngs_type": "illumina"}
  ],
  "options": {
    "threads": 8,
    "validation_level": "trust"
  }
}
```

**Result:** ALL files will use `threads=8` and `validation_level='trust'` by default.

**Validation:**
- Invalid option names → `ConfigurationError` (e.g., `"plasmid_split"` not allowed)
- Invalid threads → `ConfigurationError` (e.g., negative numbers)
- Invalid validation_level → `ConfigurationError` (e.g., `"invalid_level"`)

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
    "validation_level": "strict"
  },
  "mod_genome_filename": {
    "filename": "mod.fasta.gz"
  },
  "ref_plasmid_filename": {"filename": "plasmid_ref.gbk"},
  "mod_plasmid_filename": {"filename": "plasmid_mod.fasta"},
  "reads": [
    {
      "filename": "illumina_R1.fastq.gz",
      "ngs_type": "illumina"
    },
    {
      "filename": "illumina_R2.fastq.gz",
      "ngs_type": "illumina"
    },
    {
      "directory": "ont_reads/",
      "ngs_type": "ont",
      "validation_level": "minimal"
    }
  ],
  "ref_feature_filename": {
    "filename": "features_ref.gff3",
  },
  "mod_feature_filename": {"filename": "features_mod.bed"},
  "options": {
    "threads": 8,
    "validation_level": "trust"
  }
}
```

**Note:**
- Global `validation_level='trust'` applies to all files
- `ref_genome` overrides with `validation_level='strict'` (WARNING logged)
- `ont_reads` overrides with `validation_level='minimal'` (WARNING logged)
- All other files use global settings

### Directory-Based Reads

Process all files in a directory with the same settings:

```json
{
  "ref_genome_filename": {"filename": "genome.fasta"},
  "reads": [
    {
      "directory": "illumina_reads/",
      "ngs_type": "illumina",
      "validation_level": "trust"
    },
    {
      "directory": "ont_reads/",
      "ngs_type": "ont",
      "validation_level": "minimal"
    }
  ]
}
```

---

## Validation Rules

The package validates your configuration before processing:

1. **Required fields:** At least `ref_genome_filename`,`mod_genome_filename` AND `reads` must be present
2. **File existence:** All specified files must exist and be accessible
3. **Path resolution:** All paths are resolved relative to config directory
4. **Format detection:** File extensions must match supported formats
5. **NGS type:** Read files must specify a valid `ngs_type`
7. **Global options validation:**
   - Only `threads` and `validation_level` allowed in `options`
   - `threads` must be positive integer (or null)
   - `validation_level` must be `"strict"`, `"trust"`, or `"minimal"`

**Common validation errors:**
- Missing required fields
- File not found
- Invalid NGS type
- Invalid validator setting name
- Unsupported file format
- Invalid global option
- Invalid threads value (negative or zero)
- Invalid validation_level value
