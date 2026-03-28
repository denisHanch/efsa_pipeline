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
| `mod_genome_filename` | GenomeConfig | ❌ | Modified genome (FASTA or GenBank) |
| `ref_plasmid_filename` | GenomeConfig | ❌ | Reference plasmid (FASTA or GenBank) |
| `mod_plasmid_filename` | GenomeConfig | ❌ | Modified plasmid (FASTA or GenBank) |
| `reads` | List[ReadConfig] | ✅ | Read files (FASTQ/BAM), minimum one |
| `ref_feature_filename` | FeatureConfig | ❌ | Features for reference genome (BED, GFF, GTF) |
| `mod_feature_filename` | FeatureConfig | ❌ | Features for modified genome (BED, GFF, GTF) |
| `options` | dict | ❌ | Additional options (e.g., `{"threads": 8}`) |

**Important:**
- All file paths must be **relative to the config file directory**
- At minimum, provide `ref_genome_filename` and `reads`

---

## File Format Support

### Genome and Plasmid Files

| Format | Extensions | 
|--------|-----------|
| FASTA | `.fa`, `.fasta`, `.fna`, `.faa` |
| GenBank | `.gb`, `.gbk`, `.genbank` |

### Feature Files

| Format | Extensions | 
|--------|-----------|
| GFF | `.gff`, `.gff3` |
| GTF | `.gtf`, `.gff2` |
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

**Optional fields:**
```json
{
  "filename": "path/to/genome.fasta.gz",
  "validation_level": "trust",
  "n_sequence_limit": 10
}
```

| Field | Type | Default | Description |
|---|---|---|---|
| `filename` | string | — | Path to genome file (relative to config) |
| `validation_level` | string | global / `"trust"` | Per-file validation level override |
| `threads` | integer | global / auto | Per-file thread count override |
| `n_sequence_limit` | integer | `5` | Maximum allowed number of sequences. Applies to `ref_genome_filename` and `mod_genome_filename` only — **ignored with a warning on plasmids**. When the genome contains more sequences than this limit, the assembly is considered too fragmented: a warning is logged, the file is copied to `data/valid/` as-is, and the pipeline will not run SyRI or ref-vs-mod comparison. Set higher for highly fragmented assemblies. |

**Example:**
```json
"mod_genome_filename": {
  "filename": "data/modified.fasta",
  "n_sequence_limit": 10
}
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

`ngs_type` is case-insensitive — `"Illumina"`, `"ILLUMINA"`, and `"illumina"` are all accepted.

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

**Allowed global options:**

| Option | Type | Default | Description |
|---|---|---|---|
| `threads` | integer / `null` | auto-detect | Number of threads. `null` or omit = auto-detect from CPU cores. System warns if value exceeds available cores. |
| `validation_level` | string | `"trust"` | `"strict"`, `"trust"`, or `"minimal"`. Case-insensitive — stored as uppercase internally. |
| `logging_level` | string | `"INFO"` | `"DEBUG"`, `"INFO"`, `"WARNING"`, or `"ERROR"`. Case-insensitive. |
| `type` | string | `"prokaryote"` | `"prokaryote"` or `"eukaryote"`. Case-insensitive. Eukaryote skips inter-genome comparison. |
| `force_defragment_ref` | boolean | `false` | **Unsupported workaround — use at your own risk.** Merges all reference contigs into one sequence before validation. Use only when the reference is too fragmented for any workflow. All downstream results (alignment, variant calling, feature mapping) may be incorrect or meaningless. See warning below. |

> **⚠ Warning: `force_defragment_ref`**
>
> This option is an **unsupported workaround** and is used entirely **at your own responsibility**.
> The EFSA pipeline does not support fragmented reference genomes. Merging contigs artificially
> alters the coordinate space of the reference, which means all downstream results — inter-genome
> alignment, variant calling, feature coordinate mapping — **may be incorrect or meaningless**.
> Do **not** use these results for regulatory submissions or biological conclusions without expert review.
> Obtain a properly assembled reference genome instead.
>
> When `force_defragment_ref` is set to `true` in `config.json`, it takes priority over the
> `--force-defragment-ref` CLI flag. When set to `false` in `config.json`, the CLI flag is
> **ignored** even if provided.

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
    "logging_level": "DEBUG",
    "type": "prokaryote",
    "force_defragment_ref": false
  }
}
```

**Result:** ALL files will use the specified options. Any option not set here falls back to its default.

#### Option priority

Options can be set from three sources. The first matching source wins:

```
config.json "options"  →  CLI flags (validation.sh)  →  built-in defaults
```

A value in `config.json` always takes priority over a CLI flag, and a CLI flag takes priority over the built-in default. This applies to all options including `force_defragment_ref`: if the config sets it to `false`, passing `--force-defragment-ref` on the command line has no effect.

**Validation:**
- Invalid option names → `ConfigurationError` (e.g., `"abc"` not allowed)
- Invalid threads → `ConfigurationError` (e.g., negative numbers)
- Invalid validation_level → `ConfigurationError` (e.g., `"invalid_level"`; `"Strict"` and `"STRICT"` are both accepted)
- Invalid logging_level → `ConfigurationError` (e.g., `"verbose"`, must be DEBUG/INFO/WARNING/ERROR; case-insensitive)
- Invalid type → `ConfigurationError` (e.g., `"bacteria"`; `"Prokaryote"` and `"PROKARYOTE"` are both accepted)
- Invalid force_defragment_ref → `ConfigurationError` (must be `true` or `false`, not a string)

### File-Level Settings

Override global options for specific files by adding settings to individual file entries:
- `validation_level`: `"strict"`, `"trust"`, or `"minimal"` (overrides global; case-insensitive)
- `threads`: Number of threads for compression (int, overrides global)
- `n_sequence_limit`: Maximum number of sequences allowed in a genome file (int, default: `5`). Applies to `ref_genome_filename` and `mod_genome_filename` only; ignored with a warning on plasmids.

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
    "validation_level": "strict",
    "threads": 8,
    "n_sequence_limit": 5
  },
  "mod_genome_filename": {
    "filename": "mod.fasta.gz",
    "validation_level": "strict",
    "threads": 8,
    "n_sequence_limit": 5
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
    "logging_level": "INFO",
    "type": "prokaryote"
  }
}
```
