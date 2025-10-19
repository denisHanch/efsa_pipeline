#   How to run

download and initiate submodule validation package
```
git submodule update --init --recursive
```
update submodule package
```
git submodule update --remote --merge
```

install package and its dependencies
```
pip install -e /EFSA_workspace/modules/validation/validation_pkg
pip install -r /EFSA_workspace/modules/validation/validation_pkg/requirements.txt
```

Load your data and configuration file into `/EFSA_workspace/data/inputs/` (description of config.json bellow)

Run validation as 
```
python3 main.py /EFSA_workspace/data/inputs/config.json > validation.log
```
 - may take a lot of time - since the reads are still validated read by read.

output of the validation is on `/EFSA_workspace/data/valid/` 
logs are on `.logs/`
#   What to export

- in directory, where the python script runs, a file validation.log will be created. Please share with me.

# Configuration File — Specification & Guide

This document defines the **JSON** configuration that `config_manager.py` consumes to understand your inputs (what’s reference vs. modified; genome vs. features vs. reads) and to validate the values given by you.

---

## Quick overview

- **Name:** Configuration file is to be named `config.json`
- **Location:** Expected to be in directory together with your input files
- **Format:** JSON (UTF-8)
- **Purpose:** Tell where to find the reference/modified genomes, reads, and optional plasmids & features.
- **Who creates it:** User, Modification Applier
- **Who reads it:** `ConfigManager` (and maybe EFSA)

---

## Top-level keys

| Key | Type | Required | Description |
|---|---|:---:|---|
| `ref_genome_filename` | GenomeConfig | ✅ | Path to **reference genome** (FASTA or GenBank). |
| `mod_genome_filename` | GenomeConfig | ✅ | Path to **modified genome** (FASTA or GenBank). |
| `ref_plasmid_filename` | GenomeConfig | ❌ | Path to **reference plasmid** (FASTA or GenBank). |
| `mod_plasmid_filename` | GenomeConfig | ❌ | Path to **modified plasmid** (FASTA or GenBank). |
| `reads` | ReadConfig | ✅ | **At least one** read input. See “Reads structure” below. |
| `ref_feature_filename` | FeatureConfig | ❌ | Features for the **reference genome** (BED, GFF or GTF). |
| `mod_feature_filename` | FeatureConfig | ❌ | Features for the **modified genome** (BED, GFF or GTF). |
| `options` | dict | ❌ | Free-form extra knobs (e.g., `{"threads": 8}`) - **NOT IMPLEMENTED**. |


> ℹ️ **WARNING:** **Paths** must be relative to the **config directory** where your configuration file is.

---

## Supported file formats

- **Genomes & plasmids:** `.fa`, `.fasta`, `.fna` → FASTA; `.gb`, `.gbk`, `.genbank` → GenBank  
- **Features:** `.gff`, `.gff3` → GFF; `.gtf` → GTF; `.bed` → BED
- **Reads:** `.fastq`, `.fq` → FASTQ; `.bam` → BAM

---

## Input compression (coding)

All file-type or directory-type entries may be compressed and packed with following technologies

Allowed values:

- `"GZIP"` — gzip compression (`.gz`,`.gzip`)
- `"BZIP"` — bzip2 compression (`.bz2`,`.bzip2`)

---
## Genome Config
Each entry is **an object**: include `filename`.


Example:
```json
"ref_genome_filename": {"filename":"data/reference.gbk"}
```

---
## Read Config

#### a) file-type non-empty list
Each entry is **an object**: include `filename` and `ngs_type`.
```json
"reads": [
  { "filename": "samples/reads_R1.fastq.gz", "ngs_type": "illumina" },
  { "filename": "samples/reads_R2.fastq.gz", "ngs_type": "illumina" }
]
```

####  b) directory-type path

`reads` must be a **non-empty list** or **directory path**. Each entry is **an object**: include `directory` and `ngs_type`.
```json
"reads": [
  { "directory": "samples/", "ngs_type": "illumina" },
]
```

#### `ngs_type` allowed values

- `"illumina"` — short reads (SE or PE)
- `"ont"` — Oxford Nanopore
- `"pacbio"` — PacBio

> ℹ️ **WARNING: If the read file in format BAM is given, the ngs_type will be deprecated onto pacbio**

---
## Feature Config
Each entry is **an object**: include `filename`.


Example:
```json
"ref_feature_filename": {"filename":"data/reference.gbk"}
```

---
## Validation rules

1. `ref_genome_filename` and `mod_genome_filename` must exist.
2. `reads` must be non-empty.
3. All paths must exist or be accessible (but it will be checked).
4. File extensions determine format and coding type
