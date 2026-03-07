# Configuration File — Specification & Guide

This document describes every field in `config.json` used by the validation pipeline.

---

## Minimal example

```json
{
  "ref_genome_filename": {"filename": "reference.fasta"},
  "reads": [
    {"filename": "reads_R1.fastq.gz", "ngs_type": "illumina"},
    {"filename": "reads_R2.fastq.gz", "ngs_type": "illumina"}
  ]
}
```

---

## Top-level fields

| Field | Required | Description |
|-------|----------|-------------|
| `ref_genome_filename` | yes | Reference genome (FASTA / GenBank) |
| `reads` | yes | List of read files (FASTQ / BAM) |
| `mod_genome_filename` | no | Modified/sample genome |
| `ref_plasmid_filename` | no | Reference plasmid |
| `mod_plasmid_filename` | no | Modified plasmid |
| `ref_feature_filename` | no | Reference annotation (GFF / GTF / BED) |
| `mod_feature_filename` | no | Modified annotation |
| `options` | no | Global validation options (see below) |

---

## Genome / plasmid / feature files

Each can be a plain string or a dict with per-file options:

```json
"ref_genome_filename": "reference.fasta"

"ref_genome_filename": {
  "filename": "reference.fasta",
  "validation_level": "trust"
}
```

---

## `reads` array

Each entry must contain **either** `filename` or `directory`, plus any options.

### `filename` entry

```json
{
  "filename": "reads_R1.fastq.gz",
  "ngs_type": "illumina"
}
```

### `directory` entry — validates every file inside the folder

```json
{
  "directory": "raw_reads/",
  "ngs_type": "ont"
}
```

### Read entry fields

| Field | Required | Values | Description |
|-------|----------|--------|-------------|
| `filename` | one of these | string | Path relative to the config file |
| `directory` | one of these | string | Directory of read files |
| `ngs_type` | see below | `illumina` / `ont` / `pacbio` | Sequencing platform |
| `validation_level` | no | `strict` / `trust` / `minimal` | Per-file override |
| `threads` | no | integer | Per-file thread count override |

---

## `ngs_type` — platform specification and autodetection

`ngs_type` tells the pipeline how to handle the reads (naming conventions, expected
header format, output directory organisation when `outdir_by_ngs_type=True`).

### Explicit (always works, all validation modes)

```json
{"filename": "reads.fastq.gz", "ngs_type": "illumina"}
{"filename": "reads.fastq.gz", "ngs_type": "ont"}
{"filename": "reads.fastq.gz", "ngs_type": "pacbio"}
```

### Autodetection (strict mode only)

When `ngs_type` is **omitted**, the pipeline can detect it automatically — but only in
**strict** validation mode, because strict mode reads the entire file anyway.

```json
{
  "options": {"validation_level": "strict"},
  "reads": [
    {"filename": "unknown_reads.fastq.gz"}
  ]
}
```

**How detection works:**

The validator samples the first 20 parsed sequence headers and looks for
platform-specific signatures, using a majority vote:

| Platform | Header signatures |
|----------|------------------|
| ONT | `runid=`, `ch=`, `start_time=`, `flow_cell_id=` |
| PacBio | `/ccs`, `movie=`, `zmw` |
| Illumina | ` 1:N:`, ` 2:N:`, ` 1:Y:`, ` 2:Y:` |

If a configured `ngs_type` disagrees with what is detected in the file content, a
**warning** is logged and the content-detected value is used (content is authoritative
in strict mode).

If no signature is found and `ngs_type` was not provided, validation **fails** with an
error asking the user to set `ngs_type` explicitly.

**Trust and minimal modes require `ngs_type` to be set explicitly** — they do not read
file content and therefore cannot detect the platform. Omitting `ngs_type` with those
modes raises a `ConfigurationError` at load time.

```json
{
  "options": {"validation_level": "trust"},
  "reads": [
    {"filename": "reads.fastq.gz"}
  ]
}
```

```
ConfigurationError: 'ngs_type' is required when validation_level is 'trust'.
Must be one of: illumina, ont, pacbio.
Omit 'ngs_type' only with validation_level='strict' for autodetection.
```

### Summary

| Scenario | Result |
|----------|--------|
| `ngs_type` provided explicitly | Used as-is (validated to be a known value) |
| `ngs_type` omitted + strict mode | Autodetected from FASTQ headers |
| `ngs_type` omitted + trust/minimal | `ConfigurationError` at load time |
| `ngs_type` provided but differs from content (strict) | Warning logged, content value used |
| No signatures found + `ngs_type` omitted (strict) | `ReadValidationError` during validation |

---

## `options` — global settings

```json
{
  "options": {
    "validation_level": "strict",
    "threads": 8,
    "logging_level": "INFO"
  }
}
```

| Option | Default | Values | Description |
|--------|---------|--------|-------------|
| `validation_level` | `strict` | `strict` / `trust` / `minimal` | How thoroughly files are validated |
| `threads` | auto | integer or `null` | Parallel processing threads (`null` = auto-detect) |
| `logging_level` | `INFO` | `DEBUG` / `INFO` / `WARNING` / `ERROR` | Console log verbosity |

### Validation levels explained

| Level | File read? | Converts formats? | Use when |
|-------|-----------|-------------------|----------|
| `strict` | Fully | Yes | First-time submission; maximum accuracy |
| `trust` | First 10 reads | Yes | File was already validated before |
| `minimal` | No | No | File is already in the exact required format |

Per-file overrides are supported for `validation_level` and `threads`:

```json
{
  "options": {"validation_level": "strict"},
  "reads": [
    {"filename": "already_clean.fastq.gz", "ngs_type": "illumina", "validation_level": "minimal"},
    {"filename": "new_reads.fastq.gz"}
  ]
}
```

---

## Full example

```json
{
  "ref_genome_filename": {"filename": "GCF_ref.fasta"},
  "mod_genome_filename": {"filename": "sample_assembly.fasta"},
  "ref_plasmid_filename": {"filename": "plasmid_ref.fasta"},
  "ref_feature_filename": {"filename": "annotation_ref.gff"},
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
      "filename": "nanopore_reads.fastq.gz",
      "ngs_type": "ont"
    },
    {
      "filename": "unknown_platform_reads.fastq.gz"
    }
  ],
  "options": {
    "validation_level": "strict",
    "threads": 8,
    "logging_level": "INFO"
  }
}
```

In this example the last read entry has no `ngs_type`. Because `validation_level` is
`strict`, the pipeline will detect the platform from the file's sequence headers at
validation time.
