# Input Validation Overview


## Input Scenarios and Preprocessing Logic

The validation module not only verifies input formats, but also determines how genome assemblies are interpreted and preprocessed before entering the pipeline.

Depending on the structure of `ref.fa` and `mod.fa`, different strategies are applied for:
- chromosome and plasmid separation
- contig handling
- usage of minimap2 for sequence mapping
- preparation of files in `data/valid/`

The following table summarizes all supported scenarios:

| #     | Scenario                                                  | Mode (`config.json`)       | Input Structure                                                                                | Plasmids Handling                                                                                                                                               | `run_ref_x_mod` | **minimap2 Mapping**                                              | mod.fa Processing                                                                                          | Modules Run          | Output in `data/valid/`                                                      |
| ----- | --------------------------------------------------------- | -------------------------- | ---------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------- | ----------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------- | -------------------- | ---------------------------------------------------------------------------- |
| **1** | Single contig + plasmids                    | `prokaryote`               | **ref.fa:** 1 sequence (+ optional plasmids) <br> **mod.fa:** 1 sequence (+ optional plasmids) |  In **ref.fa**: <br> - Longest sequence → chromosome <br> - Remaining → plasmids `*ref_plasmid.fasta` <br><br> In **mod.fa**: plasmids = sequences **not mapped** to reference    | True            | **Used** to identify unmapped regions (plasmids in mod.fa)        | Reduced to **1 contig** (chromosome only)                                                                  | All modules          | `ref.fa`, `mod.fa` (1 contig) <br> `*_contig_0.fasta` <br> `*_plasmid.fasta` |
| **2** | Fragmented assembly (below limit)            | `prokaryote`               | **ref.fa:** 1 sequence <br> **mod.fa:** multiple sequences (≤ limit)                           | In reference: <br> - In **ref.fa**: <br> - Longest sequence → chromosome <br> - Remaining → plasmids `*ref_plasmid.fasta` <br><br> In **mod.fa**: <br> - Unmapped sequences → plasmids <br> - Mapped sequences → contigs | True            | **Used** to split mod.fa into mapped contigs vs unmapped plasmids | - Split into individual contigs (`*_contig.fasta`) <br> - `mod.fa` becomes multifasta **without plasmids** | All modules          | `ref.fa`, `mod.fa` + contig set <br> `*_contig.fasta` <br> `*_plasmid.fasta` |
| **3** | Fragmented assembly (above limit)             | `prokaryote`               | **ref.fa:** 1 sequence <br> **mod.fa:** multiple sequences (> limit)                           | In reference: <br> - In **ref.fa**: <br> - Longest sequence → chromosome <br> - Remaining → plasmids `*ref_plasmid.fasta` <br><br> In **mod.fa**: no plasmid detection                                       | False           | **Not used**                                                      | No processing (mod.fa copied as-is)                                                                        | Mapping-only modules | `ref.fa`, `mod.fa` (copied) <br> `*_plasmid.fasta`                           |
| **4** | Multiple sequences in reference  | `prokaryote` / `eukaryote` | **ref.fa:** multiple sequences (non-plasmid) <br> **mod.fa:** one or more sequences            | No plasmids considered                                                                                                                                          | False           | **Not used**                                                      | No processing (files copied as-is)                                                                         | Mapping-only modules | `ref.fa`, `mod.fa` (copied)                                                  |
| **5** | Fragmented reference + `force_defragment_ref` (**unsupported workaround**) | `prokaryote` / `eukaryote` | **ref.fa:** multiple sequences (fragmented, above limit) → merged to 1 before validation | In **ref.fa**: all contigs and plasmids merged into a single chromosome (`*_defragmented.fasta`) | True | Not used | Same for scenarios above | GFF annotation disabled | `ref.fa` (joined), `*_defragmented_join_order.tsv`, `mod.fa` depends on scenario |

> **Important!**
> 
> By default the limit (`n_sequence_limit`) mentioned in the table above is set to 5 for both reference and assembly fasta files. Please adjust the n_sequence_limit in `data/inputs/config.json`.
>


> **Important!**
>
> When the container is built please follow the steps to preprocess the data with a validation package.

The input validation module preprocesses and verifies all input data to ensure it meets the required format and structure before the Nextflow pipeline is executed.

## Purpose

The validation module ensures that:

- All input files are in the correct format
- Files are properly structured and can be parsed
- Data meets quality standards
- Files are converted to standardized formats for pipeline processing

## How It Works

The validation process:

1. Reads the configuration file from `data/inputs/config.json`
2. Validates each input file according to its type (genome, reads, features)
3. Converts files to standardized formats
4. Outputs validated files to a timestamped run directory `data/valid/run_YYYYMMDD_HHMMSS/`
5. Writes `data/valid/validated_params.json` (at the top level, not inside the run directory) with validated file paths and flags for Nextflow
6. Generates a log and report inside the run directory

## Running Validation

Use the `validate` wrapper script (preferred):

```bash
validate                                      # default config path
validate --config path/to/config.json         # custom config

# Global options can be set via CLI flags (config.json takes priority if the same
# option is also set there):
validate --threads 8
validate --validation-level strict
validate --logging-level DEBUG
validate --type eukaryote
validate --force-defragment-ref               # unsupported workaround — at your own responsibility
                                              # has no effect if force_defragment_ref is set in config.json
```

Priority of configurations: config.json > cli_options > defaults

## Related Documentation

- [Supported File Formats](FILE_FORMATS.md)
- [Configuration Guide](CONFIG_GUIDE.md)
- [Validation Settings](SETTINGS.md)
- [Performance Tips](PERFORMANCE.md)

## Output

After successful validation:

- Validated files are placed in `data/valid/run_YYYYMMDD_HHMMSS/` (a new timestamped directory per run; previous runs are preserved)
- `data/valid/validated_params.json` is written at the top level of `data/valid/` so Nextflow can always find it at a fixed path
- If any genome exceeds `n_sequence_limit` or `type` is `"eukaryote"`, the file is still copied but `run_ref_x_mod` will be set to `false`
- Log and report are written inside the run directory

### `validated_params.json`

This file is produced by the validation step and consumed by Nextflow via `-params-file`. It overrides the defaults in `nextflow.config` and `nextflow_schema.json`. 

#### Pipeline switches

| Parameter               | Type    | Description                                                                                           |
| ----------------------- | ------- | ----------------------------------------------------------------------------------------------------- |
| `run_ref_x_mod`         | boolean | `true` when both reference and modified genome validation succeeded and neither is fragmented; `false` when any genome exceeds `n_sequence_limit` or `type` is `"eukaryote"`. Gates all ref-vs-mod steps. |
| `run_truvari`           | boolean | Always `false` by default; can be overridden manually in `data/valid/validated_params.json`.          |
| `run_illumina`          | boolean | `true` when validated Illumina FASTQ reads are present.                                               |
| `run_nanopore`          | boolean | `true` when validated Nanopore (ONT) reads are present (FASTQ or BAM).                               |
| `run_pacbio`            | boolean | `true` when validated PacBio reads are present (FASTQ or BAM).                                       |
| `contig_file_size`      | integer | Number of contig files produced by inter-genome characterisation.                                     |
| `run_vcf_annotation`    | boolean | `true` when a validated GFF feature file is available.                                                |
| `validation_timestamp`  | string  | Timestamp of the validation run (`YYYYMMDD_HHMMSS`).                                                  |

#### File paths (null or empty list when absent)

| Parameter             | Type          | Description                                                    |
| --------------------- | ------------- | -------------------------------------------------------------- |
| `ref_fasta_validated` | string        | Path to the validated reference genome FASTA.                  |
| `mod_fasta_validated` | string        | Path to the validated modified genome FASTA.                   |
| `ref_plasmid_fasta`   | string        | Path to the validated reference plasmid FASTA (if present).    |
| `mod_plasmid_fasta`   | string        | Path to the validated modified plasmid FASTA (if present).     |
| `gff`                 | string        | Path to the validated reference GFF/GFF3 file.                 |
| `illumina_fastqs`     | string array  | Paths to all validated Illumina FASTQ files.                   |
| `ont_fastqs`          | string array  | Paths to all validated Nanopore FASTQ files.                   |
| `ont_bams`            | string array  | Paths to all validated Nanopore BAM files (copied as-is).      |
| `pacbio_fastqs`       | string array  | Paths to all validated PacBio FASTQ files.                     |
| `pacbio_bams`         | string array  | Paths to all validated PacBio BAM files (copied as-is).        |
| `contig_files`        | string array  | Paths to contig FASTA files from inter-genome characterisation. |

> **BAM files:** When a PacBio or ONT input is provided as BAM, it is copied to the run directory without conversion. Its path appears in `pacbio_bams` / `ont_bams` — separate from the FASTQ lists — so Nextflow can handle the two formats independently. Inter-read (R1/R2 pairing) validation is skipped when all reads for a sample are BAM.
