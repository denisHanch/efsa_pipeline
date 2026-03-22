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



> **Important!**
>
> When container is built please follow the steps to preprocess the data with a validation package.

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
4. Outputs validated files to `data/valid/`
5. Writes `data/valid/validated_params.json` with validated file paths and flags for Nextflow
6. Generates logs and reports in `logs/`

## Running Validation

```bash
python3 ./modules/validation/main.py ./data/inputs/config.json
```

## Related Documentation

- [Supported File Formats](FILE_FORMATS.md)
- [Configuration Guide](CONFIG_GUIDE.md)
- [Validation Settings](SETTINGS.md)
- [Performance Tips](PERFORMANCE.md)

## Output

After successful validation:

- Validated files are placed in `data/valid/`
- `data/valid/validated_params.json` is written with paths and flags for Nextflow (loaded automatically via `-params-file`)
- If any genome exceeds `n_sequence_limit` or `type` is `"eukaryote"`, the file is still copied to `data/valid/` but `run_syri` and `run_ref_x_mod` will be set to `false`
- Log file created in `./logs/validation_ID.log`
- Report file created in `./logs/report_ID.txt`

### `validated_params.json`

This file is produced by the validation step and consumed by Nextflow via `-params-file`. It overrides the defaults in `nextflow.config` and `nextflow_schema.json`.

#### Validated inputs (hidden params)

| Parameter             | Type    | Description                                                        |
| --------------------- | ------- | ------------------------------------------------------------------ |
| `ref_fasta_validated` | string  | Path to the validated reference genome FASTA.                      |
| `mod_fasta_validated` | string  | Path to the validated modified genome FASTA.                       |
| `ref_fasta_avail`     | boolean | `true` when a validated reference genome is available.             |
| `mod_fasta_avail`     | boolean | `true` when a validated modified genome is available.              |

#### Pipeline switches

| Parameter            | Type    | Description                                                                                           |
| -------------------- | ------- | ----------------------------------------------------------------------------------------------------- |
| `run_ref_x_mod`      | boolean | `true` when both reference and modified genome validation succeeded and neither is fragmented; `false` when any genome exceeds `n_sequence_limit` or `type` is `"eukaryote"`. Gates all ref-vs-mod steps. |
| `run_syri`           | boolean | `true` when the number of contig files is between 1 and `n_sequence_limit` (default `5`) of the modified genome; `false` when the assembly is too fragmented or the organism is eukaryotic. The threshold is controlled by `n_sequence_limit` in `config.json`. |
| `run_truvari`        | boolean | Always `false` by default; can be overridden manually.                                                |
| `run_illumina`       | boolean | `true` when validated Illumina reads are present.                                                     |
| `run_nanopore`       | boolean | `true` when validated Nanopore (ONT) reads are present.                                               |
| `run_pacbio`         | boolean | `true` when validated PacBio reads are present.                                                       |
| `contig_file_size`   | integer | Number of contig files produced by inter-genome characterisation.                                     |
| `run_vcf_annotation` | boolean | `true` when a validated GFF feature file is available.                                                |

#### File paths (null when absent)

| Parameter       | Type   | Description                                          |
| --------------- | ------ | ---------------------------------------------------- |
| `pacbio_fastq`  | string | Path to the first validated PacBio FASTQ file.       |
| `nanopore_fastq`| string | Path to the first validated Nanopore FASTQ file.     |
| `gff`           | string | Path to the validated reference GFF/GFF3 file.       |
