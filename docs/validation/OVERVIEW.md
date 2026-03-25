# Input Validation Overview

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
- If any genome exceeds `n_sequence_limit` or `type` is `"eukaryote"`, the file is still copied to `data/valid/` but `run_ref_x_mod` will be set to `false`
- Log file created in `./logs/validation_ID.log`
- Report file created in `./logs/report_ID.txt`

### `validated_params.json`

This file is produced by the validation step and consumed by Nextflow via `-params-file`. It overrides the defaults in `nextflow.config` and `nextflow_schema.json`. 

#### Pipeline switches

| Parameter            | Type    | Description                                                                                           |
| -------------------- | ------- | ----------------------------------------------------------------------------------------------------- |
| `run_ref_x_mod`      | boolean | `true` when both reference and modified genome validation succeeded and neither is fragmented; `false` when any genome exceeds `n_sequence_limit` or `type` is `"eukaryote"`. Gates all ref-vs-mod steps. |
| `run_truvari`        | boolean | Always `false` by default; can be overridden in `data/validation/validated_params.json`.                                                |
| `run_illumina`       | boolean | `true` when validated Illumina reads are present.                                                     |
| `run_nanopore`       | boolean | `true` when validated Nanopore (ONT) reads are present.                                               |
| `run_pacbio`         | boolean | `true` when validated PacBio reads are present.                                                       |
| `contig_file_size`   | integer | Number of contig files produced by inter-genome characterisation.                                     |
| `run_vcf_annotation` | boolean | `true` when a validated GFF feature file is available.                                                |

#### File paths (null when absent)

| Parameter             | Type   | Description                                          |
| --------------------- | ------ | ---------------------------------------------------- |
| `ref_fasta_validated` | string | Path to the validated reference genome FASTA.        |
| `mod_fasta_validated` | string | Path to the validated modified genome FASTA.         |
| `pacbio_fastq`        | string | Path to the first validated PacBio FASTQ file.       |
| `nanopore_fastq`      | string | Path to the first validated Nanopore FASTQ file.     |
| `gff`                 | string | Path to the validated reference GFF/GFF3 file.       |
