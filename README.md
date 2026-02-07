
# EFSA Pipeline

## Table of Contents
- [Quick Start](#quick-start)
- [Docker Container](#docker-container)
- [Input Validation](#input-validation)
- [Nextflow](#nextflow)
   - [Running the Pipeline](#running-the-pipeline)
   - [🔄 Pipeline Runtime Messages & Mapping Summary](#-pipeline-runtime-messages--mapping-summary)
      - [Runtime Status Messages](#runtime-status-messages)
      - [📊 Unmapped Reads Statistics](#-unmapped-reads-statistics)
      - [✅ Pipeline Execution Summary](#-pipeline-execution-summary)
      - [Removal of Nextflow Work Directory](#removal-of-nextflow-work-directory)
   - [📁 `data/valid` Directory Structure](#-datavalid-directory-structure)
   - [📁 `data/outputs` Directory Structure](#-dataoutputs-directory-structure)
      - [`fasta_ref_mod/`](#fasta_ref_mod)
      - [`illumina/`](#illumina)
      - [`pacbio/ and ont/`](#pacbio-and-ont)
      - [`truvari/`](#truvari)
      - [`unmapped_stats/`](#unmapped_stats)
  - [Graphical Representation of the Pipeline](#graphical-representation-of-the-pipeline)
  - [Generation of per structural variation (SV) type CSV tables](#generation-of-per-structral-variation-sv-type-csv-tables)
      

# Quick Start

1. **Download the repository**:

   ```bash
   git clone https://github.com/denisHanch/efsa_pipeline.git
   ```

> **Important!**
> 
> Make sure that the data for pipelines are in the folder `data/inputs`.
>

2. **Start the Docker container**:

   ```bash
   ./run_container.sh
   ```

> **Important!**
> 
> Create configuration file `config.json` on `data/inputs`.
> 

3. **Running QC** on the input data and **processing data for the Nextflow pipeline** to `data/valid` folder:

   ```bash
   python3 ./modules/validation/main.py ./data/inputs/config.json
   ```

4.  Start the pipeline with a command:

   ```bash
   nextflow run main.nf \
     -process.containerOptions "-u $(id -u):$(id -g)" \
     --max_cpu $(nproc) \
     -resume
   ```


# Docker Container

## Docker Setup for Users

This guide shows you how to run the EFSA Pipeline in a Docker container with access to input/output folders.

## Prerequisites

- Docker installed on your system
- Git (if cloning the repository)

## Quick Start

### Option 1: Using the Run Script (Recommended)

1. Make sure the script is executable:
   ```bash
   chmod +x run_container.sh
   ```

2. Run the container:
   ```bash
   ./run_container.sh
   ```

3. You'll be dropped into the container shell where you can run CLI commands
4. Type `exit` when done to return to your host system


### Option 2: Manual Docker Commands

1. Build the image:
   ```bash
   docker build -t efsa-pipeline .
   ```

2. Run interactively:
   ```bash
   docker run --privileged -d --rm \
    --network=host \
    -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -v /usr/share/ca-certificates:/usr/share/ca-certificates:ro \
    --name efsa-pipeline-container \
    -w $(pwd) \
    -v "$(pwd)/data/inputs:/EFSA_workspace/data/inputs" \
    efsa-pipeline

   docker exec -it efsa-pipeline-container /bin/sh
   ```

# Input Validation

> **Important!**
> 
> When container is build please follow the steps to preprocess the data with a validation package.
>

The input validation module preprocesses and verifies all input data to ensure it meets the required format and structure before the Nextflow pipeline is executed.

## Supported File Formats


### Genome Files
| Supported Input Formats | Final Output Format |
|------------------------|---------------------|
| FASTA: `.fasta`, `.fa`, `.fna` | `.fasta` |
| GenBank: `.gb`, `.gbk`, `.genbank` | `.fasta` |
| Compression: `.gz`, `.bz2`, `.gzip`, `.bzip2` | Uncompressed |

### Read Files
| Supported Input Formats | Final Output Format |
|------------------------|---------------------|
| FASTQ: `.fastq`, `.fq` | `.fastq` |
| BAM: `.bam` (limited support) | `.bam` |
| Compression: `.gz`, `.bz2`, `.gzip`, `.bzip2` | `.gz` |

### Feature Files
| Supported Input Formats | Final Output Format |
|------------------------|---------------------|
| GFF: `.gff`, `.gff3`, `.gtf` | `.gff3` |
| BED: `.bed` | `.gff3` |
| Compression: `.gz`, `.bz2`, `.gzip`, `.bzip2` | Uncompressed |


## Configuration

To run validation, configuration file is expected in `json` format. This config specify which file corresponds to reference genome, assembled modified genome and reads. The detailed documentation of the configuration is described [here](/docs/validation/CONFIG_GUIDE.md).

Full configuration template is prepared in `data/inputs/config.json`.

```json
{
  "ref_genome_filename": {
    "filename": "TBA",
    "validation_level": "strict/trust/minimal",
    "threads": 8
  },
  "mod_genome_filename": {
    "filename": "TBA",
    "validation_level": "strict/trust/minimal",
    "threads": 8
  },
  "ref_plasmid_filename": {
    "filename": "TBA",
    "validation_level": "strict/trust/minimal",
    "threads": 8
    },
  "mod_plasmid_filename": {
    "filename": "TBA",
    "validation_level": "strict/trust/minimal",
    "threads": 8
    },
  "reads": [
    {
      "filename": "TBA",
      "ngs_type": "illumina/pacbio/ont",
      "validation_level": "strict/trust/minimal",
      "threads": 8
    },
    {
      "directory": "TBA",
      "ngs_type": "illumina/pacbio/ont",
      "validation_level": "strict/trust/minimal",
      "threads": 8
    }
  ],
  "ref_feature_filename": {
    "filename": "TBA",
    "validation_level": "strict/trust/minimal",
    "threads": 8
  },
  "mod_feature_filename": {
    "filename": "TBA",
    "validation_level": "strict/trust/minimal",
    "threads": 8
    },
  "options": {
    "threads": 8,
    "validation_level": "strict/trust/minimal",
    "logging_level":"DEBUG/INFO/WARNING/ERROR"
  }
}
```


## Validation Levels

The package supports three validation levels to balance thoroughness and performance:

### Comparison Table

| Level | Parsing | Validation | Edits | Output | Speed | Use Case |
|-------|---------|------------|-------|--------|-------|----------|
| **strict** (default) | All data | All data | All applied | BioPython write | Slowest | Structure validation sequence by sequence, statistics gathering |
| **trust** | All data (genome)<br>First record only (reads) | First sequence only | All applied (genome)<br>None (reads) | BioPython write (genome)<br>File copy (reads) | Fast | Trust data, adapt file coding,name and location |
| **minimal** | None | None | None | File copy | Fastest | Rename and move files to meet the requirements |

## Logging

Every validation run create a new **log** file on `./logs/validation_ID.log`. With option `logging_level` in configuration file you may specify the amount of information logged.

**Supported levels:**
- `"DEBUG"`: Most verbose (all messages including debug info)
- `"INFO"`: Standard output (validation progress, file operations)
- `"WARNING"`: Warnings and errors only
- `"ERROR"`: Errors only

Every validation run create a new **report** file on `./logs/report_ID.txt`.
A validation report contains three main sections:


**1. Summary Section:**
- Overall validation status (PASSED/FAILED)
- File counts by type (genomes, reads, features)
- Inter-file validation counts (passed/failed)
- Total execution time

**2. File Validation Results**
For each validated file:
- **Input information**: Original filename, validator type, settings used
- **Output information**: Output file path, validation level
- **Metadata**: Type-specific information and statistics
  - **Genomes**: Sequence count, IDs, lengths, GC content
  - **Reads**: Read count, NGS type, pairing info, N50, lengths, total bases
  - **Features**: Feature count, types, sequence coverage
- **Timing**: Elapsed time for validation

**3. Inter-file Validation Results**
- Validation type (genome×genome, read×read, feature×genome)
- Status (PASSED/FAILED)
- Errors and warnings
- Metadata about cross-file checks

## Performance Tips

**To adapt validation performance:**

1. Install parallel compression tools: `sudo apt-get install pigz pbzip2`
2. Use `validation_level='trust'` for pre-validated data (10-15x faster)
3. Use `threads` option in config (default 8) for record-level parallelization (3-7x faster in strict mode)


# Nextflow

## Graphical Representation of the Pipeline

```mermaid
flowchart TD
    %% Inputs
    subgraph Inputs
        illumina["Illumina Short Reads"]
        pacbio["PacBio Long Reads"]
        ont["ONT Long Reads"]
        ref_fasta["Reference FASTA"]
        mod_fasta["Modified FASTA"]
    end
    style Inputs fill:#E6F4EA,stroke:#2E7D32,stroke-width:2px

    %% Short-read pipeline
    subgraph ShortReadPipeline["Short-Read Pipeline"]
        illumina --> trimgalore["TrimGalore / QC"]
        trimgalore --> bwa_index["BWA Index"]
        bwa_index --> bwa_map["BWA Mapping"]
        bwa_map --> samtools_sort["Samtools Sort & Stats"]
        samtools_sort --> picard["Picard / MultiQC"]
        samtools_sort --> sr_vcf["VCF Output (Short-Read)"]
    end
    style ShortReadPipeline fill:#D0F0C0,stroke:#388E3C,stroke-width:2px

    %% Long-read PacBio pipeline
    subgraph LongReadPacBio["Long-Read PacBio Pipeline"]
        pacbio --> minimap2_pb["Minimap2 Mapping"]
        minimap2_pb --> samtools_sort_pb["Samtools Sort & Stats"]
        samtools_sort_pb --> pb_sv["SV Calling (CuteSV / Sniffles / Debreak)"]
        pb_sv --> pb_vcf["VCF Output (PacBio)"]
    end
    style LongReadPacBio fill:#D0F0C0,stroke:#388E3C,stroke-width:2px

    %% Long-read ONT pipeline
    subgraph LongReadONT["Long-Read ONT Pipeline"]
        ont --> minimap2_ont["Minimap2 Mapping"]
        minimap2_ont --> samtools_sort_ont["Samtools Sort & Stats"]
        samtools_sort_ont --> ont_sv["SV Calling (CuteSV / Sniffles / Debreak)"]
        ont_sv --> ont_vcf["VCF Output (ONT)"]
    end
    style LongReadONT fill:#D0F0C0,stroke:#388E3C,stroke-width:2px

    %% Reference vs Modified pipeline
    subgraph RefVsMod["Reference vs Modified Pipeline"]
        ref_fasta --> nucmer["NUCmer Alignment"]
        mod_fasta --> nucmer
        nucmer --> delta_filter["Delta Filter"]
        delta_filter --> show_coords["Show Coords"]
        show_coords --> syri_vcf["VCF Output (ref_x_modsyri)"]
    end
    style RefVsMod fill:#D0F0C0,stroke:#388E3C,stroke-width:2px

    %% Truvari comparison
    subgraph Truvari["Truvari Comparison"]
        sr_vcf --> truvari["Compare SVs (Truvari)"]
        pb_vcf --> truvari
        ont_vcf --> truvari
        syri_vcf --> truvari
        truvari --> final_report["Truvari Reports / Summary"]
    end
    style Truvari fill:#D0F0C0,stroke:#2E7D32,stroke-width:2px
```

## Running the Pipeline

The main pipeline (`main.nf`) executes **all three workflows** in sequence.
Each workflow can also be executed individually if required.

### Running Main Workflow

This executes **short-read processing**, **long-read processing**, and **reference vs modified genome comparison** pipelines:

```bash
nextflow run main.nf \
  -process.containerOptions "-u $(id -u):$(id -g)" \
  --max_cpu $(nproc) \
  -resume
```

### Available Nextflow Option

| Option        | Description |
|---------------|-------------|
| `-resume`     | Resumes a pipeline run from the point where it was interrupted or previously failed. |
| `-with-report` | Generates a visual HTML report of the workflow execution, including task durations, resource usage, and statuses. |
| `-with-timeline` | Produces a timeline visualization showing when each process in the pipeline started and finished. |
| `-with-dag`     | Generates a directed acyclic graph (DAG) showing task dependencies within the workflow. |


### Available Options

| Option         | Description                                  | Default                |
| -------------- | -------------------------------------------- | ---------------------- |
| `--in_dir`     | Input directory                              | `data/valid`           |
| `--out_dir`    | Output directory                             | `data/outputs`         |
| `--registry`   | Docker/Singularity container registry        | `ghcr.io/ecomolegmo`   |
| `--max_cpu`    | Maximum CPUs per process                     | `1`                    |
| `--log`        | Enable logging                               | `true`                 |
| `--clean_work` | Remove work directory after successful run   | `true`                 |
| `--help`       | Display help message                         | –                      |

---

### Running Individual Pipelines

You can also run each of the three sub-pipelines independently.

#### Short-read Processing

For Illumina short-read data:

```bash
nextflow run workflows/short_read.nf \
  -process.containerOptions "-u $(id -u):$(id -g)" \
  --max_cpu $(nproc) \
  -resume
```

#### Long-read Processing

For Oxford Nanopore / PacBio long reads:

```bash
nextflow run workflows/long_read.nf \
  -process.containerOptions "-u $(id -u):$(id -g)" \
  --max_cpu $(nproc) \
  -resume
```

#### Reference vs Modified Genome Comparison

For comparing reference and modified FASTA assemblies:

```bash
nextflow run workflows/fasta_ref_x_mod.nf \
  -process.containerOptions "-u $(id -u):$(id -g)" \
  --max_cpu $(nproc) \
  -resume
```

 
## Generation of per structural variation (SV) type CSV tables

These utilities convert SV VCFs into compact TSV summaries and then merge available summaries into per-SV-type CSV tables.

```mermaid
flowchart LR
  subgraph TSVs
    A[vcf_to_table / vcf_to_table_long]
    B[create_empty_tbl]
  end
  A --> C[params.out_dir/tables/tsv]
  B --> C
  C --> D[restructure_sv_tbl]
  D --> E[params.out_dir/tables/csv_per_sv_sumary]
```

### Key points

- Variants are extracted into a table format with processes `vcf_to_table` and `vcf_to_table_long`
- if one of the pipelines was not running (shourt/long/assembly) an empty tsv file is generated with a process create_empty_tbl
- `restructure_sv_tbl` process: the merge step accepts any subset of (assembly, long_ont, long_pacbio, short) and ignores missing files.
- Long reads are handled as two separate sources: `long_ont` and `long_pacbio`. Output CSVs keep these in distinct `long_ont_*` and `long_pacbio_*` columns.
- Outputs overview:

```
data/outputs/tables/
├── csv_per_sv_summary
│   ├── Deletions.csv
│   ├── Insertions.csv
│   ├── Inversions.csv
│   ├── Replacements.csv
│   └── Translocations.csv
└── tsv
    ├── assembly_sv_summary.tsv
    ├── short_sv_summary.tsv
    ├── mab-pb_sv_summary.tsv
    └── map-ont_sv_summary.tsv
```


Example (local):
```bash
python3 modules/utils/create_sv_output_xlsx.py --asm assembly_sv_summary.tsv \
  --long_ont sample1_ont_sv_summary.tsv --long_pacbio sample1_pacbio_sv_summary.tsv \
  --short sample1_short_sv_summary.tsv --out csv_per_sv_sumary
```

### Explanation of csv_per_sv_summary CSV columns

The final table of each CSV file has multiple columns.
Each row represents one structural variant (SV) event, with coordinates and evidence aggregated across assembly-based, long-read, and short-read pipelines.

**Column prefixes**

- **asm_** — values reported by the assembly-based SV pipeline
- **long_** — values reported by long-read pipelines (ONT or PacBio)
- **short_** — values reported by short-read pipelines


**Common event-level and pipeline-derived columns**

These columns describe the SV event independently of any specific pipeline:

| Column name         | Description                                                                                                                                |
| ------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| **event_id**        | Unique identifier of the structural variant (SV) event.                                                                                    |
| **chrom**           | Chromosome where the SV is located (VCF `CHROM`).                                                                                          |
| **std_svtype**      | Standardized SV type harmonized across pipelines (e.g. DEL, DUP, INS, INV).                                                                |
| **event_start**     | Representative start coordinate of the SV (VCF `POS`), selected from available pipeline calls.                                             |
| **event_end**       | Representative end coordinate of the SV (VCF `END`).                                                                                       |
| **event_length_bp** | Length of the SV in base pairs, calculated as `event_end − event_start`.                                                                   |
| **support_score**   | Integrated support score summarizing evidence across all pipelines, taking into account SV agreement, confidence scores, and read support. |

| **percentage_overlap**   | Percentage of genomic overlap between two structural variants, calculated as the length of the intersection of their genomic intervals divided by the length of the smaller variant, multiplied by 100. |

**Additional pipeline-specific columns**

These columns are present only for certain pipelines:

| Column name                          | Description                                                                                                                               |
| ------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------- |
| **long_(ont|pacbio)_supporting_reads**        | Number of Oxford Nanopore / Pacbio reads supporting the structural variant (VCF `FORMAT` field `DR`) (SV).                                                                   |
| **long_(ont|pacbio)_supporting_methods**     | Number of variant callers (max 3 cute_cv, debreak, sniffles) reads supporting the structural variant (SV).                                                                            |
| **short_reads_copy_number_estimate** | Estimated copy number derived from short-read depth information (VCF `FORMAT` field `RDCN`). |


## 🔄 Pipeline Runtime Messages & Mapping Summary

During execution, the pipeline prints progress messages indicating which workflow is currently running and what type of reads are being processed.

### Runtime Status Messages

When the pipeline is running, you will see real-time messages like:

```text
ℹ️  Running pipeline: processing long-pacbio reads → mapping to the reference & modified fasta.

ℹ️  Running pipeline: processing long-ont reads → mapping to the reference & modified fasta.

ℹ️  Running pipeline: processing short reads → mapping to the reference & modified fasta.

ℹ️  Truvari: performing 3 comparisons.
```

These messages help track the execution order and confirm that all three pipelines are being executed as expected.

---

### 📊 Unmapped Reads Statistics

After mapping, the pipeline reports the number and percentage of **unmapped reads** for each analysis.
This is useful for assessing mapping efficiency and data quality.

#### Example Output

```text
📊 short-mod mapping:
    Unmapped reads: 19,880 (2.06%)
    Total input reads: 963,427

📊 short-ref mapping:
    Unmapped reads: 25,360 (2.63%)
    Total input reads: 963,362

📊 short-ref-plasmid mapping against plasmid:
    Unmapped reads: 4,677 (0.49%)
    Total input reads: 963,362

📊 ont/long-ref mapping:
    Unmapped reads: 52,745 (3.04%)
    Total input reads: 1,732,734

📊 ont/long-mod mapping:
    Unmapped reads: 48,145 (2.64%)
    Total input reads: 1,825,876

📊 pacbio/long-mod mapping:
    Unmapped reads: 41,596 (2.28%)
    Total input reads: 1,826,736

📊 pacbio/long-ref mapping:
    Unmapped reads: 47,472 (2.74%)
    Total input reads: 1,733,973
```

#### Interpretation

**Unmapped reads** represent sequences that did not align to the provided reference or modified FASTA files.

A low percentage of unmapped reads indicates:

  * High mapping quality
  * Good reference/assembly quality
  * Low contamination or sequencing errors

If the percentage of unmapped reads is unusually high, this may indicate:

* Poor read quality
* Inadequate or incomplete reference
* Contamination
* Incorrect input file selection


### ✅ Pipeline Execution Summary

The Nextflow pipelines ran successfully and produced the expected outputs. Each step completed without errors:

```text
✅ The reference to modified fasta comparision processing pipeline completed successfully.

✅ The long-ref processing pipeline completed successfully.

✅ The short-ref processing pipeline completed successfully.

✅ Truvari: the comparison of vcf files finished successfully.

✅ Execution of main.nf processing pipeline completed successfully.
```


### Removal of Nextflow Work Directory

When the pipeline is executed with the parameter:

```text
params.clean_work = true
```

Nextflow automatically deletes the temporary `work/` directory upon successful completion and logs a message to confirm this.

```text
ℹ️ Nextflow `work/` directory was removed.
```

**Notes:**

* The `work/` directory contains intermediate files and temporary outputs generated during pipeline execution.
* Removing it saves disk space while retaining all final results in the `out_dir`.
* If you want to keep intermediate files for debugging or inspection, set: `params.clean_work = false` in nextflow.config or use `--clean_work false` when running the pipeline.


## 📁 `data/valid` Directory Structure

This directory contains all input data used by the Nextflow pipeline.
```
data/valid/
├── assembled_genome.fasta
├── reference_genome.fasta
├── ref_plasmid.fa             # Reference plasmid sequences (if used)
├── mod_plasmid.fa             # Modified/assembled plasmid sequences (if used)
├── ref_feature.gff            # Genome annotation file GTF/GFF (if used)
│
├── illumina/                  
│   ├── SampleName_1.fastq.gz  
│   ├── SampleName_2.fastq.gz  
│
├── ont/                       
│   └── SampleName.fastq.gz
│
└── pacbio/                   
    └── SampleName.fastq.gz

```

| File / Folder            | Description                                           |
| ------------------------ | ----------------------------------------------------- |
| `reference_genome.fasta` | The primary reference genome sequence.                |
| `assembled_genome.fasta` | Assembled or modified genome for comparison/analysis. |
| `ref_plasmid.fa`         | Reference plasmid sequences.                          |
| `mod_plasmid.fa`         | Modified or assembled plasmid sequences.              |
| `ref_feature.gff`        | GFF feature file for annotations.                     |
| `illumina/`              | Paired-end Illumina short reads.                      |
| `ont/`                   | Oxford Nanopore long reads.                           |
| `pacbio/`                | PacBio long reads.                                    |

## 📁 `data/outputs` Directory Structure

After successful pipeline execution, the outputs are organized as follows:

```
data/outputs
├── fasta_ref_mod       → Results from reference vs modified FASTA comparison
├── illumina            → Short-read (Illumina) mapping results
├── logs                → Pipeline logs and Nextflow reports
├── ont                 → Long-read (Oxford Nanopore) mapping results
├── pacbio              → Long-read (PacBio) mapping results
├── tables              → Per-SV csv tables
├── truvari             → Variant comparison results from Truvari
└── unmapped_stats      → Summary statistics of unmapped reads for each workflow

```

A detailed description of the contents of each subfolder is provided below.

### `fasta_ref_mod/`

The pipeline overview is outlined below. 

```mermaid
%%{init: {
  "theme": "base",
  "themeVariables": {
    "primaryColor": "#B6ECE2",
    "primaryTextColor": "#160F26",
    "primaryBorderColor": "#065647",
    "lineColor": "#545555",
    "clusterBkg": "#BABCBD22",
    "clusterBorder": "#DDDEDE",
    "fontFamily": "arial"
  }
}}%%
flowchart TB

%% ===== REF_X_MOD PIPELINE =====
subgraph REF_X_MOD["Reference vs Modified Fasta Comparison Pipeline"]

    %% Inputs
    REF_FASTA["Reference FASTA"]:::input
    MOD_FASTA["Modified FASTA"]:::input

    %% Processes
    NUCMER["nucmer"]:::process
    DELTA["deltaFilter"]:::process
    SHOWCOORDS["showCoords"]:::process
    SYRI["syri"]:::process

    %% Output
    VCF_OUT["Structural Variant VCF"]:::output

    %% Connections
    REF_FASTA --> NUCMER
    MOD_FASTA --> NUCMER
    NUCMER --> DELTA --> SHOWCOORDS --> SYRI --> VCF_OUT

end

%% ===== STYLING =====
classDef input fill:#E3F2FD,stroke:#1565C0
classDef process fill:#B6ECE2,stroke:#065647
classDef output fill:#E8F5E9,stroke:#2E7D32

```

This folder contains results from the **reference vs modified FASTA comparison pipeline**:

```
fasta_ref_mod/
├── ref_x_mod.delta
├── ref_x_mod.filtered.coords
├── ref_x_mod_filtered.delta
└── ref_x_modsyri.vcf
```

**Description of files:**

* `ref_x_mod.delta`
  Raw alignment difference file between reference and modified FASTA (generated by `nucmer`/MUMmer).

* `ref_x_mod.filtered.coords`
  Filtered alignment coordinates showing high-confidence matches and structural differences.

* `ref_x_mod_filtered.delta`
  Cleaned and filtered delta file used for downstream structural comparison.

* `ref_x_modsyri.vcf`
  Structural variants and genome rearrangements detected by **SyRI**, stored in VCF format.

The table below summarises all tools used within the pipeline:

| **Tool**        | **Link for Further Information**                       |
| --------------- | ------------------------------------------------------ |
| **Nucmer**      | [MUMmer](https://mummer4.github.io/manual/manual.html) |
| **deltaFilter** | [MUMmer](https://mummer4.github.io/manual/manual.html) |
| **showCoords**  | [MUMmer](https://mummer4.github.io/manual/manual.html) |
| **Syri**        | [SyRI GitHub](https://schneebergerlab.github.io/syri/) |


---

### `illumina/`

The flowchart below summarizes the pipeline for processing short reads. VCF annotation is performed only when a GFF/GTF annotation file is provided. Delly and Freebayes are run exclusively for reference genome mapping; these steps are skipped when reads are mapped to a modified genome or a plasmid.

```mermaid
%%{init: {
  "theme": "base",
  "themeVariables": {
    "primaryColor": "#B6ECE2",
    "primaryTextColor": "#160F26",
    "primaryBorderColor": "#065647",
    "lineColor": "#545555",
    "clusterBkg": "#BABCBD22",
    "clusterBorder": "#DDDEDE",
    "fontFamily": "arial"
  }
}}%%
flowchart TB

%% ===== INPUTS =====
RAW["Raw Illumina Reads"]
REF["Reference FASTA"]
PLASMID_REF["Plasmid FASTA"]
GFF["GFF / GTF Annotation"]
CONFIG["SnpEff Config"]

style RAW fill:#E3F2FD,stroke:#1565C0,stroke-width:2px
style REF fill:#E3F2FD,stroke:#1565C0,stroke-width:2px
style PLASMID_REF fill:#E3F2FD,stroke:#1565C0,stroke-width:2px
style GFF fill:#E3F2FD,stroke:#1565C0,stroke-width:2px
style CONFIG fill:#E3F2FD,stroke:#1565C0,stroke-width:2px

%% ===== QC =====
TRIM["TrimGalore"]
FASTQC["FastQC"]
MULTIQC_QC["MultiQC (QC)"]

RAW --> TRIM --> FASTQC --> MULTIQC_QC

%% ===== QC OUTPUTS =====
TRIMMED["Trimmed FASTQ"]:::output
FASTQC_OUT["FastQC report"]:::output
MULTIQC_QC_OUT["MultiQC QC report"]:::output

TRIM --> TRIMMED
FASTQC --> FASTQC_OUT
MULTIQC_QC --> MULTIQC_QC_OUT

%% ===== SHORT REF PIPELINE =====
BWA_INDEX["BWA index"]
BWA["BWA mapping"]
SORT["Samtools sort"]
STATS["Samtools stats"]
BAM_IDX["BAM index"]
PICARD["Picard metrics"]
FREEBAYES["Freebayes variant calling"]
BCFTOOLS["BCFtools stats"]
BUILD_CFG["Build snpEff config"]
SNPEFF["snpEff annotation"]
GET_UNMAPPED["Get unmapped reads"]

REF --> BWA_INDEX --> BWA
TRIMMED --> BWA
BWA --> SORT --> SORTED_BAM["Sorted BAM"]:::output
SORT --> STATS --> STATS_OUT["Samtools stats"]:::output
SORT --> BAM_IDX --> BAM_INDEX_OUT["BAM index"]:::output

%% ===== SHORT REF SV PIPELINE (Delly) =====
DELLY["Delly (SV calling)"]
BCF2VCF["Convert BCF to VCF"]
BAM_IDX --> DELLY --> BCF2VCF --> SV_VCF["SV VCF"]:::output

BAM_IDX --> GET_UNMAPPED --> UNMAPPED_OUT["Unmapped reads FASTQ"]:::output

BAM_IDX --> FREEBAYES --> VCF_RAW["Raw VCF"]:::output
VCF_RAW --> BCFTOOLS --> BCF_STATS["BCFtools stats"]:::output
GFF --> BUILD_CFG --> SNPEFF
CONFIG --> BUILD_CFG
VCF_RAW --> SNPEFF --> VCF_ANNOT["Annotated VCF"]:::output

BAM_IDX --> PICARD --> PICARD_OUT["Picard metrics"]:::output


%% ===== PLASMID PIPELINE =====
PL_INDEX["BWA index"]
PL_BWA["BWA mapping"]
PL_SORT["Samtools sort"]
PL_STATS["Samtools stats"]
PL_BAM_IDX["BAM index"]
PL_PICARD["Picard metrics"]
GET_UNMAPPED_PL["Get unmapped reads plasmid"]

PLASMID_REF --> PL_INDEX --> PL_BWA
UNMAPPED_OUT --> PL_BWA
PL_BWA --> PL_SORT --> PL_SORTED_BAM["Sorted BAM (plasmid)"]:::output
PL_SORT --> PL_STATS --> PL_STATS_OUT["Samtools stats (plasmid)"]:::output
PL_SORT --> PL_BAM_IDX --> PL_BAM_INDEX_OUT["BAM index"]:::output
PL_BAM_IDX --> PL_PICARD --> PL_PICARD_OUT["Picard metrics (plasmid)"]:::output
PL_BAM_IDX --> GET_UNMAPPED_PL --> PL_UNMAPPED_OUT["Unmapped reads FASTQ (plasmid)"]:::output

%% ===== STYLING =====
classDef input fill:#E3F2FD,stroke:#1565C0
classDef process fill:#B6ECE2,stroke:#065647
classDef output fill:#E8F5E9,stroke:#2E7D32
class TRIMMED,FASTQC_OUT,MULTIQC_QC_OUT,SORTED_BAM,STATS_OUT,BAM_INDEX_OUT,PICARD_OUT,MULTIQC_MAP_OUT,UNMAPPED_OUT,VCF_RAW,BCF_STATS,VCF_ANNOT,SV_VCF,PL_SORTED_BAM,PL_STATS_OUT,PL_PICARD_OUT,PL_MULTIQC_OUT,PL_UNMAPPED_OUT output

%% ===== LEGEND =====
subgraph LEGEND["Legend"]
    L1["Input"]:::input
    L2["Process"]:::process
    L3["Output file"]:::output
end
```

This folder contains the full output of the **Illumina short-read processing pipeline**, including read quality control, trimming, genome mapping, and variant analysis.

```
data/outputs/illumina/
├── qc_trimming
│   ├── fastqc_out
│   ├── multiqc
│   └── trimmed_reads
├── short-mod
│   ├── bam
│   ├── bwa_index
│   ├── multiqc
│   ├── picard
│   ├── samtools_stats
│   └── unmapped_fastq
├── short-ref
│   ├── bam
│   ├── bcftools_stats
│   ├── bwa_index
│   ├── multiqc
│   ├── picard
│   ├── samtools_index_dict
│   ├── samtools_stats
│   ├── unmapped_fastq
│   └── vcf
└── short-ref-plasmid
    ├── bam
    ├── bwa_index
    ├── multiqc
    ├── picard
    ├── samtools_stats
    └── unmapped_fastq
```

#### `qc_trimming/`

This directory contains all quality control and preprocessing outputs generated from raw Illumina reads.

* `fastqc_out/`
  Raw read quality reports (per-sample) generated by **FastQC**.

* `multiqc/`
  Aggregated quality control report summarizing all FastQC results and trimming reports.

* `trimmed_reads/`
  Quality-filtered and adapter-trimmed reads used for downstream mapping.


#### `short-ref/`

This folder contains Illumina reads mapped to the **reference genome**.

Includes:

* `bam/` — Sorted and indexed BAM alignment files
* `bwa_index/` — Precomputed BWA reference genome indices
* `samtools_index_dict/` — FASTA index and sequence dictionary files
* `samtools_stats/` — Alignment and coverage statistics
* `picard/` — Alignment QC metrics
* `bcftools_stats/` — Variant calling summary statistics
* `vcf/` — Variant calls generated by Delly and (SVs) and FreeBayes (SNP and INDELs) and annotated VCF (SNP and INDELs), if gff or gtf files are present in `data/valid`
* `multiqc/` — Combined QC report from mapping and alignment metrics and variant calling metrics
* `unmapped_fastq/` — fastq file with reads that failed to align to the reference genome


#### `short-ref-plasmid/`

This folder holds the mapping results of Illumina reads aligned to the reference plasmid fasta. It is created only if a reference plasmid is present in the `data/valid` folder. A folder with a similar structure, `short-mod-plasmid/`, is created if a modified plasmid is present within the `data/valid` folder.

Includes:

* `bam/` — Aligned reads (that were not mapped to the reference) mapped to the plasmid
* `bwa_index/` — Plasmid reference index files
* `samtools_stats/` — Mapping and coverage statistics
* `picard/` — Alignment QC metrics
* `multiqc/` — Summary report of mapping and alignment metrics
* `unmapped_fastq/` — Reads not mapped to the plasmid and not mapped to the reference genome

#### `short-mod/`

This folder contains Illumina read alignments against the **modified/assembled genome**.

Includes:

* `bam/` — Sorted BAM files for modified genome mapping
* `bwa_index/` — Modified genome BWA index
* `samtools_stats/` — Mapping and coverage statistics
* `picard/` — Alignment QC metrics
* `multiqc/` — Summary report of mapping and alignment metrics
* `unmapped_fastq/` — Fastq file containing reads that failed to align to the modified genome

The table below summarises all tools used within the pipeline:

| **Tool**        | **Link for Further Information**                                     |
| --------------- | -------------------------------------------------------------------- |
| **Trim Galore** | [Trim Galore](https://github.com/FelixKrueger/TrimGalore)            |
| **FastQC**      | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) |
| **MultiQC**     | [MultiQC](https://multiqc.info/)                                     |
| **BWA**         | [BWA](http://bio-bwa.sourceforge.net/)                               |
| **Picard**      | [Picard](https://gatk.broadinstitute.org/hc/en-us/articles/360040507751-CollectAlignmentSummaryMetrics-Picard)                   |
| **Samtools**    | [Samtools](https://www.htslib.org/doc/samtools.html)                 |
| **BCFtools**    | [BCFtools](https://samtools.github.io/bcftools/)                     |
| **FreeBayes**   | [FreeBayes](https://github.com/freebayes/freebayes)                  |
| **SnpEff**      | [SnpEff](http://snpeff.sourceforge.net/)                             |
| **Delly**       | [Delly](https://github.com/dellytools/delly)                         |


---

### `pacbio/` and `ont/`

This workflow shows the processing of raw long-read sequencing data (PacBio or Nanopore) from quality control to mapping. Reads undergo NanoPlot QC, then mapped to the reference or modified genome with minimap2, followed by sorting, indexing, and calculation of unmapped reads. Structural variant calling using cute_sv, debreak, and sniffles is performed only for reads mapped to the reference genome, and results are merged with SURVIVOR and summarized with bcftools stats, producing the final long-read VCF. Reads mapped to modified or plasmid sequences skip structural variant calling.

```mermaid
%%{init: {
  "theme": "base",
  "themeVariables": {
    "primaryColor": "#B6ECE2",
    "primaryTextColor": "#160F26",
    "primaryBorderColor": "#065647",
    "lineColor": "#545555",
    "clusterBkg": "#BABCBD22",
    "clusterBorder": "#DDDEDE",
    "fontFamily": "arial"
  }
}}%%
flowchart TB

%% ===== INPUTS =====
LONG_READS["Raw Long Reads (PacBio / Nanopore)"]
REF["Reference FASTA"]
PLASMID_REF["Plasmid FASTA"]

style LONG_READS fill:#E3F2FD,stroke:#1565C0,stroke-width:2px
style REF fill:#E3F2FD,stroke:#1565C0,stroke-width:2px
style PLASMID_REF fill:#E3F2FD,stroke:#1565C0,stroke-width:2px

%% ===== QC =====
NANO_PLOT["NanoPlot QC"]
MULTIQC_QC["MultiQC (QC)"]
NANO_PLOT_OUT["NanoPlot QC report"]:::output
MULTIQC_QC_OUT["MultiQC QC report"]:::output

LONG_READS --> NANO_PLOT --> MULTIQC_QC
NANO_PLOT --> NANO_PLOT_OUT
MULTIQC_QC --> MULTIQC_QC_OUT

%% ===== LONG REF MAPPING PIPELINE =====
MINIMAP2["minimap2 mapping"]
SORT_BAM["Samtools sort"]
SORTED_BAM["Sorted BAM"]:::output
BAM_IDX["Samtools index BAM"]
BAM_IDX_OUT["BAM index"]:::output
GET_UNMAPPED["Get unmapped reads"]
UNMAPPED_OUT["Unmapped reads FASTQ"]:::output

REF --> MINIMAP2
LONG_READS --> MINIMAP2
MINIMAP2 --> SORT_BAM --> SORTED_BAM
SORT_BAM --> BAM_IDX --> BAM_IDX_OUT
BAM_IDX --> GET_UNMAPPED --> UNMAPPED_OUT

%% ===== PLASMID PIPELINE =====
MINIMAP2_PLASMID["minimap2 mapping (plasmid)"]
SORT_BAM_PLASMID["Samtools sort (plasmid)"]
SORTED_BAM_PLASMID["Sorted BAM (plasmid)"]:::output
BAM_IDX_PLASMID["Samtools index BAM (plasmid)"]
BAM_IDX_PLASMID_OUT["BAM index (plasmid)"]:::output
GET_UNMAPPED_PL["Get unmapped plasmid reads"]
UNMAPPED_PL_OUT["Unmapped plasmid reads FASTQ"]:::output

PLASMID_REF --> MINIMAP2_PLASMID
UNMAPPED_OUT --> MINIMAP2_PLASMID
MINIMAP2_PLASMID --> SORT_BAM_PLASMID --> BAM_IDX_PLASMID --> GET_UNMAPPED_PL
SORT_BAM_PLASMID --> SORTED_BAM_PLASMID
BAM_IDX_PLASMID --> BAM_IDX_PLASMID_OUT
GET_UNMAPPED_PL --> UNMAPPED_PL_OUT

%% ===== SV CALLING PIPELINE =====
CUTE_SV["cute_sv"]
DEBREAK["debreak"]
SNIFFLES["sniffles"]
SURVIVOR["survivor"]
BCFTOOLS_STATS["bcftools stats"]
LONG_VCF["Long-read VCF"]:::output

BAM_IDX --> CUTE_SV --> SURVIVOR 
BAM_IDX --> DEBREAK --> SURVIVOR
BAM_IDX --> SNIFFLES --> SURVIVOR
SURVIVOR --> BCFTOOLS_STATS --> LONG_VCF

%% ===== STYLING =====
classDef input fill:#E3F2FD,stroke:#1565C0
classDef process fill:#B6ECE2,stroke:#065647
classDef output fill:#E8F5E9,stroke:#2E7D32

%% ===== LEGEND =====
subgraph LEGEND["Legend"]
    L1["Input"]:::input
    L2["Process"]:::process
    L3["Output file"]:::output
end
```

These two folders contain the complete results from the **long-read analysis pipeline** using:

* **PacBio** reads OR
* **Oxford Nanopore Technologies (ONT)** reads

Both follow the **same folder structure** and processing logic.

```
data/outputs/ont/
data/outputs/pacbio/
├── long-mod
│   ├── bam
│   └── unmapped_fastq
├── long-ref
│   ├── bam
│   ├── bcftools_stats
│   ├── cutesv_out
│   ├── debreak_out
│   ├── sniffles_out
│   ├── survivor_out
│   └── unmapped_fastq
├── long-ref-plasmid
│   ├── bam
│   └── unmapped_fastq
└── nanoplot
    └── SampleName_report
```

#### `long-ref/`

Contains all outputs generated by mapping long reads to the **reference genome**.

Includes:

* `bam/`
  Sorted and indexed **BAM alignment files** of long reads mapped to the reference genome.

* `bcftools_stats/`
  Summary statistics of detected variants after variant calling.

* `cutesv_out/`
  Structural variants called using **cuteSV**.

* `sniffles_out/`
  Structural variants called using **Sniffles**.

* `debreak_out/`
  Structural variants detected using **DeBreak**.

* `survivor_out/`
  Merged structural variant callsets generated by **SURVIVOR**.

* `unmapped_fastq/`
  Long reads that failed to align to the reference genome.

---

#### `long-ref-plasmid/`

This folder holds the mapping results of long reads aligned to the reference plasmid sequence. It is created only if a reference plasmid is present in the `data/valid` folder. A folder with a similar structure, `long-mod-plasmid/`, is created if a modified plasmid is present within the `data/valid` folder.

Includes:

* `bam/` — Plasmid-mapped long-read alignments
* `unmapped_fastq/` — Fastq file containing reads that did not map to the plasmid


#### `long-mod/`

Contains alignments of long reads mapped to the **modified/assembled genome**.

Includes:

* `bam/` — Sorted alignment files
* `unmapped_fastq/` — Reads that failed to align to the modified genome

This enables comparison between mapping reads on reference vs modified assemblies.

#### `nanoplot/`

Contains long-read quality control and summary statistics generated using **NanoPlot**.

Example content:

* `SampleName_report/`

Inside this folder you typically find:

* Read length distributions
* N50 / N90 statistics
* Quality score profiles
* Read length vs quality plots
* Summary statistics of long-read sequencing quality

The table below summarises all tools used within the pipeline:

| **Tool**     | **Link for Further Information**                       |
| ------------ | ------------------------------------------------------ |
| **samtools** | [samtools](https://www.htslib.org/doc/samtools.html)   |
| **BCFtools** | [BCFtools](https://samtools.github.io/bcftools/)       |
| **cuteCV**   | [cuteCV](https://github.com/tjiangHIT/cuteSV?tab=readme-ov-file)             |
| **DeBreak**  | [DeBreak](https://github.com/Maggi-Chen/DeBreak)       |
| **Sniffles** | [Sniffles](https://github.com/fritzsedlazeck/Sniffles) |
| **SURVIVOR** | [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR) |
| **NanoPlot** | [NanoPlot](https://github.com/wdecoster/NanoPlot)      |


---

### `truvari/`

The flowchart illustrates the Truvari comparison pipeline for structural variant (SV) analysis. The Reference vs Modified VCF (that is outputed by  pipeline where reference and modified fasta are compared) serves as the baseline or truth-set, against which VCFs from PacBio, Nanopore, and Illumina sequencing are compared. The pipeline begins with sorting the VCF files (sortVcf), indexing them (indexVcf), and then performing the Truvari comparison to generate the final comparison results.

```mermaid
%%{init: {
  "theme": "base",
  "themeVariables": {
    "primaryColor": "#B6ECE2",
    "primaryTextColor": "#160F26",
    "primaryBorderColor": "#065647",
    "lineColor": "#545555",
    "clusterBkg": "#BABCBD22",
    "clusterBorder": "#DDDEDE",
    "fontFamily": "arial"
  }
}}%%
flowchart TB

%% ===== TRUVARI COMPARISON PIPELINE =====
subgraph Truvari_Comparision_Pipeline["Truvari Comparison Pipeline"]
    %% Inputs
    REF_MOD_VCF["Reference vs Modified VCF (Baseline / Truth-set)"]:::truthset
    PB_VCF["PacBio VCF"]:::input
    ONT_VCF["Nanopore VCF"]:::input
    IL_VCF["Illumina VCF"]:::input

    %% Processes
    SORT_VCF["sortVcf"]:::process
    INDEX_VCF["indexVcf"]:::process
    TRUVARI["truvari"]:::process

    %% Output
    TRUVARI_OUT["Truvari comparison results"]:::output

    %% Connections
    PB_VCF --> SORT_VCF
    ONT_VCF --> SORT_VCF
    IL_VCF --> SORT_VCF
    REF_MOD_VCF --> SORT_VCF

    SORT_VCF --> INDEX_VCF --> TRUVARI --> TRUVARI_OUT
end

%% ===== STYLING =====
classDef input fill:#E3F2FD,stroke:#1565C0
classDef truthset fill:#FFF3B0,stroke:#FFB300,stroke-width:2px
classDef process fill:#B6ECE2,stroke:#065647
classDef output fill:#E8F5E9,stroke:#2E7D32
```

#### Folder Structure

```
truvari
├── SampleName_sv_short_read.vcf.gz
├── SampleName_sv_short_read.vcf.gz.csi
├── SampleName.pacbio_sv_long_read.vcf.gz
├── SampleName.pacbio_sv_long_read.vcf.gz.csi
├── SampleName.ont_sv_long_read.vcf.gz
├── SampleName.ont_sv_long_read.vcf.gz.csi
├── ref_x_modsyri.vcf.gz
├── ref_x_modsyri.vcf.gz.csi
├── ref_x_modsyri_SampleName_sv_short_read_truvari       → folder comparing the short-read SV pipeline with the reference-to-modified fasta pipeline
├── ref_x_modsyri_SampleName.pacbio_sv_long_read_truvari  → folder comparing Pacbio long-read SV pipeline to reference-to-modified fasta pipeline
└── ref_x_modsyri_SampleName.ont_sv_long_read_truvari   → folder comparing Nanopore long-read SV pipeline to reference-to-modified fasta pipeline
```

#### Description

This folder contains all structural variant (SV) callsets and their **Truvari benchmarking results** comparing SVs detected from sequencing data with the structural variants derived from the **reference vs modified genome comparison (SyRI)**.


#### Reference SV Callsets

* `ref_x_modsyri.vcf.gz`
  Structural variants derived from comparing the **reference genome** and the **modified genome** using **SyRI**.

* `SampleName_sv_short_read.vcf.gz`
  Structural variants detected from **Illumina short reads**.

* `SampleName.pacbio_sv_long_read.vcf.gz`
  Structural variants detected from **PacBio long reads**.

* `SampleName.ont_sv_long_read.vcf.gz`
  Structural variants detected from **Oxford Nanopore long reads**.

All `.csi` files represent index files for fast querying of VCF contents.


### Truvari Comparison Result Folders

Each Truvari output directory contains benchmarking results comparing the **SyRI structural variants** against sequencing-based SV calls:

* `ref_x_modsyri_SampleName_sv_short_read_truvari/`
  Comparison between SyRI SVs and SVs called from **Illumina short reads**.

* `ref_x_modsyri_SampleName.pacbio_sv_long_read_truvari/`
  Comparison between SyRI SVs and SVs called from **PacBio long reads**.

* `ref_x_modsyri_SampleName.ont_sv_long_read_truvari/`
  Comparison between SyRI SVs and SVs called from **Oxford Nanopore long reads**.

Each Truvari output folder usually contains:

* Matched SV calls
* Unmatched (false negative / false positive) calls
* Precision, recall, and F1 scores
* Comparison summary statistics

The table below summarises all tools used within the pipeline:

| **Tool**     | **Link for Further Information**                      |
| ------------ | ----------------------------------------------------- |
| **Truvari**  | [Truvari GitHub](https://github.com/ACEnglish/truvari)|
| **BCFtools** | [BCFtools](https://samtools.github.io/bcftools/)      |

---

### `unmapped_stats/`

#### Folder Structure

```
unmapped_stats
├── SampleName_short_read_stats.txt
├── SampleName_pacbio_read_stats.txt
└── SampleName_ont_read_stats.txt
```

#### Description

This folder contains read mapping comparisons between the reference genome and the modified/assembled genome for short-read and long-read processing pipelines.

Each file summarizes:

* Reads mapping to both reference and modified assemblies
* Reads mapping only to reference
* Reads mapping only to modified
* Used to compare assembly quality and detect assembly-related differences

| File                               | Description                                                                                |
| ---------------------------------- | ------------------------------------------------------------------------------------------ |
| `SampleName_short_read_stats.txt`  | Mapping comparison of Illumina short reads between reference and modified assemblies       |
| `SampleName_pacbio_read_stats.txt` | Mapping comparison of PacBio long reads between reference and modified assemblies          |
| `SampleName_ont_read_stats.txt`    | Mapping comparison of Oxford Nanopore long reads between reference and modified assemblies |


Each report includes:

* Total input reads
* Number of unmapped reads
* Percentage of unmapped reads
* Mapping target (reference or modified genome)

### File Content 

Each mapping statistics (`SampleName_short_read_stats.txt`) file includes:
```
Pair IDs: <sample> reference vs <sample> modified
Intersected Reads: <number>
Unique Reads <sample> reference: <number>
Unique Reads <sample> modified: <number>
```

### Graphical Visualization of Long-Read Preprocessing

The flowchart below illustrate the preprocessing workflows applied to long reads to obtain mapping statistics.

```mermaid
%%{init: {
  "theme": "base",
  "themeVariables": {
    "primaryColor": "#B6ECE2",
    "primaryTextColor": "#160F26",
    "primaryBorderColor": "#065647",
    "lineColor": "#545555",
    "clusterBkg": "#BABCBD22",
    "clusterBorder": "#DDDEDE",
    "fontFamily": "arial"
  }
}}%%
flowchart TB

%% ===== INPUTS =====
LONG_READS["Raw Long Reads (PacBio / Nanopore)"]
REF["Reference FASTA"]
MOD_REF["Modified Genome FASTA"]
PLASMID_REF["Plasmid FASTA"]

style LONG_READS fill:#E3F2FD,stroke:#1565C0,stroke-width:2px
style REF fill:#E3F2FD,stroke:#1565C0,stroke-width:2px
style MOD_REF fill:#E3F2FD,stroke:#1565C0,stroke-width:2px
style PLASMID_REF fill:#E3F2FD,stroke:#1565C0,stroke-width:2px

%% ===== QC =====
NANO_PLOT["NanoPlot QC"]
MULTIQC_QC["MultiQC (QC)"]
NANO_PLOT_OUT["NanoPlot QC report"]:::output
MULTIQC_QC_OUT["MultiQC QC report"]:::output

LONG_READS --> NANO_PLOT --> MULTIQC_QC
NANO_PLOT --> NANO_PLOT_OUT
MULTIQC_QC --> MULTIQC_QC_OUT

%% ===============================
%% REFERENCE MAPPING PIPELINE
%% ===============================
subgraph REF_PIPE["Mapping to Reference Genome"]
    MINIMAP2_REF["minimap2 mapping (reference)"]
    SORT_REF["Samtools sort"]
    SORTED_REF["Sorted BAM (reference)"]:::output
    INDEX_REF["Samtools index BAM"]
    INDEX_REF_OUT["BAM index (reference)"]:::output
    UNMAP_REF["Get unmapped reads (reference)"]
    UNMAP_REF_OUT["Unmapped FASTQ (reference)"]:::output
end

REF --> MINIMAP2_REF
LONG_READS --> MINIMAP2_REF
MINIMAP2_REF --> SORT_REF --> SORTED_REF
SORT_REF --> INDEX_REF --> INDEX_REF_OUT
INDEX_REF --> UNMAP_REF --> UNMAP_REF_OUT


%% ==========================================
%% PLASMID RESCUE FOR REFERENCE UNMAPPED READS
%% ==========================================
subgraph PLASMID_PIPE["Map reference–unmapped reads to plasmid"]
    MINIMAP2_PLAS["minimap2 mapping (plasmid)"]
    SORT_PLAS["Samtools sort"]
    SORTED_PLAS["Sorted BAM (plasmid)"]:::output
    INDEX_PLAS["Samtools index BAM"]
    INDEX_PLAS_OUT["BAM index (plasmid)"]:::output
    UNMAP_PLAS["Get unmapped reads (plasmid)"]
    UNMAP_PLAS_OUT["Unmapped FASTQ (plasmid)"]:::output
end

UNMAP_REF_OUT --> MINIMAP2_PLAS
PLASMID_REF --> MINIMAP2_PLAS
MINIMAP2_PLAS --> SORT_PLAS --> SORTED_PLAS
SORT_PLAS --> INDEX_PLAS --> INDEX_PLAS_OUT
INDEX_PLAS --> UNMAP_PLAS --> UNMAP_PLAS_OUT


%% ===============================
%% MODIFIED GENOME MAPPING PIPELINE
%% ===============================
subgraph MOD_PIPE["Mapping to Modified Genome"]
    MINIMAP2_MOD["minimap2 mapping (modified genome)"]
    SORT_MOD["Samtools sort"]
    SORTED_MOD["Sorted BAM (modified)"]:::output
    INDEX_MOD["Samtools index BAM"]
    INDEX_MOD_OUT["BAM index (modified)"]:::output
    UNMAP_MOD["Get unmapped reads (modified genome)"]
    UNMAP_MOD_OUT["Unmapped FASTQ (modified genome)"]:::output
end

MOD_REF --> MINIMAP2_MOD
LONG_READS --> MINIMAP2_MOD
MINIMAP2_MOD --> SORT_MOD --> SORTED_MOD
SORT_MOD --> INDEX_MOD --> INDEX_MOD_OUT
INDEX_MOD --> UNMAP_MOD --> UNMAP_MOD_OUT


%% ===============================
%% FINAL UNMAPPED COMPARISON
%% ===============================
COMPARE_UNMAPPED["compare_unmapped"]
COMPARE_OUT["Comparison report"]:::output

UNMAP_PLAS_OUT --> COMPARE_UNMAPPED
UNMAP_MOD_OUT --> COMPARE_UNMAPPED
COMPARE_UNMAPPED --> COMPARE_OUT


%% ===== STYLING =====
classDef input fill:#E3F2FD,stroke:#1565C0
classDef process fill:#B6ECE2,stroke:#065647
classDef output fill:#E8F5E9,stroke:#2E7D32

%% ===== LEGEND =====
subgraph LEGEND["Legend"]
    L1["Input"]:::input
    L2["Process"]:::process
    L3["Output file"]:::output
end
```

### Graphical Visualization of Short-Read Preprocessing

The flowchart below illustrate the preprocessing workflows applied to short reads to obtain mapping statistics.

```mermaid
%%{init: {
  "theme": "base",
  "themeVariables": {
    "primaryColor": "#B6ECE2",
    "primaryTextColor": "#160F26",
    "primaryBorderColor": "#065647",
    "lineColor": "#545555",
    "clusterBkg": "#BABCBD22",
    "clusterBorder": "#DDDEDE",
    "fontFamily": "arial"
  }
}}%%

flowchart TB

%% ======================
%% INPUTS
%% ======================
RAW["Raw Illumina Reads"]
REF["Reference FASTA"]
MOD_REF["Modified Genome FASTA"]
PLASMID_REF["Plasmid FASTA"]

style RAW fill:#E3F2FD,stroke:#1565C0,stroke-width:2px
style REF fill:#E3F2FD,stroke:#1565C0,stroke-width:2px
style MOD_REF fill:#E3F2FD,stroke:#1565C0,stroke-width:2px
style PLASMID_REF fill:#E3F2FD,stroke:#1565C0,stroke-width:2px

%% ======================
%% QC + TRIMMING (COMMON)
%% ======================
subgraph QC_PIPE["Quality Control & Trimming (Common)"]
    TRIM["TrimGalore"]
    FASTQC["FastQC"]
    MULTIQC_QC["MultiQC (QC)"]

    RAW --> TRIM --> FASTQC --> MULTIQC_QC

    TRIM --> TRIMMED["Trimmed FASTQ"]:::output
    FASTQC --> FASTQC_OUT["FastQC Report"]:::output
    MULTIQC_QC --> MULTIQC_QC_OUT["MultiQC QC Report"]:::output
end

%% ======================
%% REFERENCE PIPELINE
%% ======================
subgraph REF_MAP["Mapping to Reference Genome"]
    BWA_INDEX_REF["BWA index (reference)"]
    BWA_REF["BWA mapping (reference)"]
    SORT_REF["Samtools sort (reference)"]
    BAM_IDX_REF["Samtools index (reference)"]
    UNMAP_REF["Extract unmapped reads (reference)"]

    REF --> BWA_INDEX_REF --> BWA_REF
    TRIMMED --> BWA_REF

    BWA_REF --> SORT_REF --> SORTED_REF["Sorted BAM (reference)"]:::output
    SORT_REF --> BAM_IDX_REF --> BAM_INDEX_REF["BAM index (reference)"]:::output
    BAM_IDX_REF --> UNMAP_REF --> UNMAP_REF_OUT["Unmapped reads (reference FASTQ)"]:::output
end

%% ======================
%% MODIFIED GENOME PIPELINE
%% ======================
subgraph MOD_MAP["Mapping to Modified Genome"]
    BWA_INDEX_MOD["BWA index (modified)"]
    BWA_MOD["BWA mapping (modified)"]
    SORT_MOD["Samtools sort (modified)"]
    BAM_IDX_MOD["Samtools index (modified)"]
    UNMAP_MOD["Extract unmapped reads (modified)"]

    MOD_REF --> BWA_INDEX_MOD --> BWA_MOD
    TRIMMED --> BWA_MOD

    BWA_MOD --> SORT_MOD --> SORTED_MOD["Sorted BAM (modified)"]:::output
    SORT_MOD --> BAM_IDX_MOD --> BAM_INDEX_MOD["BAM index (modified)"]:::output
    BAM_IDX_MOD --> UNMAP_MOD --> UNMAP_MOD_OUT["Unmapped reads (modified FASTQ)"]:::output
end

%% ======================
%% PLASMID PIPELINE
%% ======================
subgraph PLASMID_MAP["Mapping Unmapped-from-Reference to Plasmid"]
    PL_INDEX["BWA index (plasmid)"]
    PL_BWA["BWA mapping (plasmid)"]
    PL_SORT["Samtools sort (plasmid)"]
    PL_BAMIDX["Samtools index (plasmid)"]
    UNMAP_PL["Extract unmapped reads (plasmid)"]

    PLASMID_REF --> PL_INDEX --> PL_BWA
    UNMAP_REF_OUT --> PL_BWA

    PL_BWA --> PL_SORT --> PL_SORTED["Sorted BAM (plasmid)"]:::output
    PL_SORT --> PL_BAMIDX --> PL_IDX_OUT["BAM index (plasmid)"]:::output
    PL_BAMIDX --> UNMAP_PL --> UNMAP_PL_OUT["Unmapped reads (plasmid FASTQ)"]:::output
end

%% ======================
%% COMPARE UNMAPPED READS
%% ======================
COMPARE["compare_unmapped"]
COMPARE_OUT["Comparison Report"]:::output

UNMAP_PL_OUT --> COMPARE
UNMAP_MOD_OUT --> COMPARE
COMPARE --> COMPARE_OUT

%% ======================
%% STYLES
%% ======================
classDef input fill:#E3F2FD,stroke:#1565C0
classDef process fill:#B6ECE2,stroke:#065647
classDef output fill:#E8F5E9,stroke:#2E7D32

%% Legend
subgraph LEGEND["Legend"]
    L1["Input"]:::input
    L2["Process"]:::process
    L3["Output file"]:::output
end
```

---

### `logs/` — Nextflow Command and Log Files

#### Folder Contents

```
logs/
├── .command.begin    # Timestamp file marking the start of a process
├── .command.err      # Captures standard error output from the process
├── .command.log      # Logs process execution messages from Nextflow
├── .command.out      # Captures standard output from the process
├── .command.run      # Execution metadata (exit status, runtime, resources)
└── .command.sh       # The shell script containing the exact commands executed
```

#### Description

The `logs/` folder contains **detailed logs and command scripts** for each Nextflow process.

* **`.command.begin`** — Marks the start time of a process.
* **`.command.err`** — Captures standard error messages generated by the process.
* **`.command.log`** — General execution logs from Nextflow for the process.
* **`.command.out`** — Captures standard output of the process.
* **`.command.run`** — Metadata about process execution (e.g., exit code, runtime, resource usage).
* **`.command.sh`** — The shell script that Nextflow runs; contains the exact commands for the process.
