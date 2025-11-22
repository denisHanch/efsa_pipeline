# EFSA Pipeline

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
    -w "$WORKSPACE_PATH" \
    -v "$WORKSPACE_PATH:$WORKSPACE_PATH" \
    $INPUT_MOUNT \
    efsa-pipeline

   docker exec -it efsa-pipeline-container /bin/sh
   ```

# Input Validation

> **Important**
> When container is build please follow the steps to preprocess the data with a validation package.

The input validation module preprocesses and verifies all input data to ensure it meets the required format and structure before the Nextflow pipeline is executed.

# Nextflow

Please run the `docker login` command that is given in the efsa Slack channel - Nextflow notes.

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


| Data Type       | Supported Extensions                   |
| --------------- | -------------------------------------- |
| FASTA sequences | `.fa`, `.fna`, `.fasta`                |
| GFF annotations | `.gff`, `.gtf`                         |
| FASTQ reads     | `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz` |



## Running the Pipeline

The main pipeline (`main.nf`) executes **all three workflows** in sequence.
Each workflow can also be executed individually if required.

---

### Running Main Workflow

This executes **short-read processing**, **long-read processing**, and **reference vs modified genome comparison**:

```bash
nextflow run main.nf \
  -process.containerOptions "-u $(id -u):$(id -g)" \
  --max_cpu $(nproc) \
  -resume
```


### Available Options

| Option         | Description                                  | Default                |
| -------------- | -------------------------------------------- | ---------------------- |
| `-resume`      | Resume pipeline from last successful process | –                      |
| `--in_dir`     | Input directory                              | `data/valid`           |
| `--out_dir`    | Output directory                             | `data/outputs`         |
| `--registry`   | Docker/Singularity container registry        | `ghcr.io/kate-simonova`|
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
nextflow run workflows/short-read-ref.nf \
  -process.containerOptions "-u $(id -u):$(id -g)" \
  --max_cpu $(nproc) \
  -resume
```

#### Long-read Processing

For Oxford Nanopore / PacBio long reads:

```bash
nextflow run workflows/long-read-ref.nf \
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

Here’s a **README-ready section** that cleanly incorporates your runtime messages *and* the unmapped read statistics, with clear structure and consistent terminology.

You can paste this directly into your documentation.

## 🔄 Pipeline Runtime Messages & Mapping Summary

During execution, the pipeline prints progress messages indicating which workflow is currently running and what type of reads are being processed.

### ▶ Runtime Status Messages

When the pipeline is running, you will see real-time messages like:

```text
ℹ️  Running pipeline: processing long-pacbio reads → mapping to the reference & modified fasta.

ℹ️  Running pipeline: processing long-ont reads → mapping to the reference & modified fasta.

ℹ️  Running pipeline: processing short reads → mapping to the reference & modified fasta.

ℹ️  Truvari: performing 3 comparisons.
```

These messages help track the execution order and confirm that all three pipelines are being executed as expected.

---

## 📊 Unmapped Reads Statistics

After mapping, the pipeline reports the number and percentage of **unmapped reads** for each analysis.
This is useful for assessing mapping efficiency and data quality.

### Example Output

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
```

### Interpretation

* **Unmapped reads** represent sequences that did not align to the provided reference or modified FASTA files.
* A low percentage of unmapped reads indicates:

  * High mapping quality
  * Good reference/assembly quality
  * Low contamination or sequencing errors

If the percentage of unmapped reads is unusually high, this may indicate:

* Poor read quality
* Inadequate or incomplete reference
* Contamination
* Incorrect input file selection


## 📁 `data/outputs` Directory Structure

This directory contains all input data used by the Nextflow pipeline.
