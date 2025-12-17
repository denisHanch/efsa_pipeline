# Unmapped Reads Statistics

## Directory Structure

```
unmapped_stats/
├── SampleName_short_read_stats.txt
├── SampleName_pacbio_read_stats.txt
└── SampleName_ont_read_stats.txt
```

## Description

This folder contains read mapping comparisons between the reference genome and the modified/assembled genome for short-read and long-read processing pipelines.

Each file summarizes:

- Reads mapping to both reference and modified assemblies
- Reads mapping only to reference
- Reads mapping only to modified
- Used to compare assembly quality and detect assembly-related differences

## File Descriptions

| File                               | Description                                                                                |
| ---------------------------------- | ------------------------------------------------------------------------------------------ |
| `SampleName_short_read_stats.txt`  | Mapping comparison of Illumina short reads between reference and modified assemblies       |
| `SampleName_pacbio_read_stats.txt` | Mapping comparison of PacBio long reads between reference and modified assemblies          |
| `SampleName_ont_read_stats.txt`    | Mapping comparison of Oxford Nanopore long reads between reference and modified assemblies |

Each report includes:

- Total input reads
- Number of unmapped reads
- Percentage of unmapped reads
- Mapping target (reference or modified genome)

## File Content Format

Each mapping statistics file includes:

```
Pair IDs: <sample> reference vs <sample> modified
Intersected Reads: <number>
Unique Reads <sample> reference: <number>
Unique Reads <sample> modified: <number>
```

## Graphical Visualization of Long-Read Preprocessing Pipeline

The flowchart below illustrates the preprocessing workflows applied to long reads to obtain mapping statistics.

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

## Graphical Visualization of Short-Read Preprocessing

The flowchart below illustrates the preprocessing workflows applied to short reads to obtain mapping statistics.

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

## See Also

- [Short-Read Processing Pipeline Output](illumina.md) - Short-read mapping details
- [Long-Read Processing Pipeline Output](long-reads.md) - Long-read mapping details
- [Runtime Messages](../nextflow/runtime-messages.md) - Pipeline execution messages
