# Graphical Representation of the Pipeline

## High-Level Overview

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

## Detailed Workflow Diagrams

For detailed workflow diagrams of specific pipeline components, see the **[Outputs Documentation](../outputs/)** section:

- [Reference vs Modified FASTA Comparison](../outputs/fasta-ref-mod.md)
- [Illumina Short-Read Pipeline](../outputs/illumina.md)
- [Long-Read Pipeline (PacBio & ONT)](../outputs/long-reads.md)
- [Truvari Comparison Pipeline](../outputs/truvari.md)

## See Also

- [Running the Pipeline](running-pipeline.md) - Execution commands
- [Directory Structures](directory-structures.md) - Input/output organization
