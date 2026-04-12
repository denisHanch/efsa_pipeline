# Subworkflows

The file `workflows/subworkflows.nf` contains **reusable sub-workflows** that are shared across the main pipeline workflows. Each sub-workflow composes several processes into a logical unit with defined inputs and outputs.

## Short-Read Subworkflows

### `qc`

Quality control and adapter trimming for Illumina short reads.

| Step | Tool | Version | License |
|------|------|---------|---------|
| Adapter trimming | TrimGalore | 0.6.10 | GPL-3.0 |
| Per-read QC | FastQC | 0.11.9 | GPL-3.0 |
| Aggregated QC report | MultiQC | 1.33 | GPL-3.0 |

**Takes:** `fastqs`, `out_folder_name`
**Emits:** `trimmed` — trimmed FASTQ channel

---

### `mapping`

Read alignment and post-alignment QC for short reads.

| Step | Tool | Version | License |
|------|------|---------|---------|
| Alignment | BWA | 0.7.17 | GPL-3.0 |
| Sort | samtools sort | 1.23 | MIT |
| Alignment stats | samtools stats | 1.23 | MIT |
| BAM indexing | samtools index | 1.23 | MIT |
| Duplicate metrics | Picard MarkDuplicates | 3.4.0 | MIT |
| Aggregated report | MultiQC | 1.33 | GPL-3.0 |

**Takes:** `fasta`, `fasta_index`, `trimmed`, `out_folder_name`
**Emits:** `indexed_bam` — sorted + indexed BAM channel

---

### `sv`

Structural variant calling from short-read alignments.

| Step | Tool |
|------|------|
| FASTA indexing | samtools faidx |
| Sequence dictionary | Picard CreateSequenceDictionary |
| SV calling | Delly |
| BCF → VCF conversion | bcftools view |

**Takes:** `fasta`, `indexed_bam`, `out_folder_name`
**Emits:** `sv_vcf` — structural variant VCF channel

---

## Long-Read Subworkflows

### `mapping_long`

Read alignment for PacBio or ONT long reads.

| Step | Tool | Version | License |
|------|------|---------|---------|
| Alignment | minimap2 | 2.30 | MIT |
| Sort | samtools sort | 1.23 | MIT |
| BAM indexing | samtools index | 1.23 | MIT |

**Takes:** `fastqs`, `fasta`, `mapping_tag`, `out_folder_name`
**Emits:** `indexed_bam` — sorted + indexed BAM channel

---

### `sv_long`

Structural variant calling from long-read alignments using three callers, with consensus merging.

| Step | Tool | Version | License |
|------|------|---------|---------|
| SV calling | CuteSV | 2.1.3 | MIT |
| SV calling | Debreak | 1.2 | MIT |
| SV calling | Sniffles | 2.7.3 | MIT |
| Supporting read extraction | bcftools | 1.23 | MIT |
| Consensus merge | SURVIVOR | 1.0.7 | MIT |
| VCF stats | bcftools stats | 1.23 | MIT |

**Takes:** `fasta`, `fai`, `indexed_bam`, `mapping_tag`, `out_folder_name`
**Emits:** `merged_vcf` — SURVIVOR-merged VCF channel, `supp_reads` — supporting reads TSV channel

## How Subworkflows Are Used

The main workflow files (`short_read.nf`, `long_read.nf`) import and compose these sub-workflows:

```groovy
// In short_read.nf
qc(fastqs, out_folder)
mapping(fasta, fasta_index, qc.out.trimmed, out_folder)
sv(fasta, mapping.out.indexed_bam, out_folder)
```

```groovy
// In long_read.nf
mapping_long(fastqs, fasta, mapping_tag, out_folder)
sv_long(fasta, fai, mapping_long.out.indexed_bam, mapping_tag, out_folder)
```
