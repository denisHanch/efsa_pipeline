# Subworkflows

The file `workflows/subworkflows.nf` contains **reusable sub-workflows** that are shared across the main pipeline workflows. Each sub-workflow composes several processes into a logical unit with defined inputs and outputs.

## Short-Read Subworkflows

### `qc`

Quality control and adapter trimming for Illumina short reads.

| Step | Tool |
|------|------|
| Adapter trimming | TrimGalore |
| Per-read QC | FastQC |
| Aggregated QC report | MultiQC |

**Takes:** `fastqs`, `out_folder_name`
**Emits:** `trimmed` — trimmed FASTQ channel

---

### `mapping`

Read alignment and post-alignment QC for short reads.

| Step | Tool |
|------|------|
| Alignment | BWA |
| Sort | samtools sort |
| Alignment stats | samtools stats |
| BAM indexing | samtools index |
| Duplicate metrics | Picard MarkDuplicates |
| Aggregated report | MultiQC |

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

### `annotate_vcf`

VCF annotation with gene-level functional effects.

| Step | Tool |
|------|------|
| Build SnpEff database | build_config |
| Annotate VCF | SnpEff |

**Takes:** `fasta`, `gtf`, `vcf`, `feature_tag`, `build_setting`, `out_folder_name`
**Emits:** `qc_vcf` — SnpEff HTML report channel

---

## Long-Read Subworkflows

### `mapping_long`

Read alignment for PacBio or ONT long reads.

| Step | Tool |
|------|------|
| Alignment | minimap2 |
| Sort | samtools sort |
| BAM indexing | samtools index |

**Takes:** `fastqs`, `fasta`, `mapping_tag`, `out_folder_name`
**Emits:** `indexed_bam` — sorted + indexed BAM channel

---

### `sv_long`

Structural variant calling from long-read alignments using three callers, with consensus merging.

| Step | Tool |
|------|------|
| SV calling | CuteSV, Debreak, Sniffles |
| Supporting read extraction | extract_supp_reads (×3) |
| Consensus merge | SURVIVOR |
| VCF stats | bcftools stats |

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
