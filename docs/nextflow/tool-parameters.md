# Tool Parameter Reference

This page documents the hardcoded analysis parameters used by the pipeline's bioinformatics tools. These values have been chosen for **bacterial/prokaryotic genomes** with typical sequencing coverage and may require adjustment for other organisms or experimental designs.

---

## Short-Read Pipeline

### Freebayes (variant calling)

**Module:** `modules/variant_calling.nf`

```bash
freebayes -f <ref> --min-coverage 10 --min-base-quality 20 --min-mapping-quality 30 --min-alternate-count 3 <bam>
```

| Parameter | Value | Freebayes Default | Description |
|-----------|-------|-------------------|-------------|
| `--min-coverage` | 10 | 0 | Minimum number of reads covering a locus to call a variant. The default (0) would attempt calls at sites with a single read, producing many false positives. A value of 10 is a standard threshold for reliable variant calling with moderate Illumina coverage. |
| `--min-base-quality` | 20 | 0 | Minimum per-base Phred quality score. Q20 (99% per-base accuracy) is the widely accepted quality floor for Illumina data. The default (0) would include low-confidence bases, significantly increasing error-driven false calls. |
| `--min-mapping-quality` | 30 | 1 | Minimum read mapping Phred score. Q30 (99.9% mapping confidence) ensures only confidently placed reads contribute to variant calls. The default (1) would include multi-mapped and ambiguously placed reads, particularly problematic in repetitive bacterial regions. |
| `--min-alternate-count` | 3 | 2 | Minimum number of reads supporting the alternate allele. Slightly stricter than the default (2), providing an additional guard against sequencing-error-driven false positives in haploid genomes. |

**Rationale:** These are standard community thresholds for calling SNPs and small indels in haploid bacterial genomes with moderate Illumina coverage (30–100×). Each value departs from the freebayes default to reduce false positives — particularly important in a regulatory/safety context (GMO assessment). The freebayes defaults are intentionally permissive to support diverse use cases (e.g., polyploid organisms, low-frequency somatic variants); for haploid bacterial genomes these permissive defaults would produce excessive noise.

**Trade-offs:**

- *Increasing* thresholds improves specificity but may miss low-frequency or low-coverage variants.
- *Decreasing* thresholds (toward defaults) improves sensitivity but increases the false positive rate.

### Delly (structural variant calling)

**Module:** `modules/sv_calling.nf`

```bash
delly call -g <ref> -o <out.bcf> <bam>
```

All parameters use Delly defaults. No custom thresholds are applied.

---

## Long-Read Pipeline

### cuteSV (structural variant calling)

**Module:** `modules/sv_calling.nf`

```bash
cuteSV <bam> <ref> <out.vcf> <work_dir> -t 1
```

| Parameter | Value | Description |
|-----------|-------|-------------|
| `-t` | 1 | Number of threads. Set to 1 as a safe default for shared compute nodes. |

All other parameters use cuteSV defaults (minimum SV size 50 bp, minimum support 10 reads, etc.). These defaults are suitable for moderate-coverage PacBio/ONT data on bacterial genomes.

!!! tip
    For faster execution on dedicated compute nodes, consider increasing the thread count (e.g., `-t ${params.max_cpu}`).

### Sniffles (structural variant calling)

**Module:** `modules/sv_calling.nf`

All parameters use Sniffles defaults. No custom thresholds are applied.

### DeBreak (structural variant calling)

**Module:** `modules/sv_calling.nf`

```bash
debreak --bam <bam> -r <ref> -o <out_dir> -t ${params.max_cpu}
```

| Parameter | Value | Description |
|-----------|-------|-------------|
| `-t` | `${params.max_cpu}` | Number of threads. Uses the pipeline-wide `--max_cpu` setting. |

All other parameters use DeBreak defaults.

---

## Reference vs Modified Genome Comparison

### NUCmer (whole-genome alignment)

**Module:** `modules/assembly.nf`

```bash
nucmer --maxmatch -c 100 -b 500 -l 50 <ref> <mod> -p <prefix>
```

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--maxmatch` | *(flag)* | Use all anchor matches, not just unique ones. Required to detect duplicated and repetitive regions. |
| `-c` | 100 | Minimum alignment cluster length in bp. Filters spurious short alignments. |
| `-b` | 500 | Maximum gap (bp) between clustered matches. Allows merging across small indels or rearrangements. |
| `-l` | 50 | Minimum exact match (MEM) length in bp. Balances sensitivity with noise reduction. |

**Rationale:** These values are optimized for pairwise comparison of closely related prokaryotic genomes (reference vs. genetically modified organism). The settings are permissive enough to capture meaningful structural differences while filtering alignment noise.

**Trade-offs:**

- *Lowering* `-c` or `-l` increases sensitivity to small structural variants but may introduce spurious alignments.
- *Increasing* `-b` captures larger structural events but risks merging non-contiguous alignments.

### Delta-filter

**Module:** `modules/assembly.nf`

```bash
delta-filter -m -i 90 -l 100 <delta>
```

| Parameter | Value | Description |
|-----------|-------|-------------|
| `-m` | *(flag)* | Keep only the best set of alignments (many-to-many mapping). |
| `-i` | 90 | Minimum alignment identity (%). Filters low-quality alignments. |
| `-l` | 100 | Minimum alignment length (bp). Removes very short, unreliable alignments. |

---

## Summary Table

| Tool | Pipeline | Module | Key Non-Default Parameters |
|------|----------|--------|---------------------------|
| Freebayes | Short-read | `variant_calling.nf` | `--min-coverage 10`, `--min-base-quality 20`, `--min-mapping-quality 30`, `--min-alternate-count 3` |
| Delly | Short-read | `sv_calling.nf` | Defaults only |
| cuteSV | Long-read | `sv_calling.nf` | `-t 1` |
| Sniffles | Long-read | `sv_calling.nf` | Defaults only |
| DeBreak | Long-read | `sv_calling.nf` | `-t ${params.max_cpu}` |
| NUCmer | Ref vs Mod | `assembly.nf` | `--maxmatch`, `-c 100`, `-b 500`, `-l 50` |
| Delta-filter | Ref vs Mod | `assembly.nf` | `-m`, `-i 90`, `-l 100` |
