# Tool Parameter Reference

This page documents the hardcoded and automatically derived analysis parameters used by the pipeline's bioinformatics tools. The pipeline supports both **prokaryotic** and **eukaryotic** genomes. Validation writes the normalized organism type to `validated_params.json` as `organism_type`, and downstream Nextflow processes consume that value.

---

## Short-Read Pipeline

### Freebayes (variant calling)

**Module:** `modules/variant_calling.nf`

```bash
freebayes -f <ref> --ploidy <1|2> --min-coverage 10 --min-base-quality 20 --min-mapping-quality 30 --min-alternate-count 3 <bam>
```

| Parameter | Value | Freebayes Default | Description |
|-----------|-------|-------------------|-------------|
| `--ploidy` | Derived from `organism_type`: `1` for `prokaryote`, `2` for `eukaryote` | 2 | Sample ploidy used by FreeBayes. The value is not user-supplied directly; it is derived from the validated organism type and passed through the short-read workflow to FreeBayes. |
| `--min-coverage` | 10 | 0 | Minimum number of reads covering a locus to call a variant. The default (0) would attempt calls at sites with a single read, producing many false positives. A value of 10 is a standard threshold for reliable variant calling with moderate Illumina coverage. |
| `--min-base-quality` | 20 | 0 | Minimum per-base Phred quality score. Q20 (99% per-base accuracy) is the widely accepted quality floor for Illumina data. The default (0) would include low-confidence bases, significantly increasing error-driven false calls. |
| `--min-mapping-quality` | 30 | 1 | Minimum read mapping Phred score. Q30 (99.9% mapping confidence) ensures only confidently placed reads contribute to variant calls. The default (1) would include multi-mapped and ambiguously placed reads, which is problematic in repetitive regions of both prokaryotic and eukaryotic genomes. |
| `--min-alternate-count` | 3 | 2 | Minimum number of reads supporting the alternate allele. Slightly stricter than the default (2), providing an additional guard against sequencing-error-driven false positives. |

**Ploidy derivation:** The validated `organism_type` controls the ploidy passed to FreeBayes:

| `organism_type` | FreeBayes `--ploidy` |
|-----------------|----------------------|
| `prokaryote` | `1` |
| `eukaryote` | `2` |

**Rationale:** `--ploidy 1` matches haploid prokaryotic variant calling, while `--ploidy 2` matches diploid eukaryotic calling. The remaining thresholds are standard community thresholds for calling SNPs and small indels with moderate Illumina coverage (30–100×). Each value departs from the freebayes default to reduce false positives — particularly important in a regulatory/safety context (GMO assessment). The freebayes defaults are intentionally permissive to support diverse use cases (e.g., low-frequency somatic variants); for GMO detection pipelines these permissive defaults would produce excessive noise.

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
cuteSV <bam> <ref> <out.vcf> <work_dir> -t ${task.cpus} <read-type-specific args>
```

| Parameter | Value | Description |
|-----------|-------|-------------|
| `-t` | `${task.cpus}` | Number of threads assigned by Nextflow for this task (bounded by process-level CPU limits in `nextflow.config`). |

The pipeline also applies a read-type-specific cuteSV profile based on the validated long-read type:

| Read type | Additional cuteSV parameters |
|-----------|------------------------------|
| `pacbio-clr` | `--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5` |
| `pacbio-hifi` | `--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5` |
| `ont` | `--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3` |

Any other cuteSV parameters use tool defaults.


### Sniffles (structural variant calling)

**Module:** `modules/sv_calling.nf`

All parameters use Sniffles defaults. No custom thresholds are applied.

### DeBreak (structural variant calling)

**Module:** `modules/sv_calling.nf`

```bash
debreak --bam <bam> -o <out_dir> -t ${task.cpus} --rescue_large_ins --rescue_dup --poa --ref <ref>
```

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--bam` | `<bam>` | Sorted and indexed long-read alignment BAM from the mapping workflow. |
| `-o` | `<out_dir>` | Output directory for DeBreak results. The pipeline writes to `debreak_out`. |
| `-t` | `${task.cpus}` | Number of threads assigned by Nextflow for this task (bounded by process-level CPU limits in `nextflow.config`). |
| `--rescue_large_ins` | *(flag)* | Enables DeBreak rescue logic for large insertion calls. |
| `--rescue_dup` | *(flag)* | Enables DeBreak rescue logic for duplication calls. |
| `--poa` | *(flag)* | Enables partial order alignment refinement. |
| `--ref` | `<ref>` | Reference FASTA used for DeBreak calling and refinement. |

Any other DeBreak parameters use tool defaults.

---

## Reference vs Modified Genome Comparison

### NUCmer (whole-genome alignment)

**Module:** `modules/assembly.nf`
**Package note:** `nucmer` is a command from the **MUMmer** package. Versioning should be documented as the MUMmer package version (for this pipeline: MUMmer `4.0.1`), not as a separate `nucmer` version.

```bash
nucmer --maxmatch -c 100 -b 500 -l 50 <ref> <mod> -p <prefix>
```

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--maxmatch` | *(flag)* | Use all anchor matches, not just unique ones. Required to detect duplicated and repetitive regions. |
| `-c` | 100 | Minimum alignment cluster length in bp. Filters spurious short alignments. |
| `-b` | 500 | Maximum gap (bp) between clustered matches. Allows merging across small indels or rearrangements. |
| `-l` | 50 | Minimum exact match (MEM) length in bp. Balances sensitivity with noise reduction. |

**Rationale:** These values are suitable for pairwise comparison of closely related genomes (reference vs. genetically modified organism) for both prokaryotic and eukaryotic organisms. The settings are permissive enough to capture meaningful structural differences while filtering alignment noise.

**Trade-offs:**

- *Lowering* `-c` or `-l` increases sensitivity to small structural variants but may introduce spurious alignments.
- *Increasing* `-b` captures larger structural events but risks merging non-contiguous alignments.

### Delta-filter

**Module:** `modules/assembly.nf`
**Package note:** `delta-filter` is also a command from the **MUMmer** package and follows the same package version (MUMmer `4.0.1`).

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
| Freebayes | Short-read | `variant_calling.nf` | `--ploidy 1` for `prokaryote` or `--ploidy 2` for `eukaryote`, `--min-coverage 10`, `--min-base-quality 20`, `--min-mapping-quality 30`, `--min-alternate-count 3` |
| Delly | Short-read | `sv_calling.nf` | Defaults only |
| cuteSV | Long-read | `sv_calling.nf` | `-t ${task.cpus}` |
| Sniffles | Long-read | `sv_calling.nf` | Defaults only |
| DeBreak | Long-read | `sv_calling.nf` | `-t ${task.cpus}` |
| NUCmer | Ref vs Mod | `assembly.nf` | `--maxmatch`, `-c 100`, `-b 500`, `-l 50` |
| Delta-filter | Ref vs Mod | `assembly.nf` | `-m`, `-i 90`, `-l 100` |
