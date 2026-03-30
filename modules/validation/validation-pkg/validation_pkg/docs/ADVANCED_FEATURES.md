# Advanced Features

Documentation for advanced features in the validation_pkg package.

## Table of Contents

- [Plasmid Handling](#plasmid-handling)
- [Format Conversion](#format-conversion)
- [Sequence ID Management](#sequence-id-management)
- [Paired-End Detection](#paired-end-detection)
- [Coordinate System Conversion](#coordinate-system-conversion)
- [Parallel Compression](#parallel-compression)

---

## Plasmid Handling

Advanced plasmid detection and processing for bacterial genomes.

### Automatic Plasmid Detection

GenomeValidator automatically identifies plasmids:

```python
from validation_pkg import GenomeValidator

# Default: longest sequence = chromosome, others = plasmids
settings = GenomeValidator.Settings(main_longest=True)
result = validate_genome(config.ref_genome, settings)

print(f"Chromosome: {result.longest_sequence_id} ({result.longest_sequence_length} bp)")
print(f"Plasmids: {result.plasmid_count}")
```

**Detection Logic:**
- `main_longest=True` (default): Longest sequence is main chromosome
- `main_first=True`: First sequence is main chromosome
- All other sequences are treated as plasmids

### Split Plasmids to Separate Files

```python
settings = GenomeValidator.Settings(plasmid_split=True)
result = validate_genome(config.ref_genome, settings)

# Output files:
# - reference_genome.fasta (main chromosome)
# - ref_plasmid_1.fasta (first plasmid)
# - ref_plasmid_2.fasta (second plasmid)

print(f"Main: {result.output_file}")
print(f"Plasmids: {result.plasmid_filenames}")
```

**Use Case:** Separate analysis of chromosome and plasmids.

### Merge Plasmids to One File

```python
settings = GenomeValidator.Settings(plasmids_to_one=True)
result = validate_genome(config.ref_genome, settings)

# Output files:
# - reference_genome.fasta (chromosome only)
# - ref_plasmid.fasta (all plasmids merged)

print(f"Main: {result.output_file}")
print(f"Merged plasmids: {result.plasmid_filenames[0]}")
```

**Use Case:** Simplify file structure while separating chromosomes from plasmids.

### Plasmid-Only Files

For files containing only plasmids (no chromosome):

```python
settings = GenomeValidator.Settings(is_plasmid=True)
result = validate_genome(config.plasmid_config, settings)

# All sequences treated as plasmids
# No "main" chromosome selection
```

**Behavior:**
- If `plasmid_split=True`: Each sequence → separate file
- If `plasmids_to_one=True`: All sequences → one file
- If both False: All sequences → main output file

### Example: Complete Plasmid Workflow

```python
from validation_pkg import GenomeValidator

# Genome with chromosome + 2 plasmids
settings = GenomeValidator.Settings(
    plasmid_split=True,
    min_sequence_length=1000,  # Filter small sequences
    replace_id_with='seq'       # Rename sequences
)

result = validate_genome(config.ref_genome, settings)

# Output:
# - reference_genome.fasta: seq (chromosome)
# - ref_plasmid_1.fasta: seq1 (plasmid 1)
# - ref_plasmid_2.fasta: seq2 (plasmid 2)
```

---

## Format Conversion

Automatic format conversion between file types.

### GenBank → FASTA

GenomeValidator automatically converts GenBank to FASTA:

```python
# Input: genome.gbk (GenBank)
result = validate_genome(config.ref_genome)
# Output: genome.fasta (FASTA)
```

**Conversion Details:**
- Sequence data extracted
- Sequence IDs preserved (LOCUS name)
- Descriptions from DEFINITION line
- Annotations discarded (not in FASTA)

**Example GenBank:**
```
LOCUS       NC_000913          4641652 bp    DNA     circular BCT
DEFINITION  Escherichia coli str. K-12 substr. MG1655, complete genome.
//
```

**Output FASTA:**
```
>NC_000913 Escherichia coli str. K-12 substr. MG1655, complete genome.
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGT...
```

### BED → GFF3

FeatureValidator automatically converts BED to GFF3:

```python
# Input: features.bed
result = validate_feature(config.ref_feature)
# Output: features.gff3
```

**Coordinate Conversion:**
- BED: 0-based, half-open `[start, end)`
- GFF: 1-based, closed `[start, end]`

**Example Conversion:**
```
# Input BED (0-based):
chr1  100  200  gene1  500  +

# Output GFF3 (1-based):
chr1  .  .  101  200  500  +  .  ID=gene1
```

### BAM → FASTQ

ReadValidator supports BAM to FASTQ conversion, but **BAM files are ignored by default** (`ignore_bam=True`). To enable conversion, set `ignore_bam=False` and `keep_bam=False` (or `keep_bam=True` to also copy the original BAM to output):

```python
settings = ReadValidator.Settings(ignore_bam=False, keep_bam=False)
validator = ReadValidator(config.reads[0], settings)
result = validator.run()
# Output: reads.fastq.gz (converted from BAM)
```

**Requirements:**
- `pysam` installed or `samtools` in PATH
- BAM must contain sequence and quality scores

**Limitations:**
- Secondary/supplementary alignments are skipped
- Paired-end info may be lost
- Use with caution for production workflows

---

## Sequence ID Management

Advanced sequence ID modification and tracking.

### Replace All Sequence IDs

```python
settings = GenomeValidator.Settings(replace_id_with='chr')
result = validate_genome(config.ref_genome, settings)

# Input IDs: NC_000913.3, plasmid_pA
# Output IDs: chr, chr1
```

**Auto-Increment Behavior:**
- Single sequence: Exact value used
- Multiple sequences: Increments applied
  - First: `'chr'`
  - Second: `'chr1'`
  - Third: `'chr2'`
  - etc.

### Track Original IDs

Original IDs are stored in the description field of each sequence record:

```python
settings = GenomeValidator.Settings(replace_id_with='chromosome')
result = validate_genome(config.ref_genome, settings)

# Output FASTA:
# >chromosome NC_000913.3
# AGCTTTTCATTCTGACTGCAA...
```

The description field contains the original sequence ID. This information is available in the FASTA output but is not specially formatted with any prefix.

**Use Case:** Standardize IDs while keeping the original identifier visible in the FASTA description.

### Feature Sequence ID Replacement

Replace chromosome names in feature files (strict mode only):

```python
settings = FeatureValidator.Settings(replace_id_with='chr1')
result = validate_feature(config.ref_feature, settings)

# Input GFF:
# NC_000913.3  source  gene  100  200  .  +  .  ID=gene1
#
# Output GFF:
# chr1  source  gene  100  200  .  +  .  ID=gene1
```

**Note:** The original seqname is not preserved in the output. This operation replaces all values in column 1 and is only applied in strict mode.

---

## Paired-End Detection

Automatic detection of Illumina paired-end read patterns.

### Supported Patterns

ReadValidator detects these patterns:

| Pattern | Example Filenames |
|---------|------------------|
| Lane numbers | `sample_R1_001.fastq`, `sample_R2_001.fastq` |
| Simple numbers | `sample_1.fastq`, `sample_2.fastq` |
| Standard Illumina | `sample_R1.fastq`, `sample_R2.fastq` |
| Dot separator | `sample.R1.fastq`, `sample.R2.fastq` |
| No separator | `sampleR1.fastq`, `sampleR2.fastq` |
| With suffix | `sample_R1_suffix.fastq`, `sample_R2_suffix.fastq` |

### Access Pattern Information

```python
result = validate_read(config.reads[0])

if result.read_number:
    print(f"Paired-end detected:")
    print(f"  Base name: {result.base_name}")
    print(f"  Read: R{result.read_number}")
else:
    print("Single-end or no pattern detected")
```

**Base Name Extraction:**
```python
# File: sample_R1_001.fastq.gz
# base_name: 'sample'
# read_number: 1

# File: experiment_1_suffix.fastq.gz
# base_name: 'experiment'
# read_number: 1
```

### Validate Paired-End Completeness

Use inter-file validation to check R1↔R2 matching:

```python
from validation_pkg import readxread_validation, ReadXReadSettings

reads_results = validate_reads(config.reads)

settings = ReadXReadSettings(pair_end_basename=True)
check = readxread_validation(reads_results, settings)

if check['passed']:
    print("✓ All paired-end reads complete")
else:
    print(f"✗ Missing pairs: {check['metadata']['missing_r1']}")
```

---

## Coordinate System Conversion

BED to GFF coordinate transformation.

### Understanding Coordinate Systems

**BED Format (0-based, half-open):**
- Start: 0-based index (first position = 0)
- End: Exclusive (not included in feature)
- Example: `chr1  100  200` = positions 100-199

**GFF Format (1-based, closed):**
- Start: 1-based index (first position = 1)
- End: Inclusive (included in feature)
- Example: `chr1  source  type  101  200` = positions 101-200

### Automatic Conversion

FeatureValidator automatically converts:

```python
# Input BED:
chr1  100  200  gene1  500  +

# Conversion applied:
# BED start (0-based): 100 → GFF start (1-based): 101
# BED end (half-open): 200 → GFF end (closed): 200

# Output GFF3:
chr1  .  .  101  200  500  +  .  ID=gene1
```

### Manual Verification

```python
result = validate_feature(config.ref_feature)

# Input was BED format
# Output is GFF3 format with corrected coordinates
print(f"Converted to: {result.output_file}")
```

**Important:** Always verify coordinate systems when comparing features across tools.

---

## Parallel Compression

Automatic parallel compression for performance.

### Automatic Detection

Validators automatically use parallel compression tools if available:

```bash
# Install parallel compression tools
sudo apt-get install pigz pbzip2  # Ubuntu/Debian
brew install pigz pbzip2          # macOS
```

**Performance Gains:**
- `pigz`: 2-4x faster than `gzip`
- `pbzip2`: 3-6x faster than `bzip2`

### Usage

No code changes required:

```python
# Automatically uses pigz if available
settings = GenomeValidator.Settings(coding_type='gz')
result = validate_genome(config.ref_genome, settings)

# Falls back to standard gzip if pigz not found
```

### Thread Configuration

Control compression threads:

```python
# Via config file
{
  "options": {"threads": 16}
}

# Or programmatically
config.options['threads'] = 16

# Compression will use 16 threads
result = validate_genome(config.ref_genome)
```

### Compression Comparison

| Tool | Compression | Decompression | Best For |
|------|------------|---------------|----------|
| **gzip/pigz** | Good | Fast | General purpose, reads |
| **bzip2/pbzip2** | Better | Slower | Archival, space-critical |
| **None** | - | - | Temporary files, speed |

**Recommendation:** Use `gzip` (pigz) for reads, uncompressed for genomes/features.

---

## See Also

- [VALIDATORS.md](VALIDATORS.md) - Validator documentation
- [SETTINGS.md](SETTINGS.md) - Settings reference
- [EXAMPLES.md](EXAMPLES.md) - Code examples
- [API_REFERENCE.md](API_REFERENCE.md) - API documentation
