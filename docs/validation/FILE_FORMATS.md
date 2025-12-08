# Supported File Formats

The EFSA Pipeline validation module supports various file formats for genomic data input.

## Genome Files

| Supported Input Formats | Final Output Format |
|------------------------|---------------------|
| FASTA: `.fasta`, `.fa`, `.fna` | `.fasta` |
| GenBank: `.gb`, `.gbk`, `.genbank` | `.fasta` |
| Compression: `.gz`, `.bz2`, `.gzip`, `.bzip2` | Uncompressed |

### Details
- GenBank files are automatically converted to FASTA format
- Compressed files are decompressed during validation
- Output is always uncompressed FASTA format

## Read Files

| Supported Input Formats | Final Output Format |
|------------------------|---------------------|
| FASTQ: `.fastq`, `.fq` | `.fastq` |
| BAM: `.bam` (limited support) | `.bam` |
| Compression: `.gz`, `.bz2`, `.gzip`, `.bzip2` | `.gz` |

### Details
- FASTQ files remain in FASTQ format
- Compressed reads are kept compressed (using gzip)
- BAM files have limited support and require special handling

## Feature Files

| Supported Input Formats | Final Output Format |
|------------------------|---------------------|
| GFF: `.gff`, `.gff3`, `.gtf` | `.gff3` |
| BED: `.bed` | `.gff3` |
| Compression: `.gz`, `.bz2`, `.gzip`, `.bzip2` | Uncompressed |

### Details
- All feature formats are converted to GFF3
- BED files are automatically converted to GFF3 format
- Output is always uncompressed

## File Organization

Input files should be placed in:
```
data/inputs/
```

After validation, files are organized in:
```
data/valid/
├── assembled_genome.fasta
├── reference_genome.fasta
├── ref_plasmid.fa
├── mod_plasmid.fa
├── ref_feature.gff
├── illumina/
├── ont/
└── pacbio/
```

## See Also

- [Configuration Guide](CONFIG_GUIDE.md) - How to specify file paths
- [Validation Overview](OVERVIEW.md) - Validation process details
