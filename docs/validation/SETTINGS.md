# Validation Settings

The validation module supports three validation levels to balance thoroughness and performance.

## Validation Levels

### Comparison Table

| Level | Parsing | Validation | Edits | Output | Speed | Use Case |
|-------|---------|------------|-------|--------|-------|----------|
| **strict** | All data | All data | All applied | BioPython write | Slowest | Structure validation sequence by sequence, statistics gathering |
| **trust** (default) | All data (genome)<br>First record only (reads) | First sequence only | All applied (genome, features)<br>None (reads) | BioPython write (genome)<br>File copy (reads) | Fast | Trust data, adapt file coding, name and location |
| **minimal** | None | None | None | File copy | Fastest | Rename and move files to meet the requirements |

## Level Details

### Strict Mode
- Validates every record in the file
- Performs comprehensive quality checks
- Generates detailed statistics
- Feature files: runs coordinate validation in parallel (when `threads > 1` and file has ≥ 1000 features); falls back to direct parsing if `gffread` produces no output; fails validation if both paths return 0 features
- Recommended for first-time data processing

**Use when:**
- Processing data for the first time
- Data quality is uncertain

### Trust Mode (Default)
- Validates only the first record (reads) or a small sample (features)
- Assumes remaining data is consistent
- Faster processing for large files
- Still performs format conversions

**Use when:**
- Data has been pre-validated
- Files are from trusted sources
- You need faster processing times

### Minimal Mode
- No validation performed
- Only renames and moves files
- Fastest processing
- Use with caution

**Use when:**
- Files are already in correct format
- You're absolutely certain of data quality
- Only file organization is needed

## Configuration

Set the validation level in your `config.json`:

```json
{
  "ref_genome_filename": {
    "filename": "reference.fasta",
    "validation_level": "strict",
    "threads": 8
  }
}
```

Or set globally in options:

```json
{
  "options": {
    "validation_level": "trust",
    "threads": 8
  }
}
```

## Performance Impact

Relative speed uses `strict` as the baseline (slowest).

| Mode | Relative Speed | Resource Usage |
|------|---------------|----------------|
| strict | 1x (baseline) | High |
| trust | 10-15x faster | Moderate |
| minimal | 100x+ faster | Minimal |

## See Also

- [Performance Tips](PERFORMANCE.md) - Optimization strategies
- [Configuration Guide](CONFIG_GUIDE.md) - Full configuration options
